//! Murnaghan–Nakayama evaluator: character columns by composing sparse
//! "remove one rim hook of length ℓ" level operators (design doc
//! `02-numerics` §2).
//!
//! For a target class `ν = (l₁ ≥ l₂ ≥ … ≥ l_k)` the full character column
//! over all `ρ ⊢ n` is
//!
//! ```text
//! X[·, ν] = M_{l₁} ∘ M_{l₂} ∘ … ∘ M_{l_k} (e_∅)
//! ```
//!
//! applied smallest part first: the intermediate vector after applying the
//! ascending prefix `(l_k, …, l_j)` is the character column of that class in
//! the smaller symmetric group — well-defined and shared between targets
//! with a common ascending prefix (the suffix-sharing trie exploited by the
//! full-table generator in P2; single columns here compose the same
//! operators).
//!
//! Operators are materialized once per `(from_level, ℓ)` as CSR in *gather*
//! orientation — for each target partition, the list of `(source, sign)`
//! pairs coming from its removable ℓ-hooks — so applying an operator has no
//! partition manipulation in the loop, only indexed signed adds.
#![deny(clippy::float_arithmetic)]

use std::collections::HashMap;
use std::sync::{Arc, Mutex};

use num_traits::{One, Zero};

use crate::arith::ExactInt;
use crate::chars::rimhook::BetaSet;
use crate::partition::{canonical_cmp, partitions_in_canonical_order, Partition};

/// Partitions of one level `m ≤ n` in canonical order, with rank lookup.
struct LevelIndex {
    partitions: Vec<Partition>,
}

impl LevelIndex {
    fn build(m: u16) -> Self {
        Self {
            partitions: partitions_in_canonical_order(m),
        }
    }

    fn rank(&self, p: &Partition) -> u32 {
        self.partitions
            .binary_search_by(|probe| canonical_cmp(probe, p))
            .expect("partition belongs to this level") as u32
    }

    fn count(&self) -> usize {
        self.partitions.len()
    }
}

/// Sparse level operator `M_ℓ : level m → level m+ℓ` in gather orientation.
struct CsrOp {
    /// `row_ptr[t]..row_ptr[t+1]` indexes `entries` for target rank `t`.
    row_ptr: Vec<u32>,
    /// `(source rank at level m, MN sign ±1)`.
    entries: Vec<(u32, i8)>,
}

impl CsrOp {
    fn build(from: &LevelIndex, to: &LevelIndex, l: u8, slots: u16) -> Self {
        let mut row_ptr = Vec::with_capacity(to.count() + 1);
        let mut entries = Vec::new();
        row_ptr.push(0u32);
        for target in &to.partitions {
            BetaSet::of(target, slots).for_each_hook_removal(l, |beta, leg| {
                let source = from.rank(&beta.to_partition());
                let sign = if leg % 2 == 0 { 1 } else { -1 };
                entries.push((source, sign));
            });
            row_ptr.push(entries.len() as u32);
        }
        Self { row_ptr, entries }
    }

    fn apply_exact(&self, src: &[ExactInt]) -> Vec<ExactInt> {
        let targets = self.row_ptr.len() - 1;
        let mut out = Vec::with_capacity(targets);
        for t in 0..targets {
            let mut acc = ExactInt::zero();
            for &(source, sign) in
                &self.entries[self.row_ptr[t] as usize..self.row_ptr[t + 1] as usize]
            {
                let v = &src[source as usize];
                if sign > 0 {
                    acc += v;
                } else {
                    acc -= v;
                }
            }
            out.push(acc);
        }
        out
    }
}

/// The evaluator for one `n`: level indexes plus a lazily built, cached set
/// of CSR operators. Cheap to construct; operators materialize on first use.
pub struct MnEvaluator {
    n: u16,
    levels: Vec<LevelIndex>, // m = 0..=n
    ops: Mutex<HashMap<(u16, u8), Arc<CsrOp>>>,
}

impl MnEvaluator {
    pub fn new(n: u16) -> Self {
        assert!(n <= 255);
        let levels = (0..=n).map(LevelIndex::build).collect();
        Self {
            n,
            levels,
            ops: Mutex::new(HashMap::new()),
        }
    }

    pub fn n(&self) -> u16 {
        self.n
    }

    fn op(&self, from_level: u16, l: u8) -> Arc<CsrOp> {
        let key = (from_level, l);
        let mut ops = self.ops.lock().expect("ops mutex poisoned");
        if let Some(op) = ops.get(&key) {
            return Arc::clone(op);
        }
        let op = Arc::new(CsrOp::build(
            &self.levels[from_level as usize],
            &self.levels[(from_level + l as u16) as usize],
            l,
            self.n, // uniform slot count: enough for every partition of any level
        ));
        ops.insert(key, Arc::clone(&op));
        op
    }

    /// Exact character column: `χ^ρ(ν)` for every `ρ ⊢ n` in canonical order.
    pub fn column_exact(&self, nu: &Partition) -> Vec<ExactInt> {
        assert_eq!(nu.n(), self.n, "target class must partition n");
        let mut v = vec![ExactInt::one()]; // e_∅ at level 0
        let mut level = 0u16;
        // apply smallest part first
        for &l in nu.parts().iter().rev() {
            let op = self.op(level, l);
            v = op.apply_exact(&v);
            level += l as u16;
        }
        debug_assert_eq!(v.len(), self.levels[self.n as usize].count());
        v
    }

    /// Exact single value `χ^ρ(ν)` (tests/fixtures convenience — computes the
    /// whole column).
    pub fn value_exact(&self, rho: &Partition, nu: &Partition) -> ExactInt {
        assert_eq!(rho.n(), self.n);
        let column = self.column_exact(nu);
        let rank = self.levels[self.n as usize].rank(rho);
        column[rank as usize].clone()
    }

    /// Exact full table, `table[nu_id][rho_id]` in canonical order both ways
    /// (columns by target class). Intended for small `n` (tests, fixtures).
    pub fn full_table_exact(&self) -> Vec<Vec<ExactInt>> {
        self.levels[self.n as usize]
            .partitions
            .iter()
            .map(|nu| self.column_exact(nu))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chars::degrees::degrees;
    use crate::partition::{PartitionId, PartitionIndex};
    use num_bigint::BigInt;

    #[test]
    fn base_cases() {
        let ev = MnEvaluator::new(0);
        assert_eq!(ev.column_exact(&Partition::identity(0)), vec![ExactInt::one()]);

        let ev = MnEvaluator::new(1);
        assert_eq!(ev.column_exact(&Partition::identity(1)), vec![ExactInt::one()]);
    }

    #[test]
    fn known_table_n3() {
        // rows (canonical): [3], [2,1], [1^3]; columns: same order.
        // Standard S_3 character table (rows = irreps, columns = classes):
        //   χ^{(3)}    = (1, 1, 1)      on classes ([3], [2,1], [1^3])
        //   χ^{(2,1)}  = (-1, 0, 2)
        //   χ^{(1^3)}  = (1, -1, 1)
        let ev = MnEvaluator::new(3);
        let idx = PartitionIndex::build(3).unwrap();
        let table: Vec<Vec<i64>> = ev
            .full_table_exact()
            .iter()
            .map(|col| col.iter().map(|v| i64::try_from(v).unwrap()).collect())
            .collect();
        // table[nu][rho]
        assert_eq!(table[0], vec![1, -1, 1]); // ν = [3]
        assert_eq!(table[1], vec![1, 0, -1]); // ν = [2,1]
        assert_eq!(table[2], vec![1, 2, 1]); // ν = [1^3] (degrees)
        let _ = idx;
    }

    #[test]
    fn identity_column_equals_hook_degrees() {
        for n in [1u16, 4, 7, 10, 12] {
            let ev = MnEvaluator::new(n);
            let idx = PartitionIndex::build(n).unwrap();
            let col = ev.column_exact(&Partition::identity(n));
            let degs = degrees(&idx);
            for i in 0..idx.count() {
                assert_eq!(
                    col[i],
                    BigInt::from(degs[i].clone()),
                    "n={n}, rho={:?}",
                    idx.partition(i as PartitionId)
                );
            }
        }
    }

    #[test]
    fn trivial_and_sign_rows() {
        for n in [2u16, 5, 8, 11] {
            let ev = MnEvaluator::new(n);
            let idx = PartitionIndex::build(n).unwrap();
            for (nu_id, nu) in idx.partitions().iter().enumerate() {
                let col = ev.column_exact(nu);
                // ρ = [n] is canonical index 0: trivial character
                assert_eq!(col[0], ExactInt::one(), "trivial at {nu:?}");
                // ρ = [1^n] is the last index: sign character
                assert_eq!(
                    col[idx.count() - 1],
                    ExactInt::from(idx.sign(nu_id as PartitionId)),
                    "sign at {nu:?}"
                );
            }
        }
    }

    #[test]
    fn n_cycle_column_supported_on_hooks_with_alternating_signs() {
        for n in [3u16, 5, 6, 9] {
            let ev = MnEvaluator::new(n);
            let idx = PartitionIndex::build(n).unwrap();
            let col = ev.column_exact(&Partition::new(vec![n as u8]));
            for (i, value) in col.iter().enumerate() {
                let rho = idx.partition(i as PartitionId);
                // hook shape (n−k, 1^k): first part a, then all ones
                let parts = rho.parts();
                let is_hook =
                    parts.len() == 1 || parts[1..].iter().all(|&p| p == 1);
                if is_hook {
                    let k = parts.len() - 1; // leg of the whole hook
                    let expected = if k % 2 == 0 { 1i64 } else { -1 };
                    assert_eq!(*value, ExactInt::from(expected), "n={n}, rho={rho:?}");
                } else {
                    assert_eq!(*value, ExactInt::zero(), "n={n}, rho={rho:?}");
                }
            }
        }
    }

    #[test]
    fn transpose_sign_relation_full_table() {
        // χ^{ρ'}(ν) = sgn(ν)·χ^ρ(ν)  (spec §11.1)
        for n in [3u16, 6, 9] {
            let ev = MnEvaluator::new(n);
            let idx = PartitionIndex::build(n).unwrap();
            let table = ev.full_table_exact();
            for (nu_id, column) in table.iter().enumerate() {
                let sgn = ExactInt::from(idx.sign(nu_id as PartitionId));
                for rho_id in 0..idx.count() {
                    let t = idx.transpose_id(rho_id as PartitionId) as usize;
                    assert_eq!(
                        column[t],
                        &sgn * &column[rho_id],
                        "n={n}, nu={nu_id}, rho={rho_id}"
                    );
                }
            }
        }
    }

    #[test]
    fn self_transpose_rows_vanish_on_odd_classes() {
        for n in [4u16, 7, 10] {
            let ev = MnEvaluator::new(n);
            let idx = PartitionIndex::build(n).unwrap();
            let table = ev.full_table_exact();
            for rho_id in 0..idx.count() {
                if idx.transpose_id(rho_id as PartitionId) != rho_id as PartitionId {
                    continue;
                }
                for (nu_id, column) in table.iter().enumerate() {
                    if idx.sign(nu_id as PartitionId) == -1 {
                        assert_eq!(
                            column[rho_id],
                            ExactInt::zero(),
                            "n={n}, self-transpose rho={rho_id}, odd nu={nu_id}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn row_and_column_orthogonality() {
        for n in [2u16, 4, 6, 8] {
            let ev = MnEvaluator::new(n);
            let idx = PartitionIndex::build(n).unwrap();
            let table = ev.full_table_exact(); // [nu][rho]
            let q = idx.count();
            let factorial = ExactInt::from(idx.factorial_n().clone());
            // rows: Σ_ν |C_ν| χ^ρ χ^σ = n!·δ
            for rho in 0..q {
                for sigma in rho..q {
                    let sum: ExactInt = (0..q)
                        .map(|nu| {
                            ExactInt::from(idx.class_size(nu as PartitionId).clone())
                                * &table[nu][rho]
                                * &table[nu][sigma]
                        })
                        .sum();
                    let expected = if rho == sigma {
                        factorial.clone()
                    } else {
                        ExactInt::zero()
                    };
                    assert_eq!(sum, expected, "n={n}, rows {rho},{sigma}");
                }
            }
            // columns: Σ_ρ χ^ρ(μ) χ^ρ(ν) = z_ν·δ
            for mu in 0..q {
                for nu in mu..q {
                    let sum: ExactInt = (0..q)
                        .map(|rho| &table[mu][rho] * &table[nu][rho])
                        .sum();
                    let expected = if mu == nu {
                        ExactInt::from(idx.z_value(mu as PartitionId).clone())
                    } else {
                        ExactInt::zero()
                    };
                    assert_eq!(sum, expected, "n={n}, cols {mu},{nu}");
                }
            }
        }
    }
}
