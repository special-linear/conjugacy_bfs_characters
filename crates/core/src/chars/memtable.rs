//! Resident transpose-paired multi-prime character table (design doc 02
//! §3.3/§4.2): only representative rows (`ρ_id ≤ ρ'_id`) are stored, and
//! target columns are grouped even-sign first, then odd-sign, so every
//! transform call binds exactly one of the `W₊`/`W₋` weight vectors.
#![deny(clippy::float_arithmetic)]

use rayon::prelude::*;

use crate::arith::{ModCtx, Prime31};
use crate::chars::MnEvaluator;
use crate::partition::{PartitionId, PartitionIndex};

pub struct PairedModTable {
    primes: Vec<Prime31>,
    /// Representative rows (canonical order): `ρ_id ≤ transpose(ρ_id)`.
    rep_rows: Vec<PartitionId>,
    /// Targets: even-sign classes in canonical order, then odd-sign.
    targets: Vec<PartitionId>,
    even_count: usize,
    /// `planes[lane][target_pos * rep_count + row_pos]`, fully reduced.
    planes: Vec<Vec<u32>>,
}

impl PairedModTable {
    /// Generate by rayon-parallel MN column evaluation (all operators are
    /// prebuilt first so workers never contend on the operator cache).
    pub fn generate(index: &PartitionIndex, mn: &MnEvaluator, primes: &[Prime31]) -> Self {
        assert_eq!(index.n(), mn.n());
        let ctxs: Vec<ModCtx> = primes.iter().copied().map(ModCtx::new).collect();
        mn.prebuild_all_ops();

        let q = index.count();
        let rep_rows: Vec<PartitionId> = (0..q as u32)
            .filter(|&rho| index.transpose_id(rho) >= rho)
            .collect();
        let mut targets: Vec<PartitionId> =
            (0..q as u32).filter(|&nu| index.sign(nu) == 1).collect();
        let even_count = targets.len();
        targets.extend((0..q as u32).filter(|&nu| index.sign(nu) == -1));

        let rep_count = rep_rows.len();
        let columns: Vec<Vec<Vec<u32>>> = targets
            .par_iter()
            .map(|&nu| {
                let full = mn.column_mod(index.partition(nu), &ctxs);
                full.iter()
                    .map(|lane| {
                        rep_rows
                            .iter()
                            .map(|&rho| lane[rho as usize])
                            .collect::<Vec<u32>>()
                    })
                    .collect()
            })
            .collect();

        let mut planes: Vec<Vec<u32>> = ctxs
            .iter()
            .map(|_| Vec::with_capacity(rep_count * targets.len()))
            .collect();
        for column in &columns {
            for (lane, values) in column.iter().enumerate() {
                planes[lane].extend_from_slice(values);
            }
        }

        Self {
            primes: primes.to_vec(),
            rep_rows,
            targets,
            even_count,
            planes,
        }
    }

    pub fn primes(&self) -> &[Prime31] {
        &self.primes
    }

    pub fn rep_rows(&self) -> &[PartitionId] {
        &self.rep_rows
    }

    pub fn rep_count(&self) -> usize {
        self.rep_rows.len()
    }

    /// Targets in table order (even-sign block, then odd-sign block).
    pub fn targets(&self) -> &[PartitionId] {
        &self.targets
    }

    pub fn even_count(&self) -> usize {
        self.even_count
    }

    /// Residue column for one target position and prime lane.
    pub fn column(&self, lane: usize, target_pos: usize) -> &[u32] {
        let r = self.rep_count();
        &self.planes[lane][target_pos * r..(target_pos + 1) * r]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arith::screening_primes;
    use num_bigint::BigInt;

    /// Every stored residue equals the exact character value reduced —
    /// including the row restriction and the parity-grouped column order.
    #[test]
    fn table_matches_exact_values() {
        let n = 9u16;
        let index = PartitionIndex::build(n).unwrap();
        let mn = MnEvaluator::new(n);
        let primes = {
            let mut p = screening_primes(2);
            p.push(Prime31(11)); // small prime to exercise wraparound
            p
        };
        let table = PairedModTable::generate(&index, &mn, &primes);
        let ctxs: Vec<ModCtx> = primes.iter().copied().map(ModCtx::new).collect();

        // parity grouping
        for (pos, &nu) in table.targets().iter().enumerate() {
            let expected_sign = if pos < table.even_count() { 1 } else { -1 };
            assert_eq!(index.sign(nu), expected_sign);
        }
        // rep rows: exactly one of {rho, rho'} with rho <= rho'
        for &rho in table.rep_rows() {
            assert!(index.transpose_id(rho) >= rho);
        }
        assert_eq!(
            table.rep_count(),
            (index.count()
                + (0..index.count() as u32)
                    .filter(|&i| index.transpose_id(i) == i)
                    .count())
                / 2
        );

        for (pos, &nu) in table.targets().iter().enumerate() {
            let exact: Vec<BigInt> = mn.column_exact(index.partition(nu));
            for (lane, ctx) in ctxs.iter().enumerate() {
                let column = table.column(lane, pos);
                for (row_pos, &rho) in table.rep_rows().iter().enumerate() {
                    assert_eq!(
                        column[row_pos],
                        ctx.reduce_bigint(&exact[rho as usize]),
                        "nu={nu}, rho={rho}, p={}",
                        ctx.prime().0
                    );
                }
            }
        }
    }
}
