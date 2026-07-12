//! Central characters and union spectra (spec §4, §16).
//!
//! `ω_ρ(λ) = |C_λ|·χ^ρ(λ) / f_ρ` — an integer for symmetric groups; the
//! division is performed exactly and a failure is reported as an error
//! (it almost always means an indexing/orientation bug, spec §4).
//! A union's spectrum is the row-wise sum `θ_ρ(U) = Σ_j ω_ρ(λ_j)`.
//!
//! Active rows (spec §7) are computed on EXACT values only — never on
//! residues, which could vanish accidentally (critique-confirmed rule).
#![deny(clippy::float_arithmetic)]

use num_bigint::{BigInt, BigUint};
use num_traits::Zero;
use serde::{Deserialize, Serialize};

use crate::arith::{exact_div_checked, ExactInt};
use crate::chars::{degrees, MnEvaluator};
use crate::error::ClassdiamError;
use crate::partition::{CycleTypeTemplate, PartitionId, PartitionIndex};

/// Parity composition of a generating union (spec §17.1). Drives parity
/// filtering (valid only for single-parity unions, spec §9.6/F5).
#[derive(Copy, Clone, PartialEq, Eq, Debug, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum UnionParity {
    Even,
    Odd,
    Mixed,
}

/// A validated generating union for a concrete `n`: deduplicated full
/// conjugacy classes (identity excluded unless explicitly allowed).
#[derive(Clone, Debug)]
pub struct ResolvedUnion {
    pub class_ids: Vec<PartitionId>,
    pub union_size: BigUint,
    pub parity: UnionParity,
    pub includes_identity: bool,
}

/// Resolve cycle-type templates against a concrete `n`.
///
/// Errors: template does not fit (`TemplateDoesNotFit`), identity class
/// without `allow_identity` (spec §5.3), empty union after validation.
/// Duplicate templates are deduplicated (idempotent unions).
pub fn resolve_union(
    index: &PartitionIndex,
    templates: &[CycleTypeTemplate],
    allow_identity: bool,
) -> Result<ResolvedUnion, ClassdiamError> {
    let n = index.n();
    let mut class_ids: Vec<PartitionId> = Vec::new();
    let mut includes_identity = false;
    for t in templates {
        let padded = t.pad(n)?;
        if padded.is_identity_type() {
            if !allow_identity {
                return Err(ClassdiamError::IdentityGenerator { n });
            }
            includes_identity = true;
        }
        let id = index
            .id_of(&padded)
            .expect("padded template is a partition of n");
        if !class_ids.contains(&id) {
            class_ids.push(id);
        }
    }
    if class_ids.is_empty() {
        return Err(ClassdiamError::EmptyUnion { n });
    }
    class_ids.sort_unstable();

    let union_size: BigUint = class_ids
        .iter()
        .map(|&id| index.class_size(id).clone())
        .sum();

    let non_identity_signs: Vec<i8> = class_ids
        .iter()
        .filter(|&&id| id != index.identity_id())
        .map(|&id| index.sign(id))
        .collect();
    // The identity class is even; for parity classification of the GENERATED
    // walk it behaves like an even generator.
    let all_signs: Vec<i8> = if includes_identity {
        class_ids.iter().map(|&id| index.sign(id)).collect()
    } else {
        non_identity_signs
    };
    let parity = if all_signs.iter().all(|&s| s == 1) {
        UnionParity::Even
    } else if all_signs.iter().all(|&s| s == -1) {
        UnionParity::Odd
    } else {
        UnionParity::Mixed
    };

    Ok(ResolvedUnion {
        class_ids,
        union_size,
        parity,
        includes_identity,
    })
}

/// Exact central-eigenvalue columns `ω_ρ(λ_j)` for a set of base classes,
/// plus degrees — the shared preprocessing of spec §6.
pub struct BaseSpectra {
    base_classes: Vec<PartitionId>,
    /// `omega[j][rho]`, exact, over ALL `q` rows in canonical order.
    omega: Vec<Vec<ExactInt>>,
    degrees: Vec<BigUint>,
}

impl BaseSpectra {
    pub fn build(
        index: &PartitionIndex,
        mn: &MnEvaluator,
        base_classes: &[PartitionId],
    ) -> Result<Self, ClassdiamError> {
        assert_eq!(mn.n(), index.n());
        let degs = degrees(index);
        let mut omega = Vec::with_capacity(base_classes.len());
        for &lambda in base_classes {
            let chi_column = mn.column_exact(index.partition(lambda));
            let class_size = BigInt::from(index.class_size(lambda).clone());
            let mut col = Vec::with_capacity(chi_column.len());
            for (rho, chi) in chi_column.iter().enumerate() {
                let numerator = &class_size * chi;
                let degree = BigInt::from(degs[rho].clone());
                let value = exact_div_checked(&numerator, &degree).ok_or(
                    ClassdiamError::OmegaNotIntegral {
                        rho,
                        lambda: lambda as usize,
                    },
                )?;
                col.push(value);
            }
            omega.push(col);
        }
        Ok(Self {
            base_classes: base_classes.to_vec(),
            omega,
            degrees: degs,
        })
    }

    pub fn base_classes(&self) -> &[PartitionId] {
        &self.base_classes
    }

    pub fn degrees(&self) -> &[BigUint] {
        &self.degrees
    }

    /// `ω` column for one base class (by position in `base_classes`).
    pub fn omega_column(&self, j: usize) -> &[ExactInt] {
        &self.omega[j]
    }

    /// Exact union spectrum `θ_ρ(U) = Σ_j ω_ρ(λ_j)` over the given classes
    /// (must be a subset of the base classes).
    pub fn theta(&self, class_ids: &[PartitionId]) -> Vec<ExactInt> {
        let q = self.degrees.len();
        let mut theta = vec![ExactInt::zero(); q];
        for &id in class_ids {
            let j = self
                .base_classes
                .iter()
                .position(|&b| b == id)
                .expect("class must be among the base classes");
            for (t, w) in theta.iter_mut().zip(self.omega[j].iter()) {
                *t += w;
            }
        }
        theta
    }

    /// Rows `ρ` where some base class has `ω_ρ ≠ 0` (spec §7): the safe
    /// shared active set for every union over these base classes, valid for
    /// radii `r ≥ 1` only (radius-0 exception).
    pub fn active_rows(&self) -> Vec<PartitionId> {
        let q = self.degrees.len();
        (0..q)
            .filter(|&rho| self.omega.iter().any(|col| !col[rho].is_zero()))
            .map(|rho| rho as PartitionId)
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::partition::Partition;

    fn setup(n: u16) -> (PartitionIndex, MnEvaluator) {
        (PartitionIndex::build(n).unwrap(), MnEvaluator::new(n))
    }

    /// Sum of contents of ρ: Σ over cells (j − i), 0-based.
    fn content_sum(p: &Partition) -> i64 {
        let mut sum = 0i64;
        for (i, &row) in p.parts().iter().enumerate() {
            for j in 0..row as i64 {
                sum += j - i as i64;
            }
        }
        sum
    }

    #[test]
    fn omega_integral_for_all_classes_small_n() {
        for n in [2u16, 5, 8, 11] {
            let (idx, mn) = setup(n);
            let all: Vec<PartitionId> = (0..idx.count() as u32).collect();
            // must not error: divisibility holds for every (rho, lambda)
            BaseSpectra::build(&idx, &mn, &all).unwrap();
        }
    }

    #[test]
    fn omega_transposition_equals_content_sum() {
        for n in [2u16, 4, 6, 9, 12] {
            let (idx, mn) = setup(n);
            let transposition = idx
                .id_of(&CycleTypeTemplate::new([2]).unwrap().pad(n).unwrap())
                .unwrap();
            let spectra = BaseSpectra::build(&idx, &mn, &[transposition]).unwrap();
            for rho in 0..idx.count() {
                let expected = ExactInt::from(content_sum(idx.partition(rho as PartitionId)));
                assert_eq!(
                    spectra.omega_column(0)[rho],
                    expected,
                    "n={n}, rho={:?}",
                    idx.partition(rho as PartitionId)
                );
            }
        }
    }

    #[test]
    fn omega_transpose_sign_relation() {
        // ω_{ρ'}(λ) = sgn(λ)·ω_ρ(λ)  (spec §11.1)
        for n in [4u16, 7, 10] {
            let (idx, mn) = setup(n);
            let all: Vec<PartitionId> = (0..idx.count() as u32).collect();
            let spectra = BaseSpectra::build(&idx, &mn, &all).unwrap();
            for (j, &lambda) in spectra.base_classes().iter().enumerate() {
                let sgn = ExactInt::from(idx.sign(lambda));
                for rho in 0..idx.count() {
                    let t = idx.transpose_id(rho as PartitionId) as usize;
                    assert_eq!(
                        spectra.omega_column(j)[t],
                        &sgn * &spectra.omega_column(j)[rho],
                        "n={n}, lambda={lambda}, rho={rho}"
                    );
                }
            }
        }
    }

    #[test]
    fn theta_is_additive_over_disjoint_classes() {
        let (idx, mn) = setup(7);
        let a = idx
            .id_of(&CycleTypeTemplate::new([2]).unwrap().pad(7).unwrap())
            .unwrap();
        let b = idx
            .id_of(&CycleTypeTemplate::new([3]).unwrap().pad(7).unwrap())
            .unwrap();
        let spectra = BaseSpectra::build(&idx, &mn, &[a, b]).unwrap();
        let ta = spectra.theta(&[a]);
        let tb = spectra.theta(&[b]);
        let tu = spectra.theta(&[a, b]);
        for rho in 0..idx.count() {
            assert_eq!(tu[rho], &ta[rho] + &tb[rho]);
        }
    }

    #[test]
    fn active_rows_match_nonzero_chi() {
        // n-cycle: active rows are exactly the hooks (spec §7.1)
        let (idx, mn) = setup(8);
        let cycle = idx.id_of(&Partition::new(vec![8u8])).unwrap();
        let spectra = BaseSpectra::build(&idx, &mn, &[cycle]).unwrap();
        let active = spectra.active_rows();
        let expected: Vec<PartitionId> = (0..idx.count())
            .filter(|&i| {
                let parts = idx.partition(i as PartitionId).parts();
                parts.len() == 1 || parts[1..].iter().all(|&p| p == 1)
            })
            .map(|i| i as PartitionId)
            .collect();
        assert_eq!(active, expected);
    }

    #[test]
    fn omega_identity_class_is_all_ones() {
        let (idx, mn) = setup(6);
        let spectra = BaseSpectra::build(&idx, &mn, &[idx.identity_id()]).unwrap();
        for rho in 0..idx.count() {
            assert_eq!(spectra.omega_column(0)[rho], ExactInt::from(1));
        }
    }

    #[test]
    fn resolve_union_validation() {
        let idx = PartitionIndex::build(6).unwrap();
        let t2 = CycleTypeTemplate::new([2]).unwrap();
        let t3 = CycleTypeTemplate::new([3]).unwrap();
        let t22 = CycleTypeTemplate::new([2, 2]).unwrap();

        let u = resolve_union(&idx, std::slice::from_ref(&t2), false).unwrap();
        assert_eq!(u.parity, UnionParity::Odd);
        assert_eq!(u.union_size, BigUint::from(15u32));
        assert!(!u.includes_identity);

        let u = resolve_union(&idx, std::slice::from_ref(&t3), false).unwrap();
        assert_eq!(u.parity, UnionParity::Even);
        assert_eq!(u.union_size, BigUint::from(40u32));

        let u = resolve_union(&idx, &[t2.clone(), t3.clone()], false).unwrap();
        assert_eq!(u.parity, UnionParity::Mixed);
        assert_eq!(u.union_size, BigUint::from(55u32));
        assert_eq!(u.class_ids.len(), 2);

        // dedup
        let u = resolve_union(&idx, &[t2.clone(), t2.clone()], false).unwrap();
        assert_eq!(u.class_ids.len(), 1);

        // even x even
        let u = resolve_union(&idx, &[t3, t22], false).unwrap();
        assert_eq!(u.parity, UnionParity::Even);

        // identity rejected by default, allowed explicitly
        let id_t = CycleTypeTemplate::identity();
        assert!(matches!(
            resolve_union(&idx, std::slice::from_ref(&id_t), false),
            Err(ClassdiamError::IdentityGenerator { n: 6 })
        ));
        let u = resolve_union(&idx, &[id_t, t2], true).unwrap();
        assert!(u.includes_identity);
        assert_eq!(u.union_size, BigUint::from(16u32));

        // does not fit
        let t_big = CycleTypeTemplate::new([7]).unwrap();
        assert!(matches!(
            resolve_union(&idx, &[t_big], false),
            Err(ClassdiamError::TemplateDoesNotFit { .. })
        ));

        // empty
        assert!(matches!(
            resolve_union(&idx, &[], false),
            Err(ClassdiamError::EmptyUnion { n: 6 })
        ));
    }
}
