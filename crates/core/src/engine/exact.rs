//! Exact big-integer reference engine (spec §8.1, §12.1).
//!
//! Simplest correct implementation: full-row character table, exact `θ^r`
//! powers, exact division by `n!`, every spec-§9 invariant checked on every
//! radius. This is the oracle the modular engine (P2) is validated against;
//! it is deliberately unoptimized.
#![deny(clippy::float_arithmetic)]

use fixedbitset::FixedBitSet;
use num_bigint::BigUint;
use num_traits::{One, Zero};

use crate::arith::{exact_div_checked, ExactInt};
use crate::chars::MnEvaluator;
use crate::engine::{LayerRecord, StoppingRule, UnionRun};
use crate::error::ClassdiamError;
use crate::partition::{PartitionId, PartitionIndex};
use crate::spectra::{BaseSpectra, ResolvedUnion};

/// Exact spectral calculator: produces the full coefficient vector
/// `a_r(ν) = (1/n!)·Σ_ρ f_ρ·χ^ρ(ν)·θ_ρ^r` at the current radius.
///
/// Uses the FULL irreducible row set, so `r = 0` correctly reconstructs the
/// identity indicator (spec §9.1) — the restricted-row optimization and its
/// radius-0 exception belong to the modular engine.
pub struct ExactTransform<'a> {
    /// `table[ν][ρ]`, both in canonical order.
    table: &'a [Vec<ExactInt>],
    degrees: Vec<ExactInt>,
    theta: Vec<ExactInt>,
    /// `θ_ρ^radius`.
    power: Vec<ExactInt>,
    factorial: ExactInt,
    radius: u32,
}

impl<'a> ExactTransform<'a> {
    pub fn new(
        table: &'a [Vec<ExactInt>],
        degrees: &[BigUint],
        theta: Vec<ExactInt>,
        factorial: &BigUint,
    ) -> Self {
        let q = degrees.len();
        assert_eq!(table.len(), q);
        assert_eq!(theta.len(), q);
        Self {
            table,
            degrees: degrees.iter().map(|d| ExactInt::from(d.clone())).collect(),
            theta,
            power: vec![ExactInt::one(); q],
            factorial: ExactInt::from(factorial.clone()),
            radius: 0,
        }
    }

    pub fn radius(&self) -> u32 {
        self.radius
    }

    /// Advance to the next radius: `power *= θ`.
    pub fn advance(&mut self) {
        for (p, t) in self.power.iter_mut().zip(self.theta.iter()) {
            *p *= t;
        }
        self.radius += 1;
    }

    /// `a_radius(ν)` for every ν, with exact-divisibility and nonnegativity
    /// checks (spec §9.4, §9.5).
    pub fn coefficients(&self) -> Result<Vec<ExactInt>, ClassdiamError> {
        let q = self.degrees.len();
        let mut out = Vec::with_capacity(q);
        for nu in 0..q {
            let mut numerator = ExactInt::zero();
            for rho in 0..q {
                let chi = &self.table[nu][rho];
                if chi.is_zero() || self.power[rho].is_zero() {
                    continue;
                }
                numerator += &self.degrees[rho] * chi * &self.power[rho];
            }
            let coefficient = exact_div_checked(&numerator, &self.factorial).ok_or(
                ClassdiamError::NotDivisibleByFactorial {
                    radius: self.radius,
                    target: nu,
                },
            )?;
            if coefficient.sign() == num_bigint::Sign::Minus {
                return Err(ClassdiamError::NegativeCoefficient {
                    radius: self.radius,
                    target: nu,
                });
            }
            out.push(coefficient);
        }
        Ok(out)
    }
}

/// Run the exact reference engine for one union.
///
/// Always-on validation (this is the oracle, correctness beats speed):
/// radius-0 identity indicator, radius-1 union indicator, per-radius
/// word-count identity, divisibility, nonnegativity. Stopping: empty layer
/// (primary, spec §5.2) or parity-feasible cover (early exit); for `n ≥ 5`
/// the normal-closure prediction is asserted as a cross-check, never used
/// for termination.
pub fn run_exact(
    index: &PartitionIndex,
    mn: &MnEvaluator,
    union: &ResolvedUnion,
) -> Result<UnionRun, ClassdiamError> {
    let n = index.n();
    let q = index.count();
    let identity = index.identity_id();

    let spectra = BaseSpectra::build(index, mn, &union.class_ids)?;
    let theta = spectra.theta(&union.class_ids);
    let table = mn.full_table_exact();
    let mut transform = ExactTransform::new(&table, spectra.degrees(), theta, index.factorial_n());

    // spec §9.1: radius 0 with the full row set must be the identity indicator.
    let a0 = transform.coefficients()?;
    for (nu, coefficient) in a0.iter().enumerate() {
        let expected_one = nu as PartitionId == identity;
        debug_assert_eq!(
            coefficient.is_one(),
            expected_one,
            "radius-0 identity reconstruction failed at nu={nu}"
        );
        debug_assert_eq!(coefficient.is_zero(), !expected_one);
    }

    let mut distance = vec![-1i32; q];
    let mut first_hit = [vec![-1i32; q], vec![-1i32; q]]; // [even, odd]
    let mut visited = FixedBitSet::with_capacity(q);
    distance[identity as usize] = 0;
    first_hit[0][identity as usize] = 0;
    visited.insert(identity as usize);
    let mut layers = vec![LayerRecord {
        r: 0,
        new: vec![identity],
        support: vec![identity],
    }];

    let feasible = super::parity_feasible_set(index, union.parity);

    let union_size = ExactInt::from(union.union_size.clone());
    let mut union_size_pow = ExactInt::one();
    let radius_limit = 4 * u32::from(n).max(4);

    let (stop_radius, stopping) = loop {
        // Early exit BEFORE the next transform when nothing new can appear.
        if visited.count_ones(..) == feasible.count_ones(..) {
            break (
                layers.last().expect("layer 0 exists").r,
                StoppingRule::AllTypesVisited,
            );
        }

        transform.advance();
        let r = transform.radius();
        if r > radius_limit {
            return Err(ClassdiamError::RadiusLimitExceeded {
                n,
                limit: radius_limit,
            });
        }
        let coefficients = transform.coefficients()?;

        // spec §9.3 word-count identity, every radius.
        union_size_pow *= &union_size;
        let mut word_count = ExactInt::zero();
        for (nu, coefficient) in coefficients.iter().enumerate() {
            if !coefficient.is_zero() {
                word_count +=
                    ExactInt::from(index.class_size(nu as PartitionId).clone()) * coefficient;
            }
        }
        if word_count != union_size_pow {
            return Err(ClassdiamError::WordCountMismatch { radius: r });
        }

        // spec §9.2: radius 1 must be the indicator of the union's classes.
        if r == 1 {
            for (nu, coefficient) in coefficients.iter().enumerate() {
                let expected = union.class_ids.contains(&(nu as PartitionId));
                debug_assert_eq!(
                    coefficient.is_one(),
                    expected,
                    "radius-1 indicator failed at nu={nu}"
                );
            }
        }

        let parity_slot = (r % 2) as usize;
        let mut support = Vec::new();
        let mut new = Vec::new();
        for (nu, coefficient) in coefficients.iter().enumerate() {
            if coefficient.is_zero() {
                continue;
            }
            support.push(nu as PartitionId);
            if first_hit[parity_slot][nu] < 0 {
                first_hit[parity_slot][nu] = r as i32;
            }
            if !visited.contains(nu) {
                visited.insert(nu);
                distance[nu] = r as i32;
                new.push(nu as PartitionId);
            }
        }
        let layer_is_empty = new.is_empty();
        layers.push(LayerRecord { r, new, support });

        if layer_is_empty {
            break (r, StoppingRule::EmptyLayer); // spec §5.2, on exact supports
        }
    };

    // Cross-check (never used for termination): for n ≥ 5, a union containing
    // a non-identity class generates A_n or S_n, so the visited set must be
    // exactly the parity-feasible set.
    if n >= 5 && union.class_ids.iter().any(|&c| c != identity) {
        debug_assert_eq!(
            visited.count_ones(..),
            feasible.count_ones(..),
            "normal-closure prediction violated (n={n})"
        );
    }

    let diameter = distance.iter().copied().max().unwrap_or(0).max(0) as u32;
    let reachable_count = visited.count_ones(..);
    let [fh_even, fh_odd] = first_hit;
    Ok(UnionRun {
        n,
        distance,
        first_hit_even: fh_even,
        first_hit_odd: fh_odd,
        layers,
        diameter,
        stop_radius,
        stopping,
        reachable_count,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::partition::CycleTypeTemplate;
    use crate::spectra::resolve_union;

    fn run(n: u16, templates: &[&str]) -> UnionRun {
        let index = PartitionIndex::build(n).unwrap();
        let mn = MnEvaluator::new(n);
        let templates: Vec<CycleTypeTemplate> =
            templates.iter().map(|s| s.parse().unwrap()).collect();
        let union = resolve_union(&index, &templates, false).unwrap();
        run_exact(&index, &mn, &union).unwrap()
    }

    /// The critique-verified worked example: S_6 with transpositions
    /// (design doc 01 §9.3). d(ν) = 6 − ℓ(ν), diameter 5.
    #[test]
    fn s6_transpositions_match_worked_example() {
        let r = run(6, &["2"]);
        assert_eq!(r.distance, vec![5, 4, 4, 3, 4, 3, 2, 3, 2, 1, 0]);
        assert_eq!(r.first_hit_even, vec![-1, 4, 4, -1, 4, -1, 2, -1, 2, -1, 0]);
        assert_eq!(r.first_hit_odd, vec![5, -1, -1, 3, -1, 3, -1, 3, -1, 1, -1]);
        assert_eq!(r.diameter, 5);
        assert_eq!(r.stop_radius, 5);
        assert_eq!(r.stopping, StoppingRule::AllTypesVisited);
        assert_eq!(r.reachable_count, 11);

        let expected_layers: Vec<(u32, Vec<u32>, Vec<u32>)> = vec![
            (0, vec![10], vec![10]),
            (1, vec![9], vec![9]),
            (2, vec![6, 8], vec![6, 8, 10]),
            (3, vec![3, 5, 7], vec![3, 5, 7, 9]),
            (4, vec![1, 2, 4], vec![1, 2, 4, 6, 8, 10]),
            (5, vec![0], vec![0, 3, 5, 7, 9]),
        ];
        assert_eq!(r.layers.len(), expected_layers.len());
        for (layer, (er, enew, esupport)) in r.layers.iter().zip(expected_layers) {
            assert_eq!(layer.r, er);
            assert_eq!(layer.new, enew, "r={er}");
            assert_eq!(layer.support, esupport, "r={er}");
        }
    }

    /// Transpositions: d(ν) = n − ℓ(ν) exactly, diameter n − 1 (analytic
    /// anchor valid at every n).
    #[test]
    fn transpositions_closed_form() {
        for n in [3u16, 4, 5, 7, 8] {
            let index = PartitionIndex::build(n).unwrap();
            let r = run(n, &["2"]);
            for nu in 0..index.count() {
                let expected = n as i32 - index.partition(nu as PartitionId).len() as i32;
                assert_eq!(r.distance[nu], expected, "n={n}, nu={nu}");
            }
            assert_eq!(r.diameter, u32::from(n) - 1);
            assert_eq!(r.reachable_count, index.count());
        }
    }

    /// [2,2] in S_4 generates the Klein four-group: only {id, [2,2]}
    /// reachable, diameter 1, stops via the empty-layer rule (the parity
    /// upper bound — all even types — is NOT attained; spec §2.2's
    /// disconnected case).
    #[test]
    fn s4_double_transpositions_disconnected() {
        let index = PartitionIndex::build(4).unwrap();
        let r = run(4, &["2,2"]);
        let type_22 = index
            .id_of(&crate::partition::Partition::new(vec![2u8, 2]))
            .unwrap() as usize;
        let identity = index.identity_id() as usize;
        for nu in 0..index.count() {
            let expected = if nu == identity {
                0
            } else if nu == type_22 {
                1
            } else {
                -1
            };
            assert_eq!(r.distance[nu], expected, "nu={nu}");
        }
        assert_eq!(r.diameter, 1);
        assert_eq!(r.reachable_count, 2);
        assert_eq!(r.stopping, StoppingRule::EmptyLayer);
    }

    /// 3-cycles in S_5: reachable set = all even types (A_5), diameter 2.
    #[test]
    fn s5_three_cycles_reach_alternating_group() {
        let index = PartitionIndex::build(5).unwrap();
        let r = run(5, &["3"]);
        for nu in 0..index.count() {
            let even = index.sign(nu as PartitionId) == 1;
            assert_eq!(r.distance[nu] >= 0, even, "nu={nu}");
        }
        assert_eq!(r.diameter, 2);
        assert_eq!(r.stopping, StoppingRule::AllTypesVisited);
    }

    /// Mixed-parity union: both parities occur; everything reachable fast.
    #[test]
    fn s5_mixed_union() {
        let index = PartitionIndex::build(5).unwrap();
        let r = run(5, &["2", "3"]);
        assert_eq!(r.reachable_count, index.count());
        // mixed unions have no parity filter: some type must be hit at both
        // an even and an odd radius within the run
        let both =
            (0..index.count()).any(|nu| r.first_hit_even[nu] >= 0 && r.first_hit_odd[nu] >= 0);
        assert!(both);
    }
}
