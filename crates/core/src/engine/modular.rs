//! Production modular engine: screening on 31-bit primes with a rigorous,
//! synchronous per-radius certification gate (merged design; spec §12,
//! Failures 4/9).
//!
//! Classification of a target at radius `r`:
//! - any nonzero residue        → POSITIVE (proof: `a_r(ν)` is a nonnegative
//!   integer and `p ∤ n!`);
//! - same-parity first hit set  → POSITIVE by the `a_{r+2} ≥ a_r` lemma
//!   (counted in a metric; NEVER asserted on residues — a positive can
//!   vanish mod every prime, critique finding 2);
//! - otherwise                  → CANDIDATE, certified before layer commit:
//!   - T1 bound-zero:   `a_r(ν) ≤ ⌊|U|^r/|C_ν|⌋ = 0`;
//!   - T2 resident CRT: all residues zero and `∏ primes > ⌊|U|^r/|C_ν|⌋`
//!       determine `a_r(ν) = 0` (spec §12.3);
//!   - T3 exact:        big-integer evaluation over representative rows
//!       (cached exact character columns, incremental exact powers).
//!   Every tier terminates with a rigorous verdict; no probabilistic step.
//!
//! The stopping test runs only on fully certified supports; layer commit
//! and stopping logic mirror the exact reference engine bit for bit.
#![deny(clippy::float_arithmetic)]

use std::collections::HashMap;
use std::sync::Arc;

use fixedbitset::FixedBitSet;
use num_bigint::BigUint;
use num_traits::{One, Zero};
use serde::{Deserialize, Serialize};

use crate::arith::{exact_div_checked, ExactInt, ModCtx, Prime31};
use crate::chars::memtable::PairedModTable;
use crate::chars::MnEvaluator;
use crate::engine::{LayerRecord, StoppingRule, UnionRun};
use crate::error::ClassdiamError;
use crate::partition::{PartitionId, PartitionIndex};
use crate::spectra::{BaseSpectra, PairedSpectrum, ResolvedUnion, UnionParity};
use crate::transform::{RadiusWeights, TransformBackend};

/// Per-`n` shared state: the paired table plus precomputed modular data.
/// Built once, reused by every union at this `n` (spec §6).
pub struct ModularContext {
    pub primes: Vec<Prime31>,
    pub ctxs: Vec<ModCtx>,
    pub table: PairedModTable,
    /// `f_ρ` per representative row (exact + per-lane residues).
    rep_degrees: Vec<BigUint>,
    rep_degrees_mod: Vec<Vec<u32>>,
    /// `|C_ν|` residues in TABLE target order, per lane.
    class_size_mod: Vec<Vec<u32>>,
    factorial_mod: Vec<u32>,
    prime_product: BigUint,
}

impl ModularContext {
    pub fn build(index: &PartitionIndex, mn: &MnEvaluator, primes: &[Prime31]) -> Self {
        let ctxs: Vec<ModCtx> = primes.iter().copied().map(ModCtx::new).collect();
        let table = PairedModTable::generate(index, mn, primes);
        let rep_degrees: Vec<BigUint> = table
            .rep_rows()
            .iter()
            .map(|&rho| crate::chars::degree(index.partition(rho)))
            .collect();
        let rep_degrees_mod: Vec<Vec<u32>> = ctxs
            .iter()
            .map(|ctx| rep_degrees.iter().map(|d| ctx.reduce_biguint(d)).collect())
            .collect();
        let class_size_mod: Vec<Vec<u32>> = ctxs
            .iter()
            .map(|ctx| {
                table
                    .targets()
                    .iter()
                    .map(|&nu| ctx.reduce_biguint(index.class_size(nu)))
                    .collect()
            })
            .collect();
        let factorial_mod: Vec<u32> = ctxs
            .iter()
            .map(|ctx| ctx.reduce_biguint(index.factorial_n()))
            .collect();
        let prime_product = primes
            .iter()
            .fold(BigUint::one(), |acc, p| acc * BigUint::from(p.0));
        Self {
            primes: primes.to_vec(),
            ctxs,
            table,
            rep_degrees,
            rep_degrees_mod,
            class_size_mod,
            factorial_mod,
            prime_product,
        }
    }
}

/// Certification audit trail — the output's proof that Failures 4/9 cannot
/// have occurred silently.
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct CertificationStats {
    /// Targets with all-zero residues needing a verdict.
    pub candidates: u64,
    /// Certified zero by the `|U|^r/|C_ν|` bound (tier 1).
    pub bound_certified: u64,
    /// Certified zero by resident CRT (tier 2).
    pub crt_resident_certified: u64,
    /// Exact big-integer evaluations (tier 3).
    pub exact_evals: u64,
    /// Tier-3 evaluations that came back POSITIVE — genuinely hidden
    /// positives (expected ≈ 0 with production primes).
    pub hidden_positives: u64,
    /// Positives granted by the `a_{r+2} ≥ a_r` lemma without residues.
    pub lemma_positives: u64,
}

pub struct ModularOptions {
    /// Diagnostic radius bound multiplier (default 4·n).
    pub radius_limit_factor: u32,
}

impl Default for ModularOptions {
    fn default() -> Self {
        Self {
            radius_limit_factor: 4,
        }
    }
}

pub fn run_modular(
    index: &PartitionIndex,
    mn: &MnEvaluator,
    ctx: &ModularContext,
    spectra: &BaseSpectra,
    union: &ResolvedUnion,
    backend: &dyn TransformBackend,
    options: &ModularOptions,
) -> Result<(UnionRun, CertificationStats), ClassdiamError> {
    let n = index.n();
    let q = index.count();
    let identity = index.identity_id();
    let lanes = ctx.ctxs.len();
    let paired = PairedSpectrum::build(index, spectra, union);
    assert_eq!(paired.rep_rows, ctx.table.rep_rows(), "row order mismatch");
    let reps = paired.rep_count();

    // θ± residues per lane; modular and exact power state.
    let theta_plus_mod: Vec<Vec<u32>> = ctx
        .ctxs
        .iter()
        .map(|c| {
            paired
                .theta_plus
                .iter()
                .map(|t| c.reduce_bigint(t))
                .collect()
        })
        .collect();
    let theta_minus_mod: Vec<Vec<u32>> = ctx
        .ctxs
        .iter()
        .map(|c| {
            paired
                .theta_minus
                .iter()
                .map(|t| c.reduce_bigint(t))
                .collect()
        })
        .collect();
    let mut p_mod: Vec<Vec<u32>> = vec![vec![1u32; reps]; lanes];
    let mut p_prime_mod: Vec<Vec<u32>> = vec![vec![1u32; reps]; lanes];
    let mut p_exact: Vec<ExactInt> = vec![ExactInt::one(); reps];
    let mut p_prime_exact: Vec<ExactInt> = vec![ExactInt::one(); reps];

    // |U|^r trackers: exact (for bounds) and per-lane residues (tripwire).
    let union_size = union.union_size.clone();
    let mut union_size_pow = BigUint::one();
    let union_size_mod: Vec<u32> = ctx
        .ctxs
        .iter()
        .map(|c| c.reduce_biguint(&union_size))
        .collect();
    let mut union_size_pow_mod: Vec<u32> = vec![1u32; lanes];

    // Exact character columns over representative rows, cached per target.
    let mut column_cache: HashMap<PartitionId, Arc<Vec<ExactInt>>> = HashMap::new();
    let mut column_over_reps = |nu: PartitionId| -> Arc<Vec<ExactInt>> {
        column_cache
            .entry(nu)
            .or_insert_with(|| {
                let full = mn.column_exact(index.partition(nu));
                Arc::new(
                    ctx.table
                        .rep_rows()
                        .iter()
                        .map(|&rho| full[rho as usize].clone())
                        .collect(),
                )
            })
            .clone()
    };

    // BFS state (identical layout to the exact engine).
    let mut distance = vec![-1i32; q];
    let mut first_hit = [vec![-1i32; q], vec![-1i32; q]];
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
    let factorial_exact = ExactInt::from(index.factorial_n().clone());
    let mut stats = CertificationStats::default();
    let radius_limit = options.radius_limit_factor * u32::from(n).max(1);

    let even_count = ctx.table.even_count();
    let total_targets = ctx.table.targets().len();

    let (stop_radius, stopping) = loop {
        if visited.count_ones(..) == feasible.count_ones(..) {
            break (
                layers.last().expect("layer 0").r,
                StoppingRule::AllTypesVisited,
            );
        }

        let r = layers.last().expect("layer 0").r + 1;
        if r > radius_limit {
            return Err(ClassdiamError::RadiusLimitExceeded {
                n,
                limit: radius_limit,
            });
        }

        // Advance powers: P *= θ₊, P' *= θ₋ (modular per lane, plus exact).
        for lane in 0..lanes {
            let c = &ctx.ctxs[lane];
            for i in 0..reps {
                p_mod[lane][i] = c.mul(p_mod[lane][i], theta_plus_mod[lane][i]);
                p_prime_mod[lane][i] = c.mul(p_prime_mod[lane][i], theta_minus_mod[lane][i]);
            }
        }
        for i in 0..reps {
            p_exact[i] *= &paired.theta_plus[i];
            p_prime_exact[i] *= &paired.theta_minus[i];
        }
        union_size_pow *= &union_size;
        for lane in 0..lanes {
            union_size_pow_mod[lane] =
                ctx.ctxs[lane].mul(union_size_pow_mod[lane], union_size_mod[lane]);
        }

        // Assemble W±: paired rows f·(P ± P'), self-transpose rows f·P.
        let weights = RadiusWeights {
            w_plus: (0..lanes)
                .map(|lane| {
                    let c = &ctx.ctxs[lane];
                    (0..reps)
                        .map(|i| {
                            let combined = if paired.is_self[i] {
                                p_mod[lane][i]
                            } else {
                                c.add(p_mod[lane][i], p_prime_mod[lane][i])
                            };
                            c.mul(ctx.rep_degrees_mod[lane][i], combined)
                        })
                        .collect()
                })
                .collect(),
            w_minus: (0..lanes)
                .map(|lane| {
                    let c = &ctx.ctxs[lane];
                    (0..reps)
                        .map(|i| {
                            let combined = if paired.is_self[i] {
                                p_mod[lane][i]
                            } else {
                                c.sub(p_mod[lane][i], p_prime_mod[lane][i])
                            };
                            c.mul(ctx.rep_degrees_mod[lane][i], combined)
                        })
                        .collect()
                })
                .collect(),
        };

        // Which parity blocks can be nonzero at this radius (spec §11.3)?
        let mut ranges: Vec<std::ops::Range<usize>> = Vec::with_capacity(2);
        match union.parity {
            UnionParity::Even => ranges.push(0..even_count),
            UnionParity::Odd => {
                if r % 2 == 0 {
                    ranges.push(0..even_count);
                } else {
                    ranges.push(even_count..total_targets);
                }
            }
            UnionParity::Mixed => {
                ranges.push(0..even_count);
                ranges.push(even_count..total_targets);
            }
        }

        // Exact weight vectors for tier-3, built lazily once per radius.
        let mut exact_weights: Option<(Vec<ExactInt>, Vec<ExactInt>)> = None;

        let parity_slot = (r % 2) as usize;
        let mut support = Vec::new();
        let mut new = Vec::new();
        let mut word_count_acc = vec![0u128; lanes];

        for range in ranges {
            let numerators = backend.numerators(&ctx.table, &weights, range.clone());
            for (offset, pos) in range.clone().enumerate() {
                let nu = ctx.table.targets()[pos];
                // word-count tripwire accumulation (u128-safe: q < 2^18 terms)
                for lane in 0..lanes {
                    word_count_acc[lane] +=
                        ctx.class_size_mod[lane][pos] as u128 * numerators[lane][offset] as u128;
                }

                let any_nonzero = (0..lanes).any(|lane| numerators[lane][offset] != 0);
                let positive = if any_nonzero {
                    true
                } else if first_hit[parity_slot][nu as usize] >= 0 {
                    // a_{r} ≥ a_{r−2} > 0: positive without residues.
                    stats.lemma_positives += 1;
                    true
                } else {
                    // CANDIDATE: certify now, before anything commits.
                    stats.candidates += 1;
                    let bound = &union_size_pow / index.class_size(nu);
                    if bound.is_zero() {
                        stats.bound_certified += 1;
                        false
                    } else if ctx.prime_product > bound {
                        stats.crt_resident_certified += 1;
                        false
                    } else {
                        stats.exact_evals += 1;
                        let (w_plus_exact, w_minus_exact) =
                            exact_weights.get_or_insert_with(|| {
                                build_exact_weights(
                                    &paired,
                                    &ctx.rep_degrees,
                                    &p_exact,
                                    &p_prime_exact,
                                )
                            });
                        let w = if pos < even_count {
                            &*w_plus_exact
                        } else {
                            &*w_minus_exact
                        };
                        let column = column_over_reps(nu);
                        let mut numerator = ExactInt::zero();
                        for (chi, weight) in column.iter().zip(w.iter()) {
                            if !chi.is_zero() && !weight.is_zero() {
                                numerator += chi * weight;
                            }
                        }
                        let a = exact_div_checked(&numerator, &factorial_exact).ok_or(
                            ClassdiamError::NotDivisibleByFactorial {
                                radius: r,
                                target: nu as usize,
                            },
                        )?;
                        if a.sign() == num_bigint::Sign::Minus {
                            return Err(ClassdiamError::NegativeCoefficient {
                                radius: r,
                                target: nu as usize,
                            });
                        }
                        if a.is_zero() {
                            false
                        } else {
                            stats.hidden_positives += 1;
                            true
                        }
                    }
                };

                if positive {
                    support.push(nu);
                    if first_hit[parity_slot][nu as usize] < 0 {
                        first_hit[parity_slot][nu as usize] = r as i32;
                    }
                    if !visited.contains(nu as usize) {
                        visited.insert(nu as usize);
                        distance[nu as usize] = r as i32;
                        new.push(nu);
                    }
                }
            }
        }

        // Word-count tripwire: Σ |C_ν|·N_r(ν) ≡ n!·|U|^r (mod p), per lane.
        // (Skipped-parity targets are exactly zero by the parity theorem.)
        for lane in 0..lanes {
            let c = &ctx.ctxs[lane];
            let lhs = c.reduce_u128(word_count_acc[lane]);
            let rhs = c.mul(ctx.factorial_mod[lane], union_size_pow_mod[lane]);
            if lhs != rhs {
                return Err(ClassdiamError::WordCountMismatch { radius: r });
            }
        }

        // Sort support/new into canonical id order (table order is parity-
        // grouped; reports and the exact engine use canonical order).
        support.sort_unstable();
        new.sort_unstable();
        let layer_is_empty = new.is_empty();
        layers.push(LayerRecord { r, new, support });
        if layer_is_empty {
            break (r, StoppingRule::EmptyLayer);
        }
    };

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
    Ok((
        UnionRun {
            n,
            distance,
            first_hit_even: fh_even,
            first_hit_odd: fh_odd,
            layers,
            diameter,
            stop_radius,
            stopping,
            reachable_count,
        },
        stats,
    ))
}

/// Exact `W±` over representative rows: paired `f·(P ± P')`, self `f·P`.
fn build_exact_weights(
    paired: &PairedSpectrum,
    rep_degrees: &[BigUint],
    p_exact: &[ExactInt],
    p_prime_exact: &[ExactInt],
) -> (Vec<ExactInt>, Vec<ExactInt>) {
    let mut plus = Vec::with_capacity(paired.rep_count());
    let mut minus = Vec::with_capacity(paired.rep_count());
    for i in 0..paired.rep_count() {
        let f = ExactInt::from(rep_degrees[i].clone());
        if paired.is_self[i] {
            let w = &f * &p_exact[i];
            plus.push(w.clone());
            minus.push(w);
        } else {
            plus.push(&f * (&p_exact[i] + &p_prime_exact[i]));
            minus.push(&f * (&p_exact[i] - &p_prime_exact[i]));
        }
    }
    (plus, minus)
}
