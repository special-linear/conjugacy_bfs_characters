//! Brute-force cross-validation over the full group (spec §9.7, §22.3):
//! BFS distances, exact-length supports, and word counts, aggregated by
//! cycle type and compared with the character engine.

use classdiam_core::arith::ExactInt;
use classdiam_core::chars::MnEvaluator;
use classdiam_core::engine::exact::{run_exact, ExactTransform};
use classdiam_core::partition::{PartitionId, PartitionIndex};
use classdiam_core::spectra::BaseSpectra;
use classdiam_core::testing::bruteforce as bf;
use classdiam_core::testing::catalog::{brute_force_affordable, resolve_entry, union_catalog};
use num_traits::Zero;

/// BFS over raw permutations agrees with the character engine on distances
/// (asserted constant per class first), diameter, reachability, and the
/// per-radius NEW sets.
#[test]
fn bfs_distances_and_layers_match_engine() {
    let mut compared = 0;
    for entry in union_catalog().into_iter().filter(|e| e.n <= 8) {
        let index = PartitionIndex::build(entry.n).unwrap();
        let union = resolve_entry(&index, &entry);
        if !brute_force_affordable(&index, &union) {
            continue;
        }
        let mn = MnEvaluator::new(entry.n);
        let run = run_exact(&index, &mn, &union).unwrap();

        let generators = bf::materialize_union(&index, &union);
        let dist = bf::bfs_distances(entry.n, &generators);
        let by_type = bf::distances_by_type(&index, &dist);

        assert_eq!(run.distance, by_type, "{}: distances", entry.label);

        // layer NEW sets are exactly the distance classes
        for layer in &run.layers {
            let expected: Vec<PartitionId> = (0..index.count())
                .filter(|&nu| by_type[nu] == layer.r as i32)
                .map(|nu| nu as PartitionId)
                .collect();
            assert_eq!(layer.new, expected, "{}: new at r={}", entry.label, layer.r);
        }
        compared += 1;
    }
    assert!(
        compared >= 80,
        "too few catalog entries compared: {compared}"
    );
}

/// Exact-length supports from the set-product DP match the engine's recorded
/// layer supports (spec §5.1: supports are not frontiers), and no new type
/// appears within two radii past the stop.
#[test]
fn exact_length_supports_match_engine() {
    let mut compared = 0;
    for entry in union_catalog().into_iter().filter(|e| e.n <= 7) {
        let index = PartitionIndex::build(entry.n).unwrap();
        let union = resolve_entry(&index, &entry);
        if !brute_force_affordable(&index, &union) {
            continue;
        }
        let mn = MnEvaluator::new(entry.n);
        let run = run_exact(&index, &mn, &union).unwrap();

        let generators = bf::materialize_union(&index, &union);
        let supports = bf::exact_length_supports(&index, &generators, run.stop_radius + 2);

        for layer in &run.layers {
            assert_eq!(
                layer.support, supports[layer.r as usize],
                "{}: support at r={}",
                entry.label, layer.r
            );
        }
        for r in (run.stop_radius + 1)..=(run.stop_radius + 2) {
            for &nu in &supports[r as usize] {
                assert!(
                    run.distance[nu as usize] >= 0,
                    "{}: type {nu} in DP support at r={r} was never visited",
                    entry.label
                );
            }
        }
        compared += 1;
    }
    assert!(
        compared >= 80,
        "too few catalog entries compared: {compared}"
    );
}

/// Internal word counts `a_r(ν)` (per fixed element, spec §23 F2) match the
/// group-algebra DP for r ≤ 4 — the counting layer everything rests on, even
/// though counts never reach the user-facing output.
#[test]
fn word_counts_match_group_algebra_dp() {
    let mut compared = 0;
    for entry in union_catalog().into_iter().filter(|e| e.n <= 7) {
        let index = PartitionIndex::build(entry.n).unwrap();
        let union = resolve_entry(&index, &entry);
        if !brute_force_affordable(&index, &union) {
            continue;
        }
        let mn = MnEvaluator::new(entry.n);
        let generators = bf::materialize_union(&index, &union);
        let max_r = 4u32;
        let counts = bf::word_counts(entry.n, &generators, max_r);
        let by_type = bf::word_counts_by_type(&index, &counts);

        let spectra = BaseSpectra::build(&index, &mn, &union.class_ids).unwrap();
        let theta = spectra.theta(&union.class_ids);
        let table = mn.full_table_exact();
        let mut transform =
            ExactTransform::new(&table, spectra.degrees(), theta, index.factorial_n());

        for r in 0..=max_r {
            if r > 0 {
                transform.advance();
            }
            let coefficients = transform.coefficients().unwrap();
            for nu in 0..index.count() {
                assert_eq!(
                    coefficients[nu],
                    ExactInt::from(by_type[r as usize][nu]),
                    "{}: a_{r}({nu})",
                    entry.label
                );
            }
        }
        compared += 1;
    }
    assert!(
        compared >= 60,
        "too few catalog entries compared: {compared}"
    );
}

/// Mixed products (spec §16.3, §22.4): the spectral coefficient of
/// `K_A·K_B` matches direct enumeration of ordered pairs, and the square
/// expansion `(K_A + K_B)² = K_A² + 2·K_A·K_B + K_B²` holds at the
/// coefficient level through two independent computation routes.
#[test]
fn mixed_products_and_union_square() {
    for entry in union_catalog()
        .into_iter()
        .filter(|e| e.n <= 6 && e.templates.len() == 2)
    {
        let index = PartitionIndex::build(entry.n).unwrap();
        let union = resolve_entry(&index, &entry);
        if !brute_force_affordable(&index, &union) || union.class_ids.len() != 2 {
            continue;
        }
        let mn = MnEvaluator::new(entry.n);
        let (a, b) = (union.class_ids[0], union.class_ids[1]);
        let spectra = BaseSpectra::build(&index, &mn, &[a, b]).unwrap();
        let table = mn.full_table_exact();
        let factorial = ExactInt::from(index.factorial_n().clone());
        let q = index.count();

        // spectral coefficient vector for the product K_A·K_B:
        // (1/n!)·Σ_ρ f_ρ·χ^ρ(ν)·ω_ρ(A)·ω_ρ(B)
        let coefficient_of = |powers: &dyn Fn(usize) -> ExactInt| -> Vec<ExactInt> {
            (0..q)
                .map(|nu| {
                    let numerator: ExactInt = (0..q)
                        .map(|rho| {
                            ExactInt::from(spectra.degrees()[rho].clone())
                                * &table[nu][rho]
                                * powers(rho)
                        })
                        .sum();
                    let (quot, rem) = num_integer::Integer::div_rem(&numerator, &factorial);
                    assert!(rem.is_zero(), "mixed product not divisible by n!");
                    quot
                })
                .collect()
        };
        let om_a = spectra.omega_column(0);
        let om_b = spectra.omega_column(1);
        let ab = coefficient_of(&|rho| &om_a[rho] * &om_b[rho]);
        let aa = coefficient_of(&|rho| &om_a[rho] * &om_a[rho]);
        let bb = coefficient_of(&|rho| &om_b[rho] * &om_b[rho]);
        let uu = coefficient_of(&|rho| {
            let t = &om_a[rho] + &om_b[rho];
            &t * &t
        });

        // 1) direct enumeration of ordered pairs A × B
        let gens_a: Vec<_> = bf::materialize_union(
            &index,
            &classdiam_core::spectra::ResolvedUnion {
                class_ids: vec![a],
                union_size: index.class_size(a).clone(),
                parity: union.parity,
                includes_identity: false,
            },
        );
        let gens_b: Vec<_> = bf::materialize_union(
            &index,
            &classdiam_core::spectra::ResolvedUnion {
                class_ids: vec![b],
                union_size: index.class_size(b).clone(),
                parity: union.parity,
                includes_identity: false,
            },
        );
        let mut pair_counts: std::collections::HashMap<Vec<u8>, u64> =
            std::collections::HashMap::new();
        for ga in &gens_a {
            for gb in &gens_b {
                *pair_counts.entry(bf::compose(ga, gb)).or_insert(0) += 1;
            }
        }
        // aggregate by type with constancy assertion
        let mut by_type: Vec<Option<u64>> = vec![None; q];
        let facts: Vec<u64> = {
            let mut f = vec![1u64; entry.n as usize + 1];
            for k in 1..=entry.n as usize {
                f[k] = f[k - 1] * k as u64;
            }
            f
        };
        let total = facts[entry.n as usize];
        for r in 0..total {
            let p = bf::unrank(r, entry.n, &facts);
            let count = pair_counts.get(&p).copied().unwrap_or(0);
            let t = index.id_of(&bf::cycle_type(&p)).unwrap() as usize;
            match by_type[t] {
                None => by_type[t] = Some(count),
                Some(existing) => {
                    assert_eq!(existing, count, "{}: A·B count not constant", entry.label)
                }
            }
        }
        for nu in 0..q {
            assert_eq!(
                ab[nu],
                ExactInt::from(by_type[nu].unwrap()),
                "{}: K_A·K_B at nu={nu}",
                entry.label
            );
        }

        // 2) square expansion at the coefficient level
        for nu in 0..q {
            let rhs = &aa[nu] + ExactInt::from(2) * &ab[nu] + &bb[nu];
            assert_eq!(uu[nu], rhs, "{}: square expansion at nu={nu}", entry.label);
        }
    }
}
