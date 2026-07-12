//! Spec §9 invariants, parametrized over the deterministic union catalog.
//!
//! The exact engine also asserts most of these internally on every run; the
//! tests here exercise them explicitly (including beyond the stopping
//! radius) and check structural consistency of the reported results.

use classdiam_core::arith::ExactInt;
use classdiam_core::chars::MnEvaluator;
use classdiam_core::engine::exact::{run_exact, ExactTransform};
use classdiam_core::engine::StoppingRule;
use classdiam_core::partition::{PartitionId, PartitionIndex};
use classdiam_core::spectra::{BaseSpectra, UnionParity};
use classdiam_core::testing::catalog::{resolve_entry, union_catalog};
use num_traits::{One, Zero};

/// Every catalog entry runs to completion — which exercises the engine's
/// always-on checks (radius-0/1 indicators, per-radius word-count identity,
/// divisibility, nonnegativity) — and the reported structure is coherent.
#[test]
fn catalog_runs_with_coherent_results() {
    for entry in union_catalog() {
        let index = PartitionIndex::build(entry.n).unwrap();
        let mn = MnEvaluator::new(entry.n);
        let union = resolve_entry(&index, &entry);
        let run = run_exact(&index, &mn, &union).unwrap_or_else(|e| panic!("{}: {e}", entry.label));

        let identity = index.identity_id() as usize;
        assert_eq!(run.distance[identity], 0, "{}", entry.label);
        assert_eq!(run.layers[0].new, vec![identity as PartitionId]);

        // distance = min over the two parity chains' first hits
        for nu in 0..index.count() {
            let hits = [run.first_hit_even[nu], run.first_hit_odd[nu]];
            let min_hit = hits.iter().filter(|&&h| h >= 0).min().copied();
            assert_eq!(
                run.distance[nu],
                min_hit.unwrap_or(-1),
                "{}: distance vs first hits at nu={nu}",
                entry.label
            );
        }

        // diameter = max distance; reachable count matches
        let max_d = run.distance.iter().copied().max().unwrap();
        assert_eq!(run.diameter as i32, max_d.max(0), "{}", entry.label);
        assert_eq!(
            run.reachable_count,
            run.distance.iter().filter(|&&d| d >= 0).count(),
            "{}",
            entry.label
        );

        // recorded exact-length supports match the first-hit reconstruction
        // rule for every recorded radius (valid within the run span)
        for layer in &run.layers {
            let r = layer.r;
            let chain = if r % 2 == 0 {
                &run.first_hit_even
            } else {
                &run.first_hit_odd
            };
            let reconstructed: Vec<PartitionId> = (0..index.count())
                .filter(|&nu| chain[nu] >= 0 && chain[nu] <= r as i32)
                .map(|nu| nu as PartitionId)
                .collect();
            assert_eq!(
                layer.support, reconstructed,
                "{}: support reconstruction at r={r}",
                entry.label
            );
        }

        // parity filter (spec §9.6): single-parity unions only reach types of
        // sign ε^r at radius r
        match union.parity {
            UnionParity::Even => {
                for layer in &run.layers {
                    for &nu in &layer.support {
                        assert_eq!(index.sign(nu), 1, "{}: even union", entry.label);
                    }
                }
            }
            UnionParity::Odd => {
                for layer in &run.layers {
                    let expected = if layer.r % 2 == 0 { 1 } else { -1 };
                    for &nu in &layer.support {
                        assert_eq!(
                            index.sign(nu),
                            expected,
                            "{}: odd union at r={}",
                            entry.label,
                            layer.r
                        );
                    }
                }
            }
            UnionParity::Mixed => {}
        }
    }
}

/// Radius-0 (identity indicator, FULL rows) and radius-1 (union indicator)
/// checked explicitly through the transform (spec §9.1, §9.2).
#[test]
fn radius0_identity_and_radius1_indicator() {
    for entry in union_catalog().into_iter().filter(|e| e.n <= 8) {
        let index = PartitionIndex::build(entry.n).unwrap();
        let mn = MnEvaluator::new(entry.n);
        let union = resolve_entry(&index, &entry);
        let spectra = BaseSpectra::build(&index, &mn, &union.class_ids).unwrap();
        let theta = spectra.theta(&union.class_ids);
        let table = mn.full_table_exact();
        let mut transform =
            ExactTransform::new(&table, spectra.degrees(), theta, index.factorial_n());

        let a0 = transform.coefficients().unwrap();
        for (nu, c) in a0.iter().enumerate() {
            if nu as PartitionId == index.identity_id() {
                assert!(c.is_one(), "{}: a_0 at identity", entry.label);
            } else {
                assert!(c.is_zero(), "{}: a_0 at nu={nu}", entry.label);
            }
        }

        transform.advance();
        let a1 = transform.coefficients().unwrap();
        for (nu, c) in a1.iter().enumerate() {
            if union.class_ids.contains(&(nu as PartitionId)) {
                assert!(c.is_one(), "{}: a_1 at member {nu}", entry.label);
            } else {
                assert!(c.is_zero(), "{}: a_1 at nu={nu}", entry.label);
            }
        }
    }
}

/// Word-count identity holds beyond the stopping radius too (spec §9.3),
/// and no new type ever appears after an engine stop (spec §5.2 validated
/// spectrally two radii past the stop).
#[test]
fn word_count_and_no_new_types_beyond_stop() {
    for entry in union_catalog().into_iter().filter(|e| e.n <= 8) {
        let index = PartitionIndex::build(entry.n).unwrap();
        let mn = MnEvaluator::new(entry.n);
        let union = resolve_entry(&index, &entry);
        let run = run_exact(&index, &mn, &union).unwrap();

        let spectra = BaseSpectra::build(&index, &mn, &union.class_ids).unwrap();
        let theta = spectra.theta(&union.class_ids);
        let table = mn.full_table_exact();
        let mut transform =
            ExactTransform::new(&table, spectra.degrees(), theta, index.factorial_n());

        let union_size = ExactInt::from(union.union_size.clone());
        let mut expected_words = ExactInt::one();
        for r in 1..=run.stop_radius + 2 {
            transform.advance();
            expected_words *= &union_size;
            let coefficients = transform.coefficients().unwrap();
            let total: ExactInt = coefficients
                .iter()
                .enumerate()
                .map(|(nu, c)| ExactInt::from(index.class_size(nu as PartitionId).clone()) * c)
                .sum();
            assert_eq!(
                total, expected_words,
                "{}: word count at r={r}",
                entry.label
            );

            // beyond the stop, the support must stay inside the visited set
            if r > run.stop_radius {
                for (nu, c) in coefficients.iter().enumerate() {
                    if !c.is_zero() {
                        assert!(
                            run.distance[nu] >= 0,
                            "{}: NEW type {nu} appeared at r={r} after stop {} — stopping rule broken",
                            entry.label,
                            run.stop_radius
                        );
                    }
                }
            }
        }
    }
}

/// Mixed-parity negative test (spec §9.6/F5): some radius reaches both
/// parities at once, so no single-parity filter may be applied.
#[test]
fn mixed_unions_defeat_parity_filter() {
    let entries: Vec<_> = union_catalog().into_iter().filter(|e| e.n <= 8).collect();
    let mut saw_mixed = 0;
    for entry in entries {
        let index = PartitionIndex::build(entry.n).unwrap();
        let mn = MnEvaluator::new(entry.n);
        let union = resolve_entry(&index, &entry);
        if union.parity != UnionParity::Mixed {
            continue;
        }
        saw_mixed += 1;
        let run = run_exact(&index, &mn, &union).unwrap();
        let some_layer_has_both = run.layers.iter().any(|layer| {
            let mut even = false;
            let mut odd = false;
            for &nu in &layer.support {
                match index.sign(nu) {
                    1 => even = true,
                    _ => odd = true,
                }
            }
            even && odd
        });
        assert!(
            some_layer_has_both,
            "{}: expected mixed-parity layer",
            entry.label
        );
    }
    assert!(saw_mixed > 10, "catalog must contain mixed unions");
}

/// AllTypesVisited early exit and EmptyLayer both occur across the catalog —
/// the stopping logic's two branches are each exercised.
#[test]
fn both_stopping_rules_are_exercised() {
    let mut empty_layer = 0;
    let mut all_visited = 0;
    for entry in union_catalog().into_iter().filter(|e| e.n <= 7) {
        let index = PartitionIndex::build(entry.n).unwrap();
        let mn = MnEvaluator::new(entry.n);
        let union = resolve_entry(&index, &entry);
        match run_exact(&index, &mn, &union).unwrap().stopping {
            StoppingRule::EmptyLayer => empty_layer += 1,
            StoppingRule::AllTypesVisited => all_visited += 1,
        }
    }
    assert!(empty_layer > 0, "no EmptyLayer stop in catalog");
    assert!(all_visited > 0, "no AllTypesVisited stop in catalog");
}
