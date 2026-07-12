//! Modular-engine validation (design doc 03 §5):
//! - full differential against the exact reference engine over the catalog;
//! - the adversarial injected-prime attack on the certification gate — the
//!   one catastrophic failure mode (zero residue treated as certified zero,
//!   spec §23 F4/F9).

use classdiam_core::arith::Prime31;
use classdiam_core::chars::MnEvaluator;
use classdiam_core::engine::exact::run_exact;
use classdiam_core::engine::modular::{
    run_modular, CertificationStats, ModularContext, ModularOptions,
};
use classdiam_core::partition::{CycleTypeTemplate, PartitionIndex};
use classdiam_core::spectra::{resolve_union, BaseSpectra, ResolvedUnion};
use classdiam_core::testing::catalog::{resolve_entry, union_catalog};
use classdiam_core::transform::cpu::{CpuBlocked, CpuReference};
use classdiam_core::transform::TransformBackend;
use serde_json::Value;

fn run_both(
    n: u16,
    union: &ResolvedUnion,
    primes: &[Prime31],
    backend: &dyn TransformBackend,
) -> (
    classdiam_core::engine::UnionRun,
    classdiam_core::engine::UnionRun,
    CertificationStats,
) {
    let index = PartitionIndex::build(n).unwrap();
    let mn = MnEvaluator::new(n);
    let spectra = BaseSpectra::build(&index, &mn, &union.class_ids).unwrap();
    let exact = run_exact(&index, &mn, union).unwrap();
    let ctx = ModularContext::build(&index, &mn, primes);
    let (modular, stats) = run_modular(
        &index,
        &mn,
        &ctx,
        &spectra,
        union,
        backend,
        &ModularOptions::default(),
    )
    .unwrap();
    (exact, modular, stats)
}

/// The central differential: with production primes, the modular engine
/// reproduces the exact engine's ENTIRE output on every catalog entry —
/// distances, both first-hit arrays, every layer (new + support), stopping
/// rule, diameter.
#[test]
fn modular_engine_matches_exact_engine_full_catalog() {
    let primes = classdiam_core::arith::screening_primes(3);
    let backend = CpuBlocked;
    // one shared context per n
    let mut checked = 0;
    for entry in union_catalog() {
        let index = PartitionIndex::build(entry.n).unwrap();
        let union = resolve_entry(&index, &entry);
        let (exact, modular, stats) = run_both(entry.n, &union, &primes, &backend);
        assert_eq!(exact, modular, "{}: full run mismatch", entry.label);
        // with three 31-bit primes nothing should ever reach tier 3
        assert_eq!(
            stats.exact_evals, 0,
            "{}: unexpected exact evals",
            entry.label
        );
        assert_eq!(stats.hidden_positives, 0, "{}", entry.label);
        checked += 1;
    }
    assert!(checked > 100, "catalog too small: {checked}");
}

/// Both CPU backends drive the engine to identical results (bit-exactness
/// contract of the transform trait).
#[test]
fn reference_and_blocked_backends_agree_end_to_end() {
    let primes = classdiam_core::arith::screening_primes(2);
    for entry in union_catalog().into_iter().filter(|e| e.n <= 8) {
        let index = PartitionIndex::build(entry.n).unwrap();
        let union = resolve_entry(&index, &entry);
        let (_, run_a, _) = run_both(entry.n, &union, &primes, &CpuReference);
        let (_, run_b, _) = run_both(entry.n, &union, &primes, &CpuBlocked);
        assert_eq!(run_a, run_b, "{}", entry.label);
    }
}

/// THE adversarial test: run the modular engine on primes {11, 13} against
/// committed tuples where `a_r(ν) > 0` but `≡ 0 (mod 11·13)`. The engine
/// must certify (tier 3 fires, hidden positives found), never stop early,
/// and still reproduce the exact engine's results bit for bit.
#[test]
fn injected_primes_false_zero_does_not_corrupt_results() {
    let fixture: Value = serde_json::from_str(
        &std::fs::read_to_string(
            std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join("../../fixtures/adversarial_v1.json"),
        )
        .unwrap(),
    )
    .unwrap();

    let mut exercised_hidden = 0u64;
    let mut cases = 0;
    for tuple in fixture["tuples"].as_array().unwrap() {
        // only both-prime tuples defeat the {11,13} prime set
        let primes: Vec<u64> = serde_json::from_value(tuple["primes"].clone()).unwrap();
        if primes != vec![11, 13] {
            continue;
        }
        let n = tuple["n"].as_u64().unwrap() as u16;
        let templates: Vec<CycleTypeTemplate> = tuple["union_templates"]
            .as_array()
            .unwrap()
            .iter()
            .map(|t| {
                let parts: Vec<u8> = serde_json::from_value(t.clone()).unwrap();
                CycleTypeTemplate::new(&parts).unwrap()
            })
            .collect();

        let index = PartitionIndex::build(n).unwrap();
        let union = resolve_union(&index, &templates, false).unwrap();
        let injected = [Prime31(11), Prime31(13)];
        let (exact, modular, stats) = run_both(n, &union, &injected, &CpuBlocked);

        // results must be IDENTICAL despite the engineered false zeros
        assert_eq!(exact, modular, "adversarial n={n} run mismatch");
        // the gate must actually have worked for its living
        assert!(
            stats.candidates > 0,
            "no candidates at n={n} — tuple stale?"
        );
        exercised_hidden += stats.hidden_positives;
        cases += 1;
    }
    assert!(cases >= 10, "too few adversarial cases ran: {cases}");
    assert!(
        exercised_hidden > 0,
        "no hidden positive was ever resurrected — the adversarial path was not exercised"
    );
}

/// Certification telemetry sanity on a normal run: candidates appear (early
/// radii screen unreached types), all resolved by tiers 1/2 with production
/// primes.
#[test]
fn certification_tiers_fire_in_expected_order() {
    let primes = classdiam_core::arith::screening_primes(3);
    let index = PartitionIndex::build(9).unwrap();
    let t: CycleTypeTemplate = "2".parse().unwrap();
    let union = resolve_union(&index, &[t], false).unwrap();
    let (_, _, stats) = run_both(9, &union, &primes, &CpuBlocked);
    assert!(stats.candidates > 0);
    assert_eq!(
        stats.candidates,
        stats.bound_certified + stats.crt_resident_certified + stats.exact_evals
    );
    assert_eq!(stats.exact_evals, 0);
}
