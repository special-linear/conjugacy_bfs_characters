//! Driver-layer tests: cooperative cancellation (P5) and, as the driver
//! module grows, batch/resume orchestration and config-hash pinning.

use classdiam_core::arith::screening_primes;
use classdiam_core::chars::MnEvaluator;
use classdiam_core::checkpoint::CheckpointBody;
use classdiam_core::driver::CancelToken;
use classdiam_core::engine::modular::{
    run_modular, run_modular_resumable, ModularContext, ModularOptions, ModularOutcome,
};
use classdiam_core::partition::{CycleTypeTemplate, PartitionIndex};
use classdiam_core::spectra::{resolve_union, BaseSpectra};
use classdiam_core::testing::catalog::{resolve_entry, union_catalog};
use classdiam_core::transform::cpu::CpuBlocked;

/// Cancellation mid-run suspends with a fully-committed checkpoint, and
/// resuming completes to the exact uncancelled result (same guarantee the
/// deadline guard provides).
#[test]
fn cancel_suspends_and_resume_matches_uncancelled() {
    let primes = screening_primes(2);
    let backend = CpuBlocked;

    // transpositions at n=9: an 8-radius run, plenty of room to cancel mid-way
    let index = PartitionIndex::build(9).unwrap();
    let mn = MnEvaluator::new(9);
    let templates = vec!["2".parse::<CycleTypeTemplate>().unwrap()];
    let union = resolve_union(&index, &templates, false).unwrap();
    let spectra = BaseSpectra::build(&index, &mn, &union.class_ids).unwrap();
    let ctx = ModularContext::build(&index, &mn, &primes);

    let (reference, _) = run_modular(
        &index,
        &mn,
        &ctx,
        &spectra,
        &union,
        &backend,
        &ModularOptions::default(),
    )
    .unwrap();
    assert!(reference.stop_radius >= 3, "need a multi-radius run");

    // cancel from the layer-commit hook after two committed radii
    let cancel = CancelToken::new();
    let hook_token = cancel.clone();
    let mut on_layer = |body: &CheckpointBody| {
        if body.committed_radius >= 2 {
            hook_token.cancel();
        }
    };
    let options = ModularOptions {
        cancel: Some(cancel.clone()),
        ..Default::default()
    };
    let outcome = run_modular_resumable(
        &index,
        &mn,
        &ctx,
        &spectra,
        &union,
        &backend,
        &options,
        None,
        Some(&mut on_layer),
    )
    .unwrap();
    let ModularOutcome::Suspended(body) = outcome else {
        panic!("cancelled run must suspend");
    };
    assert_eq!(body.committed_radius, 2);
    assert_eq!(body.suspend_count, 1);
    assert!(cancel.is_cancelled());

    // resume with a fresh (untripped) token: completes and matches
    let outcome = run_modular_resumable(
        &index,
        &mn,
        &ctx,
        &spectra,
        &union,
        &backend,
        &ModularOptions {
            cancel: Some(CancelToken::new()),
            ..Default::default()
        },
        Some(body),
        None,
    )
    .unwrap();
    let ModularOutcome::Completed(resumed, _) = outcome else {
        panic!("resume must complete");
    };
    assert_eq!(resumed, reference);
}

/// A token tripped before the run starts suspends before radius 1, exactly
/// like an already-expired deadline.
#[test]
fn pre_cancelled_token_suspends_at_radius_zero() {
    let primes = screening_primes(2);
    let index = PartitionIndex::build(8).unwrap();
    let mn = MnEvaluator::new(8);
    let entry = union_catalog().into_iter().find(|e| e.n == 8).unwrap();
    let union = resolve_entry(&index, &entry);
    let spectra = BaseSpectra::build(&index, &mn, &union.class_ids).unwrap();
    let ctx = ModularContext::build(&index, &mn, &primes);

    let cancel = CancelToken::new();
    cancel.cancel();
    let options = ModularOptions {
        cancel: Some(cancel),
        ..Default::default()
    };
    let outcome = run_modular_resumable(
        &index, &mn, &ctx, &spectra, &union, &CpuBlocked, &options, None, None,
    )
    .unwrap();
    let ModularOutcome::Suspended(body) = outcome else {
        panic!("pre-cancelled run must suspend");
    };
    assert_eq!(body.committed_radius, 0);
}

/// Token clones share one flag.
#[test]
fn cancel_token_clones_share_flag() {
    let a = CancelToken::new();
    let b = a.clone();
    assert!(!a.is_cancelled() && !b.is_cancelled());
    b.cancel();
    assert!(a.is_cancelled() && b.is_cancelled());
}
