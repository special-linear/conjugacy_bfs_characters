//! Resume-equals-uninterrupted (fixed requirement 6; design doc 03 §6
//! `prop_resume_equivalence`): for every layer checkpoint a run emits,
//! resuming from it must reproduce the uninterrupted final result exactly.
//! Also: deadline suspension round-trips through the checkpoint FILE format
//! and re-enters cleanly.

use classdiam_core::arith::screening_primes;
use classdiam_core::chars::MnEvaluator;
use classdiam_core::checkpoint::{read_checkpoint, write_checkpoint, CheckpointBody};
use classdiam_core::engine::modular::{
    run_modular, run_modular_resumable, ModularContext, ModularOptions, ModularOutcome,
};
use classdiam_core::partition::PartitionIndex;
use classdiam_core::spectra::BaseSpectra;
use classdiam_core::testing::catalog::{resolve_entry, union_catalog};
use classdiam_core::transform::cpu::CpuBlocked;

/// Every mid-run checkpoint resumes to the identical final result.
#[test]
fn resume_from_every_layer_matches_uninterrupted() {
    let primes = screening_primes(2);
    let backend = CpuBlocked;
    let mut resumed_total = 0;
    for entry in union_catalog().into_iter().filter(|e| e.n <= 8) {
        let index = PartitionIndex::build(entry.n).unwrap();
        let mn = MnEvaluator::new(entry.n);
        let union = resolve_entry(&index, &entry);
        let spectra = BaseSpectra::build(&index, &mn, &union.class_ids).unwrap();
        let ctx = ModularContext::build(&index, &mn, &primes);
        let options = ModularOptions::default();

        // uninterrupted run, capturing a checkpoint at every committed layer
        let mut checkpoints: Vec<CheckpointBody> = Vec::new();
        let mut capture = |body: &CheckpointBody| checkpoints.push(body.clone());
        let outcome = run_modular_resumable(
            &index,
            &mn,
            &ctx,
            &spectra,
            &union,
            &backend,
            &options,
            None,
            Some(&mut capture),
        )
        .unwrap();
        let ModularOutcome::Completed(reference, reference_stats) = outcome else {
            panic!("no deadline, cannot suspend");
        };

        for body in checkpoints {
            let outcome = run_modular_resumable(
                &index,
                &mn,
                &ctx,
                &spectra,
                &union,
                &backend,
                &options,
                Some(body),
                None,
            )
            .unwrap();
            let ModularOutcome::Completed(resumed, resumed_stats) = outcome else {
                panic!("no deadline on resume");
            };
            assert_eq!(resumed, reference, "{}: resumed run differs", entry.label);
            // certification totals: the resumed run re-certifies only the
            // remaining radii, so totals differ — but the audit invariant
            // (every candidate resolved by some tier) must hold.
            assert_eq!(
                resumed_stats.candidates,
                resumed_stats.bound_certified
                    + resumed_stats.crt_resident_certified
                    + resumed_stats.exact_evals,
                "{}",
                entry.label
            );
            resumed_total += 1;
        }
        let _ = reference_stats;
    }
    assert!(
        resumed_total > 200,
        "too few resume points: {resumed_total}"
    );
}

/// Deadline suspension → file round trip → resume completes identically,
/// and the resumed state cannot stop on a partially-committed layer
/// (only fully certified layers are ever serialized).
#[test]
fn deadline_suspend_file_roundtrip_resume() {
    let primes = screening_primes(2);
    let backend = CpuBlocked;
    let index = PartitionIndex::build(9).unwrap();
    let mn = MnEvaluator::new(9);
    let entry = union_catalog()
        .into_iter()
        .find(|e| e.label == "n9_random0" || e.n == 9)
        .unwrap();
    let union = resolve_entry(&index, &entry);
    let spectra = BaseSpectra::build(&index, &mn, &union.class_ids).unwrap();
    let ctx = ModularContext::build(&index, &mn, &primes);

    // reference
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

    // expired deadline: suspends immediately after layer 0
    let expired = ModularOptions {
        deadline: Some(std::time::Instant::now() - std::time::Duration::from_secs(1)),
        ..Default::default()
    };
    let outcome = run_modular_resumable(
        &index, &mn, &ctx, &spectra, &union, &backend, &expired, None, None,
    )
    .unwrap();
    let ModularOutcome::Suspended(body) = outcome else {
        panic!("expired deadline must suspend");
    };
    assert_eq!(body.committed_radius, 0);
    assert_eq!(body.suspend_count, 1);

    // file round trip with hash validation
    let dir = std::env::temp_dir().join("classdiam_resume_test");
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("job.ckpt");
    let config_hash = [3u8; 32];
    write_checkpoint(&path, &config_hash, index.order_hash(), &body).unwrap();
    let (loaded, _, _) =
        read_checkpoint(&path, Some(&config_hash), Some(index.order_hash())).unwrap();
    assert_eq!(loaded, body);

    // resume without deadline: must complete and match the reference
    let outcome = run_modular_resumable(
        &index,
        &mn,
        &ctx,
        &spectra,
        &union,
        &backend,
        &ModularOptions::default(),
        Some(loaded),
        None,
    )
    .unwrap();
    let ModularOutcome::Completed(resumed, _) = outcome else {
        panic!("resume must complete");
    };
    assert_eq!(resumed, reference);

    std::fs::remove_dir_all(&dir).unwrap();
}

/// Checkpoints from a different configuration are refused by the engine.
#[test]
fn engine_refuses_mismatched_checkpoint() {
    let primes = screening_primes(2);
    let backend = CpuBlocked;
    let index = PartitionIndex::build(7).unwrap();
    let mn = MnEvaluator::new(7);
    let entries: Vec<_> = union_catalog().into_iter().filter(|e| e.n == 7).collect();
    let union_a = resolve_entry(&index, &entries[0]);
    let union_b = resolve_entry(&index, &entries[1]);
    let spectra_a = BaseSpectra::build(&index, &mn, &union_a.class_ids).unwrap();
    let spectra_b = BaseSpectra::build(&index, &mn, &union_b.class_ids).unwrap();
    let ctx = ModularContext::build(&index, &mn, &primes);

    // capture a checkpoint from union A
    let mut first: Option<CheckpointBody> = None;
    let mut capture = |b: &CheckpointBody| {
        if first.is_none() {
            first = Some(b.clone());
        }
    };
    run_modular_resumable(
        &index,
        &mn,
        &ctx,
        &spectra_a,
        &union_a,
        &backend,
        &ModularOptions::default(),
        None,
        Some(&mut capture),
    )
    .unwrap();

    // feeding it to union B must be refused
    let err = run_modular_resumable(
        &index,
        &mn,
        &ctx,
        &spectra_b,
        &union_b,
        &backend,
        &ModularOptions::default(),
        first,
        None,
    )
    .unwrap_err();
    assert!(matches!(
        err,
        classdiam_core::ClassdiamError::CheckpointMismatch { .. }
    ));
}
