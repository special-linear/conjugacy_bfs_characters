//! Driver-layer tests: cooperative cancellation, batch/resume orchestration,
//! file layout, progress events, and config-hash pinning (P5).

use std::path::PathBuf;

use classdiam_core::arith::screening_primes;
use classdiam_core::chars::MnEvaluator;
use classdiam_core::checkpoint::CheckpointBody;
use classdiam_core::driver::{
    job_config_hash, parse_n_spec, parse_union, resume_batch, run_batch, BatchConfig, CancelToken,
    DriverHooks, EngineKind, JobStatus, ProgressEvent,
};
use classdiam_core::engine::modular::{
    run_modular, run_modular_resumable, ModularContext, ModularOptions, ModularOutcome,
};
use classdiam_core::partition::{CycleTypeTemplate, PartitionIndex};
use classdiam_core::report::schema::ResultDocument;
use classdiam_core::spectra::{resolve_union, BaseSpectra};
use classdiam_core::testing::catalog::{resolve_entry, union_catalog};
use classdiam_core::transform::cpu::CpuBlocked;
use serde_json::Value;

fn temp_dir(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("classdiam_driver_{tag}"));
    let _ = std::fs::remove_dir_all(&dir);
    dir
}

fn transpositions_batch(n: u16, engine: EngineKind, out_dir: Option<PathBuf>) -> BatchConfig {
    BatchConfig {
        ns: vec![n],
        unions: vec![vec!["2".parse().unwrap()]],
        engine,
        prime_count: 2,
        deadline: None,
        allow_identity: false,
        out_dir,
        run_id: None,
    }
}

/// Strip the documented volatile fields (`run`, `timings_s`, `tool`,
/// `engine.threads`) — same rule as the golden-file test — plus
/// `config_hash_blake3`, which depends on the run configuration (the
/// golden file carries a placeholder there).
fn strip_volatile(value: &mut Value) {
    let object = value.as_object_mut().expect("document is an object");
    object.remove("run");
    object.remove("timings_s");
    object.remove("tool");
    object.remove("config_hash_blake3");
    if let Some(engine) = object.get_mut("engine").and_then(Value::as_object_mut) {
        engine.remove("threads");
    }
}

fn done_document(status: &JobStatus) -> &ResultDocument {
    match status {
        JobStatus::Done { document, .. } => document,
        _ => panic!("expected a completed job"),
    }
}

fn stripped(document: &ResultDocument) -> Value {
    let mut value = serde_json::to_value(document).unwrap();
    strip_volatile(&mut value);
    value
}

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
        &index,
        &mn,
        &ctx,
        &spectra,
        &union,
        &CpuBlocked,
        &options,
        None,
        None,
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

/// Checkpoint-compat tripwire: the `job-v1` hash derivation is frozen.
/// If this fails, existing checkpoints become unresumable — do NOT update
/// the expected value without bumping the checkpoint format version.
#[test]
fn job_config_hash_is_pinned() {
    let hash = job_config_hash(
        6,
        &[vec![2, 1, 1, 1, 1, 1]],
        false,
        &[2147483647, 2147483629],
    );
    let hex: String = hash.iter().map(|b| format!("{b:02x}")).collect();
    assert_eq!(
        hex, "f8fb02e59728b1a8931cd42095d9b1e055aafb41489e45742e7594df209be539",
        "job-v1 hash derivation changed — breaks checkpoint compatibility"
    );
}

/// `run_batch` writes the frozen file layout and the document matches the
/// committed golden file modulo volatile fields.
#[test]
fn run_batch_layout_and_golden() {
    let dir = temp_dir("batch_golden");
    let cfg = transpositions_batch(6, EngineKind::Exact, Some(dir.clone()));
    let out = run_batch(&cfg, &mut DriverHooks::default()).unwrap();
    assert!(!out.any_suspended && !out.cancelled);
    assert_eq!(out.reports.len(), 1);

    // in-memory document == on-disk file
    let document = done_document(&out.reports[0].status);
    let disk: Value =
        serde_json::from_str(&std::fs::read_to_string(dir.join("n06_g2.json")).unwrap()).unwrap();
    assert_eq!(serde_json::to_value(document).unwrap(), disk);

    // matches the committed golden file modulo volatile fields
    let golden_path =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/golden/n06_g2.json");
    let mut golden: Value =
        serde_json::from_str(&std::fs::read_to_string(golden_path).unwrap()).unwrap();
    strip_volatile(&mut golden);
    assert_eq!(stripped(document), golden);

    // manifest shape
    let manifest: Value =
        serde_json::from_str(&std::fs::read_to_string(dir.join("manifest.json")).unwrap()).unwrap();
    assert_eq!(manifest["format"], "classdiam/manifest");
    assert_eq!(manifest["format_version"], 1);
    assert_eq!(manifest["status"], "completed");
    assert_eq!(manifest["run_id"], Value::String(out.run_id.clone()));
    assert_eq!(
        manifest["config_hash_blake3"],
        Value::String(out.config_hash)
    );
    let job = &manifest["jobs"][0];
    assert_eq!(job["n"], 6);
    assert_eq!(job["union"], "g2");
    assert_eq!(job["status"], "done");
    assert_eq!(job["file"], "n06_g2.json");
    assert_eq!(job["diameter"], 5);
    assert_eq!(job["stop_radius"], 5);
    assert_eq!(job["reachable"], 11);

    std::fs::remove_dir_all(&dir).unwrap();
}

/// Interop invariant: a suspended run directory (deadline- or cancel-
/// produced) resumes to the identical result an uninterrupted batch gives.
#[test]
fn suspended_run_dir_resumes_to_uninterrupted_result() {
    // uninterrupted reference
    let ref_dir = temp_dir("roundtrip_ref");
    let out = run_batch(
        &transpositions_batch(9, EngineKind::Modular, Some(ref_dir.clone())),
        &mut DriverHooks::default(),
    )
    .unwrap();
    let reference = stripped(done_document(&out.reports[0].status));

    for (tag, deadline, cancel) in [
        (
            "deadline",
            Some(std::time::Instant::now() - std::time::Duration::from_secs(1)),
            None,
        ),
        ("cancel", None, Some(CancelToken::new())),
    ] {
        if let Some(token) = &cancel {
            token.cancel();
        }
        let dir = temp_dir(&format!("roundtrip_{tag}"));
        let mut cfg = transpositions_batch(9, EngineKind::Modular, Some(dir.clone()));
        cfg.deadline = deadline;
        let mut hooks = DriverHooks {
            cancel: cancel.clone(),
            ..Default::default()
        };
        let out = run_batch(&cfg, &mut hooks).unwrap();
        assert!(out.any_suspended, "{tag}: batch must suspend");
        assert_eq!(out.cancelled, cancel.is_some(), "{tag}");
        assert!(dir.join("checkpoints/n09_g2.ckpt").exists(), "{tag}");
        let manifest: Value =
            serde_json::from_str(&std::fs::read_to_string(dir.join("manifest.json")).unwrap())
                .unwrap();
        assert_eq!(manifest["status"], "suspended", "{tag}");
        assert_eq!(manifest["jobs"][0]["checkpoint"], "checkpoints/n09_g2.ckpt");

        // resume with fresh hooks completes and matches the reference
        let out = resume_batch(&dir, None, &mut DriverHooks::default()).unwrap();
        assert!(!out.any_suspended, "{tag}: resume must complete");
        let resumed = done_document(&out.reports[0].status);
        assert_eq!(stripped(resumed), reference, "{tag}: result differs");
        assert!(resumed.run.resumed_from_checkpoint);
        assert_eq!(resumed.run.suspend_resume_count, 1);
        // checkpoint retired, manifest patched to completed
        assert!(!dir.join("checkpoints/n09_g2.ckpt").exists(), "{tag}");
        let manifest: Value =
            serde_json::from_str(&std::fs::read_to_string(dir.join("manifest.json")).unwrap())
                .unwrap();
        assert_eq!(manifest["status"], "completed", "{tag}");
        assert_eq!(manifest["jobs"][0]["status"], "done", "{tag}");

        std::fs::remove_dir_all(&dir).unwrap();
    }
    std::fs::remove_dir_all(&ref_dir).unwrap();
}

/// In-memory batches (out_dir = None) return the document and touch no
/// files; progress fires once per committed radius with sane fields.
#[test]
fn in_memory_run_with_progress_events() {
    let mut events: Vec<ProgressEvent> = Vec::new();
    let mut on_progress = |e: &ProgressEvent| events.push(e.clone());
    let mut hooks = DriverHooks {
        progress: Some(&mut on_progress),
        ..Default::default()
    };
    let out = run_batch(
        &transpositions_batch(8, EngineKind::Modular, None),
        &mut hooks,
    )
    .unwrap();
    assert!(out.out_dir.is_none());
    let JobStatus::Done { document, file } = &out.reports[0].status else {
        panic!("expected completion");
    };
    assert!(file.is_none(), "in-memory run must not name a file");

    let stop = document.results.stopping.stop_radius;
    assert_eq!(
        events.len(),
        stop as usize,
        "one event per committed radius"
    );
    for (i, event) in events.iter().enumerate() {
        assert_eq!(event.radius, i as u32 + 1);
        assert_eq!(event.n, 8);
        assert_eq!(event.job_name, "n08_g2");
        assert_eq!((event.job_index, event.job_total), (0, 1));
    }
    // reachable counts are non-decreasing and end at the document's count
    assert!(events.windows(2).all(|w| w[0].reachable <= w[1].reachable));
    assert_eq!(
        events.last().unwrap().reachable as u64,
        document.results.reachable_count
    );
}

/// Union templates that do not fit an `n` are skipped, not fatal.
#[test]
fn oversized_template_is_skipped() {
    let cfg = BatchConfig {
        ns: vec![4],
        unions: vec![vec!["2".parse().unwrap()], vec!["7".parse().unwrap()]],
        engine: EngineKind::Exact,
        prime_count: 2,
        deadline: None,
        allow_identity: false,
        out_dir: None,
        run_id: None,
    };
    let out = run_batch(&cfg, &mut DriverHooks::default()).unwrap();
    assert!(matches!(out.reports[0].status, JobStatus::Done { .. }));
    let JobStatus::Skipped { reason } = &out.reports[1].status else {
        panic!("7-cycle cannot fit in S_4");
    };
    assert!(reason.contains("does not fit"), "{reason}");
}

#[test]
fn n_spec_and_union_parsing() {
    assert_eq!(parse_n_spec("12").unwrap(), vec![12]);
    assert_eq!(parse_n_spec("6..=8").unwrap(), vec![6, 7, 8]);
    assert_eq!(parse_n_spec(" 6, 8,10 ").unwrap(), vec![6, 8, 10]);
    assert!(parse_n_spec("8..=6").is_err());
    assert!(parse_n_spec("six").is_err());
    assert!(parse_n_spec("").is_err());

    let union = parse_union("3+2,2").unwrap();
    assert_eq!(union.len(), 2);
    assert_eq!(union[0].parts(), &[3]);
    assert_eq!(union[1].parts(), &[2, 2]);
    assert!(parse_union("0").is_err());
    assert!(parse_union("3+x").is_err());
}
