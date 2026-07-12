//! Job orchestration shared by every front end (CLI, Python).
//!
//! The math engines ([`crate::engine`]) compute one `(n, union)` result;
//! this module owns everything around them: run/job configuration hashes,
//! result/checkpoint/manifest file layout, batch and resume loops, progress
//! hooks, and cooperative cancellation. Front ends supply argument parsing
//! and presentation only, so a run directory produced by one front end is
//! always resumable by another.
//!
//! ## Hash contract (checkpoint compatibility)
//!
//! [`job_config_hash`] and [`run_config_hash`] are written into checkpoint
//! headers / manifests and re-derived on resume. Their input strings are
//! versioned (`job-v1`, `run-v1`) and MUST NOT change silently — the
//! `job_config_hash_is_pinned` test freezes the derivation.

mod cancel;
mod session;

pub use cancel::CancelToken;
pub use session::{JobOptions, Session};

use std::path::{Path, PathBuf};
use std::time::{Instant, SystemTime};

use crate::arith::{screening_primes, Prime31};
use crate::checkpoint::read_checkpoint;
use crate::error::ClassdiamError;
use crate::partition::CycleTypeTemplate;
use crate::report::schema::ResultDocument;
use crate::report::union_slug;

/// Which engine executes jobs.
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub enum EngineKind {
    /// Modular screening + certified zeros (production).
    Modular,
    /// Exact big-integer reference (small n; oracle).
    Exact,
}

/// One committed layer of progress in a modular run. The exact engine
/// emits no events (its runs are seconds).
#[derive(Clone, Debug)]
pub struct ProgressEvent {
    pub n: u16,
    /// e.g. `n30_g2`.
    pub job_name: String,
    /// Last committed radius.
    pub radius: u32,
    /// Types first reached at exactly this radius.
    pub new_count: usize,
    /// Exact-length support size of this layer.
    pub support_size: usize,
    /// Types with a known finite distance so far.
    pub reachable: usize,
    /// Position of this job in the batch (0-based) and batch size.
    pub job_index: usize,
    pub job_total: usize,
}

/// Caller-supplied hooks threaded through batch/resume/session runs.
///
/// `progress` fires once per committed layer; `cancel` is observed between
/// radii (suspending exactly like a deadline); `on_job` fires right after
/// each job finishes (the CLI prints its per-job line there).
#[derive(Default)]
pub struct DriverHooks<'a> {
    pub progress: Option<&'a mut dyn FnMut(&ProgressEvent)>,
    pub cancel: Option<CancelToken>,
    pub on_job: Option<&'a mut dyn FnMut(&JobReport)>,
}

/// How a single `(n, union)` job ended.
pub enum JobStatus {
    Done {
        document: Box<ResultDocument>,
        /// Bare result file name (e.g. `n06_g2.json`) when a result
        /// directory was configured; `None` for in-memory runs.
        file: Option<String>,
    },
    Suspended {
        /// Where the checkpoint was written; `None` when no checkpoint
        /// directory was configured (suspension state dropped).
        checkpoint_path: Option<PathBuf>,
        committed_radius: u32,
    },
    Skipped {
        reason: String,
    },
}

/// Per-job outcome plus the labels front ends print.
pub struct JobReport {
    pub n: u16,
    /// Union label used in file names and manifests, e.g. `g3+2.2`.
    pub slug: String,
    /// `n{n:02}_{slug}`, e.g. `n06_g2`.
    pub job_name: String,
    pub elapsed_s: f64,
    /// p(n) — the denominator of "reachable X/Y" displays.
    pub class_count: usize,
    pub status: JobStatus,
}

/// A batch of `(n, union)` jobs sharing per-`n` tables.
pub struct BatchConfig {
    pub ns: Vec<u16>,
    pub unions: Vec<Vec<CycleTypeTemplate>>,
    pub engine: EngineKind,
    /// Number of resident screening primes (modular engine).
    pub prime_count: usize,
    /// Wall-clock deadline shared by all jobs; on expiry remaining modular
    /// jobs suspend with checkpoints.
    pub deadline: Option<Instant>,
    pub allow_identity: bool,
    /// `None` = fully in-memory: no result files, no manifest, no
    /// checkpoints (suspension state is dropped).
    pub out_dir: Option<PathBuf>,
    /// `None` = derive via [`make_run_id`].
    pub run_id: Option<String>,
}

/// What a batch (or resume) produced.
pub struct BatchOutcome {
    pub run_id: String,
    pub config_hash: String,
    pub out_dir: Option<PathBuf>,
    pub any_suspended: bool,
    /// The cancel token fired during the batch.
    pub cancelled: bool,
    pub reports: Vec<JobReport>,
}

pub fn utc_now_rfc3339() -> String {
    humantime::format_rfc3339_seconds(SystemTime::now()).to_string()
}

/// Per-job configuration hash — written into checkpoint headers and
/// re-derived on resume from the checkpoint body itself.
pub fn job_config_hash(
    n: u16,
    resolved_classes: &[Vec<u8>],
    allow_identity: bool,
    primes: &[u32],
) -> [u8; 32] {
    let mut hasher = blake3::Hasher::new();
    hasher.update(
        format!("job-v1;n={n};classes={resolved_classes:?};allow_identity={allow_identity};primes={primes:?}")
            .as_bytes(),
    );
    *hasher.finalize().as_bytes()
}

/// Whole-run configuration hash (manifest `config_hash_blake3`, run-id
/// suffix, result-document `config_hash_blake3`).
pub fn run_config_hash(
    ns: &[u16],
    unions: &[Vec<CycleTypeTemplate>],
    allow_identity: bool,
    prime_count: usize,
) -> String {
    let mut hasher = blake3::Hasher::new();
    hasher.update(
        format!("run-v1;n={ns:?};allow_identity={allow_identity};primes={prime_count}").as_bytes(),
    );
    for u in unions {
        hasher.update(format!("{u:?};").as_bytes());
    }
    hasher.finalize().to_hex().to_string()
}

/// `{utc timestamp}-{hash prefix}`, e.g. `20260712T101500Z-8f3a2c1d`.
pub fn make_run_id(config_hash: &str) -> String {
    format!(
        "{}-{}",
        utc_now_rfc3339().replace(['-', ':'], ""),
        &config_hash[..8]
    )
}

/// Parse an `n` spec: `12`, `6..=12`, or `6,8,10`.
pub fn parse_n_spec(spec: &str) -> Result<Vec<u16>, ClassdiamError> {
    let bad = |reason: String| ClassdiamError::InvalidSpec { reason };
    let spec = spec.trim();
    if let Some((a, b)) = spec.split_once("..=") {
        let a: u16 = a
            .trim()
            .parse()
            .map_err(|_| bad(format!("bad range start in {spec:?}")))?;
        let b: u16 = b
            .trim()
            .parse()
            .map_err(|_| bad(format!("bad range end in {spec:?}")))?;
        if a > b {
            return Err(bad(format!("empty n range {spec}")));
        }
        return Ok((a..=b).collect());
    }
    spec.split(',')
        .map(|t| {
            t.trim()
                .parse::<u16>()
                .map_err(|_| bad(format!("bad n value {t:?}")))
        })
        .collect()
}

/// Parse a union spec: classes joined `+`, parts joined `,` (e.g. `3+2,2`).
pub fn parse_union(spec: &str) -> Result<Vec<CycleTypeTemplate>, ClassdiamError> {
    spec.split('+')
        .map(|part| {
            part.parse::<CycleTypeTemplate>()
                .map_err(|e| ClassdiamError::InvalidSpec {
                    reason: format!("bad class {part:?} in union {spec:?}: {e}"),
                })
        })
        .collect()
}

pub(crate) fn io_ctx(context: String, source: std::io::Error) -> ClassdiamError {
    ClassdiamError::IoContext { context, source }
}

/// Run every `(n, union)` job of the batch, sharing per-`n` tables.
/// Result/checkpoint/manifest files land under `cfg.out_dir` (if any) in
/// the layout the CLI has always produced.
pub fn run_batch(
    cfg: &BatchConfig,
    hooks: &mut DriverHooks<'_>,
) -> Result<BatchOutcome, ClassdiamError> {
    let prime_list = screening_primes(cfg.prime_count);
    let config_hash = run_config_hash(&cfg.ns, &cfg.unions, cfg.allow_identity, cfg.prime_count);
    let run_id = cfg
        .run_id
        .clone()
        .unwrap_or_else(|| make_run_id(&config_hash));
    if let Some(dir) = &cfg.out_dir {
        std::fs::create_dir_all(dir)
            .map_err(|e| io_ctx(format!("cannot create {}", dir.display()), e))?;
    }

    let job_total = cfg.ns.len() * cfg.unions.len();
    let mut job_index = 0usize;
    let mut reports: Vec<JobReport> = Vec::new();
    let mut any_suspended = false;
    for &n in &cfg.ns {
        let session = Session::with_primes(n, cfg.engine, prime_list.clone())?;
        for templates in &cfg.unions {
            let opts = JobOptions {
                allow_identity: cfg.allow_identity,
                deadline: cfg.deadline,
                result_dir: cfg.out_dir.clone(),
                checkpoint_dir: cfg.out_dir.as_ref().map(|d| d.join("checkpoints")),
                label: None,
                run_id: run_id.clone(),
                run_config_hash: config_hash.clone(),
                resume: None,
                job_index,
                job_total,
            };
            let report = match session.run_union(templates, opts, hooks) {
                Ok(report) => report,
                Err(e @ ClassdiamError::TemplateDoesNotFit { .. }) => {
                    let slug = format!("g{}", union_slug(templates));
                    JobReport {
                        n,
                        job_name: format!("n{n:02}_{slug}"),
                        slug,
                        elapsed_s: 0.0,
                        class_count: session.index().count(),
                        status: JobStatus::Skipped {
                            reason: e.to_string(),
                        },
                    }
                }
                Err(e) => return Err(e),
            };
            any_suspended |= matches!(report.status, JobStatus::Suspended { .. });
            if let Some(cb) = hooks.on_job.as_deref_mut() {
                cb(&report);
            }
            reports.push(report);
            job_index += 1;
        }
    }

    if let Some(dir) = &cfg.out_dir {
        let jobs = reports.iter().map(manifest_entry).collect();
        write_manifest(dir, &run_id, &config_hash, any_suspended, jobs)?;
    }
    Ok(BatchOutcome {
        run_id,
        config_hash,
        out_dir: cfg.out_dir.clone(),
        any_suspended,
        cancelled: hooks.cancel.as_ref().is_some_and(CancelToken::is_cancelled),
        reports,
    })
}

/// Resume every suspended job in a run directory's `checkpoints/`,
/// validating config and partition-order hashes, then patch the manifest.
pub fn resume_batch(
    run_dir: &Path,
    deadline: Option<Instant>,
    hooks: &mut DriverHooks<'_>,
) -> Result<BatchOutcome, ClassdiamError> {
    let ckpt_dir = run_dir.join("checkpoints");
    let mut paths: Vec<PathBuf> = std::fs::read_dir(&ckpt_dir)
        .map_err(|_| ClassdiamError::InvalidRunDir {
            reason: format!("no checkpoints directory in {}", run_dir.display()),
        })?
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.extension().is_some_and(|e| e == "ckpt"))
        .collect();
    paths.sort();
    if paths.is_empty() {
        return Err(ClassdiamError::InvalidRunDir {
            reason: format!("no .ckpt files in {}", ckpt_dir.display()),
        });
    }

    // reload the existing manifest to patch job entries
    let manifest_path = run_dir.join("manifest.json");
    let manifest_text =
        std::fs::read_to_string(&manifest_path).map_err(|_| ClassdiamError::InvalidRunDir {
            reason: format!("missing manifest.json in {}", run_dir.display()),
        })?;
    let manifest: serde_json::Value =
        serde_json::from_str(&manifest_text).map_err(|e| ClassdiamError::InvalidRunDir {
            reason: format!("malformed manifest.json: {e}"),
        })?;
    let run_id = manifest["run_id"].as_str().unwrap_or("resumed").to_string();
    let run_config_hash = manifest["config_hash_blake3"]
        .as_str()
        .unwrap_or("")
        .to_string();
    let mut jobs: Vec<serde_json::Value> =
        serde_json::from_value(manifest["jobs"].clone()).unwrap_or_default();

    let job_total = paths.len();
    let mut any_suspended = false;
    let mut reports: Vec<JobReport> = Vec::new();
    for (job_index, path) in paths.iter().enumerate() {
        // Parse without expected hashes first (the config is reconstructed
        // from the body), then verify the header hash against the re-derived
        // config and the order hash against the freshly built index.
        let (body, header_config, header_order) = read_checkpoint(path, None, None)?;
        let expected_config = job_config_hash(
            body.n,
            &body.resolved_classes,
            body.allow_identity_generator,
            &body.primes,
        );
        if header_config != expected_config {
            return Err(ClassdiamError::CheckpointMismatch {
                what: format!(
                    "{}: checkpoint config hash does not match its own body — refusing",
                    path.display()
                ),
            });
        }
        let prime_list: Vec<Prime31> = body.primes.iter().map(|&p| Prime31(p)).collect();
        let session = Session::with_primes(body.n, EngineKind::Modular, prime_list)?;
        if &header_order != session.index().order_hash() {
            return Err(ClassdiamError::CheckpointMismatch {
                what: format!(
                    "{}: partition-order hash mismatch — checkpoint from an incompatible version",
                    path.display()
                ),
            });
        }
        let templates: Vec<CycleTypeTemplate> = body
            .resolved_classes
            .iter()
            .map(|parts| {
                let reduced: Vec<u8> = parts.iter().copied().filter(|&x| x >= 2).collect();
                CycleTypeTemplate::new(&reduced).expect("checkpoint classes valid")
            })
            .collect();
        let opts = JobOptions {
            allow_identity: body.allow_identity_generator,
            deadline,
            result_dir: Some(run_dir.to_path_buf()),
            checkpoint_dir: Some(ckpt_dir.clone()),
            label: None,
            run_id: run_id.clone(),
            run_config_hash: run_config_hash.clone(),
            resume: Some(body),
            job_index,
            job_total,
        };
        let report = session.run_union(&templates, opts, hooks)?;
        any_suspended |= matches!(report.status, JobStatus::Suspended { .. });

        // patch the matching manifest entry (or append)
        let entry = manifest_entry(&report);
        let key = (entry["n"].clone(), entry["union"].clone());
        if let Some(existing) = jobs
            .iter_mut()
            .find(|j| (j["n"].clone(), j["union"].clone()) == key)
        {
            *existing = entry;
        } else {
            jobs.push(entry);
        }
        if let Some(cb) = hooks.on_job.as_deref_mut() {
            cb(&report);
        }
        reports.push(report);
    }

    write_manifest(run_dir, &run_id, &run_config_hash, any_suspended, jobs)?;
    Ok(BatchOutcome {
        run_id,
        config_hash: run_config_hash,
        out_dir: Some(run_dir.to_path_buf()),
        any_suspended,
        cancelled: hooks.cancel.as_ref().is_some_and(CancelToken::is_cancelled),
        reports,
    })
}

/// The manifest job entry for a report — shapes are frozen (manifest
/// format_version 1); serde_json's default map ordering keeps output
/// byte-stable.
fn manifest_entry(report: &JobReport) -> serde_json::Value {
    match &report.status {
        JobStatus::Done { document, file } => serde_json::json!({
            "n": report.n, "union": report.slug, "status": "done",
            "file": file.clone().unwrap_or_else(|| format!("{}.json", report.job_name)),
            "diameter": document.results.diameter_identity_component,
            "stop_radius": document.results.stopping.stop_radius,
            "reachable": document.results.reachable_count,
        }),
        JobStatus::Suspended {
            committed_radius, ..
        } => serde_json::json!({
            "n": report.n, "union": report.slug, "status": "suspended",
            "checkpoint": format!("checkpoints/{}.ckpt", report.job_name),
            "committed_radius": committed_radius,
        }),
        JobStatus::Skipped { reason } => serde_json::json!({
            "n": report.n, "union": report.slug, "status": "skipped", "reason": reason,
        }),
    }
}

fn write_manifest(
    out_dir: &Path,
    run_id: &str,
    config_hash: &str,
    any_suspended: bool,
    jobs: Vec<serde_json::Value>,
) -> Result<(), ClassdiamError> {
    let manifest = serde_json::json!({
        "format": "classdiam/manifest",
        "format_version": 1,
        "run_id": run_id,
        "config_hash_blake3": config_hash,
        "status": if any_suspended { "suspended" } else { "completed" },
        "jobs": jobs,
    });
    let path = out_dir.join("manifest.json");
    let text = serde_json::to_string_pretty(&manifest).expect("json value serializes");
    std::fs::write(&path, text).map_err(|e| io_ctx(format!("cannot write {}", path.display()), e))
}
