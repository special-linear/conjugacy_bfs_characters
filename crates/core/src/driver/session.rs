//! Per-`n` reusable state and single-job execution (the driver's core).
//!
//! [`Session`] holds the expensive per-`n` objects (partition index, MN
//! evaluator, modular tables) so several unions at the same `n` share them
//! — design doc 01 §12, with one deliberate divergence: the constructor
//! takes no base classes, because [`ModularContext`] is union-independent
//! and [`BaseSpectra`] is per-union and cheap.

use std::path::PathBuf;
use std::time::Instant;

use crate::arith::{screening_primes, Prime31};
use crate::chars::MnEvaluator;
use crate::checkpoint::{write_checkpoint, CheckpointBody};
use crate::engine::exact::run_exact;
use crate::engine::modular::{
    run_modular_resumable, ModularContext, ModularOptions, ModularOutcome,
};
use crate::error::ClassdiamError;
use crate::partition::{CycleTypeTemplate, PartitionIndex};
use crate::report::{build_result, union_slug, EngineDescriptor, RunMeta};
use crate::spectra::{resolve_union, BaseSpectra};
use crate::transform::cpu::CpuBlocked;
use crate::transform::TransformBackend;

use super::{
    io_ctx, job_config_hash, utc_now_rfc3339, DriverHooks, EngineKind, JobReport, JobStatus,
    ProgressEvent,
};

/// Everything one job needs beyond the union itself.
pub struct JobOptions {
    pub allow_identity: bool,
    /// Between-radii wall-clock guard (modular engine only).
    pub deadline: Option<Instant>,
    /// Where `{job_name}.json` is written; `None` = in-memory result only.
    pub result_dir: Option<PathBuf>,
    /// Where `{job_name}.ckpt` is written on suspension (and retired on
    /// completion); `None` = suspension state is dropped.
    pub checkpoint_dir: Option<PathBuf>,
    /// Generators label in the document; defaults to the union slug.
    pub label: Option<String>,
    pub run_id: String,
    pub run_config_hash: String,
    pub resume: Option<CheckpointBody>,
    /// Position within a batch, for progress display only.
    pub job_index: usize,
    pub job_total: usize,
}

/// Per-`n` reusable state: `PartitionIndex` + `MnEvaluator` (+
/// `ModularContext` for the modular engine). `Send + Sync`; build once,
/// run many unions.
pub struct Session {
    index: PartitionIndex,
    mn: MnEvaluator,
    prime_list: Vec<Prime31>,
    ctx: Option<ModularContext>,
    engine: EngineKind,
}

impl Session {
    pub fn new(n: u16, engine: EngineKind, prime_count: usize) -> Result<Self, ClassdiamError> {
        Self::with_primes(n, engine, screening_primes(prime_count))
    }

    /// Exact prime control — resume paths must reuse the checkpoint's primes.
    pub fn with_primes(
        n: u16,
        engine: EngineKind,
        prime_list: Vec<Prime31>,
    ) -> Result<Self, ClassdiamError> {
        let index = PartitionIndex::build(n)?;
        let mn = MnEvaluator::new(n);
        let ctx = (engine == EngineKind::Modular)
            .then(|| ModularContext::build(&index, &mn, &prime_list));
        Ok(Self {
            index,
            mn,
            prime_list,
            ctx,
            engine,
        })
    }

    pub fn index(&self) -> &PartitionIndex {
        &self.index
    }

    pub fn engine(&self) -> EngineKind {
        self.engine
    }

    pub fn primes(&self) -> &[Prime31] {
        &self.prime_list
    }

    /// Run one union end-to-end: resolve, compute, certify, write the
    /// result document (or a checkpoint on suspension).
    pub fn run_union(
        &self,
        templates: &[CycleTypeTemplate],
        opts: JobOptions,
        hooks: &mut DriverHooks<'_>,
    ) -> Result<JobReport, ClassdiamError> {
        let index = &self.index;
        let n = index.n();
        let slug = format!("g{}", union_slug(templates));
        let job_name = format!("n{n:02}_{slug}");
        let union = resolve_union(index, templates, opts.allow_identity)?;
        let started = utc_now_rfc3339();
        let t0 = Instant::now();
        let spectra = BaseSpectra::build(index, &self.mn, &union.class_ids)?;
        let resumed = opts.resume.is_some();
        let suspend_count_before = opts.resume.as_ref().map_or(0, |b| b.suspend_count);

        let (run, descriptor) = match self.engine {
            EngineKind::Exact => {
                let run = run_exact(index, &self.mn, &union)?;
                (run, EngineDescriptor::exact_reference())
            }
            EngineKind::Modular => {
                let ctx = self.ctx.as_ref().expect("modular context prepared");
                let options = ModularOptions {
                    deadline: opts.deadline,
                    cancel: hooks.cancel.clone(),
                    allow_identity_generator: opts.allow_identity,
                    ..Default::default()
                };

                // Adapt the caller's ProgressEvent hook to the engine's
                // per-committed-layer checkpoint callback.
                let mut layer_hook_storage;
                let mut layer_cb: Option<&mut dyn FnMut(&CheckpointBody)> = None;
                if let Some(progress) = hooks.progress.as_deref_mut() {
                    let job_name = job_name.clone();
                    let (job_index, job_total) = (opts.job_index, opts.job_total);
                    layer_hook_storage = move |body: &CheckpointBody| {
                        let last = body.layers.last();
                        let event = ProgressEvent {
                            n,
                            job_name: job_name.clone(),
                            radius: body.committed_radius,
                            new_count: last.map_or(0, |l| l.new.len()),
                            support_size: last.map_or(0, |l| l.support.len()),
                            reachable: body.distance.iter().filter(|&&d| d >= 0).count(),
                            job_index,
                            job_total,
                        };
                        progress(&event);
                    };
                    layer_cb = Some(&mut layer_hook_storage);
                }

                match run_modular_resumable(
                    index,
                    &self.mn,
                    ctx,
                    &spectra,
                    &union,
                    &CpuBlocked,
                    &options,
                    opts.resume,
                    layer_cb,
                )? {
                    ModularOutcome::Completed(run, stats) => {
                        let descriptor =
                            EngineDescriptor::modular(&self.prime_list, CpuBlocked.name(), stats);
                        (run, descriptor)
                    }
                    ModularOutcome::Suspended(body) => {
                        let mut checkpoint_path = None;
                        if let Some(dir) = &opts.checkpoint_dir {
                            std::fs::create_dir_all(dir).map_err(|e| {
                                io_ctx(format!("cannot create {}", dir.display()), e)
                            })?;
                            let path = dir.join(format!("{job_name}.ckpt"));
                            let config = job_config_hash(
                                n,
                                &body.resolved_classes,
                                body.allow_identity_generator,
                                &body.primes,
                            );
                            write_checkpoint(&path, &config, index.order_hash(), &body)?;
                            checkpoint_path = Some(path);
                        }
                        return Ok(JobReport {
                            n,
                            slug,
                            job_name,
                            elapsed_s: t0.elapsed().as_secs_f64(),
                            class_count: index.count(),
                            status: JobStatus::Suspended {
                                checkpoint_path,
                                committed_radius: body.committed_radius,
                            },
                        });
                    }
                }
            }
        };

        let elapsed = t0.elapsed().as_secs_f64();
        let document = build_result(
            index,
            templates,
            &union,
            &spectra,
            &run,
            Some(opts.label.unwrap_or_else(|| slug.clone())),
            opts.allow_identity,
            RunMeta {
                run_id: opts.run_id,
                started_utc: started,
                finished_utc: utc_now_rfc3339(),
                threads: rayon::current_num_threads() as u32,
                total_wall_s: elapsed,
                config_hash: opts.run_config_hash,
                resumed_from_checkpoint: resumed,
                suspend_resume_count: suspend_count_before,
            },
            descriptor,
        );

        let mut file = None;
        if let Some(dir) = &opts.result_dir {
            std::fs::create_dir_all(dir)
                .map_err(|e| io_ctx(format!("cannot create {}", dir.display()), e))?;
            let name = format!("{job_name}.json");
            let path = dir.join(&name);
            let text = serde_json::to_string_pretty(&document).expect("document serializes");
            std::fs::write(&path, text)
                .map_err(|e| io_ctx(format!("cannot write {}", path.display()), e))?;
            file = Some(name);
        }
        // completed: retire any leftover checkpoint
        if let Some(dir) = &opts.checkpoint_dir {
            let ckpt = dir.join(format!("{job_name}.ckpt"));
            let _ = std::fs::remove_file(&ckpt);
            let _ = std::fs::remove_file(ckpt.with_extension("ckpt.prev"));
        }

        Ok(JobReport {
            n,
            slug,
            job_name,
            elapsed_s: elapsed,
            class_count: index.count(),
            status: JobStatus::Done {
                document: Box::new(document),
                file,
            },
        })
    }
}
