//! `classdiam` CLI: `run`, `resume`, `verify`.

use std::path::{Path, PathBuf};
use std::time::{Duration, Instant, SystemTime};

use anyhow::{bail, Context, Result};
use clap::{Parser, Subcommand, ValueEnum};

use classdiam_core::arith::{screening_primes, Prime31};
use classdiam_core::chars::MnEvaluator;
use classdiam_core::checkpoint::{read_checkpoint, write_checkpoint, CheckpointBody};
use classdiam_core::engine::exact::run_exact;
use classdiam_core::engine::modular::{
    run_modular_resumable, ModularContext, ModularOptions, ModularOutcome,
};
use classdiam_core::partition::{CycleTypeTemplate, PartitionIndex};
use classdiam_core::report::{build_result, union_slug, EngineDescriptor, RunMeta};
use classdiam_core::spectra::{resolve_union, BaseSpectra, ResolvedUnion};
use classdiam_core::transform::cpu::CpuBlocked;
use classdiam_core::ClassdiamError;

#[derive(Parser)]
#[command(
    name = "classdiam",
    version,
    about = "Conjugacy-class BFS on S_n via characters"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Copy, Clone, PartialEq, Eq, ValueEnum)]
enum EngineChoice {
    /// Modular screening + certified zeros (production).
    Modular,
    /// Exact big-integer reference engine (small n; oracle).
    Exact,
}

#[derive(Subcommand)]
enum Command {
    /// Run the reduced BFS for each n in a range and each generating union.
    Run {
        /// n values: `12`, `6..=12`, or `6,8,10`.
        #[arg(short = 'n', long = "n")]
        n_spec: String,
        /// Generating union: classes joined `+`, parts joined `,`
        /// (e.g. `-u 2`, `-u "3+2,2"`). Repeat for several unions.
        #[arg(short = 'u', long = "union", required = true)]
        unions: Vec<String>,
        /// Output directory (default: results/<run_id>).
        #[arg(short = 'o', long = "out")]
        out: Option<PathBuf>,
        #[arg(long, value_enum, default_value = "modular")]
        engine: EngineChoice,
        /// Number of resident screening primes (modular engine).
        #[arg(long, default_value_t = 3)]
        primes: usize,
        /// Wall-clock budget in seconds; on expiry jobs suspend with
        /// resumable checkpoints (exit code 75).
        #[arg(long)]
        deadline: Option<u64>,
        /// Permit the identity class as a generator (spec §5.3).
        #[arg(long)]
        allow_identity: bool,
    },
    /// Resume suspended jobs from a run directory's checkpoints.
    Resume {
        run_dir: PathBuf,
        /// Wall-clock budget in seconds for this session.
        #[arg(long)]
        deadline: Option<u64>,
    },
    /// Cross-check the character engine against brute-force BFS over raw
    /// permutations for the built-in union catalog.
    Verify {
        /// Largest n to verify (group enumeration is n! — keep ≤ 9).
        #[arg(long, default_value_t = 7)]
        max_n: u16,
    },
}

fn main() -> Result<()> {
    let code = match Cli::parse().command {
        Command::Run {
            n_spec,
            unions,
            out,
            engine,
            primes,
            deadline,
            allow_identity,
        } => run(
            &n_spec,
            &unions,
            out,
            engine,
            primes,
            deadline,
            allow_identity,
        )?,
        Command::Resume { run_dir, deadline } => resume(&run_dir, deadline)?,
        Command::Verify { max_n } => {
            verify(max_n)?;
            0
        }
    };
    std::process::exit(code);
}

fn parse_n_spec(spec: &str) -> Result<Vec<u16>> {
    let spec = spec.trim();
    if let Some((a, b)) = spec.split_once("..=") {
        let a: u16 = a.trim().parse().context("bad range start")?;
        let b: u16 = b.trim().parse().context("bad range end")?;
        if a > b {
            bail!("empty n range {spec}");
        }
        return Ok((a..=b).collect());
    }
    spec.split(',')
        .map(|t| {
            t.trim()
                .parse::<u16>()
                .with_context(|| format!("bad n value {t:?}"))
        })
        .collect()
}

fn parse_union(spec: &str) -> Result<Vec<CycleTypeTemplate>> {
    spec.split('+')
        .map(|part| {
            part.parse::<CycleTypeTemplate>()
                .map_err(|e| anyhow::anyhow!("bad class {part:?} in union {spec:?}: {e}"))
        })
        .collect()
}

fn utc_now_rfc3339() -> String {
    humantime::format_rfc3339_seconds(SystemTime::now()).to_string()
}

/// Per-job configuration hash — written into checkpoint headers and
/// re-derived on resume from the checkpoint body itself.
fn job_config_hash(
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

struct JobContext<'a> {
    index: &'a PartitionIndex,
    mn: &'a MnEvaluator,
    modular_ctx: Option<&'a ModularContext>,
    prime_list: &'a [Prime31],
}

enum JobOutcome {
    Done {
        diameter: u32,
        stop_radius: u32,
        reachable: usize,
        file: String,
    },
    Suspended {
        checkpoint_file: String,
        committed_radius: u32,
    },
}

#[allow(clippy::too_many_arguments)]
fn execute_job(
    job: &JobContext<'_>,
    templates: &[CycleTypeTemplate],
    union: &ResolvedUnion,
    allow_identity: bool,
    engine: EngineChoice,
    deadline: Option<Instant>,
    resume_body: Option<CheckpointBody>,
    out_dir: &Path,
    run_id: &str,
    run_config_hash: &str,
) -> Result<JobOutcome> {
    let index = job.index;
    let mn = job.mn;
    let n = index.n();
    let slug = format!("g{}", union_slug(templates));
    let job_name = format!("n{n:02}_{slug}");
    let started = utc_now_rfc3339();
    let t0 = Instant::now();
    let spectra = BaseSpectra::build(index, mn, &union.class_ids)?;
    let resumed = resume_body.is_some();
    let suspend_count_before = resume_body.as_ref().map_or(0, |b| b.suspend_count);

    let (run, descriptor) = match engine {
        EngineChoice::Exact => {
            let run = run_exact(index, mn, union)?;
            (run, EngineDescriptor::exact_reference())
        }
        EngineChoice::Modular => {
            let ctx = job.modular_ctx.expect("modular context prepared");
            let options = ModularOptions {
                deadline,
                allow_identity_generator: allow_identity,
                ..Default::default()
            };
            match run_modular_resumable(
                index,
                mn,
                ctx,
                &spectra,
                union,
                &CpuBlocked,
                &options,
                resume_body,
                None,
            )? {
                ModularOutcome::Completed(run, stats) => {
                    let descriptor =
                        EngineDescriptor::modular(job.prime_list, CpuBlocked.name_static(), stats);
                    (run, descriptor)
                }
                ModularOutcome::Suspended(body) => {
                    let ckpt_dir = out_dir.join("checkpoints");
                    std::fs::create_dir_all(&ckpt_dir)?;
                    let path = ckpt_dir.join(format!("{job_name}.ckpt"));
                    let config = job_config_hash(
                        n,
                        &body.resolved_classes,
                        body.allow_identity_generator,
                        &body.primes,
                    );
                    write_checkpoint(&path, &config, index.order_hash(), &body)?;
                    return Ok(JobOutcome::Suspended {
                        checkpoint_file: format!("checkpoints/{job_name}.ckpt"),
                        committed_radius: body.committed_radius,
                    });
                }
            }
        }
    };

    let elapsed = t0.elapsed().as_secs_f64();
    let document = build_result(
        index,
        templates,
        union,
        &spectra,
        &run,
        Some(slug.clone()),
        allow_identity,
        RunMeta {
            run_id: run_id.to_string(),
            started_utc: started,
            finished_utc: utc_now_rfc3339(),
            threads: rayon::current_num_threads() as u32,
            total_wall_s: elapsed,
            config_hash: run_config_hash.to_string(),
            resumed_from_checkpoint: resumed,
            suspend_resume_count: suspend_count_before,
        },
        descriptor,
    );
    let file = format!("{job_name}.json");
    let path = out_dir.join(&file);
    std::fs::write(&path, serde_json::to_string_pretty(&document)?)
        .with_context(|| format!("cannot write {}", path.display()))?;
    // completed: retire any leftover checkpoint
    let ckpt = out_dir.join("checkpoints").join(format!("{job_name}.ckpt"));
    let _ = std::fs::remove_file(&ckpt);
    let _ = std::fs::remove_file(ckpt.with_extension("ckpt.prev"));
    println!(
        "{job_name}: diameter={} stop={} reachable={}/{} ({elapsed:.3}s) -> {}",
        run.diameter,
        run.stop_radius,
        run.reachable_count,
        index.count(),
        path.display()
    );
    Ok(JobOutcome::Done {
        diameter: run.diameter,
        stop_radius: run.stop_radius,
        reachable: run.reachable_count,
        file,
    })
}

trait BackendName {
    fn name_static(&self) -> &'static str;
}
impl BackendName for CpuBlocked {
    fn name_static(&self) -> &'static str {
        use classdiam_core::transform::TransformBackend;
        self.name()
    }
}

#[allow(clippy::too_many_arguments)]
fn run(
    n_spec: &str,
    union_specs: &[String],
    out: Option<PathBuf>,
    engine: EngineChoice,
    prime_count: usize,
    deadline_secs: Option<u64>,
    allow_identity: bool,
) -> Result<i32> {
    let ns = parse_n_spec(n_spec)?;
    let unions: Vec<Vec<CycleTypeTemplate>> = union_specs
        .iter()
        .map(|s| parse_union(s))
        .collect::<Result<_>>()?;
    let prime_list = screening_primes(prime_count);
    let deadline = deadline_secs.map(|s| Instant::now() + Duration::from_secs(s));

    let run_config_hash = {
        let mut hasher = blake3::Hasher::new();
        hasher.update(
            format!("run-v1;n={ns:?};allow_identity={allow_identity};primes={prime_count}")
                .as_bytes(),
        );
        for u in &unions {
            hasher.update(format!("{u:?};").as_bytes());
        }
        hasher.finalize().to_hex().to_string()
    };
    let run_id = format!(
        "{}-{}",
        utc_now_rfc3339().replace(['-', ':'], ""),
        &run_config_hash[..8]
    );
    let out_dir = out.unwrap_or_else(|| PathBuf::from("results").join(&run_id));
    std::fs::create_dir_all(&out_dir)
        .with_context(|| format!("cannot create {}", out_dir.display()))?;

    let mut jobs = Vec::new();
    let mut any_suspended = false;
    for &n in &ns {
        let index = PartitionIndex::build(n)?;
        let mn = MnEvaluator::new(n);
        let modular_ctx = (engine == EngineChoice::Modular)
            .then(|| ModularContext::build(&index, &mn, &prime_list));
        let job = JobContext {
            index: &index,
            mn: &mn,
            modular_ctx: modular_ctx.as_ref(),
            prime_list: &prime_list,
        };
        for templates in &unions {
            let slug = format!("g{}", union_slug(templates));
            let union = match resolve_union(&index, templates, allow_identity) {
                Ok(u) => u,
                Err(e @ ClassdiamError::TemplateDoesNotFit { .. }) => {
                    println!("n{n:02}_{slug}: skipped ({e})");
                    jobs.push(serde_json::json!({
                        "n": n, "union": slug, "status": "skipped", "reason": e.to_string(),
                    }));
                    continue;
                }
                Err(e) => return Err(e.into()),
            };
            match execute_job(
                &job,
                templates,
                &union,
                allow_identity,
                engine,
                deadline,
                None,
                &out_dir,
                &run_id,
                &run_config_hash,
            )? {
                JobOutcome::Done {
                    diameter,
                    stop_radius,
                    reachable,
                    file,
                } => jobs.push(serde_json::json!({
                    "n": n, "union": slug, "status": "done", "file": file,
                    "diameter": diameter, "stop_radius": stop_radius,
                    "reachable": reachable,
                })),
                JobOutcome::Suspended {
                    checkpoint_file,
                    committed_radius,
                } => {
                    any_suspended = true;
                    println!("n{n:02}_{slug}: suspended at committed radius {committed_radius}");
                    jobs.push(serde_json::json!({
                        "n": n, "union": slug, "status": "suspended",
                        "checkpoint": checkpoint_file,
                        "committed_radius": committed_radius,
                    }));
                }
            }
        }
    }

    write_manifest(&out_dir, &run_id, &run_config_hash, any_suspended, jobs)?;
    Ok(if any_suspended { 75 } else { 0 })
}

fn write_manifest(
    out_dir: &Path,
    run_id: &str,
    config_hash: &str,
    any_suspended: bool,
    jobs: Vec<serde_json::Value>,
) -> Result<()> {
    let manifest = serde_json::json!({
        "format": "classdiam/manifest",
        "format_version": 1,
        "run_id": run_id,
        "config_hash_blake3": config_hash,
        "status": if any_suspended { "suspended" } else { "completed" },
        "jobs": jobs,
    });
    std::fs::write(
        out_dir.join("manifest.json"),
        serde_json::to_string_pretty(&manifest)?,
    )?;
    println!(
        "manifest -> {} ({})",
        out_dir.join("manifest.json").display(),
        if any_suspended {
            "SUSPENDED — resume with `classdiam resume`"
        } else {
            "completed"
        }
    );
    Ok(())
}

fn resume(run_dir: &Path, deadline_secs: Option<u64>) -> Result<i32> {
    let ckpt_dir = run_dir.join("checkpoints");
    let mut paths: Vec<PathBuf> = std::fs::read_dir(&ckpt_dir)
        .with_context(|| format!("no checkpoints directory in {}", run_dir.display()))?
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.extension().is_some_and(|e| e == "ckpt"))
        .collect();
    paths.sort();
    if paths.is_empty() {
        bail!("no .ckpt files in {}", ckpt_dir.display());
    }
    let deadline = deadline_secs.map(|s| Instant::now() + Duration::from_secs(s));

    // reload the existing manifest to patch job entries
    let manifest_path = run_dir.join("manifest.json");
    let manifest: serde_json::Value = serde_json::from_str(
        &std::fs::read_to_string(&manifest_path).context("missing manifest.json")?,
    )?;
    let run_id = manifest["run_id"].as_str().unwrap_or("resumed").to_string();
    let run_config_hash = manifest["config_hash_blake3"]
        .as_str()
        .unwrap_or("")
        .to_string();
    let mut jobs: Vec<serde_json::Value> =
        serde_json::from_value(manifest["jobs"].clone()).unwrap_or_default();

    let mut any_suspended = false;
    for path in paths {
        // Parse without expected hashes first (the config is reconstructed
        // from the body), then verify the header hash against the re-derived
        // config and the order hash against the freshly built index.
        let (body, header_config, header_order) = read_checkpoint(&path, None, None)?;
        let expected_config = job_config_hash(
            body.n,
            &body.resolved_classes,
            body.allow_identity_generator,
            &body.primes,
        );
        if header_config != expected_config {
            bail!(
                "{}: checkpoint config hash does not match its own body — refusing",
                path.display()
            );
        }
        let index = PartitionIndex::build(body.n)?;
        if &header_order != index.order_hash() {
            bail!(
                "{}: partition-order hash mismatch — checkpoint from an incompatible version",
                path.display()
            );
        }
        let mn = MnEvaluator::new(body.n);
        let templates: Vec<CycleTypeTemplate> = body
            .resolved_classes
            .iter()
            .map(|parts| {
                let reduced: Vec<u8> = parts.iter().copied().filter(|&x| x >= 2).collect();
                CycleTypeTemplate::new(&reduced).expect("checkpoint classes valid")
            })
            .collect();
        let allow_identity = body.allow_identity_generator;
        let union = resolve_union(&index, &templates, allow_identity)?;
        let prime_list: Vec<Prime31> = body.primes.iter().map(|&p| Prime31(p)).collect();
        let modular_ctx = ModularContext::build(&index, &mn, &prime_list);
        let job = JobContext {
            index: &index,
            mn: &mn,
            modular_ctx: Some(&modular_ctx),
            prime_list: &prime_list,
        };
        let n = body.n;
        let slug = format!("g{}", union_slug(&templates));
        let outcome = execute_job(
            &job,
            &templates,
            &union,
            allow_identity,
            EngineChoice::Modular,
            deadline,
            Some(body),
            run_dir,
            &run_id,
            &run_config_hash,
        )?;
        let entry = match outcome {
            JobOutcome::Done {
                diameter,
                stop_radius,
                reachable,
                file,
            } => serde_json::json!({
                "n": n, "union": slug, "status": "done", "file": file,
                "diameter": diameter, "stop_radius": stop_radius, "reachable": reachable,
            }),
            JobOutcome::Suspended {
                checkpoint_file,
                committed_radius,
            } => {
                any_suspended = true;
                println!("n{n:02}_{slug}: suspended again at committed radius {committed_radius}");
                serde_json::json!({
                    "n": n, "union": slug, "status": "suspended",
                    "checkpoint": checkpoint_file, "committed_radius": committed_radius,
                })
            }
        };
        // patch the matching manifest entry (or append)
        let key = (entry["n"].clone(), entry["union"].clone());
        if let Some(existing) = jobs
            .iter_mut()
            .find(|j| (j["n"].clone(), j["union"].clone()) == key)
        {
            *existing = entry;
        } else {
            jobs.push(entry);
        }
    }

    write_manifest(run_dir, &run_id, &run_config_hash, any_suspended, jobs)?;
    Ok(if any_suspended { 75 } else { 0 })
}

fn verify(max_n: u16) -> Result<()> {
    use classdiam_core::testing::bruteforce as bf;
    use classdiam_core::testing::catalog::{brute_force_affordable, resolve_entry, union_catalog};

    let mut checked = 0usize;
    let mut skipped = 0usize;
    for entry in union_catalog().into_iter().filter(|e| e.n <= max_n) {
        let index = PartitionIndex::build(entry.n)?;
        let union = resolve_entry(&index, &entry);
        if !brute_force_affordable(&index, &union) {
            skipped += 1;
            continue;
        }
        let mn = MnEvaluator::new(entry.n);
        let run = run_exact(&index, &mn, &union)?;
        let generators = bf::materialize_union(&index, &union);
        let by_type = bf::distances_by_type(&index, &bf::bfs_distances(entry.n, &generators));
        if run.distance != by_type {
            bail!(
                "MISMATCH at {}: engine {:?} vs brute force {:?}",
                entry.label,
                run.distance,
                by_type
            );
        }
        checked += 1;
        println!("ok  {}  (diameter {})", entry.label, run.diameter);
    }
    println!(
        "verify: {checked} unions checked against brute-force BFS, {skipped} skipped by cost guard"
    );
    Ok(())
}
