//! `classdiam` CLI: `run`, `resume`, `verify`.
//!
//! Thin shell over [`classdiam_core::driver`] — argument parsing, printing,
//! and exit codes live here; all orchestration (jobs, checkpoints,
//! manifests, hashes) is in the core so Python runs behave identically.

use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

use anyhow::{bail, Result};
use clap::{Parser, Subcommand, ValueEnum};

use classdiam_core::driver::{
    self, BatchConfig, DriverHooks, EngineKind, JobReport, JobStatus,
};
use classdiam_core::partition::CycleTypeTemplate;

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

impl From<EngineChoice> for EngineKind {
    fn from(choice: EngineChoice) -> Self {
        match choice {
            EngineChoice::Modular => EngineKind::Modular,
            EngineChoice::Exact => EngineKind::Exact,
        }
    }
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

/// Per-job line, byte-identical to the historical formats. `resumed`
/// switches the suspension wording ("suspended" vs "suspended again").
fn print_job(report: &JobReport, out_dir: &Path, resumed: bool) {
    match &report.status {
        JobStatus::Done { document, file } => {
            let file = file.as_deref().expect("CLI runs always write files");
            println!(
                "{}: diameter={} stop={} reachable={}/{} ({:.3}s) -> {}",
                report.job_name,
                document.results.diameter_identity_component,
                document.results.stopping.stop_radius,
                document.results.reachable_count,
                report.class_count,
                report.elapsed_s,
                out_dir.join(file).display()
            );
        }
        JobStatus::Suspended {
            committed_radius, ..
        } => {
            let again = if resumed { " again" } else { "" };
            println!(
                "{}: suspended{again} at committed radius {committed_radius}",
                report.job_name
            );
        }
        JobStatus::Skipped { reason } => {
            println!("{}: skipped ({reason})", report.job_name);
        }
    }
}

fn print_manifest_line(out_dir: &Path, any_suspended: bool) {
    println!(
        "manifest -> {} ({})",
        out_dir.join("manifest.json").display(),
        if any_suspended {
            "SUSPENDED — resume with `classdiam resume`"
        } else {
            "completed"
        }
    );
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
    let ns = driver::parse_n_spec(n_spec)?;
    let unions: Vec<Vec<CycleTypeTemplate>> = union_specs
        .iter()
        .map(|s| driver::parse_union(s))
        .collect::<Result<_, _>>()?;
    let deadline = deadline_secs.map(|s| Instant::now() + Duration::from_secs(s));

    // Default out dir results/<run_id> is CLI policy; the driver treats
    // out_dir: None as an in-memory run.
    let config_hash = driver::run_config_hash(&ns, &unions, allow_identity, prime_count);
    let run_id = driver::make_run_id(&config_hash);
    let out_dir = out.unwrap_or_else(|| PathBuf::from("results").join(&run_id));

    let cfg = BatchConfig {
        ns,
        unions,
        engine: engine.into(),
        prime_count,
        deadline,
        allow_identity,
        out_dir: Some(out_dir.clone()),
        run_id: Some(run_id),
    };
    let mut on_job = |report: &JobReport| print_job(report, &out_dir, false);
    let mut hooks = DriverHooks {
        on_job: Some(&mut on_job),
        ..Default::default()
    };
    let outcome = driver::run_batch(&cfg, &mut hooks)?;
    print_manifest_line(&out_dir, outcome.any_suspended);
    Ok(if outcome.any_suspended { 75 } else { 0 })
}

fn resume(run_dir: &Path, deadline_secs: Option<u64>) -> Result<i32> {
    let deadline = deadline_secs.map(|s| Instant::now() + Duration::from_secs(s));
    let mut on_job = |report: &JobReport| print_job(report, run_dir, true);
    let mut hooks = DriverHooks {
        on_job: Some(&mut on_job),
        ..Default::default()
    };
    let outcome = driver::resume_batch(run_dir, deadline, &mut hooks)?;
    print_manifest_line(run_dir, outcome.any_suspended);
    Ok(if outcome.any_suspended { 75 } else { 0 })
}

fn verify(max_n: u16) -> Result<()> {
    use classdiam_core::chars::MnEvaluator;
    use classdiam_core::engine::exact::run_exact;
    use classdiam_core::partition::PartitionIndex;
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
