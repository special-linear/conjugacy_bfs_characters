//! `classdiam` CLI (P1 surface: `run`, `verify`).

use std::path::PathBuf;
use std::time::{Instant, SystemTime};

use anyhow::{bail, Context, Result};
use clap::{Parser, Subcommand};

use classdiam_core::chars::MnEvaluator;
use classdiam_core::engine::exact::run_exact;
use classdiam_core::partition::{CycleTypeTemplate, PartitionIndex};
use classdiam_core::report::{build_result, union_slug, RunMeta};
use classdiam_core::spectra::resolve_union;
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
        /// Permit the identity class as a generator (spec §5.3).
        #[arg(long)]
        allow_identity: bool,
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
    match Cli::parse().command {
        Command::Run {
            n_spec,
            unions,
            out,
            allow_identity,
        } => run(&n_spec, &unions, out, allow_identity),
        Command::Verify { max_n } => verify(max_n),
    }
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

fn run(
    n_spec: &str,
    union_specs: &[String],
    out: Option<PathBuf>,
    allow_identity: bool,
) -> Result<()> {
    let ns = parse_n_spec(n_spec)?;
    let unions: Vec<Vec<CycleTypeTemplate>> = union_specs
        .iter()
        .map(|s| parse_union(s))
        .collect::<Result<_>>()?;

    let config_hash = {
        let mut hasher = blake3::Hasher::new();
        hasher.update(format!("v1;n={ns:?};allow_identity={allow_identity};").as_bytes());
        for u in &unions {
            hasher.update(format!("{u:?};").as_bytes());
        }
        hasher.finalize().to_hex().to_string()
    };
    let run_id = format!(
        "{}-{}",
        utc_now_rfc3339().replace(['-', ':'], ""),
        &config_hash[..8]
    );
    let out_dir = out.unwrap_or_else(|| PathBuf::from("results").join(&run_id));
    std::fs::create_dir_all(&out_dir)
        .with_context(|| format!("cannot create {}", out_dir.display()))?;

    let mut jobs = Vec::new();
    for &n in &ns {
        let index = PartitionIndex::build(n)?;
        let mn = MnEvaluator::new(n);
        for templates in &unions {
            let slug = format!("g{}", union_slug(templates));
            let job_name = format!("n{n:02}_{slug}");
            let union = match resolve_union(&index, templates, allow_identity) {
                Ok(u) => u,
                Err(e @ ClassdiamError::TemplateDoesNotFit { .. }) => {
                    println!("{job_name}: skipped ({e})");
                    jobs.push(serde_json::json!({
                        "n": n, "union": slug, "status": "skipped", "reason": e.to_string(),
                    }));
                    continue;
                }
                Err(e) => return Err(e.into()),
            };

            let started = utc_now_rfc3339();
            let t0 = Instant::now();
            let spectra =
                classdiam_core::spectra::BaseSpectra::build(&index, &mn, &union.class_ids)?;
            let run = run_exact(&index, &mn, &union)?;
            let elapsed = t0.elapsed().as_secs_f64();

            let document = build_result(
                &index,
                templates,
                &union,
                &spectra,
                &run,
                Some(slug.clone()),
                allow_identity,
                RunMeta {
                    run_id: run_id.clone(),
                    started_utc: started,
                    finished_utc: utc_now_rfc3339(),
                    threads: 1,
                    total_wall_s: elapsed,
                    config_hash: config_hash.clone(),
                },
            );
            let path = out_dir.join(format!("{job_name}.json"));
            std::fs::write(&path, serde_json::to_string_pretty(&document)?)
                .with_context(|| format!("cannot write {}", path.display()))?;
            println!(
                "{job_name}: diameter={} stop={} reachable={}/{} ({elapsed:.3}s) -> {}",
                run.diameter,
                run.stop_radius,
                run.reachable_count,
                index.count(),
                path.display()
            );
            jobs.push(serde_json::json!({
                "n": n, "union": slug, "status": "done",
                "file": format!("{job_name}.json"),
                "diameter": run.diameter, "stop_radius": run.stop_radius,
                "elapsed_s": elapsed,
            }));
        }
    }

    let manifest = serde_json::json!({
        "format": "classdiam/manifest",
        "format_version": 1,
        "run_id": run_id,
        "config_hash_blake3": config_hash,
        "status": "completed",
        "jobs": jobs,
    });
    std::fs::write(
        out_dir.join("manifest.json"),
        serde_json::to_string_pretty(&manifest)?,
    )?;
    println!("manifest -> {}", out_dir.join("manifest.json").display());
    Ok(())
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
