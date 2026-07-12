//! `classdiam._core` — native half of the `classdiam` Python package.
//!
//! Design (docs/design/01-architecture.md §12): compute runs with the GIL
//! released (`Python::detach`); an internal per-committed-layer hook
//! re-attaches to (a) check signals so Ctrl-C maps to a [`CancelToken`]
//! and a clean checkpointed suspension, and (b) invoke the user's
//! `progress` callable, throttled. All payloads cross the FFI as JSON
//! strings serialized from the same serde types as the result files —
//! one schema, three consumers (file, CLI, Python). The pure-Python layer
//! (`classdiam.helpers`) parses them into dicts and friendly wrappers.

use std::path::PathBuf;
use std::time::{Duration, Instant};

use pyo3::create_exception;
use pyo3::exceptions::{PyKeyboardInterrupt, PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyDict;
use serde::Deserialize;

use classdiam_core::driver::{
    self, BatchConfig, BatchOutcome, CancelToken, DriverHooks, EngineKind, JobOptions, JobStatus,
    ProgressEvent,
};
use classdiam_core::partition::CycleTypeTemplate;
use classdiam_core::ClassdiamError as CoreError;

create_exception!(
    classdiam,
    ClassdiamError,
    PyRuntimeError,
    "The classdiam engine failed (certification, checkpoint, or I/O)."
);

/// Input-shaped core errors become ValueError; the rest are engine faults.
fn map_core_err(e: CoreError) -> PyErr {
    match e {
        CoreError::InvalidSpec { .. }
        | CoreError::TemplateDoesNotFit { .. }
        | CoreError::MalformedTemplate { .. }
        | CoreError::IdentityGenerator { .. }
        | CoreError::EmptyUnion { .. }
        | CoreError::UnsupportedN { .. }
        | CoreError::InvalidRunDir { .. } => PyValueError::new_err(e.to_string()),
        other => ClassdiamError::new_err(other.to_string()),
    }
}

fn engine_kind(name: &str) -> PyResult<EngineKind> {
    match name {
        "modular" => Ok(EngineKind::Modular),
        "exact" => Ok(EngineKind::Exact),
        other => Err(PyValueError::new_err(format!(
            "unknown engine {other:?} (expected \"modular\" or \"exact\")"
        ))),
    }
}

fn templates_from_classes(classes: &[Vec<u8>]) -> PyResult<Vec<CycleTypeTemplate>> {
    classes
        .iter()
        .map(|parts| CycleTypeTemplate::new(parts).map_err(map_core_err))
        .collect()
}

fn deadline_instant(deadline_s: Option<f64>) -> Option<Instant> {
    deadline_s.map(|s| Instant::now() + Duration::from_secs_f64(s.max(0.0)))
}

fn build_pool(threads: Option<usize>) -> PyResult<Option<rayon::ThreadPool>> {
    threads
        .map(|t| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(t)
                .build()
                .map_err(|e| PyValueError::new_err(format!("cannot build thread pool: {e}")))
        })
        .transpose()
}

fn event_dict<'py>(py: Python<'py>, event: &ProgressEvent) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    d.set_item("n", event.n)?;
    d.set_item("job_name", &event.job_name)?;
    d.set_item("radius", event.radius)?;
    d.set_item("new", event.new_count)?;
    d.set_item("support", event.support_size)?;
    d.set_item("reachable", event.reachable)?;
    d.set_item("job_index", event.job_index)?;
    d.set_item("job_total", event.job_total)?;
    Ok(d)
}

/// Run `body` with the GIL released, marshaling progress events back into
/// Python (signal check + optional user callback, throttled). Returns the
/// outcome plus whether a KeyboardInterrupt was converted into a clean
/// cancellation; any other callback exception is re-raised.
fn execute<F>(
    py: Python<'_>,
    progress: Option<Py<PyAny>>,
    pool: Option<&rayon::ThreadPool>,
    min_progress_ms: u64,
    body: F,
) -> PyResult<(BatchOutcome, bool)>
where
    F: FnOnce(&mut DriverHooks<'_>) -> Result<BatchOutcome, CoreError> + Send,
{
    let cancel = CancelToken::new();
    let hook_cancel = cancel.clone();

    let (result, pending_err) = py.detach(move || {
        let throttle = Duration::from_millis(min_progress_ms);
        let mut last: Option<Instant> = None;
        let mut pending_err: Option<PyErr> = None;
        // Installed even without a user callback: this is where Ctrl-C is
        // observed while the GIL is released. Latency is one committed
        // radius (same granularity as the deadline guard).
        let mut progress_hook = |event: &ProgressEvent| {
            if pending_err.is_some() || last.is_some_and(|t| t.elapsed() < throttle) {
                return;
            }
            last = Some(Instant::now());
            Python::attach(|py| {
                if let Err(e) = py.check_signals() {
                    hook_cancel.cancel();
                    pending_err.get_or_insert(e);
                    return;
                }
                if let Some(cb) = progress.as_ref() {
                    let call = event_dict(py, event).and_then(|d| cb.call1(py, (d,)));
                    if let Err(e) = call {
                        hook_cancel.cancel();
                        pending_err.get_or_insert(e);
                    }
                }
            });
        };
        let mut hooks = DriverHooks {
            progress: Some(&mut progress_hook),
            cancel: Some(cancel),
            on_job: None,
        };
        let result = match pool {
            Some(pool) => pool.install(|| body(&mut hooks)),
            None => body(&mut hooks),
        };
        (result, pending_err)
    });

    let outcome = result.map_err(map_core_err)?;
    let interrupted = match pending_err {
        Some(e) if e.is_instance_of::<PyKeyboardInterrupt>(py) => true,
        Some(e) => return Err(e),
        None => false,
    };
    Ok((outcome, interrupted))
}

/// The outcome JSON handed to the helper layer. Done jobs embed the full
/// result document so in-memory runs never re-read files.
fn outcome_json(outcome: &BatchOutcome, interrupted: bool) -> String {
    let jobs: Vec<serde_json::Value> = outcome
        .reports
        .iter()
        .map(|report| {
            let mut entry = serde_json::json!({
                "n": report.n,
                "slug": report.slug,
                "job_name": report.job_name,
                "elapsed_s": report.elapsed_s,
                "class_count": report.class_count,
            });
            match &report.status {
                JobStatus::Done { document, file } => {
                    entry["status"] = "done".into();
                    entry["file"] = file.clone().map_or(serde_json::Value::Null, Into::into);
                    entry["document"] =
                        serde_json::to_value(document.as_ref()).expect("document serializes");
                }
                JobStatus::Suspended {
                    checkpoint_path,
                    committed_radius,
                } => {
                    entry["status"] = "suspended".into();
                    entry["checkpoint"] = checkpoint_path
                        .as_ref()
                        .map_or(serde_json::Value::Null, |p| p.display().to_string().into());
                    entry["committed_radius"] = (*committed_radius).into();
                }
                JobStatus::Skipped { reason } => {
                    entry["status"] = "skipped".into();
                    entry["reason"] = reason.clone().into();
                }
            }
            entry
        })
        .collect();
    serde_json::json!({
        "run_id": outcome.run_id,
        "config_hash": outcome.config_hash,
        "out_dir": outcome.out_dir.as_ref().map(|p| p.display().to_string()),
        "any_suspended": outcome.any_suspended,
        "cancelled": outcome.cancelled,
        "interrupted": interrupted,
        "jobs": jobs,
    })
    .to_string()
}

/// Batch configuration as the helper layer sends it.
#[derive(Deserialize)]
struct PyBatchConfig {
    ns: Vec<u16>,
    /// Unions as lists of classes, each class a list of parts >= 2.
    unions: Vec<Vec<Vec<u8>>>,
    #[serde(default = "default_engine")]
    engine: String,
    #[serde(default = "default_primes")]
    primes: usize,
    #[serde(default)]
    deadline_s: Option<f64>,
    #[serde(default)]
    allow_identity: bool,
    #[serde(default)]
    out_dir: Option<String>,
}

fn default_engine() -> String {
    "modular".into()
}

fn default_primes() -> usize {
    3
}

impl PyBatchConfig {
    fn into_batch_config(self) -> PyResult<BatchConfig> {
        let unions = self
            .unions
            .iter()
            .map(|classes| templates_from_classes(classes))
            .collect::<PyResult<Vec<_>>>()?;
        Ok(BatchConfig {
            ns: self.ns,
            unions,
            engine: engine_kind(&self.engine)?,
            prime_count: self.primes,
            deadline: deadline_instant(self.deadline_s),
            allow_identity: self.allow_identity,
            out_dir: self.out_dir.map(PathBuf::from),
            run_id: None,
        })
    }
}

#[pyfunction]
#[pyo3(signature = (cfg_json, progress=None, threads=None, min_progress_ms=250))]
fn run_batch(
    py: Python<'_>,
    cfg_json: &str,
    progress: Option<Py<PyAny>>,
    threads: Option<usize>,
    min_progress_ms: u64,
) -> PyResult<String> {
    let cfg: PyBatchConfig = serde_json::from_str(cfg_json)
        .map_err(|e| PyValueError::new_err(format!("bad run config: {e}")))?;
    let cfg = cfg.into_batch_config()?;
    let pool = build_pool(threads)?;
    let (outcome, interrupted) = execute(py, progress, pool.as_ref(), min_progress_ms, move |h| {
        driver::run_batch(&cfg, h)
    })?;
    Ok(outcome_json(&outcome, interrupted))
}

#[pyfunction]
#[pyo3(signature = (run_dir, deadline_s=None, progress=None, threads=None, min_progress_ms=250))]
fn resume_batch(
    py: Python<'_>,
    run_dir: PathBuf,
    deadline_s: Option<f64>,
    progress: Option<Py<PyAny>>,
    threads: Option<usize>,
    min_progress_ms: u64,
) -> PyResult<String> {
    let deadline = deadline_instant(deadline_s);
    let pool = build_pool(threads)?;
    let (outcome, interrupted) = execute(py, progress, pool.as_ref(), min_progress_ms, move |h| {
        driver::resume_batch(&run_dir, deadline, h)
    })?;
    Ok(outcome_json(&outcome, interrupted))
}

#[pyfunction]
fn parse_n_spec(spec: &str) -> PyResult<Vec<u16>> {
    driver::parse_n_spec(spec).map_err(map_core_err)
}

#[pyfunction]
fn parse_union_spec(spec: &str) -> PyResult<Vec<Vec<u16>>> {
    // u16 parts: pyo3 would map Vec<u8> to Python `bytes`, not a list.
    Ok(driver::parse_union(spec)
        .map_err(map_core_err)?
        .into_iter()
        .map(|t| t.parts().iter().map(|&p| u16::from(p)).collect())
        .collect())
}

/// Options for `Session.run_union` as the helper layer sends them.
#[derive(Deserialize)]
struct PySessionJobOpts {
    #[serde(default)]
    label: Option<String>,
    #[serde(default)]
    allow_identity: bool,
    #[serde(default)]
    deadline_s: Option<f64>,
    #[serde(default)]
    checkpoint_dir: Option<String>,
    #[serde(default)]
    result_dir: Option<String>,
}

/// Per-`n` reusable state (partition index + character tables): build once,
/// run many unions. Frozen: all mutation lives in the core's interior.
#[pyclass(frozen)]
struct Session {
    inner: driver::Session,
    pool: Option<rayon::ThreadPool>,
}

#[pymethods]
impl Session {
    #[new]
    #[pyo3(signature = (n, engine="modular", primes=3, threads=None))]
    fn new(n: u16, engine: &str, primes: usize, threads: Option<usize>) -> PyResult<Self> {
        let inner = driver::Session::new(n, engine_kind(engine)?, primes).map_err(map_core_err)?;
        Ok(Self {
            inner,
            pool: build_pool(threads)?,
        })
    }

    #[pyo3(signature = (classes, opts_json, progress=None, min_progress_ms=250))]
    fn run_union(
        &self,
        py: Python<'_>,
        classes: Vec<Vec<u8>>,
        opts_json: &str,
        progress: Option<Py<PyAny>>,
        min_progress_ms: u64,
    ) -> PyResult<String> {
        let opts: PySessionJobOpts = serde_json::from_str(opts_json)
            .map_err(|e| PyValueError::new_err(format!("bad job options: {e}")))?;
        let templates = templates_from_classes(&classes)?;
        // Same hash/run-id derivation as a one-job batch, so optional
        // checkpoints stay CLI-resumable.
        let config_hash = driver::run_config_hash(
            &[self.inner.index().n()],
            std::slice::from_ref(&templates),
            opts.allow_identity,
            self.inner.primes().len(),
        );
        let run_id = driver::make_run_id(&config_hash);
        let job_opts = JobOptions {
            allow_identity: opts.allow_identity,
            deadline: deadline_instant(opts.deadline_s),
            result_dir: opts.result_dir.map(PathBuf::from),
            checkpoint_dir: opts.checkpoint_dir.map(PathBuf::from),
            label: opts.label,
            run_id: run_id.clone(),
            run_config_hash: config_hash.clone(),
            resume: None,
            job_index: 0,
            job_total: 1,
        };
        let inner = &self.inner;
        let (outcome, interrupted) = execute(
            py,
            progress,
            self.pool.as_ref(),
            min_progress_ms,
            move |hooks| {
                let report = inner.run_union(&templates, job_opts, hooks)?;
                let any_suspended = matches!(report.status, JobStatus::Suspended { .. });
                Ok(BatchOutcome {
                    run_id,
                    config_hash,
                    out_dir: None,
                    any_suspended,
                    cancelled: hooks.cancel.as_ref().is_some_and(CancelToken::is_cancelled),
                    reports: vec![report],
                })
            },
        )?;
        Ok(outcome_json(&outcome, interrupted))
    }

    /// The canonical partition order (reduced form, parts >= 2) as JSON.
    fn partition_order_json(&self) -> String {
        let index = self.inner.index();
        let reduced: Vec<Vec<u8>> = index
            .partitions()
            .iter()
            .map(|p| p.parts().iter().copied().filter(|&x| x >= 2).collect())
            .collect();
        serde_json::to_string(&reduced).expect("partitions serialize")
    }

    #[getter]
    fn n(&self) -> u16 {
        self.inner.index().n()
    }

    #[getter]
    fn class_count(&self) -> usize {
        self.inner.index().count()
    }
}

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(run_batch, m)?)?;
    m.add_function(wrap_pyfunction!(resume_batch, m)?)?;
    m.add_function(wrap_pyfunction!(parse_n_spec, m)?)?;
    m.add_function(wrap_pyfunction!(parse_union_spec, m)?)?;
    m.add_class::<Session>()?;
    m.add("ClassdiamError", m.py().get_type::<ClassdiamError>())?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
