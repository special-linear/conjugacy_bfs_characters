# Running classdiam on Kaggle

The `classdiam` Python package is a self-contained abi3 wheel (Python ≥ 3.9,
manylinux x86_64) — no dependencies, pandas optional. The workflow below
fits Kaggle's session model: bounded wall-clock budgets, checkpointed
suspension, resume across sessions.

## 1. Install

**Internet ON** (default for non-competition notebooks) — install straight
from the GitHub release:

```python
!pip install -q https://github.com/<owner>/<repo>/releases/download/v0.1.0/classdiam-0.1.0-cp39-abi3-manylinux_2_17_x86_64.whl
```

(Copy the exact wheel URL from the release page — the manylinux tag may
carry a dual `manylinux_2_17…manylinux2014` spelling.)

**Internet OFF** (competition notebooks) — upload the wheel once as a
Kaggle Dataset, attach it to the notebook, then:

```python
!pip install -q --no-deps /kaggle/input/<your-wheel-dataset>/classdiam-*.whl
```

Check:

```python
import classdiam
classdiam.__version__
```

## 2. Compute

```python
import classdiam

res = classdiam.run(n="6..=30", unions=["2"])   # transpositions in S_6..S_30
r = res[0]                                       # one RunResult per n
r.diameter                                       # 5 for S_6
r.distances                                      # per conjugacy class, canonical order
r.partitions                                     # that order (reduced cycle types)
r.distance_of([3, 2])                            # look up one class
r.support(4)                                     # classes with a length-4 factorization
r.to_dataframe()                                 # pandas is preinstalled on Kaggle
```

Union grammar (same as the CLI): classes joined `+`, parts joined `,` —
`"2"` = transpositions, `"3+2,2"` = 3-cycles ∪ double transpositions.
Nested lists work too: `[[3], [2, 2]]`.

`RunResult.raw` is the full versioned result document (schema
`classdiam/result` v1, see `docs/output_schema.md`) — identical to the JSON
files a run directory holds, so everything is dict-processable.

## 3. Long runs: deadline + checkpoints

A Kaggle session has a hard wall clock (~9 h CPU/GPU sessions, 20 min idle
cutoff; interactive sessions die on browser close). Give the run a budget
comfortably inside your session and a directory under `/kaggle/working`
(the only path persisted on Save & Run All):

```python
try:
    res = classdiam.run(
        n=35, unions=["2"],
        deadline_s=6600,                       # ~110 min budget
        out_dir="/kaggle/working/results",
    )
except classdiam.Suspended as s:
    print(s)                                    # which jobs suspended, where
    res = s.completed                           # jobs that did finish
```

On expiry every unfinished job writes a checkpoint of its last fully
certified radius and `classdiam.Suspended` is raised. The same happens on
**notebook interrupt / Ctrl-C** (`s.reason == "interrupt"`); interruption is
observed between radii, so at large `n` expect up to one radius of latency.
The exact engine (`engine="exact"`) has no checkpoints — it is a small-`n`
oracle whose runs take seconds.

Progress for long radii:

```python
res = classdiam.run(n=35, unions=["2"], deadline_s=6600,
                    out_dir="/kaggle/working/results",
                    progress=lambda e: print(e))   # one dict per committed radius
```

## 4. Resume in the next session

Save the notebook so `/kaggle/working/results` becomes its output. In the
next session, attach that output as an input dataset. Inputs are
**read-only**, so copy the run directory back into `/kaggle/working` first:

```python
!cp -r /kaggle/input/<prev-notebook>/results /kaggle/working/results

res = classdiam.resume("/kaggle/working/results", deadline_s=6600)
```

`resume` validates the checkpoint hashes (configuration + partition-order
version), continues from the committed radius, and raises `Suspended` again
if the new deadline also expires — repeat over as many sessions as needed.
Results accumulated so far:

```python
all_results = classdiam.load_results("/kaggle/working/results")
```

Run directories are fully interchangeable with the CLI: a directory
suspended on Kaggle can be resumed locally with
`classdiam resume results/…` and vice versa.

## 5. Many unions at one n: Session

`classdiam.Session` builds the per-`n` character tables once and reuses
them across unions — much cheaper than repeated `run()` calls at the same
`n`:

```python
s = classdiam.Session(30)                # engine="modular", primes=3
a = s.run_union("2")
b = s.run_union([[3], [2, 2]])
s.partition_order                        # shared canonical order
```

`Session.run_union` accepts `deadline_s` + `checkpoint_dir` (and optional
`result_dir`); a suspended union is resumable with
`classdiam.resume(<parent dir>)`.

## 6. Threads and determinism

The engine uses all vCPUs by default (`threads=` limits it). Results are
bit-identical across thread counts — the test suite proves it — so Kaggle
results reproduce locally.

## Not yet available

`estimate()` / planner-driven auto-`n` (design doc §12) arrive with the P3
planner; until then choose `n` and budgets manually.
