# classdiam

Exact distances, layer structure, and diameters of conjugacy-invariant Cayley graphs of
symmetric groups, computed via characters — without constructing the Cayley graph.

For a generating set `U ⊆ S_n` that is a union of conjugacy classes (given as cycle
types), multiplication by the class sum `K_U` is diagonal in the irreducible-character
basis, so the number of length-`r` factorizations of any permutation of cycle type `ν` is

```
a_r(ν) = (1/n!) · Σ_ρ f_ρ · χ^ρ(ν) · θ_ρ(U)^r,
```

and `ν` is reachable in exactly `r` steps iff `a_r(ν) > 0`. Distance from the identity is
constant on cycle types, so the "reduced BFS" over the `p(n)` cycle types is exact.
The full mathematical specification is in
[`notes/character_method_cayley_diameters.md`](notes/character_method_cayley_diameters.md);
the engineering design in [`docs/design/`](docs/design/).

## What it computes

Per `(n, generating union)`:

- distance from the identity for every cycle type (`-1` = unreachable);
- newly-reached cycle types and exact-length supports at every radius;
- diameter of the identity component, reachability/parity metadata;
- all recorded as versioned, self-describing JSON embedding the canonical partition
  order (`lex_desc_full_parts_v1`) and full run metadata.

All zero/positive decisions are exact: modular screening (31-bit primes) with a
rigorous per-radius certification gate (word-count bound → resident CRT → exact
evaluation). No floating point, no probabilistic step, no uncertified stopping.

## Non-goals

- **No geodesic witnesses**: the output is distance by cycle type, not explicit
  shortest factorizations (spec §18).
- **No split `A_n` classes**: generators are full `S_n` conjugacy classes by
  construction; half of a split `A_n` class cannot be expressed (spec §17.2).
- **No non-conjugacy-invariant generators**: the cycle-type reduction would be
  invalid (spec §17.3).
- **No factorization counts in output**: `a_r(ν)` values are internal only.

## Workspace

- `crates/core` — `classdiam-core`: partitions, Murnaghan–Nakayama evaluator,
  spectra, modular/exact arithmetic, transform backends, diameter engine,
  checkpointing, validation, and the shared batch/resume driver. No file
  I/O in math modules.
- `crates/cli` — `classdiam` binary: `run`, `resume`, `verify` (a thin shell
  over the core driver; `estimate` et al. arrive with the P3 planner).
- `crates/py` — `classdiam` Python package (PyO3/maturin, abi3 ≥ 3.9) for
  Kaggle notebooks; same driver, interchangeable run directories.
- `crates/gpu` (reserved) — GPU `TransformBackend` implementation.
- `fixtures/` — committed SymPy-generated ground truth; `tools/` — fixture and
  adversarial-case generators (Python).

## Quickstart

```
cargo test --workspace                 # T0+T1 test suite
cargo run -p classdiam-cli --profile ci -- run -n 6..=30 -u 2 -u "3+2,2" -o results/quick
cargo run -p classdiam-cli --profile ci -- run -n 35 -u 2 --deadline 6900   # suspends on expiry (exit 75)
cargo run -p classdiam-cli --profile ci -- resume results/<run_id>          # continues from checkpoints
cargo run -p classdiam-cli --profile ci -- verify --max-n 8                 # brute-force cross-check
```

The default engine is the modular production engine (31-bit prime screening
with a rigorous per-radius certification gate); `--engine exact` selects the
big-integer reference engine. Both produce identical results — the test
suite proves it on every catalog union and under adversarial prime choices.

Union grammar: classes joined by `+`, parts within a class by `,` — `-u 2` is the
transposition class in every `S_n`, `-u "3+2,2"` is the union of 3-cycles and double
transpositions. Cycle types are written without fixed points and padded per `n`.

## Python (Kaggle)

The same engine ships as the `classdiam` Python package (abi3 wheels for
Linux x86_64 and Windows on GitHub releases — see
[`docs/kaggle.md`](docs/kaggle.md) for the full notebook workflow):

```python
!pip install -q https://github.com/special-linear/conjugacy_bfs_characters/releases/download/v0.1.0/classdiam-0.1.0-cp39-abi3-manylinux_2_17_x86_64.manylinux2014_x86_64.whl

import classdiam
res = classdiam.run(n="6..=30", unions=["2"])       # one RunResult per n
res[0].diameter, res[0].distances                    # exact, canonical order
res[0].to_dataframe()                                # pandas view

# long runs: budget + checkpoints, resumable next session (or by the CLI)
try:
    res = classdiam.run(n=35, unions=["2"], deadline_s=6600,
                        out_dir="/kaggle/working/results")
except classdiam.Suspended:
    pass                                             # later: classdiam.resume(...)
```

Build locally with `pip install ./crates/py` (needs Rust; maturin arrives
via pip build isolation). Run directories, hashes, and the JSON schema are
identical between the CLI and Python — each can resume the other's runs.
