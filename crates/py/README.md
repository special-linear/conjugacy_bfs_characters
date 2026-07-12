# classdiam

Exact conjugacy-class BFS on symmetric groups via characters: distances,
layer structure, and diameters of Cayley graphs `Cay(S_n, U)` for unions `U`
of conjugacy classes — without enumerating `n!` permutations. The engine
works on the p(n) conjugacy classes (partitions of n), using modular
character screening with rigorous per-radius certification; results are
exact, not probabilistic.

## Install

Wheels (Linux x86_64 + Windows x86_64, Python ≥ 3.9) are attached to GitHub
releases:

```python
!pip install -q https://github.com/<owner>/<repo>/releases/download/vX.Y.Z/<wheel>.whl
```

On an internet-off Kaggle notebook, attach a dataset containing the wheel:

```python
!pip install -q --no-deps /kaggle/input/<your-wheel-dataset>/classdiam-*.whl
```

## Use

```python
import classdiam

# transpositions in S_6..S_20: one RunResult per n
res = classdiam.run(n="6..=20", unions=["2"])
res[0].diameter          # 5  (S_6, transpositions)
res[0].distances         # distance per conjugacy class, canonical order
res[0].to_dataframe()    # pandas DataFrame, one row per class

# unions of classes: "3+2,2" = 3-cycles ∪ double transpositions
classdiam.run(n=12, unions=["3+2,2", [[2]]])
```

Long computations: give a wall-clock budget and a directory; on expiry (or
notebook interrupt) the run suspends into resumable checkpoints:

```python
try:
    res = classdiam.run(n=35, unions=["2"], deadline_s=6600,
                        out_dir="/kaggle/working/results")
except classdiam.Suspended as s:
    res = s.completed            # finished jobs
# next session:
res = classdiam.resume("/kaggle/working/results", deadline_s=6600)
```

Every result is a versioned JSON document (`RunResult.raw`); run directories
are fully interchangeable with the `classdiam` CLI. See `docs/kaggle.md` in
the repository for the complete Kaggle workflow.
