# Design: `classdiam` — Conjugacy-Invariant Cayley Diameters of S_n via Characters

A Rust workspace implementing the spec in `notes/character_method_cayley_diameters.md` (spectral powers of class sums, support-only output, modular screening with certified zeros). Scope of this document: crate architecture, public API, data model, I/O, orchestration, checkpointing, CLI, and the future Python surface.

---

## 0. Design principles (derived from fixed requirements)

1. **Support-only hot path.** Since word counts `a_r(ν)` are not output, the engine never divides by `n!` in the hot loop: `a_r(ν) ≠ 0 ⇔ N_r(ν) ≠ 0 (mod p)` because `gcd(n!, p) = 1`. Division/divisibility checks live only in the validation layer.
2. **Every zero is certified.** Distances and the stopping rule (spec §5.2, Failures 4/9) require exact supports at every radius — not only at the last one (an uncertified hidden positive at radius `r` would record distance `r+2`, silently wrong). The certification gate therefore runs every radius, and is engineered to be nearly free (§7.4).
3. **One canonical partition order, owned by the library**, never by a backend (spec §3, §19.3, Failure 7). Every serialized artifact embeds the order and its hash.
4. **The transform is a trait**; the math engine never knows whether residue GEMMs run on rayon or a future GPU.
5. **Interrupt-anywhere.** The Kaggle 2 h wall clock makes suspension a first-class state, not an error.
6. **Pure Rust, Windows-first dev.** No GMP/MSYS2. `num-bigint` by default.

---

## 1. Repository layout and Cargo workspace

```
d:\Math\self\cayleypy\conjugacy_bfs_characters\
├── Cargo.toml                  # [workspace] members = ["crates/core", "crates/cli"]
├── rust-toolchain.toml         # pinned stable, edition 2024, MSRV 1.85
├── notes/                      # (existing) spec
├── fixtures/                   # committed SymPy-generated ground truth (JSON)
│   ├── characters/             #   char_n{04..14}.json: sampled (rho, nu, value)
│   ├── degrees/                #   deg_n{04..30}.json: full degree vectors
│   └── bfs/                    #   bfs_n{04..09}_{slug}.json: brute-force layer data
├── tools/
│   ├── gen_fixtures.py         # SymPy generator (run once, output committed)
│   └── check_result.py         # independent sanity-checker for result JSON
├── crates/
│   ├── core/                   # classdiam-core  (lib; no file I/O in math modules)
│   │   ├── Cargo.toml
│   │   ├── src/
│   │   │   ├── lib.rs
│   │   │   ├── error.rs        # ClassdiamError (thiserror)
│   │   │   ├── partition/
│   │   │   │   ├── mod.rs      # Partition, PartitionId, sign, z_lambda, transpose
│   │   │   │   ├── gen.rs      # canonical-order generation
│   │   │   │   ├── index.rs    # PartitionIndex, order hash
│   │   │   │   └── template.rs # CycleTypeTemplate, padding/validation
│   │   │   ├── arith/
│   │   │   │   ├── modp.rs     # Prime31, ModCtx (Barrett), ResidueMat, prime table
│   │   │   │   ├── exact.rs    # ExactInt abstraction over num-bigint / malachite
│   │   │   │   └── bounds.rs   # coefficient bounds B_r(nu), CRT bound logic
│   │   │   ├── chars/
│   │   │   │   ├── mn.rs       # Murnaghan–Nakayama trie DFS (mod-p multi-prime + exact)
│   │   │   │   ├── degrees.rs  # hook-length degrees (exact BigUint)
│   │   │   │   ├── source.rs   # CharacterSource + TableBlock + BlockSink
│   │   │   │   ├── memtable.rs # resident multi-prime table source
│   │   │   │   └── streamed.rs # per-sweep MN streaming source; eigen-group builder
│   │   │   ├── spectra.rs      # exact omega/theta, active rows, eigenvalue grouping
│   │   │   ├── transform/
│   │   │   │   ├── backend.rs  # TransformBackend trait (the GPU seam)
│   │   │   │   └── cpu.rs      # reference + blocked rayon CPU backends
│   │   │   ├── engine/
│   │   │   │   ├── mod.rs      # per-n job driver
│   │   │   │   ├── state.rs    # UnionState, phases, layer records
│   │   │   │   ├── batch.rs    # union batching + column compaction
│   │   │   │   ├── certify.rs  # tiered zero-certification gate
│   │   │   │   └── invariants.rs # cheap always-on runtime checks (word-count mod p)
│   │   │   ├── orchestrate/
│   │   │   │   ├── plan.rs     # per-n ResourcePlan, mode selection
│   │   │   │   └── estimate.rs # memory/time model + calibration
│   │   │   ├── report/
│   │   │   │   ├── schema.rs   # serde types for result JSON (versioned)
│   │   │   │   └── manifest.rs # run manifest + JSONL events
│   │   │   ├── checkpoint.rs   # postcard body + header, atomic write, resume validation
│   │   │   ├── progress.rs     # ProgressSink trait, CancelToken, DeadlineGuard
│   │   │   └── validate.rs     # spec §9 invariants (full-row r=0/r=1, divisibility, …)
│   │   ├── tests/              # integration tests (fixtures, brute-force BFS, proptest)
│   │   └── benches/            # criterion: MN gen, modular GEMM, certification
│   ├── cli/                    # classdiam-cli → binary `classdiam`
│   │   └── src/{main.rs, args.rs, config.rs, run.rs, estimate.rs, verify.rs, inspect.rs}
│   ├── py/                     # classdiam-py (PyO3/maturin) — Phase 5, dir reserved
│   └── gpu/                    # classdiam-gpu (cudarc/wgpu TransformBackend) — future
└── results/                    # default output root (gitignored)
```

### Dependencies (versions current as of 2026-07; pin exactly at implementation time)

**`classdiam-core`**

| crate | ver | why |
|---|---|---|
| `num-bigint`, `num-integer`, `num-traits` | 0.4 / 0.1 / 0.2 | default `ExactInt` (MIT/Apache-2.0, pure Rust, Windows-clean) |
| `rayon` | 1.10 | block/target parallelism, MN subtree parallelism |
| `serde` (derive) | 1 | all report/checkpoint types |
| `serde_json` | 1 | result JSON emission (also needed by PyO3 later) |
| `postcard` | 1.1 | checkpoint body encoding |
| `crc32fast` | 1.4 | checkpoint integrity |
| `blake3` | 1.5 | partition-order hash, config hash |
| `fixedbitset` | 0.5 | visited/support bitsets |
| `smallvec` | 1.13 | `Partition` parts inline storage |
| `thiserror` | 2 | error taxonomy |

Feature flags: `default = []`; `malachite` (swaps `ExactInt` backend — **LGPL-3.0-only**, off by default and documented as such; `ibig` was rejected as default due to low maintenance since ~2022, `rug` rejected per no-GMP-on-Windows requirement); `mmap-cache` (adds `memmap2` 0.9 for on-disk resident-table cache); `slow-tests`.

**`classdiam-cli`**: `clap` 4.5 (derive), `anyhow` 1, `toml` 0.8, `tracing` 0.1 + `tracing-subscriber` 0.3, `serde_json` 1, `humantime` 2; optional `indicatif` 0.17 (`progress` feature).

**dev-dependencies (core)**: `proptest` 1, `criterion` 0.5.

The math core does **no** file I/O; all persistence (checkpoints, results, table cache) goes through injected sinks so the same core runs under CLI and PyO3.

---

## 2. Canonical partition order and core data model

### 2.1 The order: `lex_desc_full_parts_v1`

Partitions of `n` are written as weakly decreasing part lists and ordered **lexicographically descending**: index 0 is `[n]`, the last index is `[1,…,1]` (the identity type). Justification:

- generated directly by the standard descending-parts recursion — no sort pass, deterministic, backend-free;
- stable under implementation changes (unlike "backend-native order", which spec §3 explicitly warns about);
- SymPy fixtures never rely on order: every fixture row carries its explicit partition, and the loader maps through `PartitionIndex`.

The order is versioned by name (`"lex_desc_full_parts_v1"`) and by hash: `order_hash = blake3( LE(n as u16) ‖ for each partition in order: len(parts) as u8 ‖ parts as u8… )` (full parts, including 1s). Every result file, checkpoint, and cache embeds both name and hash; mismatch is a hard error on load.

### 2.2 Types

```rust
/// Weakly decreasing parts, each >= 1, sum == n. n <= 255 by construction (target n <= ~50).
#[derive(Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Partition { parts: SmallVec<[u8; 16]> }

impl Partition {
    pub fn n(&self) -> u16;
    pub fn parts(&self) -> &[u8];
    pub fn num_parts(&self) -> usize;
    pub fn sign(&self) -> i8;                    // (-1)^(n - ell)
    pub fn transpose(&self) -> Partition;
    pub fn z_value(&self) -> BigUint;            // prod i^{m_i} m_i!
    pub fn reduced(&self) -> Vec<u8>;            // parts > 1 only (JSON encoding)
}

#[derive(Copy, Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
pub struct PartitionId(pub u32);                 // index into canonical order; p(50)=204226 fits easily

pub struct PartitionIndex {
    n: u16,
    partitions: Vec<Partition>,                  // canonical order
    lookup: HashMap<Partition, PartitionId>,
    pub identity: PartitionId,                   // = last index, (1^n)
    sign: Vec<i8>,
    class_size: Vec<BigUint>,                    // n!/z — exact
    class_size_mod: Vec<Vec<u32>>,               // [prime][id], filled per resident prime
    transpose: Vec<PartitionId>,
    pub order_hash: [u8; 32],
}

impl PartitionIndex {
    pub fn build(n: u16) -> Self;                // pure, deterministic
    pub fn count(&self) -> usize;                // q = p(n)
    pub fn id_of(&self, p: &Partition) -> Option<PartitionId>;
    pub fn factorial_n(&self) -> &BigUint;
}
```

### 2.3 Cycle-type templates (input convention)

```rust
/// User-facing cycle type WITHOUT fixed points, e.g. [3] or [2,2].
/// Parts must be >= 2 (a literal 1 is rejected: fixed points are implicit).
/// The empty template [] denotes the identity class.
#[derive(Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct CycleTypeTemplate(Vec<u8>);

impl CycleTypeTemplate {
    pub fn parse(s: &str) -> Result<Self, TemplateError>;   // "2,2" → [2,2]
    pub fn min_n(&self) -> u16;                             // sum of parts
    /// Pad with fixed points to a partition of n. Errors: DoesNotFit (sum > n).
    pub fn pad_to(&self, n: u16) -> Result<Partition, TemplateError>;
    pub fn slug(&self) -> String;                           // "2.2" (file names, labels)
}

#[derive(Clone, Serialize, Deserialize)]
pub struct UnionTemplate {
    pub classes: Vec<CycleTypeTemplate>,   // deduplicated; duplicates warn + collapse
    pub label: Option<String>,
    pub allow_identity_generator: bool,    // default false; [] rejected otherwise
}
```

Per-`n` resolution (`UnionTemplate::resolve(&PartitionIndex) -> Result<UnionSpec, SkipReason>`) applies the spec §17 guards:

- **DoesNotFit**: any class with `min_n() > n` ⇒ the whole `(n, union)` job is *skipped* with a recorded reason (not an error — expected in n-ranges).
- **Identity class** (`[]`, or a template that pads to `1^n`): rejected unless `allow_identity_generator`. When allowed, a result-level note records that exact-length supports become cumulative (distances unaffected), per spec §5.3.
- **Split-class guard**: inputs are by construction *full S_n classes*, so the A_n split issue (spec §17.2, Failure 6) cannot arise; nevertheless, when all classes are even the result records `generated_subgroup ∈ {A_n-or-smaller}` and `cayley_graph_on_Sn: "disconnected"` so downstream users aren't misled.
- **Empty union** (all classes skipped): job skipped.
- Mixed/odd/even parity is computed and stored (`UnionParity::{Even, Odd, Mixed}`) — it drives parity filtering and certification (§7).

---

## 3. Arithmetic layer

### 3.1 Modular

```rust
/// A prime with n < p < 2^31, drawn from a fixed descending table:
/// 2147483647, 2147483629, 2147483587, 2147483579, ... (deterministic; skip any <= n — never
/// triggered for n <= 255). Determinism of the prime sequence is part of the checkpoint contract.
#[derive(Copy, Clone, PartialEq, Eq, Debug, Serialize, Deserialize)]
pub struct Prime31(pub u32);

pub struct ModCtx { p: u64, barrett: u128 }
impl ModCtx {
    pub fn new(p: Prime31) -> Self;
    #[inline] pub fn mul(&self, a: u32, b: u32) -> u32;      // (a as u64 * b as u64) % p, Barrett
    #[inline] pub fn reduce_u128(&self, acc: u128) -> u32;
    pub fn pow(&self, a: u32, e: u64) -> u32;
    pub fn inv(&self, a: u32) -> u32;                        // Fermat
    pub fn reduce_big(&self, x: &BigInt) -> u32;             // for theta, class sizes, bounds
}

/// Column-major residue matrix (rows = active character rows R, cols = active unions m).
pub struct ResidueMat { rows: usize, cols: usize, data: Vec<u32> }
```

**Proven overflow bounds (documented in `arith::modp` and enforced by tests):** residues `< 2^31`, so each product `< 2^62`; dot products accumulate in `u128` — with `R ≤ p(50) = 204226 < 2^18` terms the accumulator stays `< 2^80 ≪ 2^128`. A future GPU backend uses the alternate documented regime: primes `< 2^30`, `u64` accumulators, mandatory reduction every ≤ 16 terms (`16·2^60 = 2^64`). Both regimes are captured as constants + `debug_assert!`s, per spec §13.1/Failure 8.

### 3.2 Exact

```rust
// arith::exact — the ONLY place that names a bigint crate.
#[cfg(not(feature = "malachite"))] pub type ExactInt = num_bigint::BigInt;
#[cfg(feature = "malachite")]      pub type ExactInt = /* malachite Integer newtype */;
```

Exact integers appear only off the hot path: degrees (hook-length formula), class sizes, `ω`/`θ` construction with the divisibility assertion of spec §4, coefficient bounds, the optional bigint certifier, and the validation suite.

---

## 4. Character layer

### 4.1 Murnaghan–Nakayama evaluator (`chars::mn`)

Column-oriented trie DP, adopted from the candidate idea: for a target class `ν = (l₁ ≥ l₂ ≥ …)` the whole column `{χ^ρ(ν)}_ρ` is built bottom-up by composing sparse "remove one rim hook of length `l`" level operators; columns are organized in a trie over sorted cycle-length multisets so shared suffixes share the DP prefix; rayon parallelizes over trie subtrees. Two additions:

- **Multi-prime emission in one traversal.** The DFS carries one residue lane per resident prime (`[u32; K]` per node). The traversal (rim-hook enumeration) dominates cost, so K primes cost ~K× only in the cheap arithmetic, not in the structure walk.
- **Exact lane on demand.** The same walker runs with `ExactInt` values for single columns (base-class `ω` construction, bigint certification, fixtures cross-check).

```rust
pub struct MnEvaluator { /* per-n trie caches, arena-allocated */ }
impl MnEvaluator {
    pub fn new(index: &PartitionIndex) -> Self;
    /// Full column over `rows` for target nu, one output lane per ModCtx.
    pub fn column_mod(&self, nu: &Partition, rows: &[PartitionId],
                      ctxs: &[ModCtx], out: &mut [Vec<u32>]);
    pub fn column_exact(&self, nu: &Partition, rows: &[PartitionId]) -> Vec<ExactInt>;
    pub fn value_exact(&self, rho: &Partition, nu: &Partition) -> ExactInt; // tests/fixtures
}
```

### 4.2 `CharacterSource` — the table abstraction

Push-based: the source drives block production (so a streamed MN source can parallelize generation internally), the engine consumes through a `Sync` sink.

```rust
/// One block of character residues for a subset of target classes, all resident primes.
pub struct TableBlock<'a> {
    pub block_id: u32,
    /// Canonical PartitionIds of the targets in this block, sorted ascending.
    /// Blocks partition the full target set exactly (each target in exactly one block per sweep).
    pub targets: &'a [PartitionId],
    pub rows: u32,                        // == R, shared active-row count
    pub planes: &'a [PrimePlane<'a>],     // one per resident prime, same layout each
}
pub struct PrimePlane<'a> {
    pub p: Prime31,
    /// Column-major: data[t * rows + r] = chi^{rows[r]}(targets[t]) mod p, fully reduced.
    pub data: &'a [u32],
}

pub enum Interrupt { Deadline, Cancelled }

pub trait BlockSink: Sync {
    /// May be called concurrently from multiple threads; each block delivered exactly once
    /// per sweep, order unspecified. Break stops the sweep promptly (suspension path).
    fn consume(&self, block: &TableBlock<'_>) -> ControlFlow<Interrupt>;
}

pub trait CharacterSource: Send + Sync {
    fn n(&self) -> u16;
    fn rows(&self) -> &[PartitionId];             // active rows, fixed for the source's lifetime
    fn primes(&self) -> &[Prime31];
    fn num_blocks(&self) -> u32;
    /// One full sweep over all target blocks. Cost model: resident = memcpy-free slice
    /// handout; streamed = full MN regeneration.
    fn sweep(&self, sink: &dyn BlockSink) -> Result<SweepOutcome, CharError>;
    /// Stable across sweeps: true if per-sweep cost is ~zero (resident) — the planner uses this.
    fn is_resident(&self) -> bool;
}
```

Implementations:

- `ResidentTable` — generation sweep runs once (via `MnEvaluator`), stores `K × R × q` u32s; subsequent sweeps hand out slices. Optional `mmap-cache` feature persists it to the run dir keyed by `(n, primes, rows_hash, order_hash)`.
- `StreamedMnSource` — regenerates per sweep; blocks follow trie-subtree structure (hence `targets` is an id list, not a contiguous range).
- `EigenGroupedSource` — see §7.6; rows are *eigenvalue groups*, not irreps, but the type shape is identical (that is why `rows()` returns ids and the engine treats row semantics opaquely: the weight vector supplied per radius is what differs).

---

## 5. Spectra (`spectra.rs`)

```rust
pub struct BaseSpectra {
    pub base_classes: Vec<PartitionId>,          // t columns
    pub omega: Vec<Vec<ExactInt>>,               // [base][rho over ALL q rows] — exact, asserts
                                                 //   class_size * chi ≡ 0 (mod degree) per spec §4
    pub active_rows: Vec<PartitionId>,           // rows with omega != 0 for SOME base class (spec §7)
}
pub struct UnionSpectrum {
    pub theta_exact: Vec<ExactInt>,              // over active_rows; zero entries allowed
    pub theta_mod: Vec<Vec<u32>>,                // [prime][row]
    pub union_size: BigUint,                     // |U| = sum of class sizes (bound engine input)
}
impl BaseSpectra {
    pub fn build(idx: &PartitionIndex, mn: &MnEvaluator, degrees: &[BigUint],
                 bases: &[PartitionId]) -> Result<Self, ClassdiamError>;
    pub fn union(&self, incidence: &[bool]) -> UnionSpectrum;   // theta = Omega·b, exact
    /// Exact eigenvalue grouping for one union (used by the planner and EigenGroupedSource).
    /// NOTE: grouping is done on EXACT theta values — never on residues (collision-unsafe).
    pub fn eigen_groups(&self, u: &UnionSpectrum) -> Vec<(ExactInt, Vec<PartitionId>)>;
}
```

Building `Ω` needs only `t` exact MN columns (one per base class) — cheap even at `n = 50`. This is also where the planner counts distinct eigenvalues *exactly* before committing to eigen-grouped mode.

**Zero-eigenvalue rows and the radius-0 exception (spec §7, §9.1).** The shared row set `active_rows` drops rows that vanish for *all* base classes. Rows where a *particular* union's `θ_ρ(U) = 0` are not dropped structurally: their `power` entries become 0 after the first multiplication and contribute nothing for `r ≥ 1`. The engine never runs a transform at `r = 0`: `distance[identity] = 0` is preset, layer 0 is recorded synthetically, and iteration starts at `r = 1`. The full-row radius-0 identity test (`K⁰ = K_{(1^n)}`) lives exclusively in `validate.rs` against a full-row source.

---

## 6. `TransformBackend` — the GPU seam (documented contract)

```rust
pub struct TransformSpec<'a> {
    pub n: u16,
    pub rows: u32,                       // R
    pub primes: &'a [Prime31],
    pub max_unions: u32,                 // m upper bound (for buffer sizing)
    pub num_blocks: u32,
    pub max_block_targets: u32,
}

pub trait TransformBackend: Send + Sync {
    fn name(&self) -> &'static str;

    /// Exactly once per (n, prime-set, row-set), before the first radius. Exclusive access.
    /// The backend may allocate device buffers here. It must NOT assume table data yet.
    fn init(&mut self, spec: &TransformSpec<'_>) -> Result<(), BackendError>;

    /// Once per radius, before any apply_block. Exclusive access. `w` is the weight matrix
    /// for this radius: w[prime] is (R x m) col-major, fully reduced. For per-irrep mode
    /// W[r,j] = degree[r] * theta[r,j]^radius (mod p); for eigen-grouped mode W[r,j] = alpha_r^radius.
    /// A device backend uploads `w` here once.
    fn begin_radius(&mut self, radius: u32, w: &[ResidueMatView<'_>]) -> Result<(), BackendError>;

    /// out[prime][t, j] = sum_r block.planes[prime].data[t*R + r] * w[prime][r, j]  (mod p).
    /// &self: called CONCURRENTLY from rayon workers, one call per (block, radius).
    fn apply_block(&self, block: &TableBlock<'_>,
                   out: &mut [ResidueMatViewMut<'_>]) -> Result<(), BackendError>;

    /// Once per radius after all blocks. Exclusive. Device backends flush/synchronize here.
    fn end_radius(&mut self) -> Result<(), BackendError> { Ok(()) }
}
```

**Contract (normative, to be reproduced as rustdoc):**

- **Memory ownership.** All inputs are borrowed for the duration of the call only; a backend must not retain references (the lifetimes make stashing impossible without copying — copying/uploading is explicitly allowed). A backend MAY cache device-side copies of table blocks keyed by `(prime, block_id)`; the engine guarantees `block_id ↔ (targets, data)` is immutable within one `(n, prime-set, row-set)` session, so caches uploaded during radius 1 are valid for all later radii. Resident sources thereby give GPU backends a "upload once, GEMM forever" pattern; streamed sources imply per-sweep re-upload (the planner accounts for this).
- **Block iteration.** The engine sweeps all blocks exactly once per radius, unordered, possibly concurrently. Output buffers are per-block disjoint slices of the radius output; no synchronization between blocks is needed or allowed to matter.
- **Determinism.** Output residues are mathematically determined (`ℤ/p` is exact); every backend must produce bit-identical, fully reduced values in `[0, p)`. There is no floating-point tolerance anywhere. This is what makes checkpoints, resumes, and CPU/GPU mixing safe. A backend using float tricks (e.g., FP64 GEMM) must *prove* exactness under the documented bound regime and is still required to match the reference backend bit-for-bit in tests.
- **Where reduction happens.** Inside the backend, entirely. Accumulator strategy is the backend's choice (CPU default: `u128` accumulate, single reduce; GPU: sub-2^30 primes with block-wise `u64` reduction every ≤ 16 MACs) but must be stated in `name()`-discoverable docs and covered by the overflow test suite (spec §13, Failure 8).
- **Errors.** Any `BackendError` fails the radius; the engine converts it into job suspension (checkpoint remains valid — checkpoints only ever contain *committed* layers).

`transform::cpu` ships two impls: `CpuReference` (naive, obviously-correct, used by tests as the oracle) and `CpuBlocked` (tiled, cache-aware; per-call single-threaded since the engine parallelizes over blocks). The future `classdiam-gpu` crate implements the same trait; nothing in `engine/` changes.

---

## 7. The diameter engine

### 7.1 Per-union state

```rust
pub enum UnionPhase {
    Running,                          // advancing radii
    Suspended { at_radius: u32 },     // checkpointed mid-run (deadline)
    Done(Completion),
}
pub enum Completion {
    EmptyLayer { stop_radius: u32 },      // spec §5.2 rule, exact supports only
    AllTypesVisited { stop_radius: u32 }, // early exit: visited == parity-feasible upper set
}

pub struct LayerRecord {
    pub r: u32,
    pub new: Vec<PartitionId>,            // first-ever hits at radius r
    pub support_size: u32,                // |{nu : a_r(nu) > 0}| (exact)
    pub cert: LayerCertStats,
}

pub struct UnionState {
    pub spec: UnionSpec,
    pub spectrum: UnionSpectrum,
    /// First radius of EVEN length with a_r(nu) > 0, and same for ODD; -1 = none yet.
    /// distance(nu) = min of the two; these two arrays are the exact-support oracle (§7.3).
    pub first_hit: [Vec<i32>; 2],
    pub visited: FixedBitSet,             // union of both parities
    pub power_mod: Vec<Vec<u32>>,         // [prime][row] = theta^r — NOT checkpointed (recomputable)
    pub layers: Vec<LayerRecord>,
    pub phase: UnionPhase,
}
```

### 7.2 Radius loop (batched over active unions, lockstep)

For the active batch (columns compacted so `W` has exactly `m_active` columns; a `Vec<UnionId>` maps columns back):

```
r += 1
1. power_mod[p][:, j] *= theta_mod[p][:, j]          (all active j, all primes)
2. Assemble W[p] (R x m_active): per-irrep mode  W = degree_mod[p] ⊙ power_mod[p]
                                 eigen-grouped   W = power_mod[p]         (degrees folded into H)
3. backend.begin_radius(r, W)
4. source.sweep(sink)  — sink, per block and rayon-parallel:
     backend.apply_block(block, out)
     for each union j, target t in block:
        residues = out[*][t, j]                      (numerator N_r mod each prime; n! never divided out)
        if any residue != 0        → POSITIVE
        else if parity-infeasible  → ZERO (proven, free)
        else if first_hit[r%2][t] set and first_hit[r%2][t] <= r-2 ... (cannot happen: r+2 lemma
             guarantees a nonzero residue is *mathematically* present; all-zero residues here would
             mean a hidden positive — recorded as candidate)                → CANDIDATE
        else                       → CANDIDATE
     accumulate per-union positive bitset + candidate list
     word-count tripwire (always on, O(q)): Σ_ν class_size_mod[p][ν]·N_r(ν) ≡ n!·|U|^r (mod p)
5. backend.end_radius()
6. certify_candidates(r, j, candidates)              → each candidate becomes POSITIVE or ZERO (§7.4)
7. Layer commit per union j:
     support_r = positives (exact now)
     for nu in support_r with first_hit[r%2][nu] < 0: first_hit[r%2][nu] = r
     new = support_r \ visited;  distance implicit (min of first_hits);  visited |= new
     push LayerRecord
     if visited ⊇ parity_feasible_upper_set(j)  → Done(AllTypesVisited)   // rigorous: parity bound only
     else if new.is_empty()                     → Done(EmptyLayer)         // spec §5.2, valid: supports exact
8. Compact finished unions out of the batch; checkpoint (§10); deadline check (§8.3).
```

Notes:

- **`r+2` lemma (free positivity).** Classes are inverse-closed, so `a_r(ν) > 0 ⇒ a_{r+2}(ν) > 0` (append `g·g⁻¹`). Consequently a type whose same-parity first hit is already recorded is *proven* in `support_r` without looking at residues — all-zero residues for such a type would expose a backend bug (asserted). This removes already-visited types from certification entirely and makes reported supports exact for free.
- **Parity filtering (spec §11.3).** Single-parity unions skip infeasible target parity outright (~½ the classification work and the certification set). Mixed-parity unions get no filter — this is precisely why `first_hit` is a *pair* of arrays: exact per-radius supports for mixed unions are `{ν : first_hit[r mod 2][ν] ∈ [0, r]}`, which the two arrays capture without storing per-radius bitsets.
- **Stopping.** Spec §5.2's single-empty-layer rule is used with *cumulative* visited (valid for mixed parity too: supports of products are monotone in support sets). The `AllTypesVisited` early exit fires first when the union generates `S_n` (or covers all even types for even-only unions) — it is justified purely by the parity upper bound, never by subgroup theory (spec §2.2's "shortcuts must not affect correctness").

### 7.3 Exact supports in the output

Layers report `new` explicitly; full exact-length supports are reported either as explicit index arrays (default for `n ≤ 30`) or via the documented reconstruction rule from `first_hit_even/odd` (always emitted). Supports are exact for every reported radius `0..=stop_radius` because every layer's candidates were certified before commit.

### 7.4 The certification gate (`engine::certify`) — divergence from the plan, with reasons

The candidate ideas prescribe "bigint fallback for candidate zeros". Adopted but **demoted to tier 4**, because rigorous certification is available much cheaper:

Let `B_r(ν) = ⌊|U|^r / |C_ν|⌋`. From the word-count identity `Σ_ν |C_ν| a_r(ν) = |U|^r` and nonnegativity, `a_r(ν) ≤ B_r(ν)` — a *much* sharper rigorous bound than `|U|^r`. Tiers, applied per candidate `(ν, r, j)`:

1. **Bound-zero:** `B_r(ν) = 0` (i.e. `|U|^r < |C_ν|`) ⇒ `a_r(ν) = 0`. Free; certifies most of the early-radius frontier where big classes are still unreached.
2. **Resident CRT:** if `∏ resident primes > B_r(ν)`, all-zero residues already *uniquely determine* `a_r(ν) = 0` (spec §12.3). With 3 × 31-bit primes this covers `B_r < 2^93` — for typical class sizes this carries certification deep into the run.
3. **Dynamic extra primes:** take fresh primes from the fixed table; for each, evaluate `N_r(ν) mod p'` directly: one multi-prime MN *column* for `ν` (`MnEvaluator::column_mod`), `θ mod p'` from exact `Ω` (`ModCtx::reduce_big`), `θ^r` by fast exponentiation (`O(R log r)`), one length-`R` dot product. Stop at the first nonzero residue (⇒ POSITIVE, a genuinely hidden positive — counted and reported) or when the accumulated prime product exceeds `B_r(ν)` (⇒ certified ZERO). **This tier always terminates with a rigorous verdict**; no probabilistic step exists anywhere.
4. **Bigint exact** (`certifier = "bigint"` or `"both"`): exact `Σ f_ρ χ^ρ(ν) θ_ρ^r` with `ExactInt`. Kept as a cross-check oracle (tests, `--certifier both` paranoia mode), not as the default — tier 3 strictly dominates it in cost while being equally exact.

All four tiers are counted per layer and reported (`certification` block in the JSON), including `hidden_positives_found` — expected to be 0 essentially always, but its presence in output is the audit trail that Failure 4/9 cannot occur silently.

### 7.5 Union batching and compaction

All unions for a given `n` run in one lockstep batch sharing the row set, the source, and each radius's sweep (the dominant cost — table access — is paid once per radius for the whole batch, per spec §8.2). When a union completes, its column is swap-removed from `Θ`/`Power`/`W` at the next radius boundary (compaction), and its result is finalized+written immediately (crash safety). Batches larger than `max_unions` (planner-derived from memory) are chunked.

### 7.6 Eigen-grouped mode (stretch-goal path, spec §15) — adopted with guardrails

For `n ≥ ~48` (or whenever the planner's resident-table estimate exceeds the memory budget), the per-irrep table is never materialized. Instead, in **one streamed MN sweep** the engine accumulates, per union `j` and per resident prime, `H_α(ν) = Σ_{ρ: θ_ρ(U_j)=α} f_ρ χ^ρ(ν) (mod p)` — grouping by **exact** `θ` values (residue-equality grouping is forbidden: collisions mod p would merge distinct eigenvalues). Afterwards each radius costs `K × E_j × q` MACs with `E_j = #distinct eigenvalues`.

Answering the "verify" question honestly: `E_j ≪ R` is **not guaranteed** in general. What is known:

- For **even-only** unions, `θ_ρ' = θ_ρ` (transpose pairing with `sgn(λ)=+1`), so `E ≤ ⌈R/2⌉ + #self-transpose` — a guaranteed ~2× reduction. For odd/mixed unions pairs give `±θ` / no relation (spec §11.1), no guarantee.
- Empirically eigenvalue multiplicities grow with `n` for small-support classes, but the design does not bet on it: **`E_j` is computed exactly and cheaply at plan time** (from exact `Ω`, which needs only `t` MN columns), so mode selection is data-driven per `(n, union)`, not hoped-for.
- **Fallback if `E ≈ R`:** streamed per-radius regeneration (cost ≈ table-gen × diameter — the honest price), and for `n ≥ 50` an optional `--distances-only` mode that shrinks the target set each radius to unvisited types (forfeiting mixed-parity support reporting, which is then emitted only via the single-parity derivation or marked unavailable). The architecture (same `CharacterSource`/`TransformBackend` shapes) supports all three without engine changes — the requirement that `n ≥ 50` not be *precluded* is met.

Cost caveat recorded in the planner: `H` construction is per-union work inside a shared sweep (all unions' `H` accumulated in the same pass), so batching still amortizes generation.

---

## 8. Orchestration across an n-range

### 8.1 Resource model (`orchestrate::estimate`)

Memory (exact formulas, u32 residues, `R ≤ q`, `K` resident primes):

| n | q = p(n) | table/prime (q²·4) | 3 primes | verdict on 300 GB |
|---|---|---|---|---|
| 30 | 5 604 | 126 MB | 0.38 GB | resident, trivial |
| 35 | 14 883 | 886 MB | 2.7 GB | resident |
| 40 | 37 338 | 5.6 GB | 16.7 GB | resident |
| 45 | 89 134 | 31.8 GB | 95.3 GB | resident (fits, checked at runtime) |
| 50 | 204 226 | 166.9 GB | 500 GB | **not resident** → eigen-grouped / streamed |

Plus per-union `O(q)` arrays and `K·R·m` power/weight matrices — negligible. Time model: per-radius transform `≈ K·R·q·m` modular MACs (e.g. `n=45`, 3 primes, 1 union: `2.4·10¹⁰` MACs ⇒ well under a second on 96 cores) — **table generation dominates**; MN cost has no clean closed form, so the estimator is calibration-based: micro-benchmark MN generation at `n ∈ {18, 22, 26}` on startup (~seconds), fit a power law, refine online from actually-measured sweeps, persist calibration in the run dir. Estimates carry stated uncertainty (×0.5–×2 band) and the planner treats them as advisory except for hard memory caps.

```rust
pub struct ResourcePlan {
    pub n: u16,
    pub mode: EngineMode,                 // ResidentTable | Streamed | EigenGrouped { e_per_union: Vec<u32> }
    pub resident_primes: u8,
    pub est_table_bytes: u64,
    pub est_prepare_s: RangeF64,
    pub est_per_radius_s: RangeF64,
    pub est_radii: u32,                   // heuristic diameter guess (n / smallest support), advisory
    pub verdict: PlanVerdict,             // Ok | NeedsMode(EngineMode) | Skip(SkipReason)
}
pub fn plan(n: u16, unions: &[UnionTemplate], budget: &Budget, cal: &Calibration) -> ResourcePlan;
```

### 8.2 Range and auto-n drivers

`n` is processed ascending; each `n` is an independent job (nothing math-reusable across `n` per spec §6 — only the calibration data carries over).

- **Explicit range** (`--n 8..=40`): jobs whose plan says `Skip(Memory)` are skipped with reason; since memory grows monotonically in `n`, the first memory skip aborts the remaining range (recorded as `not_attempted`).
- **Auto-n** (`--auto-n --min-n 10`): before each `n`, check `est_prepare + est_radii·est_per_radius` (upper band) against remaining `--max-seconds` budget and `est_table_bytes` against `--max-memory`; run while feasible, then stop gracefully with manifest status `budget_exhausted`, leaving instructions (the manifest records the next `n` and its estimate, so a follow-up Kaggle run resumes the ladder).
- Per-`(n, union)` template resolution failures (does-not-fit, identity) are per-job skips, never run-fatal.

### 8.3 Deadline and cancellation

```rust
pub struct CancelToken(Arc<AtomicBool>);           // Python/KeyboardInterrupt & Ctrl-C
pub struct DeadlineGuard { deadline: Instant, safety: Duration }
```

Checked at every block boundary (via `BlockSink` returning `ControlFlow::Break`) and radius boundary. Policy: if remaining time `< 1.5 ×` the measured last-radius time (or the sweep is interrupted mid-radius), abandon the *uncommitted* radius, checkpoint the last committed state, mark jobs `Suspended`, write manifest, exit with code 75. A killed-without-warning process loses at most one radius of work thanks to per-layer checkpoints.

---

## 9. Output: files, schemas, example

### 9.1 Decision: file layout

**One JSON per `(n, union)` + one run manifest + one append-only JSONL event log.** Rationale: per-job files are the natural unit for Kaggle artifact collection and for resuming partial batch runs (finished jobs are immutable); a single combined JSON would be rewritten constantly and lost on kill. The JSONL log is *not* the source of truth — checkpoints are — it exists for live monitoring and post-mortems; recommending JSONL *as* the primary format was considered and rejected (self-describing single-document results are easier for downstream analysis, and crash-safety is already covered by checkpoints).

```
results/<run_id>/                      # run_id = 20260712T081500Z-3f9a (UTC + 4-hex of config hash)
├── manifest.json                      # rewritten atomically (tmp+rename) after every job event
├── events.jsonl                       # append-only: {ts, event, n, union, r, ...}
├── config.resolved.toml               # exact input echo
├── n06_g2.json                        # result files: n{:02}_g{slug}.json
├── n06_g3+2.2.json                    #   union slug: parts joined ".", classes joined "+"
├── checkpoints/n12_g2.ckpt            # removed on job completion
└── cache/                             # optional mmap table cache (feature mmap-cache)
```

Windows/Unix-safe slugs (digits, `.`, `+` only). Collisions impossible (slug is a canonical encoding).

### 9.2 Result schema (`classdiam/result`, version 1) — field by field

| field | type | meaning |
|---|---|---|
| `format`, `format_version` | str, int | `"classdiam/result"`, `1`; parsers must reject unknown majors |
| `spec_version` | str | notes file identity this run implements |
| `tool` | obj | `name`, `version`, `core_commit` |
| `run` | obj | `run_id`, `started_utc`, `finished_utc`, `resumed_from_checkpoint` (bool), `suspend_resume_count` |
| `n`, `factorial_n` | int, str | `n!` as decimal string |
| `generators` | obj | `input_templates` (as given), per-class `{template, padded, index, class_size(str), sign}`, `union_size` (str), `parity` (`"even"|"odd"|"mixed"`), `allow_identity_generator`, `label` |
| `partition_order` | obj | `convention` (`"lex_desc_full_parts_v1"`), `count` (=q), `hash_blake3`, `partitions_reduced` — full ordered list, each as parts>1 only (unambiguous given `n`; identity = `[]`) |
| `class_data` | obj | parallel arrays over the order: `sign` (int8), `class_size` (dec strings) — makes files self-contained |
| `arithmetic` | obj | `resident_primes` (ints), `screening` (`"numerator-mod-p"` — documents that n! is never divided in screening), `certifier` (`"crt"|"bigint"|"both"`), overflow regime id |
| `engine` | obj | `mode`, `backend`, `active_row_count`, `zero_rows_all_bases` (ids), `threads` |
| `results.distance` | int[q] | distance per canonical index; **sentinel `-1` = unreachable** (also stated in `unreachable_value`) |
| `results.first_hit_even/odd` | int[q] | min even/odd radius with `a_r>0`; `-1` = never (up to stop radius); `distance = min` of the two; enables exact support reconstruction: `support_r = {ν : 0 ≤ fh[r mod 2][ν] ≤ r}` |
| `results.diameter_identity_component` | int | max finite distance |
| `results.reachable_count`, `generated_subgroup`, `cayley_graph_on_Sn`, `bipartite` | — | reachability/parity metadata (`generated_subgroup` derived from the *computed* visited set: `"S_n"|"A_n"|"proper_subgroup"|"trivial"`) |
| `results.stopping` | obj | `rule` (`"all_types_visited"|"empty_layer"`), `stop_radius` |
| `results.layers[]` | obj | `r`, `new` (indices, exact first hits), `support_size` (exact), `support` (indices; emitted when `emit_supports="indices"`, default for n≤30, else omitted — reconstruction rule documented in `results.support_reconstruction`) |
| `certification` | obj | totals: `candidates`, `bound_certified`, `crt_resident_certified`, `extra_prime_evals`, `extra_primes_max_used`, `bigint_evals`, `hidden_positives_found` |
| `timings_s` | obj | `prepare`, `table_generation`, `transform_total`, `certification_total`, `per_radius[]`, `total_wall` |
| `config_hash_blake3` | str | resolved-config hash (matches checkpoint header) |

### 9.3 Full example — `n = 6`, generators `[[2]]` (transpositions), file `n06_g2.json`

```json
{
  "format": "classdiam/result",
  "format_version": 1,
  "spec_version": "character_method_cayley_diameters.md@2026-07",
  "tool": { "name": "classdiam", "version": "0.1.0", "core_commit": "4a1f9c2e" },
  "run": {
    "run_id": "20260712T081500Z-3f9a",
    "started_utc": "2026-07-12T08:15:03Z",
    "finished_utc": "2026-07-12T08:15:03Z",
    "resumed_from_checkpoint": false,
    "suspend_resume_count": 0
  },
  "n": 6,
  "factorial_n": "720",
  "generators": {
    "input_templates": [[2]],
    "classes": [
      { "template": [2], "padded": [2, 1, 1, 1, 1], "index": 9, "class_size": "15", "sign": -1 }
    ],
    "union_size": "15",
    "parity": "odd",
    "allow_identity_generator": false,
    "label": "g2"
  },
  "partition_order": {
    "convention": "lex_desc_full_parts_v1",
    "count": 11,
    "hash_blake3": "7c9f1e2ab8d4460f3a5e9b0c1d2e3f4a5b6c7d8e9f0a1b2c3d4e5f6a7b8c9d0e",
    "partitions_reduced": [[6],[5],[4,2],[4],[3,3],[3,2],[3],[2,2,2],[2,2],[2],[]]
  },
  "class_data": {
    "sign":       [-1, 1, 1, -1, 1, -1, 1, -1, 1, -1, 1],
    "class_size": ["120","144","90","90","40","120","40","15","45","15","1"]
  },
  "arithmetic": {
    "resident_primes": [2147483647, 2147483629, 2147483587],
    "screening": "numerator-mod-p",
    "certifier": "crt",
    "overflow_regime": "p31-u128-accumulate-v1"
  },
  "engine": {
    "mode": "resident",
    "backend": "cpu-blocked-v1",
    "active_row_count": 10,
    "zero_rows_all_bases": [5],
    "threads": 8
  },
  "results": {
    "unreachable_value": -1,
    "distance":       [5, 4, 4, 3, 4, 3, 2, 3, 2, 1, 0],
    "first_hit_even": [-1, 4, 4, -1, 4, -1, 2, -1, 2, -1, 0],
    "first_hit_odd":  [5, -1, -1, 3, -1, 3, -1, 3, -1, 1, -1],
    "diameter_identity_component": 5,
    "reachable_count": 11,
    "generated_subgroup": "S_n",
    "cayley_graph_on_Sn": "connected",
    "bipartite": true,
    "stopping": { "rule": "all_types_visited", "stop_radius": 5 },
    "support_reconstruction": "support_r = { i : 0 <= first_hit[r mod 2][i] <= r }",
    "layers": [
      { "r": 0, "new": [10],      "support": [10],                 "support_size": 1 },
      { "r": 1, "new": [9],       "support": [9],                  "support_size": 1 },
      { "r": 2, "new": [6, 8],    "support": [6, 8, 10],           "support_size": 3 },
      { "r": 3, "new": [3, 5, 7], "support": [3, 5, 7, 9],         "support_size": 4 },
      { "r": 4, "new": [1, 2, 4], "support": [1, 2, 4, 6, 8, 10],  "support_size": 6 },
      { "r": 5, "new": [0],       "support": [0, 3, 5, 7, 9],      "support_size": 5 }
    ]
  },
  "certification": {
    "candidates": 8,
    "bound_certified": 3,
    "crt_resident_certified": 5,
    "extra_prime_evals": 0,
    "extra_primes_max_used": 0,
    "bigint_evals": 0,
    "hidden_positives_found": 0
  },
  "timings_s": {
    "prepare": 0.004, "table_generation": 0.001, "transform_total": 0.002,
    "certification_total": 0.000,
    "per_radius": [0.001, 0.000, 0.000, 0.001, 0.000],
    "total_wall": 0.009
  },
  "config_hash_blake3": "d41a0c77b2e94f01aa38c5d6e7f8091a2b3c4d5e6f708192a3b4c5d6e7f80912"
}
```

(Values above are the true mathematical answers for S₆/transpositions: `d(ν) = 6 − #cycles(ν)`, diameter 5, stop at `r = 5` via full cover; row `[3,2,1]` really has `θ = 0` and is dropped from active rows; the 8 candidates at radii 1–3 are all certified by tiers 1–2.)

### 9.4 Manifest schema (`classdiam/manifest`, version 1)

```json
{
  "format": "classdiam/manifest", "format_version": 1,
  "run_id": "20260712T081500Z-3f9a",
  "config_hash_blake3": "d41a…",
  "budget": { "deadline_seconds": 6900, "max_memory_gb": 250, "elapsed_seconds": 512.4 },
  "status": "completed",              // completed | budget_exhausted | suspended | failed
  "auto_n": { "enabled": false, "next_n": null, "next_n_estimate_s": null },
  "jobs": [
    { "n": 6, "union": "g2", "status": "done", "file": "n06_g2.json",
      "diameter": 5, "stop_radius": 5, "elapsed_s": 0.009 },
    { "n": 6, "union": "g7", "status": "skipped", "reason": "class_does_not_fit" },
    { "n": 12, "union": "g2", "status": "suspended", "checkpoint": "checkpoints/n12_g2.ckpt",
      "committed_radius": 7 }
  ]
}
```

`events.jsonl` lines: `{"ts":"…","event":"layer_committed","n":12,"union":"g2","r":7,"new":34,"support":210,"candidates":3}` etc. — append-only, fsync'd per radius.

---

## 10. Checkpoint / resume

**Key simplification vs the candidate idea:** modular power vectors are *not* serialized — they are recomputed on resume by fast exponentiation (`power[p][ρ] = θ[p][ρ]^r`, `O(K·R·log r)`, milliseconds), since the prime sequence and `θ` are deterministic functions of the config. Checkpoints therefore shrink to the irreproducible-only state: first-hit arrays, layer log, phase.

**Format:** fixed 80-byte header + `postcard` body + CRC.

```
magic "CDCK" | format_version u16 | flags u16 | config_hash blake3 [32] |
order_hash blake3 [32] | body_len u64          ...then body, then crc32(header‖body)
```

```rust
#[derive(Serialize, Deserialize)]
pub struct CheckpointBody {
    pub n: u16,
    pub union_slug: String,
    pub resolved_union: Vec<Vec<u8>>,          // padded partitions, belt-and-braces vs slug
    pub resident_primes: Vec<u32>,
    pub engine_mode: EngineMode,
    pub committed_radius: u32,                 // all layers <= this are final & certified
    pub first_hit_even: Vec<i32>,
    pub first_hit_odd: Vec<i32>,
    pub layers: Vec<LayerRecord>,
    pub cert_totals: CertStats,
    pub elapsed_before_s: f64,
    pub calibration: Calibration,              // time-model state carries across sessions
}
```

- **When written:** after every *committed* layer (post-certification), atomically (`.ckpt.tmp` + rename); previous checkpoint kept until the new rename succeeds (`keep = 2`). Never mid-radius — an interrupted sweep is simply discarded, guaranteeing the invariant "checkpoint state is exact".
- **Resume validation:** header magic/version; `config_hash` must equal the hash of the *resolved* config of the resuming run (n, union set, primes, mode, format version — CLI prints a diff on mismatch and refuses; `--force-resume` is deliberately not offered: a mismatched resume is mathematically meaningless); `order_hash` re-derived from `PartitionIndex::build(n)` and compared; CRC checked.
- **Resume procedure:** rebuild `PartitionIndex`, `Ω`, sources (deterministic; resident tables regenerated or mmap-cache-loaded), recompute `power` to `committed_radius`, continue at `committed_radius + 1`.
- **Kaggle interplay:** notebook calls `classdiam run --config … --deadline 6900` (2 h minus margin); on deadline the process suspends per §8.3 and exits 75; the notebook persists `results/<run_id>/` as a Kaggle dataset; the next notebook run calls `classdiam resume results/<run_id>` and continues. Because checkpoints are tiny (`2·q·4` bytes + layer log ≈ 1 MB at n = 45) and per-layer, at most one radius of transform work is ever repeated.

---

## 11. CLI UX (`classdiam`)

```
classdiam run       --config runs/main.toml [--resume] [--deadline 6900] [--out DIR]
classdiam run       -n 8..=20 -u 2 -u "3+2,2" -o results/quick --primes 3 --threads 0
classdiam run       --auto-n --min-n 10 -u 3 --deadline 6900 --max-memory 250GB
classdiam resume    results/20260712T081500Z-3f9a [--deadline 6900]
classdiam estimate  -n 45 -u 3 [--mode auto]        # prints ResourcePlan, no compute
classdiam verify    [--max-n 9] [--certifier both]  # invariant suite + brute-force BFS cross-check
classdiam fixtures  --check fixtures/               # MN vs committed SymPy fixtures
classdiam table     -n 12 --degrees --class-sizes   # inspection/debug dumps
classdiam inspect   results/…/n06_g2.json [--layers] [--diff OTHER.json]
```

Union grammar for `-u`: classes joined by `+`, parts within a class by `,` — `-u 2` (transpositions), `-u 2,2` (double transpositions), `-u "3+2,2"` (union). Repeat `-u` for multiple unions. `--n` accepts `A..=B`, comma lists, or is replaced by `--auto-n --min-n A`.

Config file (TOML; JSON accepted with the same schema):

```toml
[run]
n = "6..=30"                    # or [6,10,12]; ignored when auto_n = true
auto_n = false
min_n = 10                      # auto-n start
out_dir = "results/transpositions"
deadline_seconds = 6900
max_memory_gb = 250
threads = 0                     # 0 = all cores
emit_supports = "auto"          # "indices" | "none" | "auto" (indices for n<=30)

[[unions]]
classes = [[2]]

[[unions]]
classes = [[3], [2, 2]]
label = "3cyc+dbl"

[generators]
allow_identity_generator = false

[arithmetic]
resident_primes = 3
certifier = "crt"               # crt | bigint | both

[engine]
mode = "auto"                   # auto | resident | streamed | eigen-grouped
backend = "cpu-blocked"         # cpu-reference | cpu-blocked | (future: cuda, wgpu)

[checkpoint]
every_layers = 1
keep = 2
```

Exit codes: `0` completed; `75` suspended (deadline) — resumable; `2` config error; `1` internal error. Logging via `tracing` (`-v/-vv`), machine events to `events.jsonl` regardless.

---

## 12. Future PyO3 surface (signatures only; constraints on core API)

Core API constraints already honored so the wrapper is mechanical: public entry points take/return **owned, `Send + Sync`, lifetime-free, serde-serializable** types; progress is a `Sync` trait object; cancellation is an `Arc`-backed token; errors implement `std::error::Error`.

```rust
// classdiam-core — the exact surface classdiam-py wraps:
pub fn run(config: RunConfig, sink: Arc<dyn ProgressSink>, cancel: CancelToken)
    -> Result<RunOutcome, ClassdiamError>;                    // RunOutcome: Serialize
pub fn plan_only(config: RunConfig) -> Result<Vec<ResourcePlan>, ClassdiamError>;
pub struct Session { /* per-n: PartitionIndex + Omega + sources; Send + Sync */ }
impl Session {
    pub fn new(n: u16, base_classes: Vec<CycleTypeTemplate>, opts: SessionOpts) -> Result<Self, _>;
    pub fn run_union(&self, union: UnionTemplate, opts: JobOpts,
                     cancel: CancelToken) -> Result<UnionResult, _>;   // reuses tables across calls
}
```

```python
# python module `classdiam` (abi3, maturin), Kaggle usage:
def run(config: dict | str, *, progress: Callable[[dict], None] | None = None) -> dict: ...
def estimate(n: int, unions: list[list[list[int]]], **budget) -> list[dict]: ...
def resume(run_dir: str, *, deadline_s: float | None = None) -> dict: ...

class Session:
    def __init__(self, n: int, base_classes: list[list[int]],
                 *, primes: int = 3, threads: int | None = None) -> None: ...
    def run_union(self, classes: list[list[int]], *, label: str | None = None,
                  deadline_s: float | None = None,
                  checkpoint_dir: str | None = None) -> dict: ...   # parsed result JSON
    @property
    def partition_order(self) -> list[list[int]]: ...
```

GIL policy: every compute call is wrapped in `py.allow_threads(...)`; `progress` callbacks are invoked from worker threads via `Python::with_gil`, throttled to ≥ 250 ms; `KeyboardInterrupt` maps to `CancelToken.cancel()` → clean suspension with checkpoint, surfaced as a Python exception carrying the checkpoint path. Return values are `dict`s obtained by serializing `RunOutcome`/`UnionResult` through the same serde schema as the JSON files — one schema, three consumers (file, CLI, Python).

---

## 13. Divergences from the candidate ideas (summary)

| candidate idea | verdict | notes |
|---|---|---|
| Modular-first table, exact columns on demand | **Adopted** | exact lanes used only for `Ω`, degrees, certification, validation |
| MN column trie with suffix sharing | **Adopted**, extended | multi-prime lanes in one traversal; trie subtrees define the block structure (hence `TableBlock.targets` is an id list, not a range) |
| ~31-bit primes, u64 products, u128 accumulate | **Adopted** | codified as regime `p31-u128-accumulate-v1`; GPU regime (`p<2^30`, u64 acc, reduce/16) pre-specified |
| Per-radius modular GEMV/GEMM, unions as columns | **Adopted** | plus: `n!` never divided out in screening (support-equivalence mod p) |
| Bigint fallback certifies candidate zeros | **Adapted** | demoted to tier 4; primary certifier is rigorous dynamic CRT with the sharper bound `a_r(ν) ≤ |U|^r/|C_ν|` (from the word-count identity) — strictly cheaper, equally exact; bigint kept for `--certifier both` and tests. Also: certification must run **every radius** (not only at stopping), otherwise distances can be silently off by 2; the parity + `r+2`-lemma freebies keep this cheap |
| Eigenvalue grouping for `n ≥ 50` | **Adopted with guardrails** | grouping on *exact* θ only; `E` computed exactly at plan time (no wishful compression); guaranteed ≤ R/2 only for even-only unions; fallbacks = streamed sweeps, distances-only mode |
| Workspace core/cli/pyo3 | **Adopted** | + reserved `crates/gpu`; fixtures generated by `tools/gen_fixtures.py`, committed |
| Checkpoint with power vectors | **Simplified** | power vectors recomputed via modexp on resume; checkpoints shrink to first-hit arrays + layer log + hashes, written per committed layer |
| (new) per-parity first-hit arrays | **Added** | delivers the required *exact* exact-length support sets for mixed-parity unions at `O(q)` extra memory, and yields free positivity proofs that shrink the certification set |

---

## 14. Testing map (requirement 9 → suite)

- `validate.rs` + `tests/invariants.rs`: full-row radius-0 identity, radius-1 indicator (single class / union / weighted), `Σ|C_ν|a_r = |U|^r` (exact for `n ≤ 12`, and *always on* mod-p in the engine), `n! | N_r` divisibility, nonnegativity, parity, transpose symmetry `χ^{ρ'} = sgn·χ^ρ` and mixed-union even/odd θ decomposition (spec §11.1).
- `tests/bfs_crosscheck.rs`: enumerate S_n (`n ≤ 9`; `n = 10` under `slow-tests`), explicit BFS, compare distances, layers, supports, diameters against the engine for a fixed union menu including mixed-parity unions.
- `tests/fixtures.rs`: MN vs SymPy fixtures (values, degrees, class sizes), including the spec §22.2 hard cases (repeated cycle lengths, long rim hooks, self-transpose rows on odd classes, identity column = degrees).
- `tests/modular.rs`: every residue vs exact for `n ≤ 10`; **engineered prime-divisible coefficients** (search small cases where `p | a_r(ν)` for a table prime, plus synthetic dot products) proving tier-3 certification triggers and never mislabels; overflow-bound tests at both arithmetic regimes; `CpuReference` vs `CpuBlocked` bit-equality (the same harness later pins the GPU backend).
- `proptest`: partition round-trips, template padding, order-hash stability, checkpoint encode/decode/corrupt-header rejection, resume-equals-uninterrupted (kill after every layer, resume, compare final JSON byte-for-byte).

## 15. Phasing

1. **P1 (n ≤ 30, "simple and reliable"):** partitions, MN, exact `Ω`, resident tables, `CpuReference`, single-prime screening + bigint certifier, full test suite, result JSON + manifest.
2. **P2:** multi-prime residency, tiered CRT certifier, union batching/compaction, `CpuBlocked` + rayon, checkpoint/resume, deadline, CLI complete. Target n = 35–40.
3. **P3 (n = 45):** streamed source, mmap cache, planner/auto-n calibration, estimator hardening.
4. **P4 (stretch, n ≥ 50):** eigen-grouped mode, distances-only mode.
5. **P5:** `classdiam-py` (PyO3/maturin, abi3); **P6:** `classdiam-gpu` against the frozen `TransformBackend` contract.