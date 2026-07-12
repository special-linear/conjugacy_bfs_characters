# Test, Validation, Benchmarking, CI, and Packaging Design
## Conjugacy-Invariant Cayley Diameter Engine (Rust core)

Scope: this document covers the quality-engineering half of the framework design — test plan, brute-force oracles, SymPy fixture pipeline, adversarial modular tests, property tests, benchmarks, CI, and packaging. It is written against the spec at `d:\Math\self\cayleypy\conjugacy_bfs_characters\notes\character_method_cayley_diameters.md` (sections cited as §N) and the nine fixed project decisions. Crate names below (`ccd-core`, `ccd-cli`, `ccd-py`) are placeholders.

---

## 1. Testing philosophy

Three principles drive everything below:

1. **Every optimized path is differentially tested against a slower reference path that lives in the repo.** The dependency chain of trust is: brute-force BFS over raw permutations (n ≤ 10) → exact bigint spectral engine (n ≤ ~16) → modular engine with exact fallback (all n). Each layer is validated against the layer below on the overlap, plus independent closed-form anchors (hook lengths, z_λ, content sums, orthogonality, d(ν) = n − ℓ(ν) for transpositions) that work at *any* n.
2. **Orientation errors must fail loudly at the index level, not silently at the value level** (§23 Failure 7). Every fixture, checkpoint, and output file embeds the explicit partition list; every consumer asserts list equality before comparing a single character value.
3. **The stopping criterion is the single most dangerous piece of logic** (§5.2, §12, Failures 4/9). It gets its own adversarial test family in which modular screening is *forced* to produce false-zero layers via injected small primes.

### 1.1 Test tiers

| Tier | When | Budget | Contents |
|---|---|---|---|
| T0 | every commit, debug profile | < 2 min | unit tests, invariants n ≤ 8, kernel differentials, low-iteration proptest |
| T1 | every commit, release profile with `debug-assertions = true` | < 15 min | invariants n ≤ 12, brute-force BFS n ≤ 8, word-count DP n ≤ 7, modular-vs-exact differential n ≤ 10, full-table fixtures n ≤ 14, checkpoint/resume, JSON golden files |
| T2 | nightly cron + `workflow_dispatch`, `--ignored` | hours | brute force n = 9, 10; invariants n ≤ 16; fixtures n ≤ 16 full + spot to n ≈ 28; overflow stress at p(50)-scale lengths; high-iteration proptest; adversarial regeneration search |

Mechanics: T2 tests are `#[ignore]`-annotated; T1 runs via a dedicated Cargo profile (`[profile.ci]` inheriting `release` with `debug-assertions = true`, so §4's divisibility asserts and §9 debug invariants stay armed). Recommend `cargo-nextest` for per-test timeouts and for sharding T2 across CI jobs (`--partition count:k/m`).

### 1.2 Workspace layout (test-relevant view)

```
Cargo.toml                     # workspace
rust-toolchain.toml            # pinned stable
deny.toml                      # cargo-deny config
clippy.toml                    # disallowed-types: f32/f64 in ccd-core (see §23 F3)
crates/
  ccd-core/
    src/{partitions, hooks, characters(mn), spectra, transform, engine, checkpoint, output}.rs
    src/testing/               # #[cfg(any(test, feature="test-utils"))]: brute force,
                               #   naive MN reference, naive bigint kernels, union catalog
    tests/                     # integration tests (see below)
    benches/                   # criterion
  ccd-cli/
    tests/                     # CLI golden-file + checkpoint-kill tests
  ccd-py/                      # later phase (design in section 10)
fixtures/                      # committed SymPy-generated JSON (section 4)
schemas/result.v1.schema.json
scripts/gen_fixtures.py
scripts/find_adversarial.py    # one-time search, results committed as fixtures (section 5)
docs/
```

The `test-utils` feature gates all oracle code out of release builds while letting the CLI's integration tests reuse it.

---

## 2. Unit test plan (spec §22.1, §9) — concrete modules and names

All T0 unless marked. Names are `module::tests::name`.

### 2.1 `partitions::tests`
| Test | Checks | Spec |
|---|---|---|
| `partition_count_matches_oeis` | p(n) for n = 1..=50 against hardcoded OEIS A000041 values (generation only run to n ≤ 30 in T0; counts to 50 via the counting recurrence) | §3 |
| `canonical_order_is_documented_order` | first/last elements, strict ordering predicate holds pairwise, identity type `(1^n)` at its documented index, for n ≤ 12 | §3, §19.3 |
| `index_roundtrip` | `partition_to_index[partitions[i]] == i` for all i, n ≤ 20 | §3 |
| `parts_sum_to_n_and_sorted_descending` | structural invariant, n ≤ 20 | §3 |
| `cycle_type_spec_padding` | `[3]` in S₅ → `(3,1,1)`; `[2,2]` in S₄ → `(2,2)`; spec with sum > n rejected with typed error; identity spec `[]` maps to `(1^n)` and is rejected as a generator by default (§5.3) | fixed req. 2 |
| `transpose_is_involution` | λ'' = λ, all λ ⊢ n, n ≤ 20 | §11.1 |
| `self_transpose_count_matches_distinct_odd_parts` | #{λ = λ'} equals #partitions of n into distinct odd parts, n ≤ 20 | §11.2 |
| `z_lambda_known_values` | hand-computed z for a table of ~20 partitions across n = 3..8 | §3 |
| `class_sizes_sum_to_factorial` | Σ n!/z_λ = n! (bigint), n ≤ 30 | §3 |
| `sign_formula_matches_materialized_permutation` | sgn(λ) = (−1)^(n−ℓ) equals sign of an actual permutation of that type, n ≤ 10 | §3 |
| `even_class_sizes_sum_to_half_factorial` | Σ over even classes = n!/2 for n ≥ 2, n ≤ 30 | §17.1 |

### 2.2 `hooks::tests`
| Test | Checks |
|---|---|
| `degrees_match_hardcoded_tables` | full degree lists for n ≤ 6 against literature values |
| `degree_squares_sum_to_factorial` | Σ f_ρ² = n! (bigint) for n ≤ 40 — cheap, run in T0 to n ≤ 20, T1 to 40 |
| `degree_transpose_symmetry` | f_ρ = f_ρ′ for n ≤ 30 |
| `hook_degree_matches_mn_identity_column` | f_ρ = χ^ρ(1ⁿ) from the MN evaluator, n ≤ 12 (two independent formulas) |

### 2.3 `characters::mn::tests`
| Test | Checks | Spec |
|---|---|---|
| `identity_column_equals_hook_degrees` | n ≤ 14 | §22.2 |
| `trivial_and_sign_rows` | χ^{(n)} ≡ 1; χ^{(1ⁿ)}(ν) = sgn(ν), n ≤ 14 | §22.2 |
| `transposition_column_content_formula` | ω_ρ((2,1^{n−2})) = Σ contents of ρ, hence χ^ρ = f_ρ·(Σ contents)/C(n,2); closed form vs MN, n ≤ 16 | §4 |
| `n_cycle_column_hooks_only` | nonzero exactly on hook shapes, values ±1, n ≤ 16 | §7.1, §22.2 |
| `transpose_row_sign_relation` | χ^{ρ′}(ν) = sgn(ν)·χ^ρ(ν) full table, n ≤ 12 | §11.1 |
| `self_transpose_vanishes_on_odd_classes` | n ≤ 12 | §11.2 |
| `row_orthogonality_exact` / `column_orthogonality_exact` | Σ_ν |C_ν| χ^ρ χ^σ = n!·δ; Σ_ρ χ^ρ(ν)χ^ρ(μ) = z_ν·δ, bigint, n ≤ 10 (T0), n ≤ 12 (T1) | §22.2 |
| `repeated_cycle_lengths_hard_cases` | curated (ρ, ν) with many equal parts, long rim hooks, self-transpose ρ, at n = 10..14 — hardcoded expected values sourced from the fixture pipeline | §10.2, §22.2 |
| `modular_mn_matches_exact_mn` | full table mod p equals exact table reduced mod p, p ∈ {31, 2³¹−1-region primes}, n ≤ 12 | fixed req. 5 |
| `trie_shared_suffix_columns_match_naive_reference` | production suffix-sharing column DP vs a deliberately naive memoized MN in `testing::naive_mn`, random columns, n ≤ 14 | candidate idea 2 |

### 2.4 `spectra::tests`
| Test | Checks | Spec |
|---|---|---|
| `central_eigenvalue_exact_divisibility` | |C_λ|·χ^ρ(λ) ≡ 0 mod f_ρ for every (ρ, λ), n ≤ 14; failure message prints ρ, λ, both conventions (this assert also stays on in production, §4) | §4, §22.1 |
| `omega_transposition_equals_content_sum` | independent closed form | §4 |
| `omega_transpose_sign_relation` | ω_{ρ′}(λ) = sgn(λ)·ω_ρ(λ) | §11.1 |
| `theta_union_is_column_sum` | Θ = Ω·B consistency, random incidence matrices, n ≤ 12 | §16 |
| `mixed_parity_theta_even_odd_split` | θ_{ρ′} = θ_even − θ_odd reconstruction equals directly computed θ_{ρ′}, for unions with both parities — the §11.1 trap | §11.1, F5 |
| `active_row_set_definition` | ρ ∈ R iff some base χ^ρ(λ_j) ≠ 0; dropped rows all have θ = 0 for every union over the base | §7 |

### 2.5 Transform invariants — `tests/invariants.rs` (spec §9, all seven)

Parametrized over the deterministic **union catalog** (`testing::catalog`): all single classes for n ≤ 7; all 2-class unions for n ≤ 6; seeded random unions for n = 7..12; special cases: `[2]` (odd-only), `[3]` (even-only, A_n), `[2]∪[3]` (mixed parity), `[2,2]` in S₄ (non-generating, component = V₄), n-cycle (sparse column), near-identity `[2]` in S₂, empty union (rejected).

| Test | Checks | Tier | Spec |
|---|---|---|---|
| `radius0_identity_with_full_row_set` | with the FULL irreducible row set, r = 0 coefficient vector is the indicator of (1ⁿ) | T0 n ≤ 10 | §9.1 |
| `radius0_restricted_rows_documented_failure` | negative test: restricted-row transform at r = 0 does *not* reconstruct identity; engine is asserted to initialize identity manually and start at r = 1 | T0 | §7 radius-zero exception |
| `radius1_single_class_indicator` / `radius1_union_indicator` / `radius1_weighted_equals_weights` | best orientation test in the suite | T0 | §9.2 |
| `total_word_count_equals_union_size_power` | Σ_ν |C_ν|·a_r(ν) = |U|^r for r = 0..=min(diam+2, 8), exact bigint | T0 n ≤ 8, T1 n ≤ 12 | §9.3 |
| `numerators_divisible_by_factorial` / `coefficients_nonnegative` | every entry, every radius to stopping | T0/T1 | §9.4, §9.5 |
| `parity_filter_single_parity` | odd-only and even-only unions: support at radius r contained in sign ε^r classes | T0 | §9.6 |
| `parity_filter_invalid_for_mixed` | negative test: mixed union produces both parities at some radius | T0 | §9.6, F5 |
| `stopping_no_new_types_after_empty_layer` | after engine stops at radius r*, run 2 further exact radii and assert new-set stays empty | T0 n ≤ 8 | §5.2 |
| `union_square_expansion` | (K_A + K_B)² = K_A² + 2·K_A K_B + K_B² at the coefficient level via the mixed-product spectral vectors (§16.3), n ≤ 7 | T1 | §22.4 |
| `eigenvalue_grouping_matches_ungrouped` | §15 H_α-grouped path produces identical a_r to the ungrouped path, n ≤ 10, several unions | T1 | §15 |
| `batched_unions_match_single_runs` | m unions batched vs run one-by-one: identical distances, supports, diameters; also after mid-run column compaction | T1 | §8.2 |

---

## 3. Brute-force cross-validation harness (Rust, test-only)

Location: `ccd-core/src/testing/bruteforce.rs` behind `test-utils`; drivers in `tests/bruteforce.rs`.

### 3.1 Group enumeration and BFS oracle
- Permutations as `Vec<u8>`; ranking via Lehmer code for O(1) index into flat arrays. n = 10 → 3,628,800 elements ≈ 36 MB of permutation data + 3.6 MB `u8` distance array. Fine.
- Generator materialization: filter all of S_n by cycle type (cycle-type computation of 3.6M permutations is < 1 s in release; simplicity beats a direct class enumerator here).
- BFS from identity over raw permutations, frontier as ranked indices, edge relaxation by composition. Parallelize frontier expansion with rayon for n = 9, 10.
- **Cost guard**: skip (n, U) pairs with |S_n|·|U| > 10¹⁰ edges (e.g. n = 10 with n-cycles: 1.3×10¹² — excluded); the nightly catalog for n = 9, 10 uses classes with |C| ≤ ~5000 ([2], [3], [2,2], [4], small unions).

### 3.2 Comparisons (all against the character engine, exact backend)
| Test | What is compared | Tier |
|---|---|---|
| `bfs_distances_by_type_match_engine_n4_8` | per-cycle-type distances (BFS distance asserted **constant on each class** first — validates conjugacy invariance and catches F6-adjacent bugs), diameter of identity component, reachable-type set, per-radius *new* sets | T1 |
| `bfs_distances_n9`, `bfs_distances_n10` | same, restricted catalog | T2, `#[ignore]` |
| `exact_length_supports_match_setproduct_n4_8` | boolean set-product DP over the group: S_{r+1} = S_r·U (bitset over ranked perms), run to r = diam+2; per-radius **exact-length supports** classified by type — this checks support_r, not just first-seen distance (§5.1) | T1 |
| `word_counts_match_group_algebra_dp_n4_7` | full DP c_{r+1}[g·h] += c_r[g] over Z^{n!} with `u64` + `checked_add` (bound |U|^4 ≤ 2520⁴ ≈ 4×10¹³ < 2⁶³, asserted); compare a_r(ν) for r ≤ 4 **and assert count constancy on every class**. Counts are not user-facing output but are the internal quantity everything rests on — they must be validated | T1 |
| `word_counts_n8_r3` | same at n = 8, r ≤ 3 | T2 |
| `disconnected_component_v4_in_s4` | `[2,2]` in S₄: engine reports exactly {(1⁴), (2,2)} reachable, diameter 1, disconnected flag set | T0 |
| `mixed_union_bfs` | `[2]∪[3]` and other mixed-parity unions through the full comparison stack — the parity-trap case | T1 |

---

## 4. SymPy fixture pipeline

### 4.1 Honest feasibility assessment
SymPy provides exact integers, `IntegerPartition`, and p(n), but (verify at implementation time) **no ready-made S_n character table**. "SymPy fixtures" therefore means: a Python script implementing an *independent* Murnaghan–Nakayama evaluator, cross-certified inside the script by closed forms that do not share code with either MN implementation. A candidate third evaluator via the Frobenius determinant/coefficient formula was evaluated and **rejected as a general mechanism**: extracting the monomial coefficient costs ~ℓ_min!·(coefficient-DP) and the DP state space explodes beyond ℓ_min ≈ 6, so it adds little beyond n ≈ 14 where the layered certificate below is already stronger per CPU-second. (Optionally keep it as a spot-checker for ρ with min(ℓ(ρ), ρ₁) ≤ 6.)

The **layered certificate** each generated full table must pass inside the script before being written:
1. identity column = hook-length degrees (independent closed form, bigint);
2. transposition column = f_ρ·(Σ contents of ρ)/C(n,2) (independent closed form);
3. full exact orthogonality: X·diag(|C_ν|)·Xᵀ = n!·I (for n = 16: 231³ ≈ 1.2×10⁷ bigint mults, ~1 min in Python);
4. transpose-sign relation χ^{ρ′} = sgn·χ^ρ and self-transpose-vanishing-on-odd;
5. sign row and trivial row.
Orthogonality pins the table up to row permutation/negation; (1) fixes negation and pins rows up to degree ties; (2) separates transpose-paired rows (content sums are negatives); residual ties are astronomically unlikely and are additionally covered by cross-language MN agreement. Any orientation error fails check (1) or (2) loudly.

### 4.2 What `scripts/gen_fixtures.py` generates
| Artifact | Range | Est. runtime (one-time) | Est. size |
|---|---|---|---|
| Full character tables | n ≤ 14 always; n = 15, 16 with `--large` | n=14 ~minutes; n=16 ~10–30 min (memoized MN in Python) | n=14: 135² values ≈ 200 KB; n=16: 231² ≈ 0.5–1 MB; total ≲ 2 MB committed |
| Spot values | n = 17..28: ~300 seeded-random (ρ, ν) pairs per n **plus** structured columns (identity, transposition, (2,2,1^{n−4}), n-cycle) and hard cases (repeated cycle lengths, long rim hooks, self-transpose ρ) | minutes–hours per n; chunked | a few KB per n |
| Degrees + class sizes + z + signs + transpose map | all n ≤ 30 | seconds | small |
| Adversarial modular tuples (see section 5) | n ≤ 9 | one-time search | tiny |

Chunking/restartability: `--n 20 --columns 120:180` emits partial files; a `--merge` pass produces the final fixture and recomputes the checksum. Beyond n ≈ 18 full tables in Python are impractical (p(18)² ≈ 148k memoized-recursive evaluations); spot mode is the design, not a fallback.

### 4.3 Fixture file format (`fixtures/chartab_n{N}.json`, schema per spec §19.3)
```json
{
  "schema_version": "fixture.v1",
  "generator": {"script": "gen_fixtures.py", "script_version": "...",
                 "sympy_version": "...", "python_version": "...", "seed": 12345},
  "n": 14,
  "partition_order_convention": "<project canonical name, e.g. revlex-descending>",
  "partitions": [[14],[13,1], "..."],
  "row_col_convention": "rows=irreps(rho), cols=classes(nu)",
  "degrees": ["...decimal strings..."],
  "class_sizes": ["..."], "z_values": ["..."], "signs": [1,-1, "..."],
  "transpose_map": [7, 3, "..."],
  "table": [["..."]],
  "spot_values": [{"rho": 12, "nu": 30, "value": "-4084"}],
  "certificates_passed": ["hook_degrees","content_column","orthogonality","transpose_sign"],
  "payload_sha256": "..."
}
```
**Ordering caveat**: the script must not use SymPy's native partition iteration order. It implements the project's canonical enumeration itself and embeds the explicit list; there is no implicit mapping to document because there is no implicit order anywhere.

### 4.4 Rust consumers — `tests/fixtures.rs`
- `fixture_checksum_and_schema_valid` (T0)
- `fixture_partition_order_matches_rust_enumeration` — asserted **before** any value comparison; a mismatch aborts with both lists printed (T0)
- `fixture_degrees_match_hooks`, `fixture_class_sizes_match`, `fixture_transpose_map_matches` — these are *doubly* independent (Rust closed forms vs Python closed forms), so a passing run certifies orientation from both sides (T0)
- `fixture_full_table_matches_mn` — n ≤ 14 in T1, n = 15, 16 in T2
- `fixture_spot_values_match_mn` — n ≤ 20 in T1; n = 21..28 in T2 (single MN column evaluations at n = 28 are nontrivial)
- `fixture_regen_is_byte_identical` — CI job (section 8) reruns the script with pinned SymPy and diffs against the committed files

---

## 5. Modular-path adversarial tests

The modular path has exactly one catastrophic failure mode: **treating a zero residue as a certified zero** (Failures 4/9). Tests attack it at three levels.

### 5.1 Engineered prime-divisible positive coefficients (end-to-end)
Against production ~31-bit primes, a natural coefficient divisible by a given prime occurs with probability ~5×10⁻¹⁰ — random search is useless. Therefore:
- The modular backend takes its prime list as an explicit parameter (already required for testability). Tests inject **small primes** (still > n, e.g. p ∈ {11, 13} for n ≤ 10).
- `scripts/find_adversarial.py` (one-time, results committed): for n = 6..9 and the union catalog, compute exact a_r(ν) for all radii to diameter+2; search for tuples where a_r(ν) > 0 and a_r(ν) ≡ 0 mod p for *every* prime in a chosen small-prime set (i.e. divisible by the product, e.g. 143). Coefficients grow like |U|^r, so hits are plentiful by r ≈ 4–6. Commit found tuples as `fixtures/adversarial_v1.json` with provenance.
- `tests/adversarial.rs::injected_primes_false_zero_layer_does_not_stop_engine` (T1): run the full modular engine with the injected primes on a committed tuple; assert (a) the entry lands in the candidate-zero set, (b) exact fallback resurrects it as positive, (c) final distances equal the exact-backend reference, (d) the engine did **not** stop at the modularly-empty-looking layer.
- `adversarial::screened_positive_is_never_wrong` (T1): for n ≤ 9, validate every "definitely positive" screening verdict against exact coefficients — a nonzero residue must always correspond to a positive integer (this direction should never fail; it guards residue-arithmetic bugs).

### 5.2 Synthetic classifier tests (kernel plumbing, no character data)
Feed the support-classification stage hand-built numerator vectors:
- entry = n!·p₁·p₂·…·p_k (positive, ≡ 0 mod all screening primes) → must route to exact fallback → classified positive;
- entry = 0 → fallback → certified zero;
- entry with any nonzero residue → must **not** invoke fallback (asserted via an injectable/counting fallback hook on the backend trait — the trait design must expose fallback invocations for observability);
- `stopping_gate_requires_certified_layer`: the radius loop must be structurally unable to test the stopping criterion before the candidate-zero set is empty — unit-test the gate directly.

### 5.3 Overflow-bound stress (§13.1, Failure 8)
Design contract: residues < p < 2³¹ ⇒ products < 2⁶²; `u128` row accumulation is safe for row length R·(p−1)² < 2^80 ≪ 2¹²⁸ for R ≤ 2¹⁸ ≥ p(50) = 204,226. Every kernel documents its bound; tests enforce it:
- `kernels::tests::worst_case_all_entries_p_minus_1` — dot products at lengths {p(45), p(50), 10·p(50)} with all entries p−1, vs a naive `BigUint` reference kernel in `testing::naive_kernels` (T0 small, T2 full lengths);
- `kernels::tests::block_reduction_boundary` — any kernel with block-wise reduction at block size K is tested at lengths K−1, K, K+1, 2K+1;
- `kernels::tests::proptest_random_vs_bigint_reference` — random lengths/values/primes (T0, low iterations; T2, high);
- a `const_assert!` encoding the documented bound formula next to each kernel, so changing the prime width or accumulator type without re-proving the bound fails compilation.

### 5.4 Differential modular-vs-exact (randomized, seeded)
`tests/differential.rs::modular_engine_matches_exact_engine` (T1): n ∈ 5..=12, seeded random unions, both engines run to stopping; assert identical distances, per-radius new-sets, exact-length supports, diameters, reachability flags. For n ≤ 9 additionally assert every residue equals the exact coefficient reduced mod p at every radius. Seeds fixed in code; a `CCD_TEST_SEED` env override enables soak runs.

### 5.5 Production-scale self-validation (no fixtures exist at n = 40)
These run inside the engine (debug-assert level, plus an opt-in `--self-check` CLI flag) and in T2 tests at n = 30–35:
- **modular word-count identity**: Σ_ν |C_ν|·a_r(ν) ≡ |U|^r (mod p) at every radius, per prime — O(q) per radius, catches silent corruption at full scale (§9.3 modularized);
- **modular radius-1 indicator**: valid with restricted rows (dropped rows have θ = 0, which cannot happen for r ≥ 1 contributions); catches orientation errors (F1/F7) at n = 40 for pennies;
- **post-generation table self-checks** per prime: identity column ≡ hook degrees mod p (hook degrees computed exactly, reduced); transpose-sign relation over the full modular table; **randomized modular row-orthogonality**: Σ_ν |C_ν|·χ^ρ χ^σ ≡ n!·δ (mod p) for ~10⁴ sampled (ρ, σ) pairs — ~4×10⁸ ops at n = 40, seconds;
- **analytic regression anchors** at any n: for U = class of transpositions, d(ν) = n − ℓ(ν) exactly (counting fixed points as cycles) and diameter = n − 1 — a full closed-form check of the entire pipeline at production scale; odd-only unions must yield bipartite layer structure. (Curated literature results for 3-cycles and n-cycles can be added later with citations; do not hardcode uncited formulas.)

---

## 6. Property-based tests

**Choice: `proptest`** over quickcheck — better shrinking of composite structures, strategy combinators fit partition/union generation, actively maintained. Config: 64 cases in T0, 2048 in T2 via `PROPTEST_CASES`.

Strategies in `testing::strategies`: `arb_partition(n)` (random descending part sequences, not uniform — uniformity unnecessary), `arb_cycle_type_spec(max_part, max_sum)`, `arb_union(n, max_classes)`, `arb_prime_set()`, `arb_checkpoint_radius()`.

| Property | Statement |
|---|---|
| `prop_transpose_involution` | λ'' = λ |
| `prop_spec_padding_sums_to_n` | padded cycle type sums to n for every valid (spec, n) |
| `prop_class_size_sum` | Σ n!/z = n! for random n ≤ 30 (bigint) |
| `prop_sign_consistency` | sgn from formula = sgn of materialized permutation, n ≤ 10 |
| `prop_modular_mn_matches_exact` | random (ρ, ν, p), n ≤ 14 |
| `prop_batched_equals_sequential` | random union batch vs one-by-one runs, n ≤ 9 — full output equality |
| `prop_theta_linear` | θ(Ω, b₁ + b₂) = θ(b₁) + θ(b₂) for disjoint unions (weighted semantics) |
| `prop_monotone_distances_under_inclusion` | U ⊆ V ⇒ d_V(ν) ≤ d_U(ν) on common reachable set (§16.2), n ≤ 8 |
| `prop_determinism_across_thread_counts` | identical outputs with rayon pool sizes 1 vs many (modular sums are order-independent; assert it stays true) |
| `prop_resume_equivalence` | run to completion vs forced checkpoint at random radius + resume in fresh state ⇒ byte-identical semantic output (first-class requirement, fixed req. 6) |

---

## 7. Checkpoint, output-schema, and CLI tests

### 7.1 `tests/checkpoint.rs`
- `resume_after_kill_matches_uninterrupted` — CLI-level: spawn `ccd-cli run`, SIGKILL/TerminateProcess mid-run at n ~ 14, rerun with `--resume`; compare final JSON (T1);
- `resume_rejects_config_hash_mismatch` / `resume_rejects_partition_order_hash_mismatch` — mutate spec or order metadata → typed refusal, never silent reuse (§19.3);
- `checkpoint_write_is_atomic` — checkpoint protocol must be write-temp-then-rename; test simulates a truncated temp file and asserts the previous checkpoint remains loadable;
- `corrupted_checkpoint_detected` — bit-flip → checksum failure with clear error;
- `checkpoint_version_field_forward_compat` — unknown future version → refusal, not panic.

### 7.2 `tests/output.rs`
- `output_validates_against_json_schema` — every emitted result validated against `schemas/result.v1.schema.json` (via the `jsonschema` crate in tests);
- `output_embeds_partition_order_and_metadata` — explicit partition list, order-convention name, spec §19.3 metadata fields, versioned `schema_version`;
- `golden_outputs_small_n` — committed expected JSON for the n ≤ 8 catalog; comparison strips volatile fields (timestamp, duration, host), which the schema isolates under a dedicated `run_env` key;
- `output_roundtrip_parse` — serialize → parse → semantic equality;
- `no_word_counts_in_output` — negative test: the schema and serializer expose supports/distances/diameter/parity/reachability only (fixed req. 3).

---

## 8. Failure-mode coverage matrix (spec §23)

| Failure | Covering tests |
|---|---|
| F1 wrong normalization | `radius1_*`, modular radius-1 self-check at scale |
| F2 coefficient vs total words | `total_word_count_equals_union_size_power`, `word_counts_match_group_algebra_dp` |
| F3 floating-point zero tests | structural: `clippy.toml` `disallowed-types = [f32, f64]` in `ccd-core`, enforced by CI clippy `-D warnings` |
| F4/F9 stopping on uncertified layer | section 5.1/5.2 adversarial family, `stopping_gate_requires_certified_layer` |
| F5 parity on mixed unions | `parity_filter_invalid_for_mixed`, `mixed_union_bfs`, `mixed_parity_theta_even_odd_split` |
| F6 split A_n classes | unrepresentable by construction (input = S_n cycle types); documented in README + `bfs` tests include even classes lying in A_n |
| F7 lost partition order | order-first assertion in every fixture/checkpoint/output consumer; hash-mismatch refusal tests |
| F8 modular overflow | section 5.3 stress + `const_assert!` bounds |
| F10 Booleanization | word-count sum at r ≥ 2 and DP count comparison would both fail loudly |
| F11 geodesic witnesses | out of scope; README states it |

---

## 9. Benchmarks and the scaling harness

### 9.1 Criterion benches (`ccd-core/benches/`, not run in gating CI)
| Bench | Sweep | Purpose |
|---|---|---|
| `bench_mn_table` | full modular table, n ∈ {20, 25, 30} (35 behind env flag), per prime; reports cells/sec | table-generation budget; validates modular-first idea (candidate 1) |
| `bench_mn_column` | single columns: transposition-like, n-cycle, balanced ν; with vs without suffix-sharing trie | quantifies candidate idea 2 rather than assuming it |
| `bench_modular_gemv` | per-radius transform at n ∈ {25, 30, 35}; 1/8/64 union columns; rayon threads 1/2/4/8/16 | amortization + scaling curves for backend-trait doc |
| `bench_exact_fallback` | BigUint recomputation of one coefficient at radii of varying magnitude | prices the fallback so screening-prime count can be tuned |
| `bench_end_to_end` | full run per (n, class), n ∈ {12, 16, 20}, classes [2], [3], [2,2], n-cycle | regression tracking |

Nightly CI stores criterion baselines as artifacts; `critcmp` against the previous baseline, > 20 % regression posts an informational annotation (never a gate — runner noise).

### 9.2 Scaling-run harness (feeds auto-n mode)
Not criterion: a CLI subcommand `ccd-cli bench-scaling --n 20..=32 --classes '[2];[3];[2,2]' --primes 3 --emit scaling.json` measuring per (n, U): table-generation wall time and bytes per prime (model: 4·R·q bytes; n = 40 ⇒ ~5.6 GB/prime, n = 45 ⇒ ~32 GB/prime, n = 50 ⇒ ~167 GB/prime — the harness validates the model against measured RSS), per-radius transform time, radii to stopping, fallback count, peak RSS, **and eigenvalue-multiplicity statistics** (#distinct θ vs #active rows) — turning candidate idea "eigenvalue grouping" (§15) from an assumption into a measured quantity with a documented fallback (streamed ungrouped transform) if multiplicity is low. Output `scaling.json` is the resource-estimation table consumed by auto-n budget selection and is versioned like other outputs.

---

## 10. CI design (GitHub Actions)

**Step 0 (prerequisite): the repo is not yet a git repository.** `git init` with default branch `main`, `.gitignore` (`/target`, criterion output, fixture scratch dirs), commit fixtures and schemas, then add workflows.

### 10.1 `ci.yml` — push + PR
- Concurrency group with cancel-in-progress; `rust-toolchain.toml`-pinned stable; `Swatinem/rust-cache`; all cargo calls `--locked`.
- **Linux job (ubuntu-latest)**: `cargo fmt --check` → `cargo clippy --workspace --all-targets -- -D warnings` (picks up the f64 ban) → `cargo deny check` (EmbarkStudios action) → T0 `cargo test --workspace` → T1 `cargo test --workspace --profile ci`.
- **Windows job (windows-latest)**: clippy + T0 + a T1 subset (checkpoint kill/resume tests MUST run here — process-termination semantics differ on Windows, and the dev machine is Windows 11). Skip fmt/deny (redundant).
- Optional allowed-failure beta-Rust job; optional MSRV check against `rust-version`.

### 10.2 `slow.yml` — nightly cron + manual
`cargo nextest run --profile ci -- --ignored`, sharded 2–4 ways via nextest partitioning; 6 h timeout; uploads logs + criterion baselines. Linux only.

### 10.3 `fixtures.yml` — manual
Sets up Python with **pinned** SymPy version, reruns `gen_fixtures.py` for n ≤ 14 + spot ranges, `git diff --exit-code fixtures/` — guards fixture drift and keeps the script executable.

### 10.4 `cargo-deny` and dependency policy
Licenses allowlist: MIT, Apache-2.0, BSD-*, Unicode, Zlib. **Bigint decision consequence**: default `num-bigint` (MIT/Apache, pure Rust); optional `fast-bigint` feature using `ibig` (MIT/Apache) if profiling demands it; `malachite` is LGPL-3.0 — usable for a research tool but would force a deny.toml exception and wheel license notices, so it is *not* the default; `rug` excluded entirely (GMP/MSYS2 breaks the Windows-native requirement). RUSTSEC advisories deny; duplicate-version bans warn.

### 10.5 `wheels.yml` — on tags (later phase)
`PyO3/maturin-action`: manylinux_2_17 x86_64 (Kaggle target) + windows x86_64; abi3 so one wheel per platform covers Python ≥ 3.9; wheels attached to GitHub release and uploaded as artifacts for Kaggle-dataset packaging.

---

## 11. Packaging path (design only, later phase)

### 11.1 Crate layout
- `ccd-core`: math + engine, no I/O, no CLI deps; `test-utils` feature.
- `ccd-cli`: clap; subcommands `run`, `resume`, `bench-scaling`, `self-check`; owns JSON emission and checkpoint driving.
- `ccd-py`: `crate-type = ["cdylib"]`, `pyo3` with `abi3-py39`; `pyproject.toml` with maturin backend; Python module name e.g. `cayley_chars`.

### 11.2 Python API sketch
```python
import cayley_chars as cc
spec = cc.RunSpec(n=range(20, 36), generators=[[[3]], [[2],[3]]],
                  budget=cc.Budget(wall_seconds=6600, mem_bytes=250_000_000_000),
                  checkpoint_dir="/kaggle/working/ckpt")
result = cc.run(spec, resume=True, progress=print)   # plain dicts, same schema as CLI JSON
cc.load_result("/kaggle/working/out/n30_3cycles.json")
```
Design points: the Python layer is a thin driver over the identical serializer the CLI uses (one schema doc serves both); compute releases the GIL (`py.allow_threads`); `KeyboardInterrupt` and the internal budget watchdog both flush a checkpoint before returning; a `progress` callback reports (n, union, radius, unresolved-zeros count).

### 11.3 Kaggle deployment notes
- Build abi3 manylinux wheels; upload as a **Kaggle Dataset**; notebooks (often internet-off) install with `pip install /kaggle/input/<ds>/cayley_chars-*.whl --no-deps`.
- Checkpoints to `/kaggle/working` (persisted on save). The 2 h hard limit ⇒ default budget ≈ 105–110 min with margin; checkpoint every radius *and* on watchdog expiry.
- Resume across sessions: previous notebook output attached as an input dataset is **read-only** — the driver copies checkpoint files into `/kaggle/working` before resuming (documented in `docs/kaggle.md`).
- ~96 vCPUs: `RAYON_NUM_THREADS` honored; determinism-across-thread-counts is already a tested property, so Kaggle results are reproducible locally.

---

## 12. Documentation set

| Doc | Contents |
|---|---|
| `README.md` | what/why, quickstart (CLI + Python), input spec format, explicit non-goals (no geodesic witnesses F11, no split A_n classes §17.2, no non-class generators §17.3) |
| `docs/math_to_code.md` | formula → function table with spec-section anchors: ω_ρ(λ) → `spectra::central_eigenvalue` (§4), a_r formula → `transform::inverse_transform` (§4.2/§25), MN recursion → `characters::mn` (§10.2), hooks → `hooks::degree`, H_α grouping → `engine::grouped` (§15), stopping gate → `engine::stopping` (§5.2/§12.2) |
| `docs/output_schema.md` + `schemas/result.v1.schema.json` | versioned JSON schema, partition-order embedding rules, worked example, stability policy |
| `docs/gpu_backend.md` | the transform backend trait contract: residue preconditions (< p < 2³¹), accumulator bounds an implementation must prove, batching shapes, determinism requirement, candidate-zero/fallback protocol, so a cudarc/wgpu backend can be added without touching math |
| `docs/testing.md` | tier system, how to run T2, seeds, regenerating fixtures and adversarial tuples |
| `docs/checkpoints.md` | format, atomicity protocol, hash-refusal semantics |
| `docs/kaggle.md` | wheel-from-dataset install, budget/resume workflow |
| `docs/fixtures.md` | the layered-certificate argument (section 4.1) — why the fixtures are trustworthy without GAP |

---

## 13. Critical evaluation of candidate ideas (test-relevant verdicts)

- **Modular-first MN table**: adopt; testing consequences fully covered by 2.3 differentials + section 5.5 at-scale self-checks (which are *mandatory* because no fixture exists at n ≥ 30).
- **Suffix-sharing column trie**: adopt but treat as an optimization that must remain differentially tested against the committed naive MN reference forever (`trie_shared_suffix_columns_match_naive_reference`); the bench `bench_mn_column` quantifies whether it earns its complexity.
- **~31-bit primes, u128 accumulation**: adopt; bounds proven and `const_assert!`ed (5.3). The injectable prime list is not just a testing convenience — it is the mechanism that makes false-zero layers testable at all.
- **Eigenvalue grouping (§15)**: adopt as an optional path, gated on measurement: the scaling harness emits distinct-eigenvalue counts; the correctness test `eigenvalue_grouping_matches_ungrouped` exists from day one; fallback is the streamed ungrouped transform.
- **SymPy full tables for n ≤ 16 "or so"**: confirmed feasible at n ≤ 16 only; beyond that, spot values plus the layered certificate replace full tables (section 4.1); a one-off optional GAP/Sage export could add third-party tables for n ≤ 12 but is not required by the trust argument.
- **Checkpoint per-(n, union) state**: adopt; elevated to a *tested property* (`prop_resume_equivalence`) plus process-kill integration tests on both OSes, because Kaggle's 2 h wall is a hard functional requirement, not an ops nicety.

Open risks to track: (1) verify at implementation time whether current SymPy has any native S_n character support — if it does, add it as one more independent certificate layer; (2) GitHub-hosted runners may be too small for n = 16 T2 invariants — measure early, shard or self-host if needed; (3) the golden-output files freeze the schema — land `schema_version` and the volatile-field isolation (`run_env`) before the first golden file is committed.