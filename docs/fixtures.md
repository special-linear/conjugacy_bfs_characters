# Fixture trust chain (why no GAP is needed)

The runtime character backend is the Rust Murnaghan–Nakayama evaluator
(`chars::mn`, abacus/beta-set rim hooks). It is validated against committed
fixtures generated once by `tools/gen_fixtures.py` — a **second, independent
MN implementation** (Python, hook-cell/rim-path enumeration on the Young
diagram; different algorithm, different language, different author path).
SymPy provides an independent partition-enumeration cross-check; SymPy 1.14
has no S_n character tables of its own (verified at implementation time).

Agreement between two implementations is not blind trust: the Python
generator **refuses to write** any table that does not pass the layered
certificate (design doc 03 §4.1):

1. identity column = hook-length degrees (independent closed form);
2. transposition column: `χ·C(n,2) = f·(Σ contents)` (independent closed form);
3. full exact row orthogonality `X·diag(|C|)·Xᵀ = n!·I` (sampled for n = 16);
4. transpose-sign relation + self-transpose rows vanish on odd classes;
5. trivial row ≡ 1, sign row ≡ sgn.

Orthogonality pins the table up to symmetries that (1) and (2) then fix; a
row/column orientation error fails (1) or (2) loudly. The same closed forms
are asserted independently on the Rust side (`chars::mn` tests), and the Rust
consumer tests (`crates/core/tests/fixtures.rs`) additionally:

- verify the payload sha256 (rule: sha256 of compact-JSON segments joined by
  `|` — see `payload_sha256` in both the script and the tests);
- assert the **partition order matches the Rust canonical enumeration before
  any value comparison** (spec Failure 7 discipline);
- compare every fixture value against the Rust evaluator.

## Files

| file | contents | range |
|---|---|---|
| `fixtures/characters/char_nNN.json` | full character tables `table[nu][rho]` | n ≤ 14 (n = 15–16 via `--full-max 16`) |
| `fixtures/characters/spot_nNN.json` | spot values `[rho, nu, value]`: structured columns (identity, transpositions, `[2,2]`, n-cycle), self-transpose/rectangular hard cases, seeded randoms | n = 17–22 |
| `fixtures/degrees/deg_nNN.json` | degrees, class sizes, `z`, signs, transpose map (decimal strings) | n ≤ 30 |
| `fixtures/adversarial_v1.json` | `(n, union, r, ν, primes)` with `a_r(ν) > 0` divisible by every listed small prime — forces the P2 modular engine through exact certification; `masks_first_hit` marks tuples that would corrupt distances if screening were trusted (spec §23 F4/F9) | n = 6–9 |
| `fixtures/golden/n06_g2.json` | golden result document (S₆, transpositions) | — |

Regeneration: `python tools/gen_fixtures.py && python tools/find_adversarial.py`
(pinned SymPy; CI job `fixtures.yml` re-runs and diffs). Golden file:
`UPDATE_GOLDEN=1 cargo test -p classdiam-core --test output`.

Trust chain for the whole system: brute-force BFS over raw permutations
(n ≤ 10) → exact bigint spectral engine (validated against BFS and these
fixtures) → the P2 modular engine (differentially tested against the exact
engine, with the adversarial tuples attacking its certification gate).
