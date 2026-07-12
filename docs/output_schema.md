# `classdiam/result` schema, version 1

One JSON document per `(n, generating union)`. Parsers must check
`format == "classdiam/result"` and reject unknown `format_version` majors.
Serde types: `crates/core/src/report/schema.rs`. A worked, hand-verified
example is committed as `fixtures/golden/n06_g2.json` (S₆, transpositions).

All indices below refer to the **canonical partition order**
(`partition_order.convention = "lex_desc_full_parts_v1"`): partitions of `n`
as weakly decreasing part lists, lexicographically descending — index 0 is
`[n]`, the last index is the identity type `[1,…,1]`. Every document embeds
the explicit order (`partitions_reduced`, parts ≥ 2 only, identity = `[]`)
and a blake3 hash of its encoding; consumers must verify the order before
using any indexed field (spec §19.3 / Failure 7).

## Fields

| field | meaning |
|---|---|
| `format`, `format_version`, `spec_version` | schema identity; `spec_version` names the mathematical spec revision implemented |
| `tool` | producing binary and version. **Volatile** for golden comparisons |
| `run` | run id, UTC timestamps, resume counters. **Volatile** |
| `n`, `factorial_n` | degree; `n!` as a decimal string |
| `generators.input_templates` | cycle types exactly as given, WITHOUT fixed points |
| `generators.classes[]` | per class: `template`, `padded` (full cycle type), canonical `index`, `class_size` (decimal string), `sign` |
| `generators.union_size` | `\|U\| = Σ \|C_λ\|` (decimal string) |
| `generators.parity` | `even` / `odd` / `mixed` — drives the parity notes below |
| `partition_order` | convention name, `count = p(n)`, blake3 hash, explicit reduced list |
| `row_col_convention` | orientation of every character-indexed structure |
| `class_data.sign`, `class_data.class_size` | parallel arrays over the order — documents are self-contained |
| `arithmetic` | `mode` (`exact-bigint` for the reference engine; `modular+certified` later), resident primes, certifier |
| `engine` | engine mode/backend, `rows_used` (`full`/`restricted`), `active_row_count`, `zero_rows_all_bases` (rows with `ω_ρ(λ) = 0` for every generator class — droppable for `r ≥ 1`, spec §7), `threads` (**volatile**) |
| `results.distance[]` | exact Cayley distance from the identity per type; `unreachable_value` (−1) marks types outside the identity component |
| `results.first_hit_even[]` / `first_hit_odd[]` | first EVEN/ODD radius `r` with `a_r(ν) > 0`, **within the run span**; −1 = none up to `stop_radius`. `distance = min` of the two |
| `results.diameter_identity_component` | max finite distance — the graph diameter of the identity component (vertex-transitivity, spec §5) |
| `results.generated_subgroup` | `S_n` / `A_n` / `proper_subgroup` / `trivial`, derived from the **computed** reachable set, never from theory |
| `results.cayley_graph_on_sn` | `connected` iff every type is reachable |
| `results.bipartite` | true iff the union is odd-only (spec §17.1) |
| `results.stopping` | `rule` (`empty_layer` = spec §5.2 on exact supports; `all_types_visited` = parity-feasible cover early exit) and `stop_radius` |
| `results.support_reconstruction` | the documented rule below, as a string (changes when the identity class is a generator) |
| `results.layers[]` | per radius `r ≤ stop_radius`: `new` (types first reached at exactly `r`), `support` (exact-length support: `a_r(ν) > 0`), `support_size` |
| `timings_s` | wall-clock seconds. **Volatile** |
| `config_hash_blake3` | hash of the resolved run configuration |

## Semantics notes

- **Exact-length supports are not BFS frontiers** (spec §5.1): `support` may
  contain types seen at smaller radii.
- **Support reconstruction.** For unions without the identity class, within
  the run span: `support_r = { i : 0 ≤ first_hit[r mod 2][i] ≤ r }` — exact by
  the `a_{r+2}(ν) ≥ a_r(ν)` monotonicity lemma (append `u·u⁻¹`). The rule is
  stated **valid for `r ≤ stop_radius` only**: an opposite-parity first hit
  occurring after the stopping radius is not recorded (example: `[2,2]` in
  S₄ stops at r = 2 but the identity first appears at an odd radius at r = 3).
  With the identity class as a generator, supports are cumulative:
  `support_r = { i : 0 ≤ distance[i] ≤ r }`.
- **Word counts are deliberately absent** everywhere (support-only output —
  fixed project decision); `a_r(ν)` values are internal.
- **Golden-file comparisons** strip exactly: `run`, `timings_s`,
  `engine.threads`, `tool`.
