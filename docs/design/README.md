# Design documents

Produced 2026-07-12 by a three-designer + adversarial-critic panel over the
mathematical spec (`notes/character_method_cayley_diameters.md`).

| File | Contents |
|---|---|
| `00-plan.md` | **The merged, approved implementation plan — authoritative.** Where the documents below conflict, the resolution in the plan is final. |
| `01-architecture.md` | Crate layout, data model, `CharacterSource`/`TransformBackend` traits, engine state machine, orchestration, JSON schema (incl. the verified n=6 worked example), checkpointing, CLI, PyO3 surface. |
| `02-numerics.md` | MN evaluator (CSR level operators, suffix-sharing column DP), modular kernels with proven overflow bounds, transpose pairing, certification subsystem, n ≥ 50 strategies, end-to-end cost model, GPU contract. |
| `03-testing.md` | Test tiers, brute-force oracles, SymPy fixture pipeline with layered certificate, adversarial modular tests, property tests, benches, CI, packaging. |
| `04-critique.md` | Adversarial review: 18 findings (1 blocker, 5 major) with fixes, plus the list of independently verified design decisions. All findings are resolved in `00-plan.md`. |

Known conflicts resolved by the plan (critique finding 1): stopping rule and
per-radius certification gate follow `01`; transpose pairing (W±) and the cost model
follow `02`; prime regime is 31-bit/u128 on CPU (30-bit/u64 documented for GPU only);
`TransformBackend` follows `01` extended with the parity-split weight pair; eigen-grouped
mode is a separate per-union path; naming follows `01` (`classdiam-*`,
`lex_desc_full_parts_v1`, `tools/`).
