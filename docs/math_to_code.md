# Formula → code map

Spec references are to `notes/character_method_cayley_diameters.md`.

| Formula / concept | Spec | Code |
|---|---|---|
| partitions of `n`, canonical order `lex_desc_full_parts_v1` | §3, §19.3 | `partition::gen::partitions_in_canonical_order`, `partition::index::PartitionIndex` (order hash: `PartitionIndex::order_hash`) |
| `z_λ = ∏ i^{m_i} m_i!`, `\|C_λ\| = n!/z_λ`, `sgn(λ) = (−1)^{n−ℓ}` | §3 | `Partition::z_value`, `PartitionIndex::class_size`, `Partition::sign` |
| transpose partition `λ'` | §11 | `Partition::transpose`, `PartitionIndex::transpose_id` |
| cycle types w/o fixed points, padded per `n` | project input format | `partition::template::CycleTypeTemplate::pad` |
| hook lengths, `f_ρ = n!/∏h` | §22.1 | `chars::degrees::{hook_lengths, degree}` |
| rim hooks / border strips, sign `(−1)^{ht−1}` | §10.2 | `chars::rimhook::BetaSet::{for_each_hook_removal, for_each_hook_addition}` (abacus/beta-set form) |
| MN recursion; column = `M_{l₁}∘…∘M_{l_k}(e_∅)` | §10.2 | `chars::mn::MnEvaluator` (CSR level operators, gather orientation); naive reference: `testing::naive_mn::NaiveMn` |
| `ω_ρ(λ) = \|C_λ\|χ^ρ(λ)/f_ρ`, divisibility assert | §4, §25 | `spectra::BaseSpectra::build` (error `OmegaNotIntegral`) |
| `θ_ρ(U) = Σ_j ω_ρ(λ_j)` | §4.2, §25 | `spectra::BaseSpectra::theta` |
| active rows `R = ∪_j {ρ : χ^ρ(λ_j) ≠ 0}` (exact values only) | §7 | `spectra::BaseSpectra::active_rows` |
| `a_r(ν) = (1/n!)Σ_ρ f_ρ χ^ρ(ν) θ_ρ^r`, exact ÷, ≥ 0 | §4.1, §25, §9.4–9.5 | `engine::exact::ExactTransform::coefficients` |
| radius-0 identity reconstruction (full rows) | §9.1, §7 | asserted in `engine::exact::run_exact`; restricted rows start at r = 1 |
| radius-1 indicator | §9.2 | asserted in `run_exact` |
| word-count identity `Σ \|C_ν\| a_r(ν) = \|U\|^r` | §9.3 | checked every radius in `run_exact` (`WordCountMismatch`) |
| distances, `d_U(ν) = min{r : a_r > 0}`; per-parity first hits | §5 | `engine::UnionRun::{distance, first_hit_even, first_hit_odd}` |
| stopping: empty `new_r` on exact supports | §5.2 | `engine::exact::run_exact` (`StoppingRule::EmptyLayer`); parity-cover early exit = `AllTypesVisited`; normal-closure prediction is a cross-check assert only |
| coefficient bound `a_r(ν) ≤ ⌊\|U\|^r/\|C_ν\|⌋` | §9.3 corollary | `arith::bounds::coefficient_bound` (P2 certification tier 1/2) |
| screening primes `p > n`, `p < 2³¹`; overflow bound | §12–13 | `arith::modp::{screening_primes, ModCtx}`; compile-time bound in `arith` (`MAX_ACCUM_TERMS`) |
| brute-force oracles (BFS, set-product DP, word-count DP) | §9.7, §22.3 | `testing::bruteforce` |
| mixed products `K_A K_B` | §16.3, §22.4 | test `mixed_products_and_union_square` (crates/core/tests/bruteforce.rs) |
| result JSON, order embedding, metadata | §19.3 | `report::{build_result, schema}` + `docs/output_schema.md` |
