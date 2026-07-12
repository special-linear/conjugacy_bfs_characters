//! # classdiam-core
//!
//! Exact computation of distances, layer structure, and diameters of
//! conjugacy-invariant Cayley graphs of symmetric groups via characters.
//!
//! Mathematical specification: `notes/character_method_cayley_diameters.md`
//! (referenced throughout as "spec §N"). Engineering design: `docs/design/`.
//!
//! ## Module map (formula → code)
//!
//! | Concept | Module |
//! |---|---|
//! | partitions, canonical order, `z_λ`, class sizes, signs, transpose | [`partition`] |
//! | modular contexts (Barrett), exact integers, coefficient bounds | [`arith`] |
//! | hook-length degrees `f_ρ`, Murnaghan–Nakayama `χ^ρ(ν)` | [`chars`] |
//! | central eigenvalues `ω_ρ(λ)`, union spectra `θ_ρ(U)`, active rows | [`spectra`] |
//! | the per-radius transform `N_r = Xᵀ·W` behind a backend trait | [`transform`] |
//! | radius loop, certification gate, stopping, layer records | [`engine`] |
//! | spec §9 invariants (validation layer) | [`validate`] |
//!
//! The math modules ban floating point entirely (spec §23 Failure 3); all
//! reachability decisions are exact integer/residue arithmetic.

pub mod arith;
pub mod chars;
pub mod engine;
pub mod error;
pub mod partition;
pub mod spectra;

#[cfg(any(test, feature = "test-utils"))]
pub mod testing;

pub use error::ClassdiamError;
