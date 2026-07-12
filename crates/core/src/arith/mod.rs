//! Arithmetic layer: modular residues (screening) and exact integers
//! (certification, spectra, validation).
//!
//! Regime `p31-u128-accumulate-v1` (spec §13.1, critique-confirmed): screening
//! primes are `< 2³¹`, so each product of two residues is `< 2⁶²`; dot
//! products accumulate in `u128`, which is provably overflow-free for any row
//! length up to `2⁶⁶` terms — vastly beyond `p(50) = 204 226 < 2¹⁸`. The
//! compile-time assertion below pins the bound; changing prime width or
//! accumulator type without re-proving it fails compilation (Failure 8).
#![deny(clippy::float_arithmetic)]

mod bounds;
mod exact;
mod modp;

pub use bounds::coefficient_bound;
pub use exact::{exact_div_checked, ExactInt};
pub use modp::{screening_primes, ModCtx, Prime31};

/// Upper bound on transform row length (`R ≤ p(n)`; `p(50) = 204 226`).
/// Kernels may accumulate at most this many `u62` products into a `u128`
/// between reductions.
pub const MAX_ACCUM_TERMS: u128 = 1 << 18;

// Proof of the accumulation bound, checked at compile time:
// terms * (p-1)^2 <= 2^18 * (2^31 - 1)^2 < 2^80 << 2^128.
const _: () = {
    let worst = MAX_ACCUM_TERMS * (((1u128 << 31) - 1) * ((1u128 << 31) - 1));
    assert!(
        worst < (1u128 << 80),
        "p31-u128-accumulate-v1 bound violated"
    );
};
