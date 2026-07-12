//! Exact (arbitrary-precision) integers.
//!
//! This is the ONLY module that names a bigint crate; everything else uses
//! the [`ExactInt`] alias so the backend can be swapped by feature flag
//! later (design decision: `num-bigint` default — pure Rust, MIT/Apache,
//! Windows-clean; `malachite` is LGPL and deliberately not the default).
#![deny(clippy::float_arithmetic)]

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::Zero;

/// Exact signed integer used off the hot path: degrees, class sizes,
/// eigenvalues, certification, validation.
pub type ExactInt = BigInt;

/// Exact division that reports a non-zero remainder instead of truncating
/// (spec §9.4: never silently round).
pub fn exact_div_checked(numerator: &ExactInt, denominator: &ExactInt) -> Option<ExactInt> {
    let (q, r) = numerator.div_rem(denominator);
    if r.is_zero() {
        Some(q)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_division() {
        let a = ExactInt::from(720);
        assert_eq!(
            exact_div_checked(&a, &ExactInt::from(6)),
            Some(ExactInt::from(120))
        );
        assert_eq!(exact_div_checked(&a, &ExactInt::from(7)), None);
        assert_eq!(
            exact_div_checked(&ExactInt::from(-720), &ExactInt::from(6)),
            Some(ExactInt::from(-120))
        );
    }
}
