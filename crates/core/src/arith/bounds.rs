//! Rigorous coefficient bounds used by the certification tiers.
#![deny(clippy::float_arithmetic)]

use num_bigint::BigUint;

/// `B_r(ν) = ⌊|U|^r / |C_ν|⌋` — an upper bound on `a_r(ν)`, derived from the
/// word-count identity `Σ_ν |C_ν|·a_r(ν) = |U|^r` and nonnegativity
/// (spec §9.3; certification tier 1/2 in the plan).
///
/// `union_size_pow_r` is `|U|^r` (callers cache the power across targets).
pub fn coefficient_bound(union_size_pow_r: &BigUint, class_size: &BigUint) -> BigUint {
    union_size_pow_r / class_size
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_traits::Zero;

    #[test]
    fn bound_examples() {
        // |U| = 15, r = 1: bound for a class of size 45 is 0 => certified zero.
        let pow = BigUint::from(15u32);
        assert!(coefficient_bound(&pow, &BigUint::from(45u32)).is_zero());
        assert_eq!(
            coefficient_bound(&pow, &BigUint::from(15u32)),
            BigUint::from(1u32)
        );
        assert_eq!(
            coefficient_bound(&pow, &BigUint::from(1u32)),
            BigUint::from(15u32)
        );
    }
}
