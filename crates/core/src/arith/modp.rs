//! Screening primes and modular contexts.
#![deny(clippy::float_arithmetic)]

use num_bigint::{BigInt, BigUint, Sign};
use serde::{Deserialize, Serialize};

/// A screening prime with `n < p < 2³¹`.
///
/// The prime sequence is deterministic (descending from `2³¹ − 1`, primality
/// by Miller–Rabin with bases proven deterministic for `u32`), so a prime
/// list recorded in a checkpoint or result file can be re-derived and
/// verified. `p > n` is automatic for all supported `n ≤ 255`.
#[derive(Copy, Clone, PartialEq, Eq, Debug, Hash, Serialize, Deserialize)]
#[serde(transparent)]
pub struct Prime31(pub u32);

/// The first `count` screening primes, descending from `2³¹ − 1`.
///
/// The first few: 2147483647, 2147483629, 2147483587, 2147483579, …
pub fn screening_primes(count: usize) -> Vec<Prime31> {
    let mut out = Vec::with_capacity(count);
    let mut candidate: u32 = (1 << 31) - 1;
    while out.len() < count {
        if is_prime_u32(candidate) {
            out.push(Prime31(candidate));
        }
        candidate = candidate
            .checked_sub(2)
            .expect("ran out of u32 primes (unreachable)");
    }
    out
}

/// Deterministic Miller–Rabin for u32: bases {2, 7, 61} are proven
/// deterministic for all n < 4 759 123 141 > 2³².
fn is_prime_u32(n: u32) -> bool {
    if n < 2 {
        return false;
    }
    for small in [2u32, 3, 5, 7, 11, 13] {
        if n == small {
            return true;
        }
        if n % small == 0 {
            return false;
        }
    }
    let n64 = n as u64;
    let mut d = n64 - 1;
    let mut s = 0u32;
    while d % 2 == 0 {
        d /= 2;
        s += 1;
    }
    'witness: for a in [2u64, 7, 61] {
        if a % n64 == 0 {
            continue;
        }
        let mut x = pow_mod_u64(a, d, n64);
        if x == 1 || x == n64 - 1 {
            continue;
        }
        for _ in 1..s {
            x = mul_mod_u64(x, x, n64);
            if x == n64 - 1 {
                continue 'witness;
            }
        }
        return false;
    }
    true
}

fn mul_mod_u64(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 * b as u128) % m as u128) as u64
}

fn pow_mod_u64(mut base: u64, mut exp: u64, m: u64) -> u64 {
    let mut acc = 1u64;
    base %= m;
    while exp > 0 {
        if exp & 1 == 1 {
            acc = mul_mod_u64(acc, base, m);
        }
        base = mul_mod_u64(base, base, m);
        exp >>= 1;
    }
    acc
}

/// Modular arithmetic context for one prime. All residues handled by this
/// context are fully reduced values in `[0, p)` (plain representation, not
/// Montgomery — see design doc `02-numerics` §3.2).
#[derive(Clone, Debug)]
pub struct ModCtx {
    p: u64,
}

impl ModCtx {
    pub fn new(p: Prime31) -> Self {
        debug_assert!(is_prime_u32(p.0), "ModCtx requires a prime");
        Self { p: p.0 as u64 }
    }

    pub fn prime(&self) -> Prime31 {
        Prime31(self.p as u32)
    }

    #[inline]
    pub fn add(&self, a: u32, b: u32) -> u32 {
        let s = a as u64 + b as u64;
        (if s >= self.p { s - self.p } else { s }) as u32
    }

    #[inline]
    pub fn sub(&self, a: u32, b: u32) -> u32 {
        let (a, b) = (a as u64, b as u64);
        (if a >= b { a - b } else { a + self.p - b }) as u32
    }

    #[inline]
    pub fn neg(&self, a: u32) -> u32 {
        if a == 0 {
            0
        } else {
            (self.p - a as u64) as u32
        }
    }

    #[inline]
    pub fn mul(&self, a: u32, b: u32) -> u32 {
        ((a as u64 * b as u64) % self.p) as u32
    }

    /// Reduce a `u128` accumulator (see [`crate::arith::MAX_ACCUM_TERMS`] for
    /// the proven bound on what may be accumulated).
    #[inline]
    pub fn reduce_u128(&self, acc: u128) -> u32 {
        (acc % self.p as u128) as u32
    }

    pub fn pow(&self, a: u32, e: u64) -> u32 {
        pow_mod_u64(a as u64, e, self.p) as u32
    }

    /// Multiplicative inverse by Fermat (`a != 0`).
    pub fn inv(&self, a: u32) -> u32 {
        assert!(a != 0, "inverse of zero");
        self.pow(a, self.p - 2)
    }

    /// Reduce an exact unsigned integer.
    pub fn reduce_biguint(&self, x: &BigUint) -> u32 {
        let r = x % BigUint::from(self.p);
        r.to_u32_digits().first().copied().unwrap_or(0)
    }

    /// Reduce an exact signed integer into `[0, p)`.
    pub fn reduce_bigint(&self, x: &BigInt) -> u32 {
        let m = self.reduce_biguint(x.magnitude());
        match x.sign() {
            Sign::Minus => self.neg(m),
            _ => m,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigInt;

    #[test]
    fn first_primes_are_the_documented_ones() {
        let ps = screening_primes(4);
        assert_eq!(
            ps,
            vec![
                Prime31(2147483647),
                Prime31(2147483629),
                Prime31(2147483587),
                Prime31(2147483579)
            ]
        );
    }

    #[test]
    fn primes_are_prime_and_descending() {
        let ps = screening_primes(20);
        for w in ps.windows(2) {
            assert!(w[0].0 > w[1].0);
        }
        // trial-division verification, independent of Miller-Rabin
        for p in ps {
            let mut d = 3u64;
            let n = p.0 as u64;
            assert!(n % 2 != 0);
            while d * d <= n {
                assert!(n % d != 0, "{n} divisible by {d}");
                d += 2;
            }
        }
    }

    #[test]
    fn miller_rabin_agrees_with_trial_division_small() {
        fn trial(n: u32) -> bool {
            if n < 2 {
                return false;
            }
            let mut d = 2u64;
            while d * d <= n as u64 {
                if n as u64 % d == 0 {
                    return false;
                }
                d += 1;
            }
            true
        }
        for n in 0..5000u32 {
            assert_eq!(is_prime_u32(n), trial(n), "n={n}");
        }
    }

    #[test]
    fn modctx_basics() {
        let ctx = ModCtx::new(Prime31(2147483647));
        let p = 2147483647u64;
        assert_eq!(ctx.add(p as u32 - 1, 5), 4);
        assert_eq!(ctx.sub(3, 8), (p - 5) as u32);
        assert_eq!(ctx.neg(0), 0);
        assert_eq!(ctx.mul(p as u32 - 1, p as u32 - 1), 1); // (-1)^2
        assert_eq!(ctx.pow(3, 0), 1);
        let a = 123456789u32;
        assert_eq!(ctx.mul(a, ctx.inv(a)), 1);
    }

    #[test]
    fn reduce_signed() {
        let ctx = ModCtx::new(Prime31(2147483647));
        assert_eq!(ctx.reduce_bigint(&BigInt::from(-1)), 2147483646);
        assert_eq!(ctx.reduce_bigint(&BigInt::from(0)), 0);
        let big = BigInt::from(2147483647u64) * 3 + 7;
        assert_eq!(ctx.reduce_bigint(&big), 7);
        assert_eq!(ctx.reduce_bigint(&(-big)), ctx.neg(7));
    }

    #[test]
    fn reduce_u128_matches_repeated_mul() {
        let ctx = ModCtx::new(Prime31(2147483629));
        let mut acc: u128 = 0;
        let mut expect: u32 = 0;
        for i in 0..1000u64 {
            let a = (i * 2654435761 % 2147483629) as u32;
            let b = (i * 40503 % 2147483629) as u32;
            acc += a as u128 * b as u128;
            expect = ctx.add(expect, ctx.mul(a, b));
        }
        assert_eq!(ctx.reduce_u128(acc), expect);
    }
}
