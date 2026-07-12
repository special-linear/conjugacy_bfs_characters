//! Partitions of `n` and their combinatorial data (spec §3).
//!
//! A partition is stored as a weakly decreasing list of parts, each ≥ 1, with
//! parts as `u8` (so `n ≤ 255`). The **canonical order** used everywhere in
//! this crate — for conjugacy classes, irreducible characters, distance
//! arrays, and every serialized artifact — is [`ORDER_CONVENTION`]
//! (`lex_desc_full_parts_v1`): full part lists compared lexicographically
//! *descending*, so `[n]` has index 0 and `[1,…,1]` (the identity cycle type)
//! comes last. The order is versioned by name and by a blake3 hash of its
//! explicit encoding; consumers must refuse artifacts whose order hash
//! differs (spec §19.3, Failure 7).
#![deny(clippy::float_arithmetic)]

mod gen;
mod index;
mod template;

pub use gen::{partitions_in_canonical_order, PartitionCountTable};
pub use index::{factorial, PartitionId, PartitionIndex, ORDER_CONVENTION};
pub use template::CycleTypeTemplate;

use std::cmp::Ordering;
use std::fmt;

use num_bigint::BigUint;
use num_traits::One;
use serde::{Deserialize, Serialize};
use smallvec::SmallVec;

/// A partition: weakly decreasing parts, each in `1..=255`, sum ≤ 255.
///
/// The empty partition (of `n = 0`) is valid and is the DP base case for the
/// character evaluator.
#[derive(Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(transparent)]
pub struct Partition {
    parts: SmallVec<[u8; 16]>,
}

impl Partition {
    /// Build from a weakly decreasing list of positive parts.
    ///
    /// Panics if parts are unsorted or contain 0 — construction sites inside
    /// the crate are all order-preserving; external inputs go through
    /// [`Partition::from_unsorted`] or [`CycleTypeTemplate`].
    pub fn new(parts: impl Into<SmallVec<[u8; 16]>>) -> Self {
        let parts = parts.into();
        assert!(
            parts.windows(2).all(|w| w[0] >= w[1]),
            "partition parts must be weakly decreasing: {parts:?}"
        );
        assert!(!parts.contains(&0), "partition parts must be positive");
        assert!(
            parts.iter().map(|&p| p as u32).sum::<u32>() <= 255,
            "partition sum exceeds 255"
        );
        Self { parts }
    }

    /// Build from parts in any order (0 parts rejected by panic as in `new`).
    pub fn from_unsorted(mut parts: Vec<u8>) -> Self {
        parts.sort_unstable_by(|a, b| b.cmp(a));
        Self::new(SmallVec::from_vec(parts))
    }

    /// The identity cycle type `(1^n)`.
    pub fn identity(n: u16) -> Self {
        assert!(n <= 255);
        Self {
            parts: std::iter::repeat_n(1u8, n as usize).collect(),
        }
    }

    /// The sum of parts (the `n` this partitions).
    pub fn n(&self) -> u16 {
        self.parts.iter().map(|&p| p as u16).sum()
    }

    /// Number of parts `ℓ(λ)`.
    pub fn len(&self) -> usize {
        self.parts.len()
    }

    pub fn is_empty(&self) -> bool {
        self.parts.is_empty()
    }

    pub fn parts(&self) -> &[u8] {
        &self.parts
    }

    /// `true` iff every part equals 1 (includes the empty partition of n=0).
    pub fn is_identity_type(&self) -> bool {
        self.parts.iter().all(|&p| p == 1)
    }

    /// Multiplicity `m_i` of part `i`.
    pub fn multiplicity(&self, part: u8) -> u32 {
        self.parts.iter().filter(|&&p| p == part).count() as u32
    }

    /// `(part, multiplicity)` pairs, parts strictly decreasing.
    pub fn multiplicities(&self) -> Vec<(u8, u32)> {
        let mut out: Vec<(u8, u32)> = Vec::new();
        for &p in self.parts.iter() {
            match out.last_mut() {
                Some((q, m)) if *q == p => *m += 1,
                _ => out.push((p, 1)),
            }
        }
        out
    }

    /// Sign of the conjugacy class: `sgn(λ) = (−1)^{n − ℓ(λ)}` (spec §3).
    pub fn sign(&self) -> i8 {
        if (self.n() as usize - self.len()) % 2 == 0 {
            1
        } else {
            -1
        }
    }

    /// `z_λ = ∏ i^{m_i} · m_i!` — the centralizer order (spec §3).
    pub fn z_value(&self) -> BigUint {
        let mut z = BigUint::one();
        for (part, mult) in self.multiplicities() {
            for _ in 0..mult {
                z *= BigUint::from(part);
            }
            for k in 2..=mult {
                z *= BigUint::from(k);
            }
        }
        z
    }

    /// The transpose (conjugate) partition `λ'`: `λ'_j = #{i : λ_i ≥ j}`.
    pub fn transpose(&self) -> Self {
        let mut parts: SmallVec<[u8; 16]> = SmallVec::new();
        if let Some(&first) = self.parts.first() {
            for j in 1..=first {
                let count = self.parts.iter().filter(|&&p| p >= j).count();
                parts.push(count as u8);
            }
        }
        Self { parts }
    }

    /// `true` iff `λ = λ'`.
    pub fn is_self_transpose(&self) -> bool {
        *self == self.transpose()
    }
}

impl fmt::Debug for Partition {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self.parts.as_slice())
    }
}

impl fmt::Display for Partition {
    /// Parts joined by `.` — the slug form used in file names (`docs/design/01` §9.1).
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut first = true;
        for p in self.parts.iter() {
            if !first {
                write!(f, ".")?;
            }
            write!(f, "{p}")?;
            first = false;
        }
        Ok(())
    }
}

/// Comparator realizing the canonical order `lex_desc_full_parts_v1`:
/// `Ordering::Less` means "comes earlier in canonical order", i.e. is
/// lexicographically **greater** as a full part list.
///
/// Only meaningful for partitions of the same `n`; for different sums it
/// still yields a total order (plain lex-desc) but nothing relies on that.
pub fn canonical_cmp(a: &Partition, b: &Partition) -> Ordering {
    b.parts().cmp(a.parts())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn multiplicities_and_z() {
        // λ = (3,3,2,1,1,1): z = 3^2·2! · 2^1·1! · 1^3·3! = 18·2·6 = 216
        let p = Partition::new(vec![3u8, 3, 2, 1, 1, 1]);
        assert_eq!(p.n(), 11);
        assert_eq!(p.multiplicities(), vec![(3, 2), (2, 1), (1, 3)]);
        assert_eq!(p.z_value(), BigUint::from(216u32));
    }

    #[test]
    fn z_hand_values() {
        // (n): z = n; (1^n): z = n!; (2,1^{n-2}) in S_4: 2·1·2! = 4
        assert_eq!(Partition::new(vec![5u8]).z_value(), BigUint::from(5u32));
        assert_eq!(Partition::identity(4).z_value(), BigUint::from(24u32));
        assert_eq!(
            Partition::new(vec![2u8, 1, 1]).z_value(),
            BigUint::from(4u32)
        );
        // (2,2) in S_4: 2^2·2! = 8; class size 24/8 = 3 ✓
        assert_eq!(Partition::new(vec![2u8, 2]).z_value(), BigUint::from(8u32));
        // (3,2) in S_5: 3·2 = 6; class size 120/6 = 20 ✓
        assert_eq!(Partition::new(vec![3u8, 2]).z_value(), BigUint::from(6u32));
    }

    #[test]
    fn signs() {
        assert_eq!(Partition::new(vec![2u8, 1, 1]).sign(), -1); // transposition: odd
        assert_eq!(Partition::new(vec![3u8, 1]).sign(), 1); // 3-cycle: even
        assert_eq!(Partition::new(vec![2u8, 2]).sign(), 1); // double transposition: even
        assert_eq!(Partition::identity(7).sign(), 1);
        assert_eq!(Partition::new(vec![6u8]).sign(), -1); // 6-cycle: odd
    }

    #[test]
    fn transpose_examples_and_involution() {
        let p = Partition::new(vec![4u8, 2, 1]);
        assert_eq!(p.transpose().parts(), &[3, 2, 1, 1]);
        assert_eq!(p.transpose().transpose(), p);
        assert!(Partition::new(vec![3u8, 2, 1]).is_self_transpose());
        assert!(!Partition::new(vec![3u8, 1]).is_self_transpose());
        assert_eq!(
            Partition::new(vec![5u8]).transpose(),
            Partition::identity(5)
        );
        // empty partition
        let e = Partition::new(SmallVec::<[u8; 16]>::new());
        assert_eq!(e.transpose(), e);
    }

    #[test]
    fn canonical_cmp_order() {
        let a = Partition::new(vec![4u8]);
        let b = Partition::new(vec![3u8, 1]);
        let c = Partition::new(vec![2u8, 2]);
        let d = Partition::new(vec![2u8, 1, 1]);
        let e = Partition::identity(4);
        let mut v = vec![e.clone(), c.clone(), a.clone(), d.clone(), b.clone()];
        v.sort_by(canonical_cmp);
        assert_eq!(v, vec![a, b, c, d, e]);
    }

    #[test]
    fn display_slug() {
        assert_eq!(Partition::new(vec![3u8, 2, 2]).to_string(), "3.2.2");
        assert_eq!(Partition::identity(0).to_string(), "");
    }
}
