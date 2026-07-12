//! Beta-sets (abacus displays) and rim-hook enumeration (spec §10.2).
//!
//! A partition `λ = (λ_1 ≥ … ≥ λ_k)` displayed on `s ≥ k` slots has beta
//! numbers `b_i = λ_i + (s − i)` for `i = 1..=s` (with `λ_i = 0` beyond the
//! parts) — `s` distinct values. Removing an `ℓ`-rim-hook corresponds to
//! moving a bead `b → b − ℓ` onto a free position; the leg length of the
//! removed hook equals the number of beads strictly between `b − ℓ` and `b`,
//! and the Murnaghan–Nakayama sign is `(−1)^leg`. Both the removal count and
//! the sign are independent of the slot count `s` (standard abacus facts;
//! exercised by tests).
#![deny(clippy::float_arithmetic)]

use smallvec::SmallVec;

use crate::partition::Partition;

/// A beta-set on a fixed number of slots; values sorted ascending.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct BetaSet {
    /// Distinct values, ascending. `values.len()` is the slot count.
    values: SmallVec<[u16; 32]>,
}

impl BetaSet {
    /// Display `p` on `slots ≥ p.len()` slots.
    pub fn of(p: &Partition, slots: u16) -> Self {
        let k = p.len();
        assert!(slots as usize >= k, "not enough slots for {p:?}");
        let mut values: SmallVec<[u16; 32]> = SmallVec::with_capacity(slots as usize);
        // ascending: slot i (1-based from the top) has b_i = λ_i + s − i;
        // iterate i = s..1 so values come out ascending.
        for i in (1..=slots).rev() {
            let part = if (i as usize) <= k {
                p.parts()[i as usize - 1] as u16
            } else {
                0
            };
            values.push(part + slots - i);
        }
        debug_assert!(values.windows(2).all(|w| w[0] < w[1]));
        Self { values }
    }

    pub fn slots(&self) -> u16 {
        self.values.len() as u16
    }

    fn contains(&self, v: u16) -> bool {
        self.values.binary_search(&v).is_ok()
    }

    /// Number of values strictly less than `v`.
    fn count_below(&self, v: u16) -> usize {
        self.values.partition_point(|&x| x < v)
    }

    /// Recover the partition (zero parts stripped).
    pub fn to_partition(&self) -> Partition {
        let mut parts: SmallVec<[u8; 16]> = SmallVec::new();
        // ascending index j: λ_{s−j} = values[j] − j
        for (j, &v) in self.values.iter().enumerate() {
            let part = v as usize - j;
            if part > 0 {
                parts.push(u8::try_from(part).expect("part fits u8"));
            }
        }
        parts.reverse();
        Partition::new(parts)
    }

    /// Enumerate `ℓ`-rim-hook removals: calls `f(result, leg)` for every bead
    /// `b` with `b ≥ ℓ` and `b − ℓ` free. `leg` is the hook's leg length
    /// (rows spanned − 1); the MN sign is `(−1)^leg`.
    pub fn for_each_hook_removal(&self, l: u8, mut f: impl FnMut(BetaSet, u32)) {
        let l = l as u16;
        for (i, &b) in self.values.iter().enumerate() {
            if b < l || self.contains(b - l) {
                continue;
            }
            // beads strictly between b−ℓ and b: those below b (index i counts
            // them) minus those below or equal b−ℓ (= those below, since b−ℓ free)
            let leg = (i - self.count_below(b - l)) as u32;
            let mut values = self.values.clone();
            values.remove(i);
            let insert_at = values.partition_point(|&x| x < b - l);
            values.insert(insert_at, b - l);
            f(BetaSet { values }, leg);
        }
    }

    /// Enumerate `ℓ`-rim-hook additions (inverse of removal): calls
    /// `f(result, leg)` for every bead `b` with `b + ℓ` free.
    pub fn for_each_hook_addition(&self, l: u8, mut f: impl FnMut(BetaSet, u32)) {
        let l = l as u16;
        for (i, &b) in self.values.iter().enumerate() {
            let target = b + l;
            if self.contains(target) {
                continue;
            }
            // beads strictly between b and b+ℓ
            let leg = (self.count_below(target) - i - 1) as u32;
            let mut values = self.values.clone();
            values.remove(i);
            let insert_at = values.partition_point(|&x| x < target);
            values.insert(insert_at, target);
            f(BetaSet { values }, leg);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn part(parts: &[u8]) -> Partition {
        Partition::new(parts.to_vec())
    }

    #[test]
    fn roundtrip() {
        for parts in [vec![], vec![3u8], vec![4, 2, 1], vec![2, 2, 2], vec![1, 1, 1]] {
            let p = Partition::new(parts);
            for slots in [p.len() as u16, p.len() as u16 + 1, p.len() as u16 + 5, 10] {
                if (slots as usize) < p.len() {
                    continue;
                }
                assert_eq!(BetaSet::of(&p, slots).to_partition(), p, "{p:?} @ {slots}");
            }
        }
    }

    /// Removals: count and (multiset of) results are slot-independent.
    #[test]
    fn removals_slot_independent() {
        let p = part(&[4, 3, 1]);
        for l in 1..=8u8 {
            let mut base: Vec<(Partition, u32)> = Vec::new();
            BetaSet::of(&p, 3).for_each_hook_removal(l, |b, leg| {
                base.push((b.to_partition(), leg));
            });
            base.sort_by(|a, b| a.0.parts().cmp(b.0.parts()));
            for slots in [4u16, 6, 9] {
                let mut got: Vec<(Partition, u32)> = Vec::new();
                BetaSet::of(&p, slots).for_each_hook_removal(l, |b, leg| {
                    got.push((b.to_partition(), leg));
                });
                got.sort_by(|a, b| a.0.parts().cmp(b.0.parts()));
                assert_eq!(got, base, "l={l}, slots={slots}");
            }
        }
    }

    #[test]
    fn known_removals() {
        // (2), remove 2-hook: horizontal domino, leg 0, result empty.
        let mut got = Vec::new();
        BetaSet::of(&part(&[2]), 3).for_each_hook_removal(2, |b, leg| {
            got.push((b.to_partition(), leg));
        });
        assert_eq!(got, vec![(part(&[]), 0)]);

        // (1,1), remove 2-hook: vertical domino, leg 1.
        let mut got = Vec::new();
        BetaSet::of(&part(&[1, 1]), 2).for_each_hook_removal(2, |b, leg| {
            got.push((b.to_partition(), leg));
        });
        assert_eq!(got, vec![(part(&[]), 1)]);

        // (3,1): 3-hooks — the L-shaped border strip {(1,2),(1,3)}+... :
        // removals: beads of (3,1) @2 slots: {4,1}. l=3: 4->1 occupied; 1->-2 no.
        // Hmm: (3,1) has two 3-rim-hooks? Cells: row1: 3, row2: 1.
        // Border strips of size 3: {(1,1),(2,1)}∪? — strip must be connected
        // along the rim: {(2,1),(1,1),(1,2)} (leg 1) and {(1,2),(1,3)} is size 2 — no.
        // Actually with beads {4,1}: 4->1 blocked, 1 can't. So ZERO removals?
        // Check: removing a 3-strip from (3,1) leaves a partition of 1 = (1).
        // Possible strips: cells {(1,2),(1,3),(2,1)}? not connected. {(1,1),(1,2),(1,3)}
        // leaves (0,1) — not a partition shape. So indeed no valid 3-hook. ✓
        let mut count = 0;
        BetaSet::of(&part(&[3, 1]), 2).for_each_hook_removal(3, |_, _| count += 1);
        assert_eq!(count, 0);

        // (3,1): 4-hook = whole shape, ht = 2, leg 1.
        let mut got = Vec::new();
        BetaSet::of(&part(&[3, 1]), 2).for_each_hook_removal(4, |b, leg| {
            got.push((b.to_partition(), leg));
        });
        assert_eq!(got, vec![(part(&[]), 1)]);
    }

    #[test]
    fn addition_inverts_removal() {
        let p = part(&[4, 2, 2, 1]);
        let slots = 9u16;
        for l in 1..=9u8 {
            BetaSet::of(&p, slots).for_each_hook_removal(l, |smaller, leg| {
                // adding the same hook back must reproduce p with the same leg
                let mut found = false;
                smaller.for_each_hook_addition(l, |bigger, leg2| {
                    if bigger.to_partition() == p && leg2 == leg {
                        found = true;
                    }
                });
                assert!(found, "l={l}");
            });
        }
    }

    /// nnz formula from the numerics design (§2.1): the number of addable
    /// ℓ-hooks over ALL partitions μ of m equals ℓ·(p(m) + p(m−ℓ) + p(m−2ℓ) + …)
    /// — equivalently, total (source, target) pairs for M_ℓ : m → m+ℓ.
    #[test]
    fn nnz_formula() {
        use crate::partition::partitions_in_canonical_order;
        fn p_of(m: i32) -> u64 {
            if m < 0 {
                return 0;
            }
            partitions_in_canonical_order(m as u16).len() as u64
        }
        for m in 0..=9u16 {
            for l in 1..=6u8 {
                let slots = m + l as u16; // enough for any target of m+l
                let mut nnz = 0u64;
                for mu in partitions_in_canonical_order(m) {
                    BetaSet::of(&mu, slots).for_each_hook_addition(l, |_, _| nnz += 1);
                }
                let mut expected = 0u64;
                let mut j = 0i32;
                loop {
                    let remaining = m as i32 - j * (l as i32);
                    if remaining < 0 {
                        break;
                    }
                    expected += p_of(remaining);
                    j += 1;
                }
                expected *= l as u64;
                assert_eq!(nnz, expected, "m={m}, l={l}");
            }
        }
    }
}
