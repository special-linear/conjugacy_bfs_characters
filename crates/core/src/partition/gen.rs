//! Generation of partitions in the canonical order, and partition counting.
#![deny(clippy::float_arithmetic)]

use smallvec::SmallVec;

use super::Partition;

/// `counts[m][l]` = number of partitions of `m` with every part ≤ `l`,
/// for `0 ≤ m, l ≤ n`. `counts[m][n]` is `p(m)`.
///
/// Values fit `u64` comfortably for all supported `n ≤ 255`
/// (`p(255) ≈ 1.3·10¹⁴`). Used for dense ranking at every level of the
/// character DP and for pre-sizing allocations.
pub struct PartitionCountTable {
    n: usize,
    counts: Vec<u64>, // (n+1) x (n+1), row-major by m
}

impl PartitionCountTable {
    pub fn build(n: u16) -> Self {
        let n = n as usize;
        let w = n + 1;
        let mut counts = vec![0u64; w * w];
        counts[..w].fill(1); // p(0 | parts ≤ l) = 1 (empty partition)
        for m in 1..=n {
            for l in 1..=n {
                // partitions of m with parts ≤ l: either no part equals... standard
                // recurrence: p(m | ≤l) = p(m | ≤l−1) + p(m−l | ≤l)
                let without = counts[m * w + (l - 1)];
                let with = if m >= l { counts[(m - l) * w + l] } else { 0 };
                counts[m * w + l] = without
                    .checked_add(with)
                    .expect("partition count overflow (unreachable for n <= 255)");
            }
        }
        Self { n, counts }
    }

    /// Number of partitions of `m` with all parts ≤ `l`.
    pub fn count_max_part(&self, m: u16, l: u8) -> u64 {
        let l = (l as usize).min(self.n);
        self.counts[m as usize * (self.n + 1) + l]
    }

    /// `p(m)` — the number of partitions of `m`.
    pub fn p(&self, m: u16) -> u64 {
        self.counts[m as usize * (self.n + 1) + self.n]
    }
}

/// All partitions of `n` in the canonical order `lex_desc_full_parts_v1`:
/// full part lists in lexicographically descending order — `[n]` first,
/// `[1,…,1]` last. Generated directly by the descending-parts recursion,
/// which emits exactly this order.
pub fn partitions_in_canonical_order(n: u16) -> Vec<Partition> {
    assert!(n <= 255, "n must be <= 255");
    let table = PartitionCountTable::build(n);
    let q = usize::try_from(table.p(n)).expect("p(n) exceeds usize");
    let mut out = Vec::with_capacity(q);
    let mut current: SmallVec<[u8; 16]> = SmallVec::new();
    emit(n, n.min(255) as u8, &mut current, &mut out);
    debug_assert_eq!(out.len(), q);
    out
}

fn emit(remaining: u16, max_part: u8, current: &mut SmallVec<[u8; 16]>, out: &mut Vec<Partition>) {
    if remaining == 0 {
        out.push(Partition::new(current.clone()));
        return;
    }
    let hi = (max_part as u16).min(remaining) as u8;
    for k in (1..=hi).rev() {
        current.push(k);
        emit(remaining - k as u16, k, current, out);
        current.pop();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::partition::canonical_cmp;

    /// OEIS A000041: p(0) .. p(50).
    const A000041: [u64; 51] = [
        1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42, 56, 77, 101, 135, 176, 231, 297, 385, 490, 627, 792,
        1002, 1255, 1575, 1958, 2436, 3010, 3718, 4565, 5604, 6842, 8349, 10143, 12310, 14883,
        17977, 21637, 26015, 31185, 37338, 44583, 53174, 63261, 75175, 89134, 105558, 124754,
        147273, 173525, 204226,
    ];

    #[test]
    fn partition_counts_match_oeis() {
        let table = PartitionCountTable::build(50);
        for (m, &expected) in A000041.iter().enumerate() {
            assert_eq!(table.p(m as u16), expected, "p({m})");
        }
    }

    #[test]
    fn generated_count_matches_table_small_n() {
        for n in 0..=30u16 {
            let parts = partitions_in_canonical_order(n);
            assert_eq!(parts.len() as u64, A000041[n as usize], "n={n}");
        }
    }

    #[test]
    fn order_is_lex_descending_and_endpoints_correct() {
        for n in 1..=12u16 {
            let parts = partitions_in_canonical_order(n);
            assert_eq!(parts[0].parts(), &[n as u8], "first is [n], n={n}");
            assert!(
                parts[parts.len() - 1].is_identity_type(),
                "last is identity, n={n}"
            );
            for w in parts.windows(2) {
                assert_eq!(
                    canonical_cmp(&w[0], &w[1]),
                    std::cmp::Ordering::Less,
                    "strictly ordered: {:?} before {:?}",
                    w[0],
                    w[1]
                );
            }
            // every element sums to n, sorted descending (Partition::new asserts)
            for p in &parts {
                assert_eq!(p.n(), n);
            }
        }
    }

    #[test]
    fn known_order_n5() {
        let parts = partitions_in_canonical_order(5);
        let expected: Vec<Vec<u8>> = vec![
            vec![5],
            vec![4, 1],
            vec![3, 2],
            vec![3, 1, 1],
            vec![2, 2, 1],
            vec![2, 1, 1, 1],
            vec![1, 1, 1, 1, 1],
        ];
        let got: Vec<Vec<u8>> = parts.iter().map(|p| p.parts().to_vec()).collect();
        assert_eq!(got, expected);
    }

    #[test]
    fn count_max_part_consistency() {
        let table = PartitionCountTable::build(20);
        // partitions of 6 with parts ≤ 2: [2,2,2],[2,2,1,1],[2,1,1,1,1],[1^6] = 4
        assert_eq!(table.count_max_part(6, 2), 4);
        // parts ≤ 1: exactly one
        for m in 0..=20u16 {
            assert_eq!(table.count_max_part(m, 1), 1);
        }
    }

    #[test]
    fn n_zero() {
        let parts = partitions_in_canonical_order(0);
        assert_eq!(parts.len(), 1);
        assert!(parts[0].is_empty());
    }
}
