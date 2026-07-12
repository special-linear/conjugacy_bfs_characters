//! Irreducible degrees `f_ρ = χ^ρ(1)` via the hook-length formula.
//!
//! Independent of the Murnaghan–Nakayama evaluator; the identity column of
//! the MN table must reproduce these values (a permanent cross-check between
//! two formulas, spec §22.2).
#![deny(clippy::float_arithmetic)]

use num_bigint::BigUint;
use num_traits::{One, Zero};

use crate::partition::{factorial, Partition, PartitionId, PartitionIndex};

/// Hook lengths of all cells, row by row: `h(i,j) = λ_i − j + λ'_j − i − 1`
/// in 0-based coordinates.
pub fn hook_lengths(p: &Partition) -> Vec<u32> {
    let transpose = p.transpose();
    let parts = p.parts();
    let cols = transpose.parts();
    let mut hooks = Vec::with_capacity(p.n() as usize);
    for (i, &row_len) in parts.iter().enumerate() {
        for (j, &col_len) in cols[..row_len as usize].iter().enumerate() {
            let h = row_len as u32 - j as u32 + col_len as u32 - i as u32 - 1;
            hooks.push(h);
        }
    }
    hooks
}

/// `f_ρ = n! / ∏ hooks` — exact division, asserted.
pub fn degree(p: &Partition) -> BigUint {
    let mut hook_product = BigUint::one();
    for h in hook_lengths(p) {
        hook_product *= BigUint::from(h);
    }
    let n_factorial = factorial(p.n());
    let (q, r) = num_integer::Integer::div_rem(&n_factorial, &hook_product);
    assert!(r.is_zero(), "hook product must divide n! for {p:?}");
    q
}

/// Degrees for every irreducible, in canonical order.
pub fn degrees(index: &PartitionIndex) -> Vec<BigUint> {
    (0..index.count())
        .map(|i| degree(index.partition(i as PartitionId)))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::partition::PartitionIndex;

    #[test]
    fn known_degree_tables() {
        // n=3, canonical order [3],[2,1],[1^3]
        let idx = PartitionIndex::build(3).unwrap();
        let d: Vec<u32> = degrees(&idx)
            .iter()
            .map(|d| u32::try_from(d.clone()).unwrap())
            .collect();
        assert_eq!(d, vec![1, 2, 1]);

        // n=4: [4],[3,1],[2,2],[2,1,1],[1^4]
        let idx = PartitionIndex::build(4).unwrap();
        let d: Vec<u32> = degrees(&idx)
            .iter()
            .map(|d| u32::try_from(d.clone()).unwrap())
            .collect();
        assert_eq!(d, vec![1, 3, 2, 3, 1]);

        // n=5: [5],[4,1],[3,2],[3,1,1],[2,2,1],[2,1,1,1],[1^5]
        let idx = PartitionIndex::build(5).unwrap();
        let d: Vec<u32> = degrees(&idx)
            .iter()
            .map(|d| u32::try_from(d.clone()).unwrap())
            .collect();
        assert_eq!(d, vec![1, 4, 5, 6, 5, 4, 1]);

        // n=6: [6],[5,1],[4,2],[4,1,1],[3,3],[3,2,1],[3,1^3],[2,2,2],[2,2,1,1],[2,1^4],[1^6]
        let idx = PartitionIndex::build(6).unwrap();
        let d: Vec<u32> = degrees(&idx)
            .iter()
            .map(|d| u32::try_from(d.clone()).unwrap())
            .collect();
        assert_eq!(d, vec![1, 5, 9, 10, 5, 16, 10, 5, 9, 5, 1]);
    }

    #[test]
    fn degree_squares_sum_to_factorial() {
        for n in [1u16, 2, 5, 8, 12, 16, 20] {
            let idx = PartitionIndex::build(n).unwrap();
            let sum: BigUint = degrees(&idx).iter().map(|d| d * d).sum();
            assert_eq!(&sum, idx.factorial_n(), "n={n}");
        }
    }

    #[test]
    fn degree_transpose_symmetry() {
        for n in [4u16, 9, 15, 20] {
            let idx = PartitionIndex::build(n).unwrap();
            let d = degrees(&idx);
            for i in 0..idx.count() {
                let t = idx.transpose_id(i as PartitionId) as usize;
                assert_eq!(d[i], d[t], "n={n}, i={i}");
            }
        }
    }

    #[test]
    fn hook_lengths_example() {
        // λ = (3,2): hooks are [4,3,1 / 2,1]
        let mut h = hook_lengths(&Partition::new(vec![3u8, 2]));
        h.sort_unstable();
        assert_eq!(h, vec![1, 1, 2, 3, 4]);
        // degree = 5!/(4·3·1·2·1) = 120/24 = 5
        assert_eq!(degree(&Partition::new(vec![3u8, 2])), BigUint::from(5u32));
    }

    #[test]
    fn trivial_edge_cases() {
        assert_eq!(degree(&Partition::identity(0)), BigUint::one());
        assert_eq!(degree(&Partition::new(vec![7u8])), BigUint::one());
        assert_eq!(degree(&Partition::identity(7)), BigUint::one());
    }
}
