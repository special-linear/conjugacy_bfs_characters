//! `PartitionIndex`: the canonical, versioned index of all partitions of `n`.
#![deny(clippy::float_arithmetic)]

use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::{One, Zero};

use super::{canonical_cmp, gen::partitions_in_canonical_order, Partition};
use crate::error::ClassdiamError;

/// Dense index of a partition within the canonical order for a fixed `n`.
pub type PartitionId = u32;

/// Name of the canonical partition order. Recorded in every serialized
/// artifact together with [`PartitionIndex::order_hash`]; the pair versions
/// the convention (spec §19.3).
pub const ORDER_CONVENTION: &str = "lex_desc_full_parts_v1";

/// Everything indexed by the canonical order for one `n`: the partitions
/// themselves, signs, `z_λ`, class sizes, and the transpose involution.
///
/// Built once per `n`; all hot paths use dense `PartitionId`s into the
/// parallel arrays here. No hash maps: lookup is binary search in the
/// canonical order.
pub struct PartitionIndex {
    n: u16,
    partitions: Vec<Partition>,
    signs: Vec<i8>,
    z_values: Vec<BigUint>,
    class_sizes: Vec<BigUint>,
    transpose_map: Vec<PartitionId>,
    order_hash: [u8; 32],
    factorial_n: BigUint,
}

impl PartitionIndex {
    pub fn build(n: u16) -> Result<Self, ClassdiamError> {
        if n > 255 {
            return Err(ClassdiamError::UnsupportedN { n: n as u32 });
        }
        let partitions = partitions_in_canonical_order(n);
        if u32::try_from(partitions.len()).is_err() {
            return Err(ClassdiamError::UnsupportedN { n: n as u32 });
        }

        let factorial_n = factorial(n);
        let mut signs = Vec::with_capacity(partitions.len());
        let mut z_values = Vec::with_capacity(partitions.len());
        let mut class_sizes = Vec::with_capacity(partitions.len());
        for p in &partitions {
            signs.push(p.sign());
            let z = p.z_value();
            let (size, rem) = factorial_n.div_rem(&z);
            debug_assert!(rem.is_zero(), "z_lambda must divide n!");
            z_values.push(z);
            class_sizes.push(size);
        }

        let order_hash = hash_order(n, &partitions);

        let mut index = Self {
            n,
            partitions,
            signs,
            z_values,
            class_sizes,
            transpose_map: Vec::new(),
            order_hash,
            factorial_n,
        };
        index.transpose_map = index
            .partitions
            .iter()
            .map(|p| {
                index
                    .id_of(&p.transpose())
                    .expect("transpose of a partition of n is a partition of n")
            })
            .collect();
        Ok(index)
    }

    pub fn n(&self) -> u16 {
        self.n
    }

    /// `q = p(n)`, the number of partitions / conjugacy classes / irreducibles.
    pub fn count(&self) -> usize {
        self.partitions.len()
    }

    pub fn partition(&self, id: PartitionId) -> &Partition {
        &self.partitions[id as usize]
    }

    pub fn partitions(&self) -> &[Partition] {
        &self.partitions
    }

    /// Canonical index of `p`, if `p` is a partition of this `n`.
    pub fn id_of(&self, p: &Partition) -> Option<PartitionId> {
        self.partitions
            .binary_search_by(|probe| canonical_cmp(probe, p))
            .ok()
            .map(|i| i as PartitionId)
    }

    /// The identity cycle type `(1^n)` — always the last index.
    pub fn identity_id(&self) -> PartitionId {
        (self.partitions.len() - 1) as PartitionId
    }

    pub fn sign(&self, id: PartitionId) -> i8 {
        self.signs[id as usize]
    }

    pub fn z_value(&self, id: PartitionId) -> &BigUint {
        &self.z_values[id as usize]
    }

    /// `|C_λ| = n! / z_λ`.
    pub fn class_size(&self, id: PartitionId) -> &BigUint {
        &self.class_sizes[id as usize]
    }

    /// Index of the transpose (conjugate) partition.
    pub fn transpose_id(&self, id: PartitionId) -> PartitionId {
        self.transpose_map[id as usize]
    }

    pub fn factorial_n(&self) -> &BigUint {
        &self.factorial_n
    }

    /// blake3 hash of the explicit order encoding:
    /// `LE(n as u16) ‖ (len(parts) as u8 ‖ parts…) per partition in order`.
    pub fn order_hash(&self) -> &[u8; 32] {
        &self.order_hash
    }

    pub fn order_hash_hex(&self) -> String {
        self.order_hash.iter().map(|b| format!("{b:02x}")).collect()
    }
}

fn hash_order(n: u16, partitions: &[Partition]) -> [u8; 32] {
    let mut hasher = blake3::Hasher::new();
    hasher.update(&n.to_le_bytes());
    for p in partitions {
        hasher.update(&[p.len() as u8]);
        hasher.update(p.parts());
    }
    *hasher.finalize().as_bytes()
}

/// `n!` as an exact integer.
pub fn factorial(n: u16) -> BigUint {
    let mut f = BigUint::one();
    for k in 2..=u64::from(n) {
        f *= BigUint::from(k);
    }
    f
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn id_roundtrip_all_partitions() {
        for n in [0u16, 1, 4, 7, 12, 20] {
            let idx = PartitionIndex::build(n).unwrap();
            for i in 0..idx.count() {
                let id = i as PartitionId;
                assert_eq!(idx.id_of(idx.partition(id)), Some(id), "n={n}, i={i}");
            }
            assert!(idx.partition(idx.identity_id()).is_identity_type());
        }
    }

    #[test]
    fn class_sizes_sum_to_factorial() {
        for n in [1u16, 2, 5, 9, 14, 20, 30] {
            let idx = PartitionIndex::build(n).unwrap();
            let total: BigUint = (0..idx.count())
                .map(|i| idx.class_size(i as PartitionId).clone())
                .sum();
            assert_eq!(&total, idx.factorial_n(), "n={n}");
        }
    }

    #[test]
    fn even_class_sizes_sum_to_half_factorial() {
        for n in [2u16, 3, 6, 10, 15, 25] {
            let idx = PartitionIndex::build(n).unwrap();
            let even_total: BigUint = (0..idx.count())
                .filter(|&i| idx.sign(i as PartitionId) == 1)
                .map(|i| idx.class_size(i as PartitionId).clone())
                .sum();
            assert_eq!(
                even_total * BigUint::from(2u32),
                idx.factorial_n().clone(),
                "n={n}"
            );
        }
    }

    #[test]
    fn z_times_class_size_is_factorial() {
        let idx = PartitionIndex::build(11).unwrap();
        for i in 0..idx.count() {
            let id = i as PartitionId;
            assert_eq!(&(idx.z_value(id) * idx.class_size(id)), idx.factorial_n());
        }
    }

    #[test]
    fn identity_class_has_size_one() {
        let idx = PartitionIndex::build(9).unwrap();
        assert_eq!(idx.class_size(idx.identity_id()), &BigUint::one());
        assert_eq!(idx.sign(idx.identity_id()), 1);
    }

    #[test]
    fn transpose_map_is_involution() {
        for n in [3u16, 8, 13, 20] {
            let idx = PartitionIndex::build(n).unwrap();
            let mut self_transpose = 0usize;
            for i in 0..idx.count() {
                let id = i as PartitionId;
                let t = idx.transpose_id(id);
                assert_eq!(idx.transpose_id(t), id);
                if t == id {
                    self_transpose += 1;
                }
            }
            assert_eq!(
                self_transpose,
                count_distinct_odd_part_partitions(n),
                "self-transpose count = #partitions into distinct odd parts, n={n}"
            );
        }
    }

    /// Independent DP: number of partitions of n into DISTINCT ODD parts.
    fn count_distinct_odd_part_partitions(n: u16) -> usize {
        let n = n as usize;
        let mut ways = vec![0u64; n + 1];
        ways[0] = 1;
        let mut part = 1usize;
        while part <= n {
            for m in (part..=n).rev() {
                ways[m] += ways[m - part];
            }
            part += 2;
        }
        ways[n] as usize
    }

    #[test]
    fn factorials() {
        assert_eq!(factorial(0), BigUint::one());
        assert_eq!(factorial(1), BigUint::one());
        assert_eq!(factorial(6), BigUint::from(720u32));
        assert_eq!(factorial(10), BigUint::from(3628800u64));
    }

    #[test]
    fn order_hash_matches_documented_encoding_and_distinguishes_n() {
        let idx6 = PartitionIndex::build(6).unwrap();
        let idx7 = PartitionIndex::build(7).unwrap();
        assert_ne!(idx6.order_hash(), idx7.order_hash());

        // Re-derive the n=6 hash from the documented encoding, independently.
        let mut hasher = blake3::Hasher::new();
        hasher.update(&6u16.to_le_bytes());
        for p in idx6.partitions() {
            hasher.update(&[p.len() as u8]);
            hasher.update(p.parts());
        }
        assert_eq!(idx6.order_hash(), hasher.finalize().as_bytes());
        assert_eq!(idx6.order_hash_hex().len(), 64);
    }

    #[test]
    fn known_index_n6() {
        // Canonical order for n=6 (lex desc), as used in the design docs' worked example.
        let idx = PartitionIndex::build(6).unwrap();
        let expected: Vec<Vec<u8>> = vec![
            vec![6],
            vec![5, 1],
            vec![4, 2],
            vec![4, 1, 1],
            vec![3, 3],
            vec![3, 2, 1],
            vec![3, 1, 1, 1],
            vec![2, 2, 2],
            vec![2, 2, 1, 1],
            vec![2, 1, 1, 1, 1],
            vec![1, 1, 1, 1, 1, 1],
        ];
        let got: Vec<Vec<u8>> = idx
            .partitions()
            .iter()
            .map(|p| p.parts().to_vec())
            .collect();
        assert_eq!(got, expected);
        // class sizes for n=6 in this order (see design doc worked example)
        let sizes: Vec<u32> = (0..idx.count())
            .map(|i| u32::try_from(idx.class_size(i as PartitionId).clone()).unwrap())
            .collect();
        assert_eq!(sizes, vec![120, 144, 90, 90, 40, 120, 40, 15, 45, 15, 1]);
        let signs: Vec<i8> = (0..idx.count())
            .map(|i| idx.sign(i as PartitionId))
            .collect();
        assert_eq!(signs, vec![-1, 1, 1, -1, 1, -1, 1, -1, 1, -1, 1]);
    }
}
