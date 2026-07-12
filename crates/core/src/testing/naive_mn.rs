//! Naive memoized Murnaghan–Nakayama evaluator (reference oracle).
//!
//! Direct recursion on the definition: remove the largest cycle of `ν` as a
//! rim hook of `ρ` in every possible way, recurse on the rest. Memoized on
//! `(ρ, remaining cycles)` exactly as spec §10.2 warns is necessary. No
//! shared level operators, no cleverness — this is the code the production
//! evaluator must always agree with.

use std::collections::HashMap;

use num_traits::{One, Zero};

use crate::arith::ExactInt;
use crate::chars::rimhook::BetaSet;
use crate::partition::Partition;

#[derive(Default)]
pub struct NaiveMn {
    memo: HashMap<(Vec<u8>, Vec<u8>), ExactInt>,
}

impl NaiveMn {
    pub fn new() -> Self {
        Self::default()
    }

    /// `χ^ρ(ν)` for `ρ, ν ⊢ n`.
    pub fn chi(&mut self, rho: &Partition, nu: &Partition) -> ExactInt {
        assert_eq!(rho.n(), nu.n(), "rho and nu must partition the same n");
        self.chi_inner(rho, nu.parts())
    }

    /// `cycles` is weakly decreasing (a suffix view of ν's parts).
    fn chi_inner(&mut self, rho: &Partition, cycles: &[u8]) -> ExactInt {
        if cycles.is_empty() {
            debug_assert_eq!(rho.n(), 0);
            return ExactInt::one();
        }
        let key = (rho.parts().to_vec(), cycles.to_vec());
        if let Some(v) = self.memo.get(&key) {
            return v.clone();
        }
        let l = cycles[0]; // largest remaining cycle
        let rest = &cycles[1..];
        let mut sum = ExactInt::zero();
        let slots = rho.len().max(1) as u16;
        let mut removals: Vec<(Partition, u32)> = Vec::new();
        BetaSet::of(rho, slots).for_each_hook_removal(l, |beta, leg| {
            removals.push((beta.to_partition(), leg));
        });
        for (sub, leg) in removals {
            let term = self.chi_inner(&sub, rest);
            if leg % 2 == 0 {
                sum += term;
            } else {
                sum -= term;
            }
        }
        self.memo.insert(key, sum.clone());
        sum
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chars::MnEvaluator;
    use crate::partition::{PartitionId, PartitionIndex};

    #[test]
    fn matches_hand_values() {
        let mut mn = NaiveMn::new();
        // S_3: χ^{(2,1)}([3]) = -1, χ^{(2,1)}([2,1]) = 0, χ^{(2,1)}([1^3]) = 2
        let rho = Partition::new(vec![2u8, 1]);
        assert_eq!(mn.chi(&rho, &Partition::new(vec![3u8])), ExactInt::from(-1));
        assert_eq!(mn.chi(&rho, &Partition::new(vec![2u8, 1])), ExactInt::from(0));
        assert_eq!(mn.chi(&rho, &Partition::identity(3)), ExactInt::from(2));
        // S_4: χ^{(2,2)}: degree 2; on [2,1,1]: 0; on [2,2]: 2; on [3,1]: -1; on [4]: 0
        let rho = Partition::new(vec![2u8, 2]);
        assert_eq!(mn.chi(&rho, &Partition::identity(4)), ExactInt::from(2));
        assert_eq!(
            mn.chi(&rho, &Partition::new(vec![2u8, 1, 1])),
            ExactInt::from(0)
        );
        assert_eq!(mn.chi(&rho, &Partition::new(vec![2u8, 2])), ExactInt::from(2));
        assert_eq!(mn.chi(&rho, &Partition::new(vec![3u8, 1])), ExactInt::from(-1));
        assert_eq!(mn.chi(&rho, &Partition::new(vec![4u8])), ExactInt::from(0));
    }

    /// The production trie-DP evaluator must agree with this oracle on full
    /// tables (the permanent differential test, design doc 03 §2.3).
    #[test]
    fn production_evaluator_matches_naive_full_tables() {
        for n in 1..=10u16 {
            let idx = PartitionIndex::build(n).unwrap();
            let ev = MnEvaluator::new(n);
            let mut naive = NaiveMn::new();
            for nu_id in 0..idx.count() {
                let nu = idx.partition(nu_id as PartitionId);
                let col = ev.column_exact(nu);
                for rho_id in 0..idx.count() {
                    let rho = idx.partition(rho_id as PartitionId);
                    assert_eq!(
                        col[rho_id],
                        naive.chi(rho, nu),
                        "n={n}, rho={rho:?}, nu={nu:?}"
                    );
                }
            }
        }
    }
}
