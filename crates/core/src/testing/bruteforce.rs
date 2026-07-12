//! Brute-force oracles over raw permutations (spec §9.7, §22.3; design doc
//! 03 §3): plain BFS, exact-length set-product DP, and group-algebra word
//! counts. Everything here enumerates all of `S_n` — intended for `n ≤ 8`
//! in regular runs, `n = 9, 10` behind `#[ignore]`.
//!
//! Composition convention (internal only; distances by cycle type are
//! convention-independent because classes are inverse-closed):
//! `compose(a, b)[x] = a[b[x]]` — apply `b` first.

use std::collections::VecDeque;

use fixedbitset::FixedBitSet;

use crate::partition::{Partition, PartitionId, PartitionIndex};
use crate::spectra::ResolvedUnion;

/// A permutation of `0..n` as its image list.
pub type Perm = Vec<u8>;

pub fn compose(a: &Perm, b: &Perm) -> Perm {
    b.iter().map(|&x| a[x as usize]).collect()
}

pub fn identity_perm(n: u16) -> Perm {
    (0..n as u8).collect()
}

/// Cycle type of a permutation.
pub fn cycle_type(perm: &Perm) -> Partition {
    let n = perm.len();
    let mut seen = vec![false; n];
    let mut parts: Vec<u8> = Vec::new();
    for start in 0..n {
        if seen[start] {
            continue;
        }
        let mut len = 0u8;
        let mut x = start;
        while !seen[x] {
            seen[x] = true;
            x = perm[x] as usize;
            len += 1;
        }
        parts.push(len);
    }
    Partition::from_unsorted(parts)
}

/// Factorials `0! ..= n!` as `u64` (`n ≤ 20`).
fn factorials(n: u16) -> Vec<u64> {
    let mut f = vec![1u64; n as usize + 1];
    for k in 1..=n as usize {
        f[k] = f[k - 1] * k as u64;
    }
    f
}

/// Lehmer rank of a permutation (0 = identity ... n!−1).
pub fn rank(perm: &Perm, facts: &[u64]) -> u32 {
    let n = perm.len();
    let mut r = 0u64;
    for i in 0..n {
        let smaller_after = perm[i + 1..].iter().filter(|&&x| x < perm[i]).count() as u64;
        r += smaller_after * facts[n - 1 - i];
    }
    u32::try_from(r).expect("rank fits u32 for n <= 12")
}

/// Inverse of [`rank`].
pub fn unrank(mut r: u64, n: u16, facts: &[u64]) -> Perm {
    let n = n as usize;
    let mut available: Vec<u8> = (0..n as u8).collect();
    let mut perm = Vec::with_capacity(n);
    for i in 0..n {
        let f = facts[n - 1 - i];
        let idx = (r / f) as usize;
        r %= f;
        perm.push(available.remove(idx));
    }
    perm
}

/// All permutations of the union's classes, as image lists.
pub fn materialize_union(index: &PartitionIndex, union: &ResolvedUnion) -> Vec<Perm> {
    let n = index.n();
    let facts = factorials(n);
    let total = facts[n as usize];
    let mut out = Vec::new();
    for r in 0..total {
        let p = unrank(r, n, &facts);
        let t = cycle_type(&p);
        let id = index.id_of(&t).expect("cycle type is a partition of n");
        if union.class_ids.contains(&id) {
            out.push(p);
        }
    }
    out
}

/// Plain BFS over `S_n`; returns the distance per permutation rank
/// (−1 unreachable).
pub fn bfs_distances(n: u16, generators: &[Perm]) -> Vec<i32> {
    let facts = factorials(n);
    let total = facts[n as usize] as usize;
    let mut dist = vec![-1i32; total];
    let id = identity_perm(n);
    let id_rank = rank(&id, &facts) as usize;
    dist[id_rank] = 0;
    let mut queue = VecDeque::new();
    queue.push_back(id);
    while let Some(v) = queue.pop_front() {
        let dv = dist[rank(&v, &facts) as usize];
        for g in generators {
            let w = compose(&v, g);
            let wr = rank(&w, &facts) as usize;
            if dist[wr] < 0 {
                dist[wr] = dv + 1;
                queue.push_back(w);
            }
        }
    }
    dist
}

/// Aggregate per-permutation distances by cycle type, asserting the distance
/// is constant on every conjugacy class (conjugacy invariance — its failure
/// would mean the generating set is not a union of full classes, spec §17.3).
pub fn distances_by_type(index: &PartitionIndex, dist: &[i32]) -> Vec<i32> {
    let n = index.n();
    let facts = factorials(n);
    let mut by_type = vec![i32::MIN; index.count()];
    for (r, &d) in dist.iter().enumerate() {
        let t = cycle_type(&unrank(r as u64, n, &facts));
        let id = index.id_of(&t).unwrap() as usize;
        if by_type[id] == i32::MIN {
            by_type[id] = d;
        } else {
            assert_eq!(
                by_type[id], d,
                "BFS distance not constant on class {t:?} — conjugacy invariance broken"
            );
        }
    }
    by_type.iter_mut().for_each(|d| {
        if *d == i32::MIN {
            *d = -1;
        }
    });
    by_type
}

/// Exact-length supports by set-product DP: `S_0 = {id}`, `S_{r+1} = S_r·U`.
/// Returns, for each `r = 0..=max_radius`, the sorted set of cycle-type ids
/// with at least one element in `S_r` (spec §5.1: exact-length supports are
/// not BFS frontiers).
pub fn exact_length_supports(
    index: &PartitionIndex,
    generators: &[Perm],
    max_radius: u32,
) -> Vec<Vec<PartitionId>> {
    let n = index.n();
    let facts = factorials(n);
    let total = facts[n as usize] as usize;
    let mut current = FixedBitSet::with_capacity(total);
    current.insert(rank(&identity_perm(n), &facts) as usize);
    let mut out = Vec::with_capacity(max_radius as usize + 1);
    for _ in 0..=max_radius {
        // record support types of `current`
        let mut types = FixedBitSet::with_capacity(index.count());
        for v in current.ones() {
            let t = cycle_type(&unrank(v as u64, n, &facts));
            types.insert(index.id_of(&t).unwrap() as usize);
        }
        out.push(types.ones().map(|t| t as PartitionId).collect());
        // advance
        let mut next = FixedBitSet::with_capacity(total);
        for v in current.ones() {
            let vp = unrank(v as u64, n, &facts);
            for g in generators {
                next.insert(rank(&compose(&vp, g), &facts) as usize);
            }
        }
        current = next;
    }
    out
}

/// Group-algebra word counts: `counts[r][w] = #(g_1,…,g_r) ∈ U^r with
/// g_1⋯g_r = w`, for `r = 0..=max_radius`. `u64` with checked adds.
pub fn word_counts(n: u16, generators: &[Perm], max_radius: u32) -> Vec<Vec<u64>> {
    let facts = factorials(n);
    let total = facts[n as usize] as usize;
    let mut current = vec![0u64; total];
    current[rank(&identity_perm(n), &facts) as usize] = 1;
    let mut out = vec![current.clone()];
    for _ in 0..max_radius {
        let mut next = vec![0u64; total];
        for (v, &c) in current.iter().enumerate() {
            if c == 0 {
                continue;
            }
            let vp = unrank(v as u64, n, &facts);
            for g in generators {
                let w = rank(&compose(&vp, g), &facts) as usize;
                next[w] = next[w].checked_add(c).expect("word count overflow");
            }
        }
        out.push(next.clone());
        current = next;
    }
    out
}

/// Word counts aggregated by cycle type: asserts the count is constant on
/// every class and returns `counts_by_type[r][type_id]` — the per-element
/// count `a_r(ν)` (spec §23 Failure 2: NOT multiplied by the class size).
pub fn word_counts_by_type(index: &PartitionIndex, counts: &[Vec<u64>]) -> Vec<Vec<u64>> {
    let n = index.n();
    let facts = factorials(n);
    let mut out = Vec::with_capacity(counts.len());
    for row in counts {
        let mut by_type = vec![u64::MAX; index.count()];
        for (v, &c) in row.iter().enumerate() {
            let t = cycle_type(&unrank(v as u64, n, &facts));
            let id = index.id_of(&t).unwrap() as usize;
            if by_type[id] == u64::MAX {
                by_type[id] = c;
            } else {
                assert_eq!(by_type[id], c, "word count not constant on class {t:?}");
            }
        }
        out.push(by_type);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rank_unrank_roundtrip() {
        for n in [1u16, 3, 5, 6] {
            let facts = factorials(n);
            let total = facts[n as usize];
            for r in 0..total {
                let p = unrank(r, n, &facts);
                assert_eq!(rank(&p, &facts) as u64, r, "n={n}, r={r}");
            }
        }
    }

    #[test]
    fn cycle_types() {
        assert_eq!(cycle_type(&vec![0, 1, 2]), Partition::identity(3));
        assert_eq!(cycle_type(&vec![1, 0, 2]).parts(), &[2, 1]);
        assert_eq!(cycle_type(&vec![1, 2, 0]).parts(), &[3]);
        assert_eq!(cycle_type(&vec![1, 0, 3, 2]).parts(), &[2, 2]);
    }

    #[test]
    fn composition_convention() {
        // a = (0 1), b = (1 2): compose(a,b) applies b first: 0->0->1? No:
        // compose(a,b)[x] = a[b[x]]: x=0: b[0]=0, a[0]=1 => 1; x=1: b[1]=2, a[2]=2;
        // x=2: b[2]=1, a[1]=0. Result [1,2,0] = 3-cycle.
        let a = vec![1u8, 0, 2];
        let b = vec![0u8, 2, 1];
        assert_eq!(compose(&a, &b), vec![1, 2, 0]);
    }

    #[test]
    fn class_sizes_via_materialization() {
        use crate::partition::CycleTypeTemplate;
        use crate::spectra::resolve_union;
        let index = PartitionIndex::build(6).unwrap();
        for (template, expected) in [("2", 15usize), ("3", 40), ("2,2", 45), ("6", 120)] {
            let t: CycleTypeTemplate = template.parse().unwrap();
            let union = resolve_union(&index, &[t], false).unwrap();
            let gens = materialize_union(&index, &union);
            assert_eq!(gens.len(), expected, "{template}");
        }
    }
}
