//! Deterministic union catalog for parametrized invariant and oracle tests
//! (design doc 03 §2.5).

use num_bigint::BigUint;

use crate::partition::{CycleTypeTemplate, Partition, PartitionId, PartitionIndex};
use crate::spectra::{resolve_union, ResolvedUnion};

pub struct CatalogEntry {
    pub n: u16,
    pub templates: Vec<CycleTypeTemplate>,
    pub label: String,
}

fn template_of(p: &Partition) -> CycleTypeTemplate {
    let non_fixed: Vec<u8> = p.parts().iter().copied().filter(|&x| x >= 2).collect();
    CycleTypeTemplate::new(&non_fixed).expect("parts >= 2")
}

/// Tiny deterministic PRNG (xorshift64*) — no external dependency, fixed
/// seeds in the catalog, reproducible everywhere.
struct XorShift(u64);

impl XorShift {
    fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x.wrapping_mul(0x2545F4914F6CDD1D)
    }

    fn below(&mut self, bound: usize) -> usize {
        (self.next() % bound as u64) as usize
    }
}

/// The full catalog:
/// - every single non-identity class for `n = 2..=7`;
/// - every 2-class union for `n = 4..=6`;
/// - named special cases (mixed parity, disconnected V₄, n-cycles);
/// - seeded random 1–3-class unions for `n = 8..=12`.
pub fn union_catalog() -> Vec<CatalogEntry> {
    let mut out = Vec::new();

    for n in 2..=7u16 {
        let index = PartitionIndex::build(n).unwrap();
        for id in 0..index.count() as u32 {
            if id == index.identity_id() {
                continue;
            }
            let t = template_of(index.partition(id));
            out.push(CatalogEntry {
                n,
                label: format!("n{n}_single_{}", t.slug()),
                templates: vec![t],
            });
        }
    }

    for n in 4..=6u16 {
        let index = PartitionIndex::build(n).unwrap();
        let classes: Vec<PartitionId> = (0..index.count() as u32)
            .filter(|&id| id != index.identity_id())
            .collect();
        for (i, &a) in classes.iter().enumerate() {
            for &b in &classes[i + 1..] {
                let ta = template_of(index.partition(a));
                let tb = template_of(index.partition(b));
                out.push(CatalogEntry {
                    n,
                    label: format!("n{n}_pair_{}+{}", ta.slug(), tb.slug()),
                    templates: vec![ta, tb],
                });
            }
        }
    }

    // Named specials beyond the systematic ranges.
    for (n, specs, tag) in [
        (8u16, vec!["2"], "transpositions"),
        (8, vec!["3"], "three_cycles"),
        (8, vec!["2", "3"], "mixed"),
        (8, vec!["8"], "full_cycle"),
        (10, vec!["2"], "transpositions"),
        (10, vec!["2,2"], "double_transpositions"),
    ] {
        out.push(CatalogEntry {
            n,
            label: format!("n{n}_special_{tag}"),
            templates: specs.iter().map(|s| s.parse().unwrap()).collect(),
        });
    }

    // Seeded random unions.
    let mut rng = XorShift(0x5EED_C1A5_5D1A_2026);
    for n in 8..=12u16 {
        let index = PartitionIndex::build(n).unwrap();
        for k in 0..3 {
            let count = 1 + rng.below(3);
            let mut templates = Vec::new();
            let mut ids = Vec::new();
            while ids.len() < count {
                let id = rng.below(index.count()) as u32;
                if id != index.identity_id() && !ids.contains(&id) {
                    ids.push(id);
                    templates.push(template_of(index.partition(id)));
                }
            }
            out.push(CatalogEntry {
                n,
                label: format!("n{n}_random{k}"),
                templates,
            });
        }
    }

    out
}

/// Resolve an entry (never identity-including, always valid by construction).
pub fn resolve_entry(index: &PartitionIndex, entry: &CatalogEntry) -> ResolvedUnion {
    resolve_union(index, &entry.templates, false).expect("catalog entries are valid")
}

/// Brute-force cost guard: `|U| · n!` compositions must stay affordable
/// (design doc 03 §3.1). Debug builds get a tighter budget.
pub fn brute_force_affordable(index: &PartitionIndex, union: &ResolvedUnion) -> bool {
    let budget: u64 = if cfg!(debug_assertions) {
        20_000_000
    } else {
        300_000_000
    };
    let cost = union.union_size.clone() * index.factorial_n();
    cost <= BigUint::from(budget)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn catalog_is_deterministic_and_valid() {
        let a = union_catalog();
        let b = union_catalog();
        assert_eq!(a.len(), b.len());
        for (x, y) in a.iter().zip(b.iter()) {
            assert_eq!(x.label, y.label);
            assert_eq!(x.templates, y.templates);
        }
        // every entry resolves
        for entry in &a {
            let index = PartitionIndex::build(entry.n).unwrap();
            let union = resolve_entry(&index, entry);
            assert!(!union.class_ids.is_empty(), "{}", entry.label);
        }
        // sanity: systematic part sizes — singles: Σ_{n=2..7} (p(n)−1) = 1+2+4+6+10+14 = 37
        let singles = a.iter().filter(|e| e.label.contains("single")).count();
        assert_eq!(singles, 37);
        // pairs: C(4,2)+C(6,2)+C(10,2) = 6+15+45 = 66
        let pairs = a.iter().filter(|e| e.label.contains("pair")).count();
        assert_eq!(pairs, 66);
    }
}
