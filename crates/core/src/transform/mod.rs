//! The per-radius modular transform behind a backend trait — the GPU seam
//! (design doc 01 §6; contract summary in `docs/gpu_backend.md`).
//!
//! One call computes numerators `N_r(ν) = Σ_{ρ ∈ reps} X[ρ,ν]·W_{sgn(ν)}[ρ]
//! (mod p)` for a contiguous, parity-aligned range of table targets. The
//! division by `n!` never happens on this path: `gcd(n!, p) = 1`, so
//! `a_r(ν) ≡ 0 ⇔ N_r(ν) ≡ 0`.
//!
//! ## Backend contract (normative)
//!
//! - All residues are fully reduced `u32` in `[0, p)`, plain representation.
//! - Outputs are **bit-exact**: every backend must produce identical values
//!   (no floating point, no tolerance) — this is what makes checkpoints and
//!   CPU/GPU mixing safe, and gives a free cross-backend test.
//! - Accumulation strategy is the backend's choice but must carry a proven
//!   overflow bound. CPU regime `p31-u128-accumulate-v1`: products `< 2⁶²`,
//!   `u128` accumulator, row length ≤ [`crate::arith::MAX_ACCUM_TERMS`].
//!   The documented GPU regime uses primes `< 2³⁰` with `u64` accumulators
//!   reduced at least every 16 products (`16·(2³⁰−1)² + p < 2⁶⁴`, strict).
//! - `targets` ranges never straddle the even/odd column boundary; the
//!   caller binds `w_plus` for even-sign blocks and `w_minus` for odd.
//! - Backends may cache device-side table planes keyed by (prime, table
//!   identity); table contents are immutable for the life of the table.
#![deny(clippy::float_arithmetic)]

pub mod cpu;

use std::ops::Range;

use crate::chars::memtable::PairedModTable;

/// Per-radius weight vectors over representative rows, one pair per prime
/// lane: `W₊ = f·(P + P′)`, `W₋ = f·(P − P′)` for paired rows and
/// `W₊ = W₋ = f·P` for self-transpose rows (assembled by the engine).
pub struct RadiusWeights {
    /// `w_plus[lane][rep]`
    pub w_plus: Vec<Vec<u32>>,
    /// `w_minus[lane][rep]`
    pub w_minus: Vec<Vec<u32>>,
}

impl RadiusWeights {
    pub fn for_parity(&self, even_targets: bool) -> &[Vec<u32>] {
        if even_targets {
            &self.w_plus
        } else {
            &self.w_minus
        }
    }
}

pub trait TransformBackend: Send + Sync {
    fn name(&self) -> &'static str;

    /// `out[lane][t − targets.start] = Σ_rep column(lane, t)[rep] · w[lane][rep] (mod p)`
    /// where `w` is `weights.for_parity(...)` for the (single) parity of the
    /// range. `targets` must lie entirely inside one parity block.
    fn numerators(
        &self,
        table: &PairedModTable,
        weights: &RadiusWeights,
        targets: Range<usize>,
    ) -> Vec<Vec<u32>>;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arith::{screening_primes, ModCtx, Prime31};
    use crate::chars::MnEvaluator;
    use crate::partition::PartitionIndex;

    /// CpuReference and CpuBlocked agree bit-for-bit on random weights —
    /// the harness that later pins a GPU backend.
    #[test]
    fn backends_agree_bitwise() {
        let n = 8u16;
        let index = PartitionIndex::build(n).unwrap();
        let mn = MnEvaluator::new(n);
        let mut primes = screening_primes(2);
        primes.push(Prime31(13));
        let table = PairedModTable::generate(&index, &mn, &primes);
        let ctxs: Vec<ModCtx> = primes.iter().copied().map(ModCtx::new).collect();

        // deterministic pseudo-random weights
        let mut state = 0x1234_5678_9abc_def0u64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let weights = RadiusWeights {
            w_plus: ctxs
                .iter()
                .map(|c| {
                    (0..table.rep_count())
                        .map(|_| (next() % c.prime().0 as u64) as u32)
                        .collect()
                })
                .collect(),
            w_minus: ctxs
                .iter()
                .map(|c| {
                    (0..table.rep_count())
                        .map(|_| (next() % c.prime().0 as u64) as u32)
                        .collect()
                })
                .collect(),
        };

        let reference = cpu::CpuReference;
        let blocked = cpu::CpuBlocked;
        for range in [
            0..table.even_count(),
            table.even_count()..table.targets().len(),
        ] {
            let a = reference.numerators(&table, &weights, range.clone());
            let b = blocked.numerators(&table, &weights, range.clone());
            assert_eq!(a, b, "backends disagree on {range:?}");
        }
    }
}
