//! CPU transform backends: an obviously-correct reference and the blocked
//! rayon production kernel. Both use the proven `p31-u128-accumulate-v1`
//! regime (see the compile-time bound in [`crate::arith`]).
#![deny(clippy::float_arithmetic)]

use std::ops::Range;

use rayon::prelude::*;

use crate::arith::ModCtx;
use crate::chars::memtable::PairedModTable;
use crate::transform::{RadiusWeights, TransformBackend};

/// Naive single-threaded oracle: modular multiply-add per term.
pub struct CpuReference;

impl TransformBackend for CpuReference {
    fn name(&self) -> &'static str {
        "cpu-reference-v1"
    }

    fn numerators(
        &self,
        table: &PairedModTable,
        weights: &RadiusWeights,
        targets: Range<usize>,
    ) -> Vec<Vec<u32>> {
        assert_parity_aligned(table, &targets);
        let w = weights.for_parity(targets.start < table.even_count());
        let ctxs: Vec<ModCtx> = table.primes().iter().copied().map(ModCtx::new).collect();
        ctxs.iter()
            .enumerate()
            .map(|(lane, ctx)| {
                targets
                    .clone()
                    .map(|t| {
                        let column = table.column(lane, t);
                        let mut acc = 0u32;
                        for (x, wv) in column.iter().zip(&w[lane]) {
                            acc = ctx.add(acc, ctx.mul(*x, *wv));
                        }
                        acc
                    })
                    .collect()
            })
            .collect()
    }
}

/// Production kernel: u128 accumulation (one reduction per dot product),
/// rayon-parallel over targets.
pub struct CpuBlocked;

impl TransformBackend for CpuBlocked {
    fn name(&self) -> &'static str {
        "cpu-blocked-v1"
    }

    fn numerators(
        &self,
        table: &PairedModTable,
        weights: &RadiusWeights,
        targets: Range<usize>,
    ) -> Vec<Vec<u32>> {
        assert_parity_aligned(table, &targets);
        let rep_count = table.rep_count();
        debug_assert!((rep_count as u128) <= crate::arith::MAX_ACCUM_TERMS);
        let w = weights.for_parity(targets.start < table.even_count());
        let ctxs: Vec<ModCtx> = table.primes().iter().copied().map(ModCtx::new).collect();
        ctxs.iter()
            .enumerate()
            .map(|(lane, ctx)| {
                let lane_w = &w[lane];
                targets
                    .clone()
                    .into_par_iter()
                    .map(|t| {
                        let column = table.column(lane, t);
                        let mut acc: u128 = 0;
                        for (x, wv) in column.iter().zip(lane_w) {
                            acc += *x as u128 * *wv as u128;
                        }
                        ctx.reduce_u128(acc)
                    })
                    .collect()
            })
            .collect()
    }
}

fn assert_parity_aligned(table: &PairedModTable, targets: &Range<usize>) {
    let e = table.even_count();
    assert!(
        targets.end <= e || targets.start >= e,
        "target range {targets:?} straddles the parity boundary {e}"
    );
    assert!(targets.end <= table.targets().len());
}
