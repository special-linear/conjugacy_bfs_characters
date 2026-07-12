//! The diameter engine: radius iteration, support classification, distance
//! updates, stopping logic (spec §5, §8).
//!
//! P1 ships the exact big-integer reference engine ([`exact`]); the modular
//! production engine with the tiered certification gate arrives in P2 and is
//! differentially tested against this one.
#![deny(clippy::float_arithmetic)]

pub mod exact;
pub mod modular;

use fixedbitset::FixedBitSet;
use serde::{Deserialize, Serialize};

use crate::partition::{PartitionId, PartitionIndex};
use crate::spectra::UnionParity;

/// Parity upper bound on the reachable set (pure sign argument, never
/// subgroup theory): even-only unions stay inside even-sign types; odd or
/// mixed unions have no restriction.
pub fn parity_feasible_set(index: &PartitionIndex, parity: UnionParity) -> FixedBitSet {
    let q = index.count();
    let mut feasible = FixedBitSet::with_capacity(q);
    match parity {
        UnionParity::Even => {
            for nu in 0..q {
                if index.sign(nu as PartitionId) == 1 {
                    feasible.insert(nu);
                }
            }
        }
        UnionParity::Odd | UnionParity::Mixed => feasible.insert_range(..),
    }
    feasible
}

/// Why a run terminated (merged design: the spec §5.2 empty-layer rule is
/// primary; the parity-bound cover check is a sound early exit).
#[derive(Copy, Clone, PartialEq, Eq, Debug, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum StoppingRule {
    /// No new cycle type appeared at this radius (spec §5.2; valid only on
    /// exact supports).
    EmptyLayer,
    /// Every parity-feasible cycle type has been visited — no new type can
    /// ever appear (pure parity upper bound, no subgroup theory).
    AllTypesVisited,
}

/// One committed BFS layer.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct LayerRecord {
    pub r: u32,
    /// Cycle types first reached at exactly this radius.
    pub new: Vec<PartitionId>,
    /// Exact-length support: types with `a_r(ν) > 0`.
    pub support: Vec<PartitionId>,
}

/// Core result of one `(n, union)` run — everything downstream serialization
/// needs, independent of the arithmetic backend that produced it.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct UnionRun {
    pub n: u16,
    /// Distance from the identity per canonical cycle type; `-1` unreachable.
    pub distance: Vec<i32>,
    /// First EVEN radius with `a_r(ν) > 0`; `-1` = none up to stop radius.
    pub first_hit_even: Vec<i32>,
    /// First ODD radius with `a_r(ν) > 0`; `-1` = none up to stop radius.
    pub first_hit_odd: Vec<i32>,
    pub layers: Vec<LayerRecord>,
    pub diameter: u32,
    pub stop_radius: u32,
    pub stopping: StoppingRule,
    pub reachable_count: usize,
}
