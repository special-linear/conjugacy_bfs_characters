//! Serde types for the `classdiam/result` JSON schema, version 1.
//!
//! Field-by-field documentation lives in `docs/output_schema.md`; parsers
//! must reject unknown `format_version` majors. Every document is
//! self-describing: it embeds the explicit partition order (reduced form:
//! parts ≥ 2 only, identity = `[]`), its hash, class data, and all run
//! metadata (spec §19.3).

use serde::{Deserialize, Serialize};

use crate::engine::StoppingRule;
use crate::partition::PartitionId;
use crate::spectra::UnionParity;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ResultDocument {
    pub format: String,
    pub format_version: u32,
    pub spec_version: String,
    pub tool: Tool,
    /// VOLATILE: stripped in golden-file comparisons.
    pub run: RunInfo,
    pub n: u16,
    pub factorial_n: String,
    pub generators: Generators,
    pub partition_order: PartitionOrder,
    pub row_col_convention: String,
    pub class_data: ClassData,
    pub arithmetic: Arithmetic,
    pub engine: EngineInfo,
    pub results: Results,
    /// VOLATILE: stripped in golden-file comparisons.
    pub timings_s: Timings,
    pub config_hash_blake3: String,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Tool {
    pub name: String,
    pub version: String,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RunInfo {
    pub run_id: String,
    pub started_utc: String,
    pub finished_utc: String,
    pub resumed_from_checkpoint: bool,
    pub suspend_resume_count: u32,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Generators {
    /// Cycle-type templates exactly as given (no fixed points).
    pub input_templates: Vec<Vec<u8>>,
    pub classes: Vec<GeneratorClass>,
    /// `|U|` as a decimal string.
    pub union_size: String,
    pub parity: UnionParity,
    pub allow_identity_generator: bool,
    pub label: String,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct GeneratorClass {
    pub template: Vec<u8>,
    /// Full cycle type after padding with fixed points.
    pub padded: Vec<u8>,
    /// Canonical index of the class.
    pub index: PartitionId,
    pub class_size: String,
    pub sign: i8,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PartitionOrder {
    pub convention: String,
    pub count: u64,
    pub hash_blake3: String,
    /// Full ordered list, reduced form (parts ≥ 2; identity type = []).
    pub partitions_reduced: Vec<Vec<u8>>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ClassData {
    pub sign: Vec<i8>,
    pub class_size: Vec<String>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Arithmetic {
    /// `"exact-bigint"` (reference engine) or `"modular+certified"` (P2+).
    pub mode: String,
    pub resident_primes: Vec<u32>,
    pub certifier: String,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct EngineInfo {
    pub mode: String,
    pub backend: String,
    /// `"full"` or `"restricted"` — which row set the transform ran on.
    pub rows_used: String,
    /// Rows with some `ω ≠ 0` over the union's classes.
    pub active_row_count: u64,
    /// Rows where `ω_ρ(λ) = 0` for EVERY generator class (dropped by a
    /// restricted-row engine for r ≥ 1; spec §7).
    pub zero_rows_all_bases: Vec<PartitionId>,
    /// VOLATILE: stripped in golden-file comparisons.
    pub threads: u32,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum GeneratedSubgroup {
    #[serde(rename = "S_n")]
    SN,
    #[serde(rename = "A_n")]
    AN,
    #[serde(rename = "proper_subgroup")]
    ProperSubgroup,
    #[serde(rename = "trivial")]
    Trivial,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Results {
    pub unreachable_value: i32,
    /// Per canonical cycle type; `-1` = unreachable.
    pub distance: Vec<i32>,
    /// First EVEN/ODD radius with `a_r(ν) > 0` within the run span;
    /// `-1` = none up to `stop_radius`.
    pub first_hit_even: Vec<i32>,
    pub first_hit_odd: Vec<i32>,
    pub diameter_identity_component: u32,
    pub reachable_count: u64,
    /// Derived from the COMPUTED visited set, never from subgroup theory.
    pub generated_subgroup: GeneratedSubgroup,
    pub cayley_graph_on_sn: String,
    pub bipartite: bool,
    pub stopping: Stopping,
    pub support_reconstruction: String,
    pub layers: Vec<LayerOut>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Stopping {
    pub rule: StoppingRule,
    pub stop_radius: u32,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct LayerOut {
    pub r: u32,
    pub new: Vec<PartitionId>,
    pub support: Vec<PartitionId>,
    pub support_size: u64,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Timings {
    pub total_wall: f64,
}
