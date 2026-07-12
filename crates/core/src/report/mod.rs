//! Result document assembly (`classdiam/result` schema, version 1).
//!
//! The core builds the fully serde-serializable document; all file I/O and
//! wall-clock measurement stay in the caller (CLI / future PyO3). One schema
//! serves files, CLI, and Python (design doc 01 §9.2 with critique findings
//! 13/15 applied: `row_col_convention`, active-row information, and a
//! `support_reconstruction` rule conditional on the identity-generator flag).
//!
//! Volatile fields are confined to `run`, `timings_s`, and `engine.threads`
//! — golden-file comparisons strip exactly those (finding 14).

pub mod schema;

use num_traits::Zero;

use crate::engine::UnionRun;
use crate::partition::{CycleTypeTemplate, PartitionId, PartitionIndex, ORDER_CONVENTION};
use crate::spectra::{BaseSpectra, ResolvedUnion, UnionParity};
use schema::*;

/// Data only the caller knows (wall clock, identifiers, thread count).
#[derive(Clone, Debug, Default)]
pub struct RunMeta {
    pub run_id: String,
    pub started_utc: String,
    pub finished_utc: String,
    pub threads: u32,
    pub total_wall_s: f64,
    pub config_hash: String,
}

/// Union slug used in file names: classes joined `+`, parts joined `.`
/// (e.g. `3+2.2`); prefixed `g` by the CLI when forming file names.
pub fn union_slug(templates: &[CycleTypeTemplate]) -> String {
    templates
        .iter()
        .map(|t| t.slug())
        .collect::<Vec<_>>()
        .join("+")
}

#[allow(clippy::too_many_arguments)]
pub fn build_result(
    index: &PartitionIndex,
    templates: &[CycleTypeTemplate],
    union: &ResolvedUnion,
    spectra: &BaseSpectra,
    run: &UnionRun,
    label: Option<String>,
    allow_identity_generator: bool,
    meta: RunMeta,
) -> ResultDocument {
    let q = index.count();

    let classes: Vec<GeneratorClass> = union
        .class_ids
        .iter()
        .map(|&id| {
            let p = index.partition(id);
            GeneratorClass {
                template: p.parts().iter().copied().filter(|&x| x >= 2).collect(),
                padded: p.parts().to_vec(),
                index: id,
                class_size: index.class_size(id).to_string(),
                sign: index.sign(id),
            }
        })
        .collect();

    let partitions_reduced: Vec<Vec<u8>> = index
        .partitions()
        .iter()
        .map(|p| p.parts().iter().copied().filter(|&x| x >= 2).collect())
        .collect();

    let zero_rows_all_bases: Vec<PartitionId> = (0..q as u32)
        .filter(|&rho| {
            (0..union.class_ids.len()).all(|j| spectra.omega_column(j)[rho as usize].is_zero())
        })
        .collect();

    // Subgroup classification from the COMPUTED visited set (never theory):
    let even_ids: Vec<usize> = (0..q)
        .filter(|&i| index.sign(i as PartitionId) == 1)
        .collect();
    let reachable_ids: Vec<usize> = (0..q).filter(|&i| run.distance[i] >= 0).collect();
    let generated_subgroup = if reachable_ids.len() == q {
        GeneratedSubgroup::SN
    } else if reachable_ids == even_ids {
        GeneratedSubgroup::AN
    } else if reachable_ids.len() == 1 {
        GeneratedSubgroup::Trivial
    } else {
        GeneratedSubgroup::ProperSubgroup
    };

    let support_reconstruction = if union.includes_identity {
        "supports are cumulative (identity generator): support_r = { i : 0 <= distance[i] <= r }"
    } else {
        "support_r = { i : 0 <= first_hit[r mod 2][i] <= r }, valid for r <= stop_radius"
    };

    ResultDocument {
        format: "classdiam/result".into(),
        format_version: 1,
        spec_version: "character_method_cayley_diameters.md@2026-07".into(),
        tool: Tool {
            name: "classdiam".into(),
            version: env!("CARGO_PKG_VERSION").into(),
        },
        run: RunInfo {
            run_id: meta.run_id,
            started_utc: meta.started_utc,
            finished_utc: meta.finished_utc,
            resumed_from_checkpoint: false,
            suspend_resume_count: 0,
        },
        n: index.n(),
        factorial_n: index.factorial_n().to_string(),
        generators: Generators {
            input_templates: templates.iter().map(|t| t.parts().to_vec()).collect(),
            classes,
            union_size: union.union_size.to_string(),
            parity: union.parity,
            allow_identity_generator,
            label: label.unwrap_or_else(|| format!("g{}", union_slug(templates))),
        },
        partition_order: PartitionOrder {
            convention: ORDER_CONVENTION.into(),
            count: q as u64,
            hash_blake3: index.order_hash_hex(),
            partitions_reduced,
        },
        row_col_convention: "rows=irreps(rho), cols=classes(nu); both in the canonical order"
            .into(),
        class_data: ClassData {
            sign: (0..q).map(|i| index.sign(i as PartitionId)).collect(),
            class_size: (0..q)
                .map(|i| index.class_size(i as PartitionId).to_string())
                .collect(),
        },
        arithmetic: Arithmetic {
            mode: "exact-bigint".into(),
            resident_primes: vec![],
            certifier: "exact".into(),
        },
        engine: EngineInfo {
            mode: "exact_reference".into(),
            backend: "exact-bigint-v1".into(),
            rows_used: "full".into(),
            active_row_count: q as u64 - zero_rows_all_bases.len() as u64,
            zero_rows_all_bases,
            threads: meta.threads,
        },
        results: Results {
            unreachable_value: -1,
            distance: run.distance.clone(),
            first_hit_even: run.first_hit_even.clone(),
            first_hit_odd: run.first_hit_odd.clone(),
            diameter_identity_component: run.diameter,
            reachable_count: run.reachable_count as u64,
            generated_subgroup,
            cayley_graph_on_sn: if reachable_ids.len() == q {
                "connected".into()
            } else {
                "disconnected".into()
            },
            bipartite: union.parity == UnionParity::Odd,
            stopping: Stopping {
                rule: run.stopping,
                stop_radius: run.stop_radius,
            },
            support_reconstruction: support_reconstruction.into(),
            layers: run
                .layers
                .iter()
                .map(|layer| LayerOut {
                    r: layer.r,
                    new: layer.new.clone(),
                    support: layer.support.clone(),
                    support_size: layer.support.len() as u64,
                })
                .collect(),
        },
        timings_s: Timings {
            total_wall: meta.total_wall_s,
        },
        config_hash_blake3: meta.config_hash,
    }
}
