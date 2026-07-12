//! Result-document schema tests: the committed golden file for the
//! critique-verified n=6 transpositions example, roundtrip stability, and
//! the support-only guarantee (no factorization counts anywhere).
//!
//! Regenerate the golden file with:
//! `UPDATE_GOLDEN=1 cargo test -p classdiam-core --test output`

use std::path::PathBuf;

use classdiam_core::chars::MnEvaluator;
use classdiam_core::engine::exact::run_exact;
use classdiam_core::partition::{CycleTypeTemplate, PartitionIndex};
use classdiam_core::report::schema::ResultDocument;
use classdiam_core::report::{build_result, RunMeta};
use classdiam_core::spectra::{resolve_union, BaseSpectra};
use serde_json::Value;

fn golden_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/golden")
}

/// Deterministic document for S_6 / transpositions with placeholder
/// volatile fields.
fn n6_transpositions_document() -> ResultDocument {
    let index = PartitionIndex::build(6).unwrap();
    let mn = MnEvaluator::new(6);
    let templates = vec![CycleTypeTemplate::new([2]).unwrap()];
    let union = resolve_union(&index, &templates, false).unwrap();
    let spectra = BaseSpectra::build(&index, &mn, &union.class_ids).unwrap();
    let run = run_exact(&index, &mn, &union).unwrap();
    build_result(
        &index,
        &templates,
        &union,
        &spectra,
        &run,
        Some("g2".into()),
        false,
        RunMeta {
            run_id: "golden".into(),
            started_utc: "1970-01-01T00:00:00Z".into(),
            finished_utc: "1970-01-01T00:00:00Z".into(),
            threads: 1,
            total_wall_s: 0.0,
            config_hash: "golden".into(),
        },
    )
}

/// Strip the documented volatile fields (design finding 14): `run`,
/// `timings_s`, `engine.threads`, and `tool` (version churn).
fn strip_volatile(value: &mut Value) {
    let object = value.as_object_mut().expect("document is an object");
    object.remove("run");
    object.remove("timings_s");
    object.remove("tool");
    if let Some(engine) = object.get_mut("engine").and_then(Value::as_object_mut) {
        engine.remove("threads");
    }
}

#[test]
fn golden_n06_g2() {
    let document = n6_transpositions_document();
    let mut current = serde_json::to_value(&document).unwrap();
    strip_volatile(&mut current);

    let path = golden_dir().join("n06_g2.json");
    if std::env::var("UPDATE_GOLDEN").is_ok() {
        std::fs::create_dir_all(golden_dir()).unwrap();
        std::fs::write(&path, serde_json::to_string_pretty(&document).unwrap()).unwrap();
        return;
    }
    let committed = std::fs::read_to_string(&path).unwrap_or_else(|e| {
        panic!(
            "missing golden file {} ({e}); run with UPDATE_GOLDEN=1",
            path.display()
        )
    });
    let mut expected: Value = serde_json::from_str(&committed).unwrap();
    strip_volatile(&mut expected);
    assert_eq!(
        current, expected,
        "golden mismatch — schema or math changed"
    );
}

#[test]
fn golden_values_match_worked_example() {
    // Independent of the file: the numbers the critique verified by hand.
    let document = n6_transpositions_document();
    assert_eq!(
        document.results.distance,
        vec![5, 4, 4, 3, 4, 3, 2, 3, 2, 1, 0]
    );
    assert_eq!(document.results.diameter_identity_component, 5);
    assert_eq!(document.results.reachable_count, 11);
    assert!(document.results.bipartite);
    assert_eq!(document.generators.union_size, "15");
    assert_eq!(document.partition_order.count, 11);
    assert_eq!(
        document.partition_order.partitions_reduced,
        vec![
            vec![6u8],
            vec![5],
            vec![4, 2],
            vec![4],
            vec![3, 3],
            vec![3, 2],
            vec![3],
            vec![2, 2, 2],
            vec![2, 2],
            vec![2],
            vec![],
        ]
    );
    assert_eq!(
        document.class_data.class_size,
        vec!["120", "144", "90", "90", "40", "120", "40", "15", "45", "15", "1"]
    );
    // ρ = [3,2,1] (index 5) is the unique zero row for transpositions
    // (content sum 0) — the row the design docs call out.
    assert_eq!(document.engine.zero_rows_all_bases, vec![5]);
    assert_eq!(document.engine.active_row_count, 10);
}

#[test]
fn document_roundtrips_through_json() {
    let document = n6_transpositions_document();
    let json = serde_json::to_string(&document).unwrap();
    let back: ResultDocument = serde_json::from_str(&json).unwrap();
    let a = serde_json::to_value(&document).unwrap();
    let b = serde_json::to_value(&back).unwrap();
    assert_eq!(a, b);
}

/// Support-only guarantee: no factorization counts of any kind appear in
/// the document (fixed requirement 2; spec §23 F2 is unrepresentable).
#[test]
fn no_word_counts_in_output() {
    let document = n6_transpositions_document();
    let value = serde_json::to_value(&document).unwrap();
    let mut keys = Vec::new();
    collect_keys(&value, &mut keys);
    for forbidden in ["word_count", "coefficient", "a_r", "counts"] {
        assert!(
            !keys.iter().any(|k| k.contains(forbidden)),
            "forbidden key fragment {forbidden:?} present"
        );
    }
}

fn collect_keys(value: &Value, out: &mut Vec<String>) {
    match value {
        Value::Object(map) => {
            for (k, v) in map {
                out.push(k.clone());
                collect_keys(v, out);
            }
        }
        Value::Array(items) => items.iter().for_each(|v| collect_keys(v, out)),
        _ => {}
    }
}
