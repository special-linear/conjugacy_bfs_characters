//! Fixture consumers: the committed Python-generated character data must
//! agree with the Rust implementations (design doc 03 §4.4).
//!
//! Discipline (spec §23 Failure 7): the partition ORDER is asserted before
//! any value comparison, and the payload hash is verified first of all.

use std::path::PathBuf;

use classdiam_core::chars::{degree, MnEvaluator};
use classdiam_core::partition::{PartitionId, PartitionIndex, ORDER_CONVENTION};
use num_bigint::BigInt;
use serde_json::Value;
use sha2::{Digest, Sha256};

fn fixtures_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures")
}

fn load(path: &PathBuf) -> Value {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| {
        panic!(
            "cannot read {} ({e}); run tools/gen_fixtures.py",
            path.display()
        )
    });
    serde_json::from_str(&text).unwrap()
}

/// sha256 over compact-serialized segments joined by `|` — the documented
/// payload-hash rule shared with tools/gen_fixtures.py.
fn payload_sha256(segments: &[&Value]) -> String {
    let joined = segments
        .iter()
        .map(|v| serde_json::to_string(v).unwrap())
        .collect::<Vec<_>>()
        .join("|");
    format!("{:x}", Sha256::digest(joined.as_bytes()))
}

/// Order-first assertion: fixture partitions must equal the Rust canonical
/// enumeration exactly; abort with both lists on mismatch.
fn assert_order_matches(fixture: &Value, index: &PartitionIndex) {
    assert_eq!(
        fixture["partition_order_convention"].as_str().unwrap(),
        ORDER_CONVENTION
    );
    let fixture_parts: Vec<Vec<u8>> =
        serde_json::from_value(fixture["partitions"].clone()).unwrap();
    let ours: Vec<Vec<u8>> = index
        .partitions()
        .iter()
        .map(|p| p.parts().to_vec())
        .collect();
    assert_eq!(
        fixture_parts,
        ours,
        "PARTITION ORDER MISMATCH at n={} — refusing all value comparisons",
        index.n()
    );
}

fn files_matching(subdir: &str, prefix: &str) -> Vec<PathBuf> {
    let dir = fixtures_dir().join(subdir);
    let mut out: Vec<PathBuf> = std::fs::read_dir(&dir)
        .unwrap_or_else(|e| panic!("missing {} ({e}); run tools/gen_fixtures.py", dir.display()))
        .map(|e| e.unwrap().path())
        .filter(|p| {
            p.file_name()
                .and_then(|f| f.to_str())
                .is_some_and(|f| f.starts_with(prefix) && f.ends_with(".json"))
        })
        .collect();
    out.sort();
    assert!(!out.is_empty(), "no fixtures {subdir}/{prefix}*");
    out
}

#[test]
fn full_character_tables_match_rust_evaluator() {
    for path in files_matching("characters", "char_n") {
        let fixture = load(&path);
        let n = fixture["n"].as_u64().unwrap() as u16;
        assert_eq!(fixture["schema_version"].as_str().unwrap(), "fixture.v1");
        assert_eq!(
            fixture["payload_sha256"].as_str().unwrap(),
            payload_sha256(&[&fixture["partitions"], &fixture["table"]]),
            "payload hash mismatch in {}",
            path.display()
        );

        let index = PartitionIndex::build(n).unwrap();
        assert_order_matches(&fixture, &index);

        // degrees vs Rust hook-length formula
        let degrees: Vec<String> = serde_json::from_value(fixture["degrees"].clone()).unwrap();
        for (i, d) in degrees.iter().enumerate() {
            assert_eq!(
                d,
                &degree(index.partition(i as PartitionId)).to_string(),
                "degree mismatch at n={n}, rho={i}"
            );
        }

        // full table vs the production evaluator
        let table: Vec<Vec<i64>> = serde_json::from_value(fixture["table"].clone()).unwrap();
        let ev = MnEvaluator::new(n);
        let ours = ev.full_table_exact();
        for (nu, row) in table.iter().enumerate() {
            for (rho, &value) in row.iter().enumerate() {
                assert_eq!(
                    ours[nu][rho],
                    BigInt::from(value),
                    "chi mismatch at n={n}, nu={nu}, rho={rho}"
                );
            }
        }
    }
}

#[test]
fn spot_values_match_rust_evaluator() {
    for path in files_matching("characters", "spot_n") {
        let fixture = load(&path);
        let n = fixture["n"].as_u64().unwrap() as u16;
        assert_eq!(
            fixture["payload_sha256"].as_str().unwrap(),
            payload_sha256(&[&fixture["partitions"], &fixture["spots"]]),
            "payload hash mismatch in {}",
            path.display()
        );
        let index = PartitionIndex::build(n).unwrap();
        assert_order_matches(&fixture, &index);

        let spots: Vec<(usize, usize, String)> =
            serde_json::from_value(fixture["spots"].clone()).unwrap();
        // group by target column to reuse column computations
        let mut by_nu: std::collections::BTreeMap<usize, Vec<(usize, String)>> = Default::default();
        for (rho, nu, value) in spots {
            by_nu.entry(nu).or_default().push((rho, value));
        }
        let ev = MnEvaluator::new(n);
        for (nu, entries) in by_nu {
            let column = ev.column_exact(index.partition(nu as PartitionId));
            for (rho, value) in entries {
                assert_eq!(
                    column[rho].to_string(),
                    value,
                    "spot mismatch at n={n}, rho={rho}, nu={nu}"
                );
            }
        }
    }
}

/// The committed adversarial tuples (P2 modular-engine test inputs) are
/// verified against the exact engine NOW: `a_r(ν)` really is positive and
/// really vanishes modulo every listed prime.
#[test]
fn adversarial_tuples_are_genuine() {
    use classdiam_core::engine::exact::ExactTransform;
    use classdiam_core::partition::CycleTypeTemplate;
    use classdiam_core::spectra::{resolve_union, BaseSpectra};
    use num_traits::Zero;

    let fixture = load(&fixtures_dir().join("adversarial_v1.json"));
    let tuples = fixture["tuples"].as_array().unwrap();
    assert!(tuples.len() >= 50, "suspiciously few adversarial tuples");

    let mut masking = 0;
    for tuple in tuples {
        let n = tuple["n"].as_u64().unwrap() as u16;
        let r = tuple["r"].as_u64().unwrap() as u32;
        let nu = tuple["nu_index"].as_u64().unwrap() as usize;
        let templates: Vec<CycleTypeTemplate> = tuple["union_templates"]
            .as_array()
            .unwrap()
            .iter()
            .map(|t| {
                let parts: Vec<u8> = serde_json::from_value(t.clone()).unwrap();
                CycleTypeTemplate::new(&parts).unwrap()
            })
            .collect();
        let primes: Vec<u64> = serde_json::from_value(tuple["primes"].clone()).unwrap();
        let expected: String = tuple["a_r"].as_str().unwrap().to_string();

        let index = PartitionIndex::build(n).unwrap();
        let ev = MnEvaluator::new(n);
        let union = resolve_union(&index, &templates, false).unwrap();
        let spectra = BaseSpectra::build(&index, &ev, &union.class_ids).unwrap();
        let theta = spectra.theta(&union.class_ids);
        let table = ev.full_table_exact();
        let mut transform =
            ExactTransform::new(&table, spectra.degrees(), theta, index.factorial_n());
        for _ in 0..r {
            transform.advance();
        }
        let a = &transform.coefficients().unwrap()[nu];
        assert_eq!(a.to_string(), expected, "tuple value mismatch");
        assert!(!a.is_zero(), "adversarial coefficient must be positive");
        for p in primes {
            assert!(
                (a % BigInt::from(p)).is_zero(),
                "coefficient not divisible by claimed prime {p}"
            );
        }
        if tuple["masks_first_hit"].as_bool().unwrap() {
            masking += 1;
        }
    }
    assert!(masking >= 20, "need first-hit-masking tuples for P2 tests");
}

#[test]
fn degrees_and_class_data_match() {
    for path in files_matching("degrees", "deg_n") {
        let fixture = load(&path);
        let n = fixture["n"].as_u64().unwrap() as u16;
        assert_eq!(
            fixture["payload_sha256"].as_str().unwrap(),
            payload_sha256(&[
                &fixture["partitions"],
                &fixture["degrees"],
                &fixture["class_sizes"],
                &fixture["z_values"],
                &fixture["signs"],
                &fixture["transpose_map"],
            ]),
            "payload hash mismatch in {}",
            path.display()
        );
        let index = PartitionIndex::build(n).unwrap();
        assert_order_matches(&fixture, &index);

        let degrees: Vec<String> = serde_json::from_value(fixture["degrees"].clone()).unwrap();
        let class_sizes: Vec<String> =
            serde_json::from_value(fixture["class_sizes"].clone()).unwrap();
        let z_values: Vec<String> = serde_json::from_value(fixture["z_values"].clone()).unwrap();
        let signs: Vec<i8> = serde_json::from_value(fixture["signs"].clone()).unwrap();
        let transpose_map: Vec<u32> =
            serde_json::from_value(fixture["transpose_map"].clone()).unwrap();

        for i in 0..index.count() {
            let id = i as PartitionId;
            assert_eq!(
                degrees[i],
                degree(index.partition(id)).to_string(),
                "n={n} deg {i}"
            );
            assert_eq!(
                class_sizes[i],
                index.class_size(id).to_string(),
                "n={n} size {i}"
            );
            assert_eq!(z_values[i], index.z_value(id).to_string(), "n={n} z {i}");
            assert_eq!(signs[i], index.sign(id), "n={n} sign {i}");
            assert_eq!(
                transpose_map[i],
                index.transpose_id(id),
                "n={n} transpose {i}"
            );
        }
    }
}
