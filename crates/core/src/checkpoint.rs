//! Checkpoint files: fixed header + postcard body + CRC32 (design doc 01
//! §10). Only committed, fully certified state is ever serialized — a
//! checkpoint can never encode an uncertified stop. Powers and tables are
//! NOT stored; they are recomputed on resume (recompute-don't-reread).
//!
//! Header (80 bytes):
//! `magic "CDCK" | format_version u16 LE | flags u16 LE | config_hash [32] |
//!  order_hash [32] | body_len u64 LE` — then the postcard body, then
//! `crc32(header ‖ body)` as u32 LE.
//!
//! Loading refuses on any mismatch: magic, version, CRC, config hash, order
//! hash. There is deliberately no `--force` override — a mismatched resume
//! is mathematically meaningless (spec §19.3 / Failure 7).

use std::io::Write as _;
use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::engine::modular::CertificationStats;
use crate::engine::LayerRecord;
use crate::error::ClassdiamError;

pub const MAGIC: &[u8; 4] = b"CDCK";
pub const FORMAT_VERSION: u16 = 1;

/// Everything needed to resume one `(n, union)` job at
/// `committed_radius + 1`. All layers `≤ committed_radius` are final and
/// certified.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct CheckpointBody {
    pub n: u16,
    /// Padded cycle types of the union's classes (canonical parts lists).
    pub resolved_classes: Vec<Vec<u8>>,
    pub allow_identity_generator: bool,
    pub primes: Vec<u32>,
    pub committed_radius: u32,
    pub distance: Vec<i32>,
    pub first_hit_even: Vec<i32>,
    pub first_hit_odd: Vec<i32>,
    pub layers: Vec<LayerRecord>,
    pub cert_stats: CertificationStats,
    pub suspend_count: u32,
}

fn header_bytes(config_hash: &[u8; 32], order_hash: &[u8; 32], body_len: u64) -> [u8; 80] {
    let mut h = [0u8; 80];
    h[0..4].copy_from_slice(MAGIC);
    h[4..6].copy_from_slice(&FORMAT_VERSION.to_le_bytes());
    // flags (6..8) reserved, zero
    h[8..40].copy_from_slice(config_hash);
    h[40..72].copy_from_slice(order_hash);
    h[72..80].copy_from_slice(&body_len.to_le_bytes());
    h
}

/// Atomic write: `.tmp` + rename; the previous checkpoint (if any) is kept
/// as `.prev` until the rename succeeds.
pub fn write_checkpoint(
    path: &Path,
    config_hash: &[u8; 32],
    order_hash: &[u8; 32],
    body: &CheckpointBody,
) -> Result<(), ClassdiamError> {
    let body_bytes = postcard::to_allocvec(body).map_err(|e| ClassdiamError::CheckpointFormat {
        reason: format!("encode: {e}"),
    })?;
    let header = header_bytes(config_hash, order_hash, body_bytes.len() as u64);
    let mut hasher = crc32fast::Hasher::new();
    hasher.update(&header);
    hasher.update(&body_bytes);
    let crc = hasher.finalize();

    let tmp = path.with_extension("ckpt.tmp");
    {
        let mut f = std::fs::File::create(&tmp)?;
        f.write_all(&header)?;
        f.write_all(&body_bytes)?;
        f.write_all(&crc.to_le_bytes())?;
        f.sync_all()?;
    }
    if path.exists() {
        let prev = path.with_extension("ckpt.prev");
        let _ = std::fs::remove_file(&prev);
        let _ = std::fs::rename(path, &prev);
    }
    std::fs::rename(&tmp, path)?;
    Ok(())
}

/// Read and fully validate a checkpoint. `expected_config`/`expected_order`
/// of `None` skip that comparison (the caller re-derives and compares the
/// config hash from the body when the config is reconstructed from the
/// checkpoint itself).
pub fn read_checkpoint(
    path: &Path,
    expected_config: Option<&[u8; 32]>,
    expected_order: Option<&[u8; 32]>,
) -> Result<(CheckpointBody, [u8; 32], [u8; 32]), ClassdiamError> {
    let data = std::fs::read(path)?;
    if data.len() < 84 {
        return Err(ClassdiamError::CheckpointFormat {
            reason: "truncated (shorter than header + crc)".into(),
        });
    }
    if &data[0..4] != MAGIC {
        return Err(ClassdiamError::CheckpointFormat {
            reason: "bad magic".into(),
        });
    }
    let version = u16::from_le_bytes([data[4], data[5]]);
    if version != FORMAT_VERSION {
        return Err(ClassdiamError::CheckpointFormat {
            reason: format!("unsupported format version {version}"),
        });
    }
    let mut config_hash = [0u8; 32];
    config_hash.copy_from_slice(&data[8..40]);
    let mut order_hash = [0u8; 32];
    order_hash.copy_from_slice(&data[40..72]);
    let body_len = u64::from_le_bytes(data[72..80].try_into().expect("8 bytes")) as usize;
    if data.len() != 80 + body_len + 4 {
        return Err(ClassdiamError::CheckpointFormat {
            reason: format!(
                "length mismatch: header says {} body bytes, file has {}",
                body_len,
                data.len().saturating_sub(84)
            ),
        });
    }
    let (payload, crc_bytes) = data.split_at(80 + body_len);
    let expected_crc = u32::from_le_bytes(crc_bytes.try_into().expect("4 bytes"));
    let mut hasher = crc32fast::Hasher::new();
    hasher.update(payload);
    if hasher.finalize() != expected_crc {
        return Err(ClassdiamError::CheckpointFormat {
            reason: "CRC mismatch (corrupted file)".into(),
        });
    }
    if let Some(expected) = expected_config {
        if expected != &config_hash {
            return Err(ClassdiamError::CheckpointMismatch {
                what:
                    "config hash differs — this checkpoint belongs to a different run configuration"
                        .into(),
            });
        }
    }
    if let Some(expected) = expected_order {
        if expected != &order_hash {
            return Err(ClassdiamError::CheckpointMismatch {
                what: "partition-order hash differs — refusing to mix order conventions".into(),
            });
        }
    }
    let body: CheckpointBody =
        postcard::from_bytes(&payload[80..]).map_err(|e| ClassdiamError::CheckpointFormat {
            reason: format!("decode: {e}"),
        })?;
    Ok((body, config_hash, order_hash))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_body() -> CheckpointBody {
        CheckpointBody {
            n: 9,
            resolved_classes: vec![vec![2, 1, 1, 1, 1, 1, 1, 1]],
            allow_identity_generator: false,
            primes: vec![2147483647, 2147483629],
            committed_radius: 3,
            distance: vec![0, 1, -1, 2],
            first_hit_even: vec![0, -1, -1, 2],
            first_hit_odd: vec![-1, 1, -1, -1],
            layers: vec![LayerRecord {
                r: 0,
                new: vec![0],
                support: vec![0],
            }],
            cert_stats: CertificationStats {
                candidates: 5,
                bound_certified: 3,
                crt_resident_certified: 2,
                ..Default::default()
            },
            suspend_count: 1,
        }
    }

    #[test]
    fn roundtrip_and_refusals() {
        let dir = std::env::temp_dir().join("classdiam_ckpt_test");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("job.ckpt");
        let config = [7u8; 32];
        let order = [9u8; 32];
        let body = sample_body();
        write_checkpoint(&path, &config, &order, &body).unwrap();

        let (back, ch, oh) = read_checkpoint(&path, Some(&config), Some(&order)).unwrap();
        assert_eq!(back, body);
        assert_eq!(ch, config);
        assert_eq!(oh, order);

        // config-hash refusal
        let wrong = [8u8; 32];
        assert!(matches!(
            read_checkpoint(&path, Some(&wrong), Some(&order)),
            Err(ClassdiamError::CheckpointMismatch { .. })
        ));
        // order-hash refusal
        assert!(matches!(
            read_checkpoint(&path, Some(&config), Some(&wrong)),
            Err(ClassdiamError::CheckpointMismatch { .. })
        ));

        // corruption: flip one body byte -> CRC failure
        let mut raw = std::fs::read(&path).unwrap();
        let mid = 80 + 3;
        raw[mid] ^= 0xFF;
        let bad = dir.join("bad.ckpt");
        std::fs::write(&bad, &raw).unwrap();
        assert!(matches!(
            read_checkpoint(&bad, None, None),
            Err(ClassdiamError::CheckpointFormat { .. })
        ));

        // truncation
        raw.truncate(50);
        std::fs::write(&bad, &raw).unwrap();
        assert!(matches!(
            read_checkpoint(&bad, None, None),
            Err(ClassdiamError::CheckpointFormat { .. })
        ));

        // unknown future version
        let mut raw = std::fs::read(&path).unwrap();
        raw[4] = 99;
        std::fs::write(&bad, &raw).unwrap();
        assert!(matches!(
            read_checkpoint(&bad, None, None),
            Err(ClassdiamError::CheckpointFormat { .. })
        ));

        // overwrite keeps .prev
        let mut body2 = sample_body();
        body2.committed_radius = 4;
        write_checkpoint(&path, &config, &order, &body2).unwrap();
        assert!(path.with_extension("ckpt.prev").exists());
        let (latest, _, _) = read_checkpoint(&path, Some(&config), Some(&order)).unwrap();
        assert_eq!(latest.committed_radius, 4);

        std::fs::remove_dir_all(&dir).unwrap();
    }
}
