//! Cycle-type templates: cycle types written WITHOUT fixed points, padded
//! per `n`, so one generator spec runs across a whole range of `n`
//! (project decision; e.g. `[3]` = the 3-cycle class in every `S_n`, n ≥ 3).
#![deny(clippy::float_arithmetic)]

use std::fmt;
use std::str::FromStr;

use serde::{Deserialize, Serialize};
use smallvec::SmallVec;

use super::Partition;
use crate::error::ClassdiamError;

/// A cycle type with fixed points omitted: weakly decreasing parts, all ≥ 2.
/// The empty template denotes the identity class (rejected as a generator by
/// default at engine level, spec §5.3).
#[derive(Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(try_from = "Vec<u8>", into = "Vec<u8>")]
pub struct CycleTypeTemplate {
    parts: SmallVec<[u8; 8]>,
}

impl CycleTypeTemplate {
    /// Build from parts in any order. Parts of 0 or 1 are rejected: fixed
    /// points must be omitted (they are implied by padding), so a part of 1
    /// is always a spec error rather than a convention choice.
    pub fn new(parts: impl AsRef<[u8]>) -> Result<Self, ClassdiamError> {
        let raw = parts.as_ref();
        if let Some(&bad) = raw.iter().find(|&&p| p < 2) {
            return Err(ClassdiamError::MalformedTemplate {
                template: raw.to_vec(),
                reason: format!(
                    "part {bad} is not allowed; write cycle types without fixed points (all parts >= 2)"
                ),
            });
        }
        let mut parts: SmallVec<[u8; 8]> = SmallVec::from_slice(raw);
        parts.sort_unstable_by(|a, b| b.cmp(a));
        if parts.iter().map(|&p| p as u32).sum::<u32>() > 255 {
            return Err(ClassdiamError::MalformedTemplate {
                template: raw.to_vec(),
                reason: "sum of parts exceeds 255".into(),
            });
        }
        Ok(Self { parts })
    }

    /// The identity-class template (empty).
    pub fn identity() -> Self {
        Self {
            parts: SmallVec::new(),
        }
    }

    pub fn parts(&self) -> &[u8] {
        &self.parts
    }

    /// `true` for the empty template (identity class).
    pub fn is_identity(&self) -> bool {
        self.parts.is_empty()
    }

    /// Smallest `n` this template fits in. The identity template needs n ≥ 0
    /// but is only meaningful for n ≥ 1.
    pub fn min_n(&self) -> u16 {
        self.parts.iter().map(|&p| p as u16).sum()
    }

    /// Pad with fixed points to a full cycle type of `S_n`.
    pub fn pad(&self, n: u16) -> Result<Partition, ClassdiamError> {
        let sum = self.min_n();
        if sum > n {
            return Err(ClassdiamError::TemplateDoesNotFit {
                template: self.parts.to_vec(),
                n,
                min_n: sum,
            });
        }
        let mut full: SmallVec<[u8; 16]> = SmallVec::from_slice(&self.parts);
        full.extend(std::iter::repeat(1u8).take((n - sum) as usize));
        Ok(Partition::new(full))
    }

    /// Slug form for file names: parts joined by `.`; identity = `id`.
    pub fn slug(&self) -> String {
        if self.is_identity() {
            "id".to_string()
        } else {
            self.parts
                .iter()
                .map(|p| p.to_string())
                .collect::<Vec<_>>()
                .join(".")
        }
    }
}

impl fmt::Debug for CycleTypeTemplate {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self.parts.as_slice())
    }
}

impl fmt::Display for CycleTypeTemplate {
    /// Human/CLI form: parts joined by `,` (e.g. `2,2`); identity = empty string.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}",
            self.parts
                .iter()
                .map(|p| p.to_string())
                .collect::<Vec<_>>()
                .join(",")
        )
    }
}

impl FromStr for CycleTypeTemplate {
    type Err = ClassdiamError;

    /// Parse the CLI form: comma-separated parts, whitespace tolerated.
    /// Empty string = identity template.
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let trimmed = s.trim();
        if trimmed.is_empty() {
            return Ok(Self::identity());
        }
        let mut parts = Vec::new();
        for token in trimmed.split(',') {
            let token = token.trim();
            let value: u8 = token.parse().map_err(|_| ClassdiamError::MalformedTemplate {
                template: Vec::new(),
                reason: format!("cannot parse part {token:?} in {s:?}"),
            })?;
            parts.push(value);
        }
        Self::new(&parts)
    }
}

impl TryFrom<Vec<u8>> for CycleTypeTemplate {
    type Error = ClassdiamError;
    fn try_from(parts: Vec<u8>) -> Result<Self, Self::Error> {
        Self::new(&parts)
    }
}

impl From<CycleTypeTemplate> for Vec<u8> {
    fn from(t: CycleTypeTemplate) -> Self {
        t.parts.to_vec()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn padding() {
        let t = CycleTypeTemplate::new([3]).unwrap();
        assert_eq!(t.pad(5).unwrap().parts(), &[3, 1, 1]);
        assert_eq!(t.pad(3).unwrap().parts(), &[3]);
        assert!(matches!(
            t.pad(2),
            Err(ClassdiamError::TemplateDoesNotFit { min_n: 3, n: 2, .. })
        ));

        let t = CycleTypeTemplate::new([2, 2]).unwrap();
        assert_eq!(t.pad(4).unwrap().parts(), &[2, 2]);
        assert_eq!(t.pad(6).unwrap().parts(), &[2, 2, 1, 1]);
    }

    #[test]
    fn identity_template() {
        let t = CycleTypeTemplate::identity();
        assert!(t.is_identity());
        assert_eq!(t.min_n(), 0);
        assert!(t.pad(4).unwrap().is_identity_type());
        assert_eq!(t.slug(), "id");
    }

    #[test]
    fn rejects_fixed_points_and_zeros() {
        assert!(CycleTypeTemplate::new([2, 1]).is_err());
        assert!(CycleTypeTemplate::new([1]).is_err());
        assert!(CycleTypeTemplate::new([0]).is_err());
    }

    #[test]
    fn sorts_input() {
        let t = CycleTypeTemplate::new([2, 4, 3]).unwrap();
        assert_eq!(t.parts(), &[4, 3, 2]);
    }

    #[test]
    fn parse_and_display() {
        let t: CycleTypeTemplate = "2,2".parse().unwrap();
        assert_eq!(t.parts(), &[2, 2]);
        assert_eq!(t.to_string(), "2,2");
        assert_eq!(t.slug(), "2.2");

        let t: CycleTypeTemplate = " 3 , 2 ".parse().unwrap();
        assert_eq!(t.parts(), &[3, 2]);

        let t: CycleTypeTemplate = "".parse().unwrap();
        assert!(t.is_identity());

        assert!("2,x".parse::<CycleTypeTemplate>().is_err());
        assert!("1".parse::<CycleTypeTemplate>().is_err());
    }

    #[test]
    fn serde_roundtrip() {
        let t = CycleTypeTemplate::new([3, 2]).unwrap();
        let json = serde_json::to_string(&t).unwrap();
        assert_eq!(json, "[3,2]");
        let back: CycleTypeTemplate = serde_json::from_str(&json).unwrap();
        assert_eq!(back, t);
        // rejects invalid on deserialize
        assert!(serde_json::from_str::<CycleTypeTemplate>("[1]").is_err());
    }
}
