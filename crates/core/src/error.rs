//! Error taxonomy for classdiam-core.

use thiserror::Error;

/// All fallible operations in the core return this error type.
#[derive(Debug, Error)]
pub enum ClassdiamError {
    /// A cycle-type template does not fit in the requested `n`
    /// (sum of non-fixed-point parts exceeds `n`).
    #[error("cycle type {template:?} does not fit in S_{n} (needs n >= {min_n})")]
    TemplateDoesNotFit {
        template: Vec<u8>,
        n: u16,
        min_n: u16,
    },

    /// The identity class was supplied as a generator without
    /// `allow_identity_generator` (spec §5.3).
    #[error(
        "identity class rejected as generator for S_{n}; set allow_identity_generator to permit it"
    )]
    IdentityGenerator { n: u16 },

    /// An empty generating set (no classes after validation).
    #[error("empty generating set for S_{n}")]
    EmptyUnion { n: u16 },

    /// A cycle-type template contains a zero part or is otherwise malformed.
    #[error("malformed cycle-type template {template:?}: {reason}")]
    MalformedTemplate { template: Vec<u8>, reason: String },

    /// `n` outside the supported range (2 ..= 255; partition parts are `u8`).
    #[error("n = {n} outside supported range 1..=255")]
    UnsupportedN { n: u32 },

    /// Central-eigenvalue integrality violated: `f_ρ ∤ |C_λ|·χ^ρ(λ)`.
    /// Almost always an indexing/orientation bug (spec §4).
    #[error("divisibility assertion failed: degree of rho #{rho} does not divide |C|*chi for lambda #{lambda}")]
    OmegaNotIntegral { rho: usize, lambda: usize },

    /// Exact division by `n!` left a remainder (spec §9.4).
    #[error("numerator not divisible by n! at radius {radius}, target #{target}")]
    NotDivisibleByFactorial { radius: u32, target: usize },

    /// A coefficient that must be a nonnegative integer came out negative
    /// (spec §9.5).
    #[error("negative coefficient at radius {radius}, target #{target}")]
    NegativeCoefficient { radius: u32, target: usize },

    /// The always-on word-count identity `Σ_ν |C_ν|·a_r(ν) = |U|^r` failed
    /// (spec §9.3) — silent arithmetic corruption.
    #[error("word-count identity violated at radius {radius}")]
    WordCountMismatch { radius: u32 },

    /// Diagnostic safety abort: the radius loop exceeded its configured
    /// bound without terminating (spec's stopping rule should always fire
    /// far earlier).
    #[error("radius limit {limit} exceeded without termination (n = {n})")]
    RadiusLimitExceeded { n: u16, limit: u32 },

    /// A checkpoint file is malformed (bad magic/version/CRC/truncation).
    #[error("invalid checkpoint: {reason}")]
    CheckpointFormat { reason: String },

    /// A checkpoint belongs to a different configuration or partition-order
    /// version; resuming would be mathematically meaningless (spec §19.3).
    #[error("checkpoint mismatch: {what}")]
    CheckpointMismatch { what: String },

    #[error("i/o error: {0}")]
    Io(#[from] std::io::Error),
}
