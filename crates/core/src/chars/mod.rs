//! Character data for `S_n`: hook-length degrees and the Murnaghan–Nakayama
//! evaluator (spec §10).
#![deny(clippy::float_arithmetic)]

pub mod degrees;
pub mod mn;
pub mod rimhook;

pub use degrees::{degree, degrees, hook_lengths};
pub use mn::MnEvaluator;
