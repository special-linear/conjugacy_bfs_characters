//! Test-only oracles: deliberately naive, obviously-correct reference
//! implementations that the optimized production code is differentially
//! tested against — forever (design doc `03-testing` §1).
//!
//! Available to this crate's tests always, and to other crates' tests via
//! the `test-utils` feature. Never part of a release build's public surface
//! in spirit; nothing here is performance-tuned.

pub mod bruteforce;
pub mod catalog;
pub mod naive_mn;
