//! Job orchestration shared by every front end (CLI, Python).
//!
//! The math engines ([`crate::engine`]) compute one `(n, union)` result;
//! this module owns everything around them: run/job configuration hashes,
//! result/checkpoint/manifest file layout, batch and resume loops, progress
//! hooks, and cooperative cancellation. Front ends supply argument parsing
//! and presentation only, so a run directory produced by one front end is
//! always resumable by another.

mod cancel;

pub use cancel::CancelToken;
