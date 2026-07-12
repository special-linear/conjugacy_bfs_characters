//! Cooperative cancellation for long runs (design doc 01 §12).
//!
//! A [`CancelToken`] is a shared flag checked by the modular engine at the
//! same between-radii point as the wall-clock deadline; cancellation
//! therefore produces the same clean, fully-certified checkpointed
//! suspension as an expired deadline. The exact engine does not observe
//! tokens (it is a small-n oracle whose runs take seconds).

use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;

/// Shared cancellation flag; clones observe the same flag.
#[derive(Clone, Debug, Default)]
pub struct CancelToken(Arc<AtomicBool>);

impl CancelToken {
    pub fn new() -> Self {
        Self::default()
    }

    /// Request cancellation. Idempotent; safe from any thread (e.g. a
    /// signal handler or a Python callback).
    pub fn cancel(&self) {
        self.0.store(true, Ordering::Release);
    }

    pub fn is_cancelled(&self) -> bool {
        self.0.load(Ordering::Acquire)
    }
}
