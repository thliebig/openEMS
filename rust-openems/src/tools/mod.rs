//! Utility functions and helper modules.
//!
//! This module provides various utility functions used throughout the codebase:
//! - Signal handling for graceful shutdown
//! - CLI option parsing and configuration
//! - Denormal number handling for performance
//! - Address operations for array indexing

pub mod signal;
pub mod signal_handler;
pub mod cli;
pub mod denormal;

// Re-exports for convenience
pub use signal_handler::{SignalHandler, SignalGuard, InterruptGuard};
pub use cli::{GlobalOptions, CliParser, EnvConfig};
pub use denormal::{DenormalGuard, DenormalMode, DenormalMonitor, DenormalStats};

/// Address operation helpers (replacement for C++ AdrOp)
pub mod address {
    use crate::arrays::Dimensions;

    /// Calculate linear index from 3D coordinates.
    #[inline]
    pub fn to_linear(dims: &Dimensions, i: usize, j: usize, k: usize) -> usize {
        dims.to_linear(i, j, k)
    }

    /// Calculate 3D coordinates from linear index.
    #[inline]
    pub fn from_linear(dims: &Dimensions, idx: usize) -> (usize, usize, usize) {
        dims.from_linear(idx)
    }
}
