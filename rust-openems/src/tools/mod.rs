//! Utility functions and helper modules.
//!
//! This module provides various utility functions used throughout the codebase.

pub mod signal;

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
