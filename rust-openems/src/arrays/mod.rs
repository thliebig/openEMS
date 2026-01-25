//! High-performance array types for FDTD field storage.
//!
//! This module provides memory-efficient, SIMD-friendly array types
//! for storing electromagnetic field components on the Yee grid.
//!
//! ## Array Types
//!
//! - [`Field3D`]: 3D field storage for a single component
//! - [`VectorField3D`]: 3D vector field (Ex, Ey, Ez or Hx, Hy, Hz)
//! - [`FieldComponent`]: Enum for field component selection
//!
//! ## Memory Layout
//!
//! Arrays use a contiguous memory layout with optional SIMD alignment
//! for optimal vectorization performance.

mod field;
pub mod simd;

pub use field::{Field3D, FieldComponent, VectorField3D};
pub use simd::SimdOps;

/// Alignment for SIMD operations (64 bytes = AVX-512 width)
pub const SIMD_ALIGN: usize = 64;

/// Index type for grid coordinates
pub type GridIndex = usize;

/// Grid dimensions (nx, ny, nz)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Dimensions {
    /// Number of cells in x direction
    pub nx: usize,
    /// Number of cells in y direction
    pub ny: usize,
    /// Number of cells in z direction
    pub nz: usize,
}

impl Dimensions {
    /// Create new dimensions
    #[inline]
    pub fn new(nx: usize, ny: usize, nz: usize) -> Self {
        Self { nx, ny, nz }
    }

    /// Total number of cells
    #[inline]
    pub fn total(&self) -> usize {
        self.nx * self.ny * self.nz
    }

    /// Check if index is valid
    #[inline]
    pub fn is_valid(&self, i: usize, j: usize, k: usize) -> bool {
        i < self.nx && j < self.ny && k < self.nz
    }

    /// Convert 3D index to linear index (row-major order: z varies fastest)
    #[inline]
    pub fn to_linear(&self, i: usize, j: usize, k: usize) -> usize {
        debug_assert!(self.is_valid(i, j, k));
        i * self.ny * self.nz + j * self.nz + k
    }

    /// Convert linear index to 3D index
    #[inline]
    pub fn from_linear(&self, idx: usize) -> (usize, usize, usize) {
        let k = idx % self.nz;
        let remaining = idx / self.nz;
        let j = remaining % self.ny;
        let i = remaining / self.ny;
        (i, j, k)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dimensions() {
        let dims = Dimensions::new(10, 20, 30);
        assert_eq!(dims.total(), 6000);

        // Test index conversion round-trip
        for i in 0..dims.nx {
            for j in 0..dims.ny {
                for k in 0..dims.nz {
                    let linear = dims.to_linear(i, j, k);
                    let (i2, j2, k2) = dims.from_linear(linear);
                    assert_eq!((i, j, k), (i2, j2, k2));
                }
            }
        }
    }
}
