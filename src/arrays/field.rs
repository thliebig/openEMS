//! Field storage types for FDTD simulations.

use super::{Dimensions, SIMD_ALIGN};
use aligned_vec::{AVec, ConstAlign};
use rayon::prelude::*;
use std::ops::{Index, IndexMut};

/// Field component identifier
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum FieldComponent {
    /// Electric field x-component
    Ex,
    /// Electric field y-component
    Ey,
    /// Electric field z-component
    Ez,
    /// Magnetic field x-component
    Hx,
    /// Magnetic field y-component
    Hy,
    /// Magnetic field z-component
    Hz,
}

impl FieldComponent {
    /// Check if this is an electric field component
    #[inline]
    pub fn is_electric(&self) -> bool {
        matches!(self, Self::Ex | Self::Ey | Self::Ez)
    }

    /// Check if this is a magnetic field component
    #[inline]
    pub fn is_magnetic(&self) -> bool {
        matches!(self, Self::Hx | Self::Hy | Self::Hz)
    }

    /// Get the direction index (0=x, 1=y, 2=z)
    #[inline]
    pub fn direction(&self) -> usize {
        match self {
            Self::Ex | Self::Hx => 0,
            Self::Ey | Self::Hy => 1,
            Self::Ez | Self::Hz => 2,
        }
    }
}

/// 3D scalar field storage with SIMD-aligned memory.
///
/// Memory is laid out contiguously with z-index varying fastest,
/// which is optimal for vectorization along the z-axis.
#[derive(Clone)]
pub struct Field3D {
    /// Field data (SIMD-aligned)
    data: AVec<f32, ConstAlign<SIMD_ALIGN>>,
    /// Grid dimensions
    dims: Dimensions,
}

impl std::fmt::Debug for Field3D {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Field3D")
            .field("dims", &self.dims)
            .field("data_len", &self.data.len())
            .finish()
    }
}

// Safety: f32 is Pod and we're using aligned allocation
unsafe impl Send for Field3D {}
unsafe impl Sync for Field3D {}

impl Field3D {
    /// Create a new field initialized to zero.
    pub fn new(dims: Dimensions) -> Self {
        let total = dims.total();
        let mut data = AVec::new(SIMD_ALIGN);
        data.resize(total, 0.0f32);
        Self { data, dims }
    }

    /// Create a new field from existing data.
    ///
    /// # Panics
    /// Panics if data length doesn't match dimensions.
    pub fn from_vec(dims: Dimensions, vec: Vec<f32>) -> Self {
        assert_eq!(vec.len(), dims.total());
        let mut data = AVec::new(SIMD_ALIGN);
        data.resize(vec.len(), 0.0);
        data.copy_from_slice(&vec);
        Self { data, dims }
    }

    /// Get the field dimensions.
    #[inline]
    pub fn dims(&self) -> Dimensions {
        self.dims
    }

    /// Get the total number of cells.
    #[inline]
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Check if field is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Get immutable access to the underlying data.
    #[inline]
    pub fn as_slice(&self) -> &[f32] {
        &self.data
    }

    /// Get mutable access to the underlying data.
    #[inline]
    pub fn as_mut_slice(&mut self) -> &mut [f32] {
        &mut self.data
    }

    /// Get a pointer to the data (for SIMD operations).
    #[inline]
    pub fn as_ptr(&self) -> *const f32 {
        self.data.as_ptr()
    }

    /// Get a mutable pointer to the data.
    #[inline]
    pub fn as_mut_ptr(&mut self) -> *mut f32 {
        self.data.as_mut_ptr()
    }

    /// Get value at 3D index.
    #[inline]
    pub fn get(&self, i: usize, j: usize, k: usize) -> f32 {
        self.data[self.dims.to_linear(i, j, k)]
    }

    /// Set value at 3D index.
    #[inline]
    pub fn set(&mut self, i: usize, j: usize, k: usize, value: f32) {
        let idx = self.dims.to_linear(i, j, k);
        self.data[idx] = value;
    }

    /// Add value at 3D index.
    #[inline]
    pub fn add(&mut self, i: usize, j: usize, k: usize, value: f32) {
        let idx = self.dims.to_linear(i, j, k);
        self.data[idx] += value;
    }

    /// Fill entire field with a value.
    pub fn fill(&mut self, value: f32) {
        self.data.par_iter_mut().for_each(|v| *v = value);
    }

    /// Set all values to zero.
    pub fn clear(&mut self) {
        self.fill(0.0);
    }

    /// Get a slice along the z-axis at position (i, j).
    #[inline]
    pub fn z_slice(&self, i: usize, j: usize) -> &[f32] {
        let start = self.dims.to_linear(i, j, 0);
        let end = start + self.dims.nz;
        &self.data[start..end]
    }

    /// Get a mutable slice along the z-axis at position (i, j).
    #[inline]
    pub fn z_slice_mut(&mut self, i: usize, j: usize) -> &mut [f32] {
        let start = self.dims.to_linear(i, j, 0);
        let end = start + self.dims.nz;
        &mut self.data[start..end]
    }

    /// Apply a function to all elements in parallel.
    pub fn par_map_inplace<F>(&mut self, f: F)
    where
        F: Fn(f32) -> f32 + Sync + Send,
    {
        self.data.par_iter_mut().for_each(|v| *v = f(*v));
    }

    /// Element-wise addition with another field.
    pub fn add_field(&mut self, other: &Field3D) {
        assert_eq!(self.dims, other.dims);
        self.data
            .par_iter_mut()
            .zip(other.data.par_iter())
            .for_each(|(a, b)| *a += *b);
    }

    /// Element-wise multiplication by scalar.
    pub fn scale(&mut self, factor: f32) {
        self.data.par_iter_mut().for_each(|v| *v *= factor);
    }

    /// Compute the sum of all elements.
    pub fn sum(&self) -> f64 {
        self.data.par_iter().map(|&v| v as f64).sum()
    }

    /// Compute the maximum absolute value.
    pub fn max_abs(&self) -> f32 {
        self.data
            .par_iter()
            .map(|&v| v.abs())
            .reduce(|| 0.0f32, |a, b| a.max(b))
    }

    /// Compute the energy (sum of squares).
    pub fn energy(&self) -> f64 {
        self.data.par_iter().map(|&v| (v as f64) * (v as f64)).sum()
    }
}

impl Index<(usize, usize, usize)> for Field3D {
    type Output = f32;

    #[inline]
    fn index(&self, (i, j, k): (usize, usize, usize)) -> &Self::Output {
        &self.data[self.dims.to_linear(i, j, k)]
    }
}

impl IndexMut<(usize, usize, usize)> for Field3D {
    #[inline]
    fn index_mut(&mut self, (i, j, k): (usize, usize, usize)) -> &mut Self::Output {
        let idx = self.dims.to_linear(i, j, k);
        &mut self.data[idx]
    }
}

/// 3D vector field (stores all three components).
pub struct VectorField3D {
    /// X component
    pub x: Field3D,
    /// Y component
    pub y: Field3D,
    /// Z component
    pub z: Field3D,
}

impl std::fmt::Debug for VectorField3D {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("VectorField3D")
            .field("dims", &self.x.dims())
            .field("x_len", &self.x.len())
            .field("y_len", &self.y.len())
            .field("z_len", &self.z.len())
            .finish()
    }
}

impl VectorField3D {
    /// Create a new zero-initialized vector field.
    pub fn new(dims: Dimensions) -> Self {
        Self {
            x: Field3D::new(dims),
            y: Field3D::new(dims),
            z: Field3D::new(dims),
        }
    }

    /// Get the field dimensions.
    #[inline]
    pub fn dims(&self) -> Dimensions {
        self.x.dims()
    }

    /// Get a reference to a component by direction index.
    #[inline]
    pub fn component(&self, dir: usize) -> &Field3D {
        match dir {
            0 => &self.x,
            1 => &self.y,
            2 => &self.z,
            _ => panic!("Invalid direction index"),
        }
    }

    /// Get a mutable reference to a component by direction index.
    #[inline]
    pub fn component_mut(&mut self, dir: usize) -> &mut Field3D {
        match dir {
            0 => &mut self.x,
            1 => &mut self.y,
            2 => &mut self.z,
            _ => panic!("Invalid direction index"),
        }
    }

    /// Set all components to zero.
    pub fn clear(&mut self) {
        self.x.clear();
        self.y.clear();
        self.z.clear();
    }

    /// Set all components to zero (alias for clear).
    pub fn zero(&mut self) {
        self.clear();
    }

    /// Compute the total energy in the vector field.
    pub fn energy(&self) -> f64 {
        self.x.energy() + self.y.energy() + self.z.energy()
    }

    /// Get the maximum absolute value across all components.
    pub fn max_abs(&self) -> f32 {
        self.x.max_abs().max(self.y.max_abs()).max(self.z.max_abs())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_field3d_basic() {
        let dims = Dimensions::new(10, 20, 30);
        let mut field = Field3D::new(dims);

        assert_eq!(field.len(), 6000);
        assert_eq!(field.get(0, 0, 0), 0.0);

        field.set(5, 10, 15, 42.0);
        assert_eq!(field.get(5, 10, 15), 42.0);
        assert_eq!(field[(5, 10, 15)], 42.0);
    }

    #[test]
    fn test_field3d_slice() {
        let dims = Dimensions::new(3, 4, 5);
        let mut field = Field3D::new(dims);

        // Fill z-slice with sequential values
        let slice = field.z_slice_mut(1, 2);
        for (k, v) in slice.iter_mut().enumerate() {
            *v = k as f32;
        }

        // Verify
        for k in 0..5 {
            assert_eq!(field.get(1, 2, k), k as f32);
        }
    }

    #[test]
    fn test_vector_field() {
        let dims = Dimensions::new(5, 5, 5);
        let mut vf = VectorField3D::new(dims);

        vf.x.set(1, 1, 1, 1.0);
        vf.y.set(1, 1, 1, 2.0);
        vf.z.set(1, 1, 1, 3.0);

        // Energy = 1^2 + 2^2 + 3^2 = 14
        let e = vf.energy();
        assert!((e - 14.0).abs() < 1e-6);
    }
}
