//! Compressed FDTD Engine - Memory-bandwidth optimized implementation.
//!
//! This engine uses compressed coefficient storage to reduce memory bandwidth
//! during FDTD updates. Instead of loading full coefficient arrays, it uses
//! index-based lookup into small coefficient tables that fit in cache.
//!
//! This is a Rust port of the C++ `engine_sse_compressed` implementation.

use crate::arrays::{Dimensions, VectorField3D};
use crate::Result;
use rayon::prelude::*;
use wide::f32x8;

use super::compressed::{CompressedECoefficients, CompressedHCoefficients};
use super::operator::Operator;

/// Wrapper for raw pointer to make it Send + Sync for parallel iteration.
#[derive(Copy, Clone)]
struct SendPtr<T>(*const T);

unsafe impl<T> Send for SendPtr<T> {}
unsafe impl<T> Sync for SendPtr<T> {}

impl<T> SendPtr<T> {
    #[inline]
    fn new(ptr: *const T) -> Self {
        Self(ptr)
    }

    #[inline]
    unsafe fn add(&self, offset: usize) -> *const T {
        self.0.add(offset)
    }
}

/// Mutable version of SendPtr.
#[derive(Copy, Clone)]
struct SendPtrMut<T>(*mut T);

unsafe impl<T> Send for SendPtrMut<T> {}
unsafe impl<T> Sync for SendPtrMut<T> {}

impl<T> SendPtrMut<T> {
    #[inline]
    fn new(ptr: *mut T) -> Self {
        Self(ptr)
    }

    #[inline]
    unsafe fn add(&self, offset: usize) -> *mut T {
        self.0.add(offset)
    }
}

/// Load 8 f32 values from a slice into a SIMD vector.
#[inline]
fn load_f32x8(slice: &[f32]) -> f32x8 {
    let arr: [f32; 8] = slice[..8].try_into().unwrap();
    f32x8::from(arr)
}

/// Store a SIMD vector to a slice.
#[inline]
fn store_f32x8(slice: &mut [f32], v: f32x8) {
    slice[..8].copy_from_slice(v.as_array_ref());
}

/// Compressed FDTD Engine using index-based coefficient lookup.
///
/// This engine stores compressed coefficient tables and index arrays,
/// which reduces memory bandwidth by:
/// 1. Keeping coefficient tables small enough to fit in L1/L2 cache
/// 2. Loading only indices during field updates
/// 3. Benefiting from cache reuse when many cells share the same coefficients
pub struct EngineCompressed {
    /// Electric field components (Ex, Ey, Ez)
    e_field: VectorField3D,
    /// Magnetic field components (Hx, Hy, Hz)
    h_field: VectorField3D,
    /// Current timestep
    timestep: u64,
    /// Grid dimensions
    dims: Dimensions,
    /// Compressed E-field coefficients
    e_coeffs: CompressedECoefficients,
    /// Compressed H-field coefficients
    h_coeffs: CompressedHCoefficients,
}

impl EngineCompressed {
    /// Create a new compressed engine from an operator.
    ///
    /// This compresses the operator's coefficient arrays into lookup tables,
    /// potentially reducing memory usage for simulations with uniform regions.
    pub fn new(operator: &Operator) -> Self {
        let dims = operator.dimensions();
        let e_coeff = operator.e_coefficients();
        let h_coeff = operator.h_coefficients();

        // Compress coefficients
        let e_coeffs = CompressedECoefficients::from_fields(
            &[&e_coeff.ca[0], &e_coeff.ca[1], &e_coeff.ca[2]],
            &[&e_coeff.cb[0], &e_coeff.cb[1], &e_coeff.cb[2]],
        );
        let h_coeffs = CompressedHCoefficients::from_fields(
            &[&h_coeff.da[0], &h_coeff.da[1], &h_coeff.da[2]],
            &[&h_coeff.db[0], &h_coeff.db[1], &h_coeff.db[2]],
        );

        Self {
            e_field: VectorField3D::new(dims),
            h_field: VectorField3D::new(dims),
            timestep: 0,
            dims,
            e_coeffs,
            h_coeffs,
        }
    }

    /// Create compressed engine with pre-compressed coefficients.
    pub fn with_compressed(
        dims: Dimensions,
        e_coeffs: CompressedECoefficients,
        h_coeffs: CompressedHCoefficients,
    ) -> Self {
        Self {
            e_field: VectorField3D::new(dims),
            h_field: VectorField3D::new(dims),
            timestep: 0,
            dims,
            e_coeffs,
            h_coeffs,
        }
    }

    /// Get current timestep number.
    #[inline]
    pub fn timestep(&self) -> u64 {
        self.timestep
    }

    /// Get reference to E-field.
    #[inline]
    pub fn e_field(&self) -> &VectorField3D {
        &self.e_field
    }

    /// Get mutable reference to E-field.
    #[inline]
    pub fn e_field_mut(&mut self) -> &mut VectorField3D {
        &mut self.e_field
    }

    /// Get reference to H-field.
    #[inline]
    pub fn h_field(&self) -> &VectorField3D {
        &self.h_field
    }

    /// Get mutable reference to H-field.
    #[inline]
    pub fn h_field_mut(&mut self) -> &mut VectorField3D {
        &mut self.h_field
    }

    /// Compute total electromagnetic energy.
    pub fn total_energy(&self) -> f64 {
        let e_energy = self.e_field.energy();
        let h_energy = self.h_field.energy();
        e_energy + h_energy
    }

    /// Get compression statistics.
    pub fn compression_stats(&self) -> (f64, f64) {
        (self.e_coeffs.compression_ratio(), self.h_coeffs.compression_ratio())
    }

    /// Get number of unique coefficient sets.
    pub fn num_unique_coefficients(&self) -> ([usize; 3], [usize; 3]) {
        (self.e_coeffs.num_unique(), self.h_coeffs.num_unique())
    }

    /// Perform one complete FDTD timestep using compressed coefficients.
    pub fn step(&mut self, _operator: &Operator) -> Result<()> {
        // Use the parallel SIMD implementation for best performance
        self.step_parallel_simd();
        self.timestep += 1;
        Ok(())
    }

    /// Parallel SIMD implementation with z-line processing.
    fn step_parallel_simd(&mut self) {
        // H-field update
        self.update_h_parallel_simd();
        // E-field update
        self.update_e_parallel_simd();
    }

    /// SIMD-accelerated z-line field update.
    #[inline]
    fn update_line_simd(field: &mut [f32], curl: &[f32], coeff_a: &[f32], coeff_b: &[f32]) {
        let n = field.len();
        let chunks = n / 8;

        // SIMD processing (8 elements at a time)
        for c in 0..chunks {
            let k = c * 8;
            let f = load_f32x8(&field[k..]);
            let curl_v = load_f32x8(&curl[k..]);
            let a = load_f32x8(&coeff_a[k..]);
            let b = load_f32x8(&coeff_b[k..]);

            let result = a * f + b * curl_v;
            store_f32x8(&mut field[k..], result);
        }

        // Scalar tail
        for k in (chunks * 8)..n {
            field[k] = coeff_a[k] * field[k] + coeff_b[k] * curl[k];
        }
    }

    /// Parallel H-field update with SIMD z-line processing.
    fn update_h_parallel_simd(&mut self) {
        let dims = self.dims;
        let nx = dims.nx;
        let ny = dims.ny;
        let nz = dims.nz;

        // Get raw pointers wrapped for thread safety
        let hx_ptr = SendPtrMut::new(self.h_field.x.as_mut_ptr());
        let hy_ptr = SendPtrMut::new(self.h_field.y.as_mut_ptr());
        let hz_ptr = SendPtrMut::new(self.h_field.z.as_mut_ptr());

        let ex_ptr = SendPtr::new(self.e_field.x.as_ptr());
        let ey_ptr = SendPtr::new(self.e_field.y.as_ptr());
        let ez_ptr = SendPtr::new(self.e_field.z.as_ptr());

        // Coefficient table pointers
        let h_table_x_ptr = SendPtr::new(self.h_coeffs.table(0).as_ptr());
        let h_table_y_ptr = SendPtr::new(self.h_coeffs.table(1).as_ptr());
        let h_table_z_ptr = SendPtr::new(self.h_coeffs.table(2).as_ptr());

        let h_indices_x_ptr = SendPtr::new(self.h_coeffs.indices(0).as_ptr());
        let h_indices_y_ptr = SendPtr::new(self.h_coeffs.indices(1).as_ptr());
        let h_indices_z_ptr = SendPtr::new(self.h_coeffs.indices(2).as_ptr());

        // Parallel iteration over i-slices
        (0..nx).into_par_iter().for_each(|i| {
            // Each thread gets its own temporary buffers
            let mut da_buf = vec![0.0f32; nz];
            let mut db_buf = vec![0.0f32; nz];
            let mut curl_buf = vec![0.0f32; nz];

            for j in 0..ny {
                // Process Hx: curl_x = dEz/dy - dEy/dz
                // Expand coefficients for this z-line
                for k in 0..nz {
                    let idx = dims.to_linear(i, j, k);
                    let table_idx = unsafe { *h_indices_x_ptr.add(idx) } as usize;
                    let coeff = unsafe { *h_table_x_ptr.add(table_idx) };
                    da_buf[k] = coeff.da;
                    db_buf[k] = coeff.db;

                    // Compute curl
                    let ez_curr = unsafe { *ez_ptr.add(idx) };
                    let ez_jp1 = if j + 1 < ny {
                        unsafe { *ez_ptr.add(dims.to_linear(i, j + 1, k)) }
                    } else {
                        ez_curr
                    };
                    let ey_curr = unsafe { *ey_ptr.add(idx) };
                    let ey_kp1 = if k + 1 < nz {
                        unsafe { *ey_ptr.add(dims.to_linear(i, j, k + 1)) }
                    } else {
                        ey_curr
                    };
                    curl_buf[k] = (ez_jp1 - ez_curr) - (ey_kp1 - ey_curr);
                }

                // Apply SIMD update for Hx
                let hx_start = dims.to_linear(i, j, 0);
                let hx_line = unsafe { std::slice::from_raw_parts_mut(hx_ptr.add(hx_start), nz) };
                Self::update_line_simd(hx_line, &curl_buf, &da_buf, &db_buf);

                // Process Hy: curl_y = dEx/dz - dEz/dx
                for k in 0..nz {
                    let idx = dims.to_linear(i, j, k);
                    let table_idx = unsafe { *h_indices_y_ptr.add(idx) } as usize;
                    let coeff = unsafe { *h_table_y_ptr.add(table_idx) };
                    da_buf[k] = coeff.da;
                    db_buf[k] = coeff.db;

                    let ex_curr = unsafe { *ex_ptr.add(idx) };
                    let ex_kp1 = if k + 1 < nz {
                        unsafe { *ex_ptr.add(dims.to_linear(i, j, k + 1)) }
                    } else {
                        ex_curr
                    };
                    let ez_curr = unsafe { *ez_ptr.add(idx) };
                    let ez_ip1 = if i + 1 < nx {
                        unsafe { *ez_ptr.add(dims.to_linear(i + 1, j, k)) }
                    } else {
                        ez_curr
                    };
                    curl_buf[k] = (ex_kp1 - ex_curr) - (ez_ip1 - ez_curr);
                }

                let hy_start = dims.to_linear(i, j, 0);
                let hy_line = unsafe { std::slice::from_raw_parts_mut(hy_ptr.add(hy_start), nz) };
                Self::update_line_simd(hy_line, &curl_buf, &da_buf, &db_buf);

                // Process Hz: curl_z = dEy/dx - dEx/dy
                for k in 0..nz {
                    let idx = dims.to_linear(i, j, k);
                    let table_idx = unsafe { *h_indices_z_ptr.add(idx) } as usize;
                    let coeff = unsafe { *h_table_z_ptr.add(table_idx) };
                    da_buf[k] = coeff.da;
                    db_buf[k] = coeff.db;

                    let ey_curr = unsafe { *ey_ptr.add(idx) };
                    let ey_ip1 = if i + 1 < nx {
                        unsafe { *ey_ptr.add(dims.to_linear(i + 1, j, k)) }
                    } else {
                        ey_curr
                    };
                    let ex_curr = unsafe { *ex_ptr.add(idx) };
                    let ex_jp1 = if j + 1 < ny {
                        unsafe { *ex_ptr.add(dims.to_linear(i, j + 1, k)) }
                    } else {
                        ex_curr
                    };
                    curl_buf[k] = (ey_ip1 - ey_curr) - (ex_jp1 - ex_curr);
                }

                let hz_start = dims.to_linear(i, j, 0);
                let hz_line = unsafe { std::slice::from_raw_parts_mut(hz_ptr.add(hz_start), nz) };
                Self::update_line_simd(hz_line, &curl_buf, &da_buf, &db_buf);
            }
        });
    }

    /// Parallel E-field update with SIMD z-line processing.
    fn update_e_parallel_simd(&mut self) {
        let dims = self.dims;
        let nx = dims.nx;
        let ny = dims.ny;
        let nz = dims.nz;

        // Get raw pointers wrapped for thread safety
        let ex_ptr = SendPtrMut::new(self.e_field.x.as_mut_ptr());
        let ey_ptr = SendPtrMut::new(self.e_field.y.as_mut_ptr());
        let ez_ptr = SendPtrMut::new(self.e_field.z.as_mut_ptr());

        let hx_ptr = SendPtr::new(self.h_field.x.as_ptr());
        let hy_ptr = SendPtr::new(self.h_field.y.as_ptr());
        let hz_ptr = SendPtr::new(self.h_field.z.as_ptr());

        // Coefficient table pointers
        let e_table_x_ptr = SendPtr::new(self.e_coeffs.table(0).as_ptr());
        let e_table_y_ptr = SendPtr::new(self.e_coeffs.table(1).as_ptr());
        let e_table_z_ptr = SendPtr::new(self.e_coeffs.table(2).as_ptr());

        let e_indices_x_ptr = SendPtr::new(self.e_coeffs.indices(0).as_ptr());
        let e_indices_y_ptr = SendPtr::new(self.e_coeffs.indices(1).as_ptr());
        let e_indices_z_ptr = SendPtr::new(self.e_coeffs.indices(2).as_ptr());

        // Parallel iteration starting from i=1
        (1..nx).into_par_iter().for_each(|i| {
            let mut ca_buf = vec![0.0f32; nz];
            let mut cb_buf = vec![0.0f32; nz];
            let mut curl_buf = vec![0.0f32; nz];

            for j in 1..ny {
                // Process Ex: curl = dHz/dy - dHy/dz
                for k in 1..nz {
                    let idx = dims.to_linear(i, j, k);
                    let table_idx = unsafe { *e_indices_x_ptr.add(idx) } as usize;
                    let coeff = unsafe { *e_table_x_ptr.add(table_idx) };
                    ca_buf[k] = coeff.ca;
                    cb_buf[k] = coeff.cb;

                    let hz_curr = unsafe { *hz_ptr.add(idx) };
                    let hz_jm1 = unsafe { *hz_ptr.add(dims.to_linear(i, j - 1, k)) };
                    let hy_curr = unsafe { *hy_ptr.add(idx) };
                    let hy_km1 = unsafe { *hy_ptr.add(dims.to_linear(i, j, k - 1)) };
                    curl_buf[k] = (hz_curr - hz_jm1) - (hy_curr - hy_km1);
                }
                ca_buf[0] = 0.0;
                cb_buf[0] = 0.0;
                curl_buf[0] = 0.0;

                let ex_start = dims.to_linear(i, j, 0);
                let ex_line = unsafe { std::slice::from_raw_parts_mut(ex_ptr.add(ex_start), nz) };
                Self::update_line_simd(ex_line, &curl_buf, &ca_buf, &cb_buf);

                // Process Ey: curl = dHx/dz - dHz/dx
                for k in 1..nz {
                    let idx = dims.to_linear(i, j, k);
                    let table_idx = unsafe { *e_indices_y_ptr.add(idx) } as usize;
                    let coeff = unsafe { *e_table_y_ptr.add(table_idx) };
                    ca_buf[k] = coeff.ca;
                    cb_buf[k] = coeff.cb;

                    let hx_curr = unsafe { *hx_ptr.add(idx) };
                    let hx_km1 = unsafe { *hx_ptr.add(dims.to_linear(i, j, k - 1)) };
                    let hz_curr = unsafe { *hz_ptr.add(idx) };
                    let hz_im1 = unsafe { *hz_ptr.add(dims.to_linear(i - 1, j, k)) };
                    curl_buf[k] = (hx_curr - hx_km1) - (hz_curr - hz_im1);
                }
                ca_buf[0] = 0.0;
                cb_buf[0] = 0.0;
                curl_buf[0] = 0.0;

                let ey_start = dims.to_linear(i, j, 0);
                let ey_line = unsafe { std::slice::from_raw_parts_mut(ey_ptr.add(ey_start), nz) };
                Self::update_line_simd(ey_line, &curl_buf, &ca_buf, &cb_buf);

                // Process Ez: curl = dHy/dx - dHx/dy
                for k in 1..nz {
                    let idx = dims.to_linear(i, j, k);
                    let table_idx = unsafe { *e_indices_z_ptr.add(idx) } as usize;
                    let coeff = unsafe { *e_table_z_ptr.add(table_idx) };
                    ca_buf[k] = coeff.ca;
                    cb_buf[k] = coeff.cb;

                    let hy_curr = unsafe { *hy_ptr.add(idx) };
                    let hy_im1 = unsafe { *hy_ptr.add(dims.to_linear(i - 1, j, k)) };
                    let hx_curr = unsafe { *hx_ptr.add(idx) };
                    let hx_jm1 = unsafe { *hx_ptr.add(dims.to_linear(i, j - 1, k)) };
                    curl_buf[k] = (hy_curr - hy_im1) - (hx_curr - hx_jm1);
                }
                ca_buf[0] = 0.0;
                cb_buf[0] = 0.0;
                curl_buf[0] = 0.0;

                let ez_start = dims.to_linear(i, j, 0);
                let ez_line = unsafe { std::slice::from_raw_parts_mut(ez_ptr.add(ez_start), nz) };
                Self::update_line_simd(ez_line, &curl_buf, &ca_buf, &cb_buf);
            }
        });
    }

    /// Reset fields to zero.
    pub fn reset(&mut self) {
        self.e_field.clear();
        self.h_field.clear();
        self.timestep = 0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fdtd::BoundaryConditions;
    use crate::geometry::{CoordinateSystem, Grid};

    #[test]
    fn test_compressed_engine_creation() {
        let grid = Grid::new(
            CoordinateSystem::Cartesian,
            vec![0.0, 0.001, 0.002, 0.003],
            vec![0.0, 0.001, 0.002, 0.003],
            vec![0.0, 0.001, 0.002, 0.003],
        );
        let op = Operator::new(grid, BoundaryConditions::default()).unwrap();
        let engine = EngineCompressed::new(&op);

        assert_eq!(engine.timestep(), 0);
    }

    #[test]
    fn test_compressed_engine_step() {
        let grid = Grid::new(
            CoordinateSystem::Cartesian,
            (0..11).map(|i| i as f64 * 0.001).collect(),
            (0..11).map(|i| i as f64 * 0.001).collect(),
            (0..11).map(|i| i as f64 * 0.001).collect(),
        );
        let op = Operator::new(grid, BoundaryConditions::all_pec()).unwrap();
        let mut engine = EngineCompressed::new(&op);

        // Add initial energy
        engine.e_field.z.set(5, 5, 5, 1.0);

        let initial_energy = engine.total_energy();
        assert!(initial_energy > 0.0);

        // Run some steps
        for _ in 0..100 {
            engine.step(&op).unwrap();
        }

        let final_energy = engine.total_energy();

        // Energy should be approximately conserved
        let relative_diff = (final_energy - initial_energy).abs() / initial_energy;
        assert!(
            relative_diff < 0.1,
            "Energy not conserved: initial={}, final={}, diff={}",
            initial_energy,
            final_energy,
            relative_diff
        );
    }

    #[test]
    fn test_compression_benefit() {
        // Create a uniform grid - should compress well
        let grid = Grid::uniform(50, 50, 50, 0.001);
        let op = Operator::new(grid, BoundaryConditions::all_pec()).unwrap();
        let engine = EngineCompressed::new(&op);

        let (e_ratio, h_ratio) = engine.compression_stats();
        let (e_unique, h_unique) = engine.num_unique_coefficients();

        println!("E compression ratio: {:.2}, unique: {:?}", e_ratio, e_unique);
        println!("H compression ratio: {:.2}, unique: {:?}", h_ratio, h_unique);

        assert!(e_ratio > 0.0);
        assert!(h_ratio > 0.0);
    }
}
