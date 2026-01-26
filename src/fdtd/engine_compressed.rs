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

        // Print compression stats
        let e_unique = e_coeffs.num_unique();
        let h_unique = h_coeffs.num_unique();
        let total_cells = dims.total();
        eprintln!(
            "Compressed engine: {} cells, E-coeffs unique: {:?}, H-coeffs unique: {:?}",
            total_cells, e_unique, h_unique
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

    /// Perform one complete FDTD timestep using compressed coefficients.
    pub fn step(&mut self) -> Result<()> {
        // Use the parallel implementation for best performance
        self.step_parallel();
        self.timestep += 1;
        Ok(())
    }

    /// Single-threaded implementation for reference/debugging.
    #[allow(dead_code)]
    fn step_basic(&mut self) {
        self.update_h_basic();
        self.update_e_basic();
    }

    /// Update H-field using compressed coefficients (basic implementation).
    fn update_h_basic(&mut self) {
        let dims = self.dims;

        for i in 0..dims.nx {
            for j in 0..dims.ny {
                for k in 0..dims.nz {
                    // Get compressed coefficients via index lookup
                    let coeff_x = self.h_coeffs.get(0, i, j, k);
                    let coeff_y = self.h_coeffs.get(1, i, j, k);
                    let coeff_z = self.h_coeffs.get(2, i, j, k);

                    // Curl of E for Hx: dEz/dy - dEy/dz
                    let dez_dy = if j + 1 < dims.ny {
                        self.e_field.z.get(i, j + 1, k) - self.e_field.z.get(i, j, k)
                    } else {
                        0.0
                    };
                    let dey_dz = if k + 1 < dims.nz {
                        self.e_field.y.get(i, j, k + 1) - self.e_field.y.get(i, j, k)
                    } else {
                        0.0
                    };
                    let curl_x = dez_dy - dey_dz;
                    let hx_old = self.h_field.x.get(i, j, k);
                    self.h_field.x.set(i, j, k, coeff_x.da * hx_old + coeff_x.db * curl_x);

                    // Curl of E for Hy: dEx/dz - dEz/dx
                    let dex_dz = if k + 1 < dims.nz {
                        self.e_field.x.get(i, j, k + 1) - self.e_field.x.get(i, j, k)
                    } else {
                        0.0
                    };
                    let dez_dx = if i + 1 < dims.nx {
                        self.e_field.z.get(i + 1, j, k) - self.e_field.z.get(i, j, k)
                    } else {
                        0.0
                    };
                    let curl_y = dex_dz - dez_dx;
                    let hy_old = self.h_field.y.get(i, j, k);
                    self.h_field.y.set(i, j, k, coeff_y.da * hy_old + coeff_y.db * curl_y);

                    // Curl of E for Hz: dEy/dx - dEx/dy
                    let dey_dx = if i + 1 < dims.nx {
                        self.e_field.y.get(i + 1, j, k) - self.e_field.y.get(i, j, k)
                    } else {
                        0.0
                    };
                    let dex_dy = if j + 1 < dims.ny {
                        self.e_field.x.get(i, j + 1, k) - self.e_field.x.get(i, j, k)
                    } else {
                        0.0
                    };
                    let curl_z = dey_dx - dex_dy;
                    let hz_old = self.h_field.z.get(i, j, k);
                    self.h_field.z.set(i, j, k, coeff_z.da * hz_old + coeff_z.db * curl_z);
                }
            }
        }
    }

    /// Update E-field using compressed coefficients (basic implementation).
    fn update_e_basic(&mut self) {
        let dims = self.dims;

        for i in 1..dims.nx {
            for j in 1..dims.ny {
                for k in 1..dims.nz {
                    // Get compressed coefficients via index lookup
                    let coeff_x = self.e_coeffs.get(0, i, j, k);
                    let coeff_y = self.e_coeffs.get(1, i, j, k);
                    let coeff_z = self.e_coeffs.get(2, i, j, k);

                    // Curl of H for Ex: dHz/dy - dHy/dz
                    let dhz_dy = self.h_field.z.get(i, j, k) - self.h_field.z.get(i, j - 1, k);
                    let dhy_dz = self.h_field.y.get(i, j, k) - self.h_field.y.get(i, j, k - 1);
                    let curl_x = dhz_dy - dhy_dz;
                    let ex_old = self.e_field.x.get(i, j, k);
                    self.e_field.x.set(i, j, k, coeff_x.ca * ex_old + coeff_x.cb * curl_x);

                    // Curl of H for Ey: dHx/dz - dHz/dx
                    let dhx_dz = self.h_field.x.get(i, j, k) - self.h_field.x.get(i, j, k - 1);
                    let dhz_dx = self.h_field.z.get(i, j, k) - self.h_field.z.get(i - 1, j, k);
                    let curl_y = dhx_dz - dhz_dx;
                    let ey_old = self.e_field.y.get(i, j, k);
                    self.e_field.y.set(i, j, k, coeff_y.ca * ey_old + coeff_y.cb * curl_y);

                    // Curl of H for Ez: dHy/dx - dHx/dy
                    let dhy_dx = self.h_field.y.get(i, j, k) - self.h_field.y.get(i - 1, j, k);
                    let dhx_dy = self.h_field.x.get(i, j, k) - self.h_field.x.get(i, j - 1, k);
                    let curl_z = dhy_dx - dhx_dy;
                    let ez_old = self.e_field.z.get(i, j, k);
                    self.e_field.z.set(i, j, k, coeff_z.ca * ez_old + coeff_z.cb * curl_z);
                }
            }
        }
    }

    /// Parallel implementation using Rayon and compressed coefficients.
    fn step_parallel(&mut self) {
        self.update_h_parallel();
        self.update_e_parallel();
    }

    /// Parallel H-field update using compressed coefficients.
    fn update_h_parallel(&mut self) {
        let dims = self.dims;
        let nx = dims.nx;
        let ny = dims.ny;
        let nz = dims.nz;

        // Get raw pointers for parallel access
        let hx_ptr = SendPtrMut::new(self.h_field.x.as_mut_ptr());
        let hy_ptr = SendPtrMut::new(self.h_field.y.as_mut_ptr());
        let hz_ptr = SendPtrMut::new(self.h_field.z.as_mut_ptr());

        let ex_ptr = SendPtr::new(self.e_field.x.as_ptr());
        let ey_ptr = SendPtr::new(self.e_field.y.as_ptr());
        let ez_ptr = SendPtr::new(self.e_field.z.as_ptr());

        // Get coefficient table and index pointers
        let h_table_x = self.h_coeffs.table(0);
        let h_table_y = self.h_coeffs.table(1);
        let h_table_z = self.h_coeffs.table(2);

        let h_indices_x = SendPtr::new(self.h_coeffs.indices(0).as_ptr());
        let h_indices_y = SendPtr::new(self.h_coeffs.indices(1).as_ptr());
        let h_indices_z = SendPtr::new(self.h_coeffs.indices(2).as_ptr());

        let h_table_x_ptr = SendPtr::new(h_table_x.as_ptr());
        let h_table_y_ptr = SendPtr::new(h_table_y.as_ptr());
        let h_table_z_ptr = SendPtr::new(h_table_z.as_ptr());

        // Parallel iteration over i-slices
        (0..nx).into_par_iter().for_each(|i| {
            for j in 0..ny {
                for k in 0..nz {
                    let idx = dims.to_linear(i, j, k);

                    unsafe {
                        // Load coefficient indices and look up coefficients
                        let coeff_idx_x = *h_indices_x.add(idx) as usize;
                        let coeff_idx_y = *h_indices_y.add(idx) as usize;
                        let coeff_idx_z = *h_indices_z.add(idx) as usize;

                        let coeff_x = *h_table_x_ptr.add(coeff_idx_x);
                        let coeff_y = *h_table_y_ptr.add(coeff_idx_y);
                        let coeff_z = *h_table_z_ptr.add(coeff_idx_z);

                        // Get E-field values
                        let ez_curr = *ez_ptr.add(idx);
                        let ez_jp1 = if j + 1 < ny {
                            *ez_ptr.add(dims.to_linear(i, j + 1, k))
                        } else {
                            ez_curr
                        };
                        let ey_curr = *ey_ptr.add(idx);
                        let ey_kp1 = if k + 1 < nz {
                            *ey_ptr.add(dims.to_linear(i, j, k + 1))
                        } else {
                            ey_curr
                        };
                        let ex_curr = *ex_ptr.add(idx);
                        let ex_kp1 = if k + 1 < nz {
                            *ex_ptr.add(dims.to_linear(i, j, k + 1))
                        } else {
                            ex_curr
                        };
                        let ex_jp1 = if j + 1 < ny {
                            *ex_ptr.add(dims.to_linear(i, j + 1, k))
                        } else {
                            ex_curr
                        };
                        let ez_ip1 = if i + 1 < nx {
                            *ez_ptr.add(dims.to_linear(i + 1, j, k))
                        } else {
                            ez_curr
                        };
                        let ey_ip1 = if i + 1 < nx {
                            *ey_ptr.add(dims.to_linear(i + 1, j, k))
                        } else {
                            ey_curr
                        };

                        // Hx update: curl_x = dEz/dy - dEy/dz
                        let curl_x = (ez_jp1 - ez_curr) - (ey_kp1 - ey_curr);
                        let hx_old = *hx_ptr.add(idx);
                        *hx_ptr.add(idx) = coeff_x.da * hx_old + coeff_x.db * curl_x;

                        // Hy update: curl_y = dEx/dz - dEz/dx
                        let curl_y = (ex_kp1 - ex_curr) - (ez_ip1 - ez_curr);
                        let hy_old = *hy_ptr.add(idx);
                        *hy_ptr.add(idx) = coeff_y.da * hy_old + coeff_y.db * curl_y;

                        // Hz update: curl_z = dEy/dx - dEx/dy
                        let curl_z = (ey_ip1 - ey_curr) - (ex_jp1 - ex_curr);
                        let hz_old = *hz_ptr.add(idx);
                        *hz_ptr.add(idx) = coeff_z.da * hz_old + coeff_z.db * curl_z;
                    }
                }
            }
        });
    }

    /// Parallel E-field update using compressed coefficients.
    fn update_e_parallel(&mut self) {
        let dims = self.dims;
        let nx = dims.nx;
        let ny = dims.ny;
        let nz = dims.nz;

        let ex_ptr = SendPtrMut::new(self.e_field.x.as_mut_ptr());
        let ey_ptr = SendPtrMut::new(self.e_field.y.as_mut_ptr());
        let ez_ptr = SendPtrMut::new(self.e_field.z.as_mut_ptr());

        let hx_ptr = SendPtr::new(self.h_field.x.as_ptr());
        let hy_ptr = SendPtr::new(self.h_field.y.as_ptr());
        let hz_ptr = SendPtr::new(self.h_field.z.as_ptr());

        // Get coefficient table and index pointers
        let e_table_x = self.e_coeffs.table(0);
        let e_table_y = self.e_coeffs.table(1);
        let e_table_z = self.e_coeffs.table(2);

        let e_indices_x = SendPtr::new(self.e_coeffs.indices(0).as_ptr());
        let e_indices_y = SendPtr::new(self.e_coeffs.indices(1).as_ptr());
        let e_indices_z = SendPtr::new(self.e_coeffs.indices(2).as_ptr());

        let e_table_x_ptr = SendPtr::new(e_table_x.as_ptr());
        let e_table_y_ptr = SendPtr::new(e_table_y.as_ptr());
        let e_table_z_ptr = SendPtr::new(e_table_z.as_ptr());

        (1..nx).into_par_iter().for_each(|i| {
            for j in 1..ny {
                for k in 1..nz {
                    let idx = dims.to_linear(i, j, k);

                    unsafe {
                        // Load coefficient indices and look up coefficients
                        let coeff_idx_x = *e_indices_x.add(idx) as usize;
                        let coeff_idx_y = *e_indices_y.add(idx) as usize;
                        let coeff_idx_z = *e_indices_z.add(idx) as usize;

                        let coeff_x = *e_table_x_ptr.add(coeff_idx_x);
                        let coeff_y = *e_table_y_ptr.add(coeff_idx_y);
                        let coeff_z = *e_table_z_ptr.add(coeff_idx_z);

                        // Get H-field values
                        let hz_curr = *hz_ptr.add(idx);
                        let hz_jm1 = *hz_ptr.add(dims.to_linear(i, j - 1, k));
                        let hz_im1 = *hz_ptr.add(dims.to_linear(i - 1, j, k));
                        let hy_curr = *hy_ptr.add(idx);
                        let hy_km1 = *hy_ptr.add(dims.to_linear(i, j, k - 1));
                        let hy_im1 = *hy_ptr.add(dims.to_linear(i - 1, j, k));
                        let hx_curr = *hx_ptr.add(idx);
                        let hx_km1 = *hx_ptr.add(dims.to_linear(i, j, k - 1));
                        let hx_jm1 = *hx_ptr.add(dims.to_linear(i, j - 1, k));

                        // Ex update: curl = dHz/dy - dHy/dz
                        let curl_x = (hz_curr - hz_jm1) - (hy_curr - hy_km1);
                        let ex_old = *ex_ptr.add(idx);
                        *ex_ptr.add(idx) = coeff_x.ca * ex_old + coeff_x.cb * curl_x;

                        // Ey update: curl = dHx/dz - dHz/dx
                        let curl_y = (hx_curr - hx_km1) - (hz_curr - hz_im1);
                        let ey_old = *ey_ptr.add(idx);
                        *ey_ptr.add(idx) = coeff_y.ca * ey_old + coeff_y.cb * curl_y;

                        // Ez update: curl = dHy/dx - dHx/dy
                        let curl_z = (hy_curr - hy_im1) - (hx_curr - hx_jm1);
                        let ez_old = *ez_ptr.add(idx);
                        *ez_ptr.add(idx) = coeff_z.ca * ez_old + coeff_z.cb * curl_z;
                    }
                }
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
            engine.step().unwrap();
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

        // For a uniform vacuum grid, we should have very few unique coefficients
        // The compression ratio might be < 1 due to index overhead, but
        // the memory bandwidth benefit comes from cache utilization
        println!("E compression ratio: {:.2}, H compression ratio: {:.2}", e_ratio, h_ratio);

        // Just verify it runs without panic
        assert!(e_ratio > 0.0);
        assert!(h_ratio > 0.0);
    }
}
