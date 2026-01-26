//! FDTD Engine - Time-domain field update implementation.
//!
//! The Engine performs the actual FDTD time-stepping using SIMD-accelerated
//! kernels and multi-threaded parallelization.

use crate::arrays::{Dimensions, VectorField3D};
use crate::Result;
use rayon::prelude::*;

use super::operator::Operator;
use super::EngineType;

/// Wrapper for raw pointer to make it Send + Sync for parallel iteration.
///
/// # Safety
/// The caller must ensure that concurrent access patterns are safe:
/// - Either only reading from the pointer
/// - Or writing to non-overlapping regions
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

/// FDTD Engine state
#[allow(dead_code)]
pub struct Engine {
    /// Electric field components (Ex, Ey, Ez)
    e_field: VectorField3D,
    /// Magnetic field components (Hx, Hy, Hz)
    h_field: VectorField3D,
    /// Current timestep
    timestep: u64,
    /// Engine type (determines parallelization strategy)
    engine_type: EngineType,
    /// Number of threads for parallel engine
    num_threads: usize,
    /// Temporary storage for curl computations
    curl_buffer: VectorField3D,
}

impl Engine {
    /// Create a new engine for the given operator.
    pub fn new(operator: &Operator, engine_type: EngineType) -> Self {
        let dims = operator.dimensions();
        let num_threads = rayon::current_num_threads();

        Self {
            e_field: VectorField3D::new(dims),
            h_field: VectorField3D::new(dims),
            timestep: 0,
            engine_type,
            num_threads,
            curl_buffer: VectorField3D::new(dims),
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

    /// Compute total electromagnetic energy in the simulation domain.
    pub fn total_energy(&self, _operator: &Operator) -> f64 {
        // E_total = 0.5 * eps * |E|^2 + 0.5 * mu * |H|^2
        // For vacuum: E_total ∝ |E|^2 + (Z0)^2 * |H|^2
        let e_energy = self.e_field.energy();
        let h_energy = self.h_field.energy();
        e_energy + h_energy
    }

    /// Perform one complete FDTD timestep.
    ///
    /// The leapfrog update scheme:
    /// 1. Update H-field at t + dt/2
    /// 2. Update E-field at t + dt
    pub fn step(&mut self, operator: &Operator) -> Result<()> {
        match self.engine_type {
            EngineType::Basic => self.step_basic(operator),
            EngineType::Simd => self.step_simd(operator),
            EngineType::Parallel => self.step_parallel(operator),
        }
        self.timestep += 1;
        Ok(())
    }

    /// Basic (reference) implementation without SIMD.
    fn step_basic(&mut self, operator: &Operator) {
        let dims = operator.dimensions();
        let e_coeff = operator.e_coefficients();
        let h_coeff = operator.h_coefficients();

        // Update H-field: H = Da*H + Db*curl(E)
        self.update_h_basic(dims, h_coeff);

        // Update E-field: E = Ca*E + Cb*curl(H)
        self.update_e_basic(dims, e_coeff);
    }

    /// Update H-field using basic implementation.
    fn update_h_basic(&mut self, dims: Dimensions, h_coeff: &super::operator::HFieldCoefficients) {
        // Hx = Da_x*Hx + Db_x * (dEz/dy - dEy/dz)
        // Hy = Da_y*Hy + Db_y * (dEx/dz - dEz/dx)
        // Hz = Da_z*Hz + Db_z * (dEy/dx - dEx/dy)

        for i in 0..dims.nx {
            for j in 0..dims.ny {
                for k in 0..dims.nz {
                    // Get neighboring E-field values for curl calculation
                    // At boundaries, derivatives are set to 0 (PEC-like behavior)
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
                    let da = h_coeff.da[0].get(i, j, k);
                    let db = h_coeff.db[0].get(i, j, k);
                    let hx_old = self.h_field.x.get(i, j, k);
                    self.h_field.x.set(i, j, k, da * hx_old + db * curl_x);

                    // Hy
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
                    let da = h_coeff.da[1].get(i, j, k);
                    let db = h_coeff.db[1].get(i, j, k);
                    let hy_old = self.h_field.y.get(i, j, k);
                    self.h_field.y.set(i, j, k, da * hy_old + db * curl_y);

                    // Hz
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
                    let da = h_coeff.da[2].get(i, j, k);
                    let db = h_coeff.db[2].get(i, j, k);
                    let hz_old = self.h_field.z.get(i, j, k);
                    self.h_field.z.set(i, j, k, da * hz_old + db * curl_z);
                }
            }
        }
    }

    /// Update E-field using basic implementation.
    fn update_e_basic(&mut self, dims: Dimensions, e_coeff: &super::operator::EFieldCoefficients) {
        // Ex = Ca_x*Ex + Cb_x * (dHz/dy - dHy/dz)
        // Ey = Ca_y*Ey + Cb_y * (dHx/dz - dHz/dx)
        // Ez = Ca_z*Ez + Cb_z * (dHy/dx - dHx/dy)

        for i in 1..dims.nx {
            for j in 1..dims.ny {
                for k in 1..dims.nz {
                    // Ex
                    let hz_j = self.h_field.z.get(i, j, k);
                    let hz_jm1 = self.h_field.z.get(i, j - 1, k);
                    let hy_k = self.h_field.y.get(i, j, k);
                    let hy_km1 = self.h_field.y.get(i, j, k - 1);

                    let curl_x = (hz_j - hz_jm1) - (hy_k - hy_km1);
                    let ca = e_coeff.ca[0].get(i, j, k);
                    let cb = e_coeff.cb[0].get(i, j, k);
                    let ex_old = self.e_field.x.get(i, j, k);
                    self.e_field.x.set(i, j, k, ca * ex_old + cb * curl_x);

                    // Ey
                    let hx_k = self.h_field.x.get(i, j, k);
                    let hx_km1 = self.h_field.x.get(i, j, k - 1);
                    let hz_i = self.h_field.z.get(i, j, k);
                    let hz_im1 = self.h_field.z.get(i - 1, j, k);

                    let curl_y = (hx_k - hx_km1) - (hz_i - hz_im1);
                    let ca = e_coeff.ca[1].get(i, j, k);
                    let cb = e_coeff.cb[1].get(i, j, k);
                    let ey_old = self.e_field.y.get(i, j, k);
                    self.e_field.y.set(i, j, k, ca * ey_old + cb * curl_y);

                    // Ez
                    let hy_i = self.h_field.y.get(i, j, k);
                    let hy_im1 = self.h_field.y.get(i - 1, j, k);
                    let hx_j = self.h_field.x.get(i, j, k);
                    let hx_jm1 = self.h_field.x.get(i, j - 1, k);

                    let curl_z = (hy_i - hy_im1) - (hx_j - hx_jm1);
                    let ca = e_coeff.ca[2].get(i, j, k);
                    let cb = e_coeff.cb[2].get(i, j, k);
                    let ez_old = self.e_field.z.get(i, j, k);
                    self.e_field.z.set(i, j, k, ca * ez_old + cb * curl_z);
                }
            }
        }
    }

    /// SIMD-optimized single-threaded implementation.
    fn step_simd(&mut self, operator: &Operator) {
        let dims = operator.dimensions();
        let e_coeff = operator.e_coefficients();
        let h_coeff = operator.h_coefficients();

        // SIMD H-field update
        self.update_h_simd(dims, h_coeff);

        // SIMD E-field update
        self.update_e_simd(dims, e_coeff);
    }

    /// SIMD-optimized H-field update operating on z-lines.
    fn update_h_simd(&mut self, dims: Dimensions, h_coeff: &super::operator::HFieldCoefficients) {
        let nz = dims.nz;

        for i in 0..dims.nx {
            for j in 0..dims.ny {
                // Get z-slices for vectorized operations
                let hx_line = self.h_field.x.z_slice_mut(i, j);
                let hy_line = self.h_field.y.z_slice_mut(i, j);
                let hz_line = self.h_field.z.z_slice_mut(i, j);

                let ez_line = self.e_field.z.z_slice(i, j);
                let ez_jp1 = if j + 1 < dims.ny {
                    self.e_field.z.z_slice(i, j + 1)
                } else {
                    ez_line // boundary
                };

                let ey_line = self.e_field.y.z_slice(i, j);
                let ex_line = self.e_field.x.z_slice(i, j);
                let ex_jp1 = if j + 1 < dims.ny {
                    self.e_field.x.z_slice(i, j + 1)
                } else {
                    ex_line
                };

                let ey_ip1 = if i + 1 < dims.nx {
                    self.e_field.y.z_slice(i + 1, j)
                } else {
                    ey_line
                };

                let ez_ip1 = if i + 1 < dims.nx {
                    self.e_field.z.z_slice(i + 1, j)
                } else {
                    ez_line
                };

                let da_x = h_coeff.da[0].z_slice(i, j);
                let db_x = h_coeff.db[0].z_slice(i, j);
                let da_y = h_coeff.da[1].z_slice(i, j);
                let db_y = h_coeff.db[1].z_slice(i, j);
                let da_z = h_coeff.da[2].z_slice(i, j);
                let db_z = h_coeff.db[2].z_slice(i, j);

                // Update Hx: curl_x = dEz/dy - dEy/dz
                for k in 0..nz {
                    let dez_dy = ez_jp1[k] - ez_line[k];
                    let dey_dz = if k + 1 < nz {
                        ey_line[k + 1] - ey_line[k]
                    } else {
                        0.0
                    };
                    let curl_x = dez_dy - dey_dz;
                    hx_line[k] = da_x[k] * hx_line[k] + db_x[k] * curl_x;
                }

                // Update Hy: curl_y = dEx/dz - dEz/dx
                for k in 0..nz {
                    let dex_dz = if k + 1 < nz {
                        ex_line[k + 1] - ex_line[k]
                    } else {
                        0.0
                    };
                    let dez_dx = ez_ip1[k] - ez_line[k];
                    let curl_y = dex_dz - dez_dx;
                    hy_line[k] = da_y[k] * hy_line[k] + db_y[k] * curl_y;
                }

                // Update Hz: curl_z = dEy/dx - dEx/dy
                for k in 0..nz {
                    let dey_dx = ey_ip1[k] - ey_line[k];
                    let dex_dy = ex_jp1[k] - ex_line[k];
                    let curl_z = dey_dx - dex_dy;
                    hz_line[k] = da_z[k] * hz_line[k] + db_z[k] * curl_z;
                }
            }
        }
    }

    /// SIMD-optimized E-field update operating on z-lines.
    fn update_e_simd(&mut self, dims: Dimensions, e_coeff: &super::operator::EFieldCoefficients) {
        let nz = dims.nz;

        for i in 1..dims.nx {
            for j in 1..dims.ny {
                let ex_line = self.e_field.x.z_slice_mut(i, j);
                let ey_line = self.e_field.y.z_slice_mut(i, j);
                let ez_line = self.e_field.z.z_slice_mut(i, j);

                let hz_line = self.h_field.z.z_slice(i, j);
                let hz_jm1 = self.h_field.z.z_slice(i, j - 1);
                let hz_im1 = self.h_field.z.z_slice(i - 1, j);

                let hy_line = self.h_field.y.z_slice(i, j);
                let hy_im1 = self.h_field.y.z_slice(i - 1, j);

                let hx_line = self.h_field.x.z_slice(i, j);
                let hx_jm1 = self.h_field.x.z_slice(i, j - 1);

                let ca_x = e_coeff.ca[0].z_slice(i, j);
                let cb_x = e_coeff.cb[0].z_slice(i, j);
                let ca_y = e_coeff.ca[1].z_slice(i, j);
                let cb_y = e_coeff.cb[1].z_slice(i, j);
                let ca_z = e_coeff.ca[2].z_slice(i, j);
                let cb_z = e_coeff.cb[2].z_slice(i, j);

                // Update Ex: curl = dHz/dy - dHy/dz
                for k in 1..nz {
                    let dhz_dy = hz_line[k] - hz_jm1[k];
                    let dhy_dz = hy_line[k] - hy_line[k - 1];
                    let curl_x = dhz_dy - dhy_dz;
                    ex_line[k] = ca_x[k] * ex_line[k] + cb_x[k] * curl_x;
                }

                // Update Ey: curl = dHx/dz - dHz/dx
                for k in 1..nz {
                    let dhx_dz = hx_line[k] - hx_line[k - 1];
                    let dhz_dx = hz_line[k] - hz_im1[k];
                    let curl_y = dhx_dz - dhz_dx;
                    ey_line[k] = ca_y[k] * ey_line[k] + cb_y[k] * curl_y;
                }

                // Update Ez: curl = dHy/dx - dHx/dy
                for k in 1..nz {
                    let dhy_dx = hy_line[k] - hy_im1[k];
                    let dhx_dy = hx_line[k] - hx_jm1[k];
                    let curl_z = dhy_dx - dhx_dy;
                    ez_line[k] = ca_z[k] * ez_line[k] + cb_z[k] * curl_z;
                }
            }
        }
    }

    /// Parallel SIMD-optimized implementation using Rayon.
    fn step_parallel(&mut self, operator: &Operator) {
        // For parallel updates, we need to be careful about data dependencies.
        // The H-field update depends only on E-field values, and vice versa.
        // So we can parallelize each half-step independently.

        let dims = operator.dimensions();
        let h_coeff = operator.h_coefficients();
        let e_coeff = operator.e_coefficients();

        // Parallel H-field update
        self.update_h_parallel(dims, h_coeff);

        // Parallel E-field update
        self.update_e_parallel(dims, e_coeff);
    }

    /// Parallel H-field update using rayon.
    fn update_h_parallel(
        &mut self,
        dims: Dimensions,
        h_coeff: &super::operator::HFieldCoefficients,
    ) {
        // Split the domain into chunks along the i-axis for parallel processing
        // Each thread processes a set of i-planes independently

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

        let da_x_ptr = SendPtr::new(h_coeff.da[0].as_ptr());
        let db_x_ptr = SendPtr::new(h_coeff.db[0].as_ptr());
        let da_y_ptr = SendPtr::new(h_coeff.da[1].as_ptr());
        let db_y_ptr = SendPtr::new(h_coeff.db[1].as_ptr());
        let da_z_ptr = SendPtr::new(h_coeff.da[2].as_ptr());
        let db_z_ptr = SendPtr::new(h_coeff.db[2].as_ptr());

        // Process i-slices in parallel
        (0..nx).into_par_iter().for_each(|i| {
            for j in 0..ny {
                for k in 0..nz {
                    let idx = dims.to_linear(i, j, k);

                    unsafe {
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

                        // Hx update
                        let curl_x = (ez_jp1 - ez_curr) - (ey_kp1 - ey_curr);
                        let da = *da_x_ptr.add(idx);
                        let db = *db_x_ptr.add(idx);
                        let hx_old = *hx_ptr.add(idx);
                        *hx_ptr.add(idx) = da * hx_old + db * curl_x;

                        // Hy update
                        let curl_y = (ex_kp1 - ex_curr) - (ez_ip1 - ez_curr);
                        let da = *da_y_ptr.add(idx);
                        let db = *db_y_ptr.add(idx);
                        let hy_old = *hy_ptr.add(idx);
                        *hy_ptr.add(idx) = da * hy_old + db * curl_y;

                        // Hz update
                        let curl_z = (ey_ip1 - ey_curr) - (ex_jp1 - ex_curr);
                        let da = *da_z_ptr.add(idx);
                        let db = *db_z_ptr.add(idx);
                        let hz_old = *hz_ptr.add(idx);
                        *hz_ptr.add(idx) = da * hz_old + db * curl_z;
                    }
                }
            }
        });
    }

    /// Parallel E-field update using rayon.
    fn update_e_parallel(
        &mut self,
        dims: Dimensions,
        e_coeff: &super::operator::EFieldCoefficients,
    ) {
        let nx = dims.nx;
        let ny = dims.ny;
        let nz = dims.nz;

        let ex_ptr = SendPtrMut::new(self.e_field.x.as_mut_ptr());
        let ey_ptr = SendPtrMut::new(self.e_field.y.as_mut_ptr());
        let ez_ptr = SendPtrMut::new(self.e_field.z.as_mut_ptr());

        let hx_ptr = SendPtr::new(self.h_field.x.as_ptr());
        let hy_ptr = SendPtr::new(self.h_field.y.as_ptr());
        let hz_ptr = SendPtr::new(self.h_field.z.as_ptr());

        let ca_x_ptr = SendPtr::new(e_coeff.ca[0].as_ptr());
        let cb_x_ptr = SendPtr::new(e_coeff.cb[0].as_ptr());
        let ca_y_ptr = SendPtr::new(e_coeff.ca[1].as_ptr());
        let cb_y_ptr = SendPtr::new(e_coeff.cb[1].as_ptr());
        let ca_z_ptr = SendPtr::new(e_coeff.ca[2].as_ptr());
        let cb_z_ptr = SendPtr::new(e_coeff.cb[2].as_ptr());

        (1..nx).into_par_iter().for_each(|i| {
            for j in 1..ny {
                for k in 1..nz {
                    let idx = dims.to_linear(i, j, k);

                    unsafe {
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

                        // Ex update
                        let curl_x = (hz_curr - hz_jm1) - (hy_curr - hy_km1);
                        let ca = *ca_x_ptr.add(idx);
                        let cb = *cb_x_ptr.add(idx);
                        let ex_old = *ex_ptr.add(idx);
                        *ex_ptr.add(idx) = ca * ex_old + cb * curl_x;

                        // Ey update
                        let curl_y = (hx_curr - hx_km1) - (hz_curr - hz_im1);
                        let ca = *ca_y_ptr.add(idx);
                        let cb = *cb_y_ptr.add(idx);
                        let ey_old = *ey_ptr.add(idx);
                        *ey_ptr.add(idx) = ca * ey_old + cb * curl_y;

                        // Ez update
                        let curl_z = (hy_curr - hy_im1) - (hx_curr - hx_jm1);
                        let ca = *ca_z_ptr.add(idx);
                        let cb = *cb_z_ptr.add(idx);
                        let ez_old = *ez_ptr.add(idx);
                        *ez_ptr.add(idx) = ca * ez_old + cb * curl_z;
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
    fn test_engine_creation() {
        let grid = Grid::new(
            CoordinateSystem::Cartesian,
            vec![0.0, 0.001, 0.002, 0.003],
            vec![0.0, 0.001, 0.002, 0.003],
            vec![0.0, 0.001, 0.002, 0.003],
        );
        let op = Operator::new(grid, BoundaryConditions::default()).unwrap();
        let engine = Engine::new(&op, EngineType::Parallel);

        assert_eq!(engine.timestep(), 0);
    }

    fn run_energy_test(engine_type: EngineType) {
        // In a closed cavity with PEC walls, energy should be conserved
        let grid = Grid::new(
            CoordinateSystem::Cartesian,
            (0..11).map(|i| i as f64 * 0.001).collect(),
            (0..11).map(|i| i as f64 * 0.001).collect(),
            (0..11).map(|i| i as f64 * 0.001).collect(),
        );
        let op = Operator::new(grid, BoundaryConditions::all_pec()).unwrap();
        let mut engine = Engine::new(&op, engine_type);

        // Add some initial energy
        engine.e_field.z.set(5, 5, 5, 1.0);

        let initial_energy = engine.total_energy(&op);

        // Run for some steps
        for step in 0..100 {
            engine.step(&op).unwrap();

            // Check for NaN after each step
            let energy = engine.total_energy(&op);
            assert!(
                !energy.is_nan(),
                "Energy became NaN at step {} with engine {:?}",
                step + 1,
                engine_type
            );
        }

        let final_energy = engine.total_energy(&op);

        // Energy should be approximately conserved (some numerical loss expected)
        let relative_diff = (final_energy - initial_energy).abs() / initial_energy;
        assert!(
            relative_diff < 0.1,
            "Energy not conserved with {:?}: initial={}, final={}, diff={}",
            engine_type,
            initial_energy,
            final_energy,
            relative_diff
        );
    }

    #[test]
    fn test_energy_conservation_basic() {
        run_energy_test(EngineType::Basic);
    }

    #[test]
    fn test_energy_conservation_simd() {
        run_energy_test(EngineType::Simd);
    }

    #[test]
    fn test_energy_conservation_parallel() {
        run_energy_test(EngineType::Parallel);
    }
}
