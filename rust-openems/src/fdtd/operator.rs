//! FDTD Operator - Material coefficient computation.
//!
//! The Operator pre-computes all material-dependent coefficients needed
//! for the FDTD time-stepping. This includes:
//!
//! - E-field update coefficients (Ca, Cb)
//! - H-field update coefficients (Da, Db)
//! - PML coefficients (if applicable)
//! - Dispersive material poles (Lorentz, Debye, etc.)

use crate::arrays::{Dimensions, Field3D};
use crate::constants::{C0, EPS0, MU0};
use crate::geometry::Grid;
use crate::Result;

use super::BoundaryConditions;

/// Material update coefficients for E-field.
///
/// E_new = Ca * E_old + Cb * curl_H
#[derive(Clone)]
pub struct EFieldCoefficients {
    /// Ca coefficient (3 components)
    pub ca: [Field3D; 3],
    /// Cb coefficient (3 components)
    pub cb: [Field3D; 3],
}

/// Material update coefficients for H-field.
///
/// H_new = Da * H_old + Db * curl_E
#[derive(Clone)]
pub struct HFieldCoefficients {
    /// Da coefficient (3 components)
    pub da: [Field3D; 3],
    /// Db coefficient (3 components)
    pub db: [Field3D; 3],
}

/// FDTD Operator containing all pre-computed coefficients.
#[allow(dead_code)]
pub struct Operator {
    /// Grid information
    grid: Grid,
    /// E-field coefficients
    e_coeff: EFieldCoefficients,
    /// H-field coefficients
    h_coeff: HFieldCoefficients,
    /// Timestep (seconds)
    dt: f64,
    /// Number of timesteps for excitation signal
    excitation_timesteps: usize,
    /// Boundary conditions
    boundaries: BoundaryConditions,
}

impl Operator {
    /// Create a new operator for the given grid.
    pub fn new(grid: Grid, boundaries: BoundaryConditions) -> Result<Self> {
        let dims = grid.dimensions();

        // Calculate stable timestep using CFL condition
        let dt = Self::calculate_timestep(&grid);

        // Initialize coefficients
        let e_coeff = Self::init_e_coefficients(dims, &grid, dt);
        let h_coeff = Self::init_h_coefficients(dims, &grid, dt);

        Ok(Self {
            grid,
            e_coeff,
            h_coeff,
            dt,
            excitation_timesteps: 0,
            boundaries,
        })
    }

    /// Calculate stable timestep using CFL condition.
    ///
    /// For a 3D Cartesian grid: dt <= 1 / (c * sqrt(1/dx^2 + 1/dy^2 + 1/dz^2))
    fn calculate_timestep(grid: &Grid) -> f64 {
        let (dx, dy, dz) = grid.cell_size();
        let inv_sum = 1.0 / (dx * dx) + 1.0 / (dy * dy) + 1.0 / (dz * dz);
        let dt_max = 1.0 / (C0 * inv_sum.sqrt());

        // Use 99% of maximum stable timestep for safety margin
        0.99 * dt_max
    }

    /// Initialize E-field update coefficients.
    fn init_e_coefficients(dims: Dimensions, grid: &Grid, dt: f64) -> EFieldCoefficients {
        let mut ca = [
            Field3D::new(dims),
            Field3D::new(dims),
            Field3D::new(dims),
        ];
        let mut cb = [
            Field3D::new(dims),
            Field3D::new(dims),
            Field3D::new(dims),
        ];

        let (dx, dy, dz) = grid.cell_size();
        let _dt_f32 = dt as f32;

        // For vacuum: Ca = 1, Cb = dt / (eps0 * delta)
        // E_new = E_old + (dt/eps0) * curl_H

        // Cb values for each direction
        let cb_x = (dt / (EPS0 * dx)) as f32;
        let cb_y = (dt / (EPS0 * dy)) as f32;
        let cb_z = (dt / (EPS0 * dz)) as f32;

        // Fill with vacuum coefficients
        for ca_dir in &mut ca {
            ca_dir.fill(1.0);
        }
        cb[0].fill(cb_x);
        cb[1].fill(cb_y);
        cb[2].fill(cb_z);

        EFieldCoefficients { ca, cb }
    }

    /// Initialize H-field update coefficients.
    fn init_h_coefficients(dims: Dimensions, grid: &Grid, dt: f64) -> HFieldCoefficients {
        let mut da = [
            Field3D::new(dims),
            Field3D::new(dims),
            Field3D::new(dims),
        ];
        let mut db = [
            Field3D::new(dims),
            Field3D::new(dims),
            Field3D::new(dims),
        ];

        let (dx, dy, dz) = grid.cell_size();

        // For vacuum: Da = 1, Db = -dt / (mu0 * delta)
        // From Maxwell: dH/dt = -(1/mu) * curl(E)
        // So: H_new = H_old - (dt/mu) * curl_E

        let db_x = -(dt / (MU0 * dx)) as f32;
        let db_y = -(dt / (MU0 * dy)) as f32;
        let db_z = -(dt / (MU0 * dz)) as f32;

        for da_dir in &mut da {
            da_dir.fill(1.0);
        }
        db[0].fill(db_x);
        db[1].fill(db_y);
        db[2].fill(db_z);

        HFieldCoefficients { da, db }
    }

    /// Get the timestep in seconds.
    #[inline]
    pub fn timestep(&self) -> f64 {
        self.dt
    }

    /// Get the grid.
    #[inline]
    pub fn grid(&self) -> &Grid {
        &self.grid
    }

    /// Get the grid dimensions.
    #[inline]
    pub fn dimensions(&self) -> Dimensions {
        self.grid.dimensions()
    }

    /// Get E-field coefficients.
    #[inline]
    pub fn e_coefficients(&self) -> &EFieldCoefficients {
        &self.e_coeff
    }

    /// Get H-field coefficients.
    #[inline]
    pub fn h_coefficients(&self) -> &HFieldCoefficients {
        &self.h_coeff
    }

    /// Get boundary conditions.
    #[inline]
    pub fn boundaries(&self) -> &BoundaryConditions {
        &self.boundaries
    }

    /// Apply material properties at a specific region.
    ///
    /// # Arguments
    /// * `eps_r` - Relative permittivity
    /// * `mu_r` - Relative permeability
    /// * `sigma_e` - Electric conductivity (S/m)
    /// * `sigma_m` - Magnetic conductivity (Ohm/m)
    /// * `region` - (i_start, i_end, j_start, j_end, k_start, k_end)
    pub fn set_material(
        &mut self,
        eps_r: f64,
        mu_r: f64,
        sigma_e: f64,
        sigma_m: f64,
        region: (usize, usize, usize, usize, usize, usize),
    ) {
        let (i0, i1, j0, j1, k0, k1) = region;
        let dt = self.dt;

        // E-field coefficients with material properties
        // Ca = (1 - sigma_e*dt/(2*eps)) / (1 + sigma_e*dt/(2*eps))
        // Cb = (dt/(eps*delta)) / (1 + sigma_e*dt/(2*eps))
        let eps = eps_r * EPS0;
        let loss_e = sigma_e * dt / (2.0 * eps);
        let ca_mat = ((1.0 - loss_e) / (1.0 + loss_e)) as f32;

        let (dx, dy, dz) = self.grid.cell_size();
        let cb_factor = (dt / eps) / (1.0 + loss_e);
        let cb_x = (cb_factor / dx) as f32;
        let cb_y = (cb_factor / dy) as f32;
        let cb_z = (cb_factor / dz) as f32;

        // H-field coefficients with material properties
        let mu = mu_r * MU0;
        let loss_m = sigma_m * dt / (2.0 * mu);
        let da_mat = ((1.0 - loss_m) / (1.0 + loss_m)) as f32;

        let db_factor = (dt / mu) / (1.0 + loss_m);
        let db_x = (db_factor / dx) as f32;
        let db_y = (db_factor / dy) as f32;
        let db_z = (db_factor / dz) as f32;

        // Apply to the region
        for i in i0..i1 {
            for j in j0..j1 {
                for k in k0..k1 {
                    // E-field coefficients
                    for dir in 0..3 {
                        self.e_coeff.ca[dir].set(i, j, k, ca_mat);
                    }
                    self.e_coeff.cb[0].set(i, j, k, cb_x);
                    self.e_coeff.cb[1].set(i, j, k, cb_y);
                    self.e_coeff.cb[2].set(i, j, k, cb_z);

                    // H-field coefficients
                    for dir in 0..3 {
                        self.h_coeff.da[dir].set(i, j, k, da_mat);
                    }
                    self.h_coeff.db[0].set(i, j, k, db_x);
                    self.h_coeff.db[1].set(i, j, k, db_y);
                    self.h_coeff.db[2].set(i, j, k, db_z);
                }
            }
        }
    }

    /// Mark a region as Perfect Electric Conductor (PEC).
    ///
    /// Sets E-field to zero in the region.
    pub fn set_pec(&mut self, region: (usize, usize, usize, usize, usize, usize)) {
        let (i0, i1, j0, j1, k0, k1) = region;

        for i in i0..i1 {
            for j in j0..j1 {
                for k in k0..k1 {
                    // PEC: Ca = 0, Cb = 0 (E-field forced to zero)
                    for dir in 0..3 {
                        self.e_coeff.ca[dir].set(i, j, k, 0.0);
                        self.e_coeff.cb[dir].set(i, j, k, 0.0);
                    }
                }
            }
        }
    }

    /// Get Nyquist rate for given frequency
    pub fn nyquist_timesteps(&self, freq: f64) -> usize {
        let period = 1.0 / freq;
        let steps = period / self.dt;
        steps.ceil() as usize
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::CoordinateSystem;

    #[test]
    fn test_timestep_calculation() {
        let grid = Grid::new(
            CoordinateSystem::Cartesian,
            vec![0.0, 0.001, 0.002], // 1mm cells
            vec![0.0, 0.001, 0.002],
            vec![0.0, 0.001, 0.002],
        );
        let op = Operator::new(grid, BoundaryConditions::default()).unwrap();

        // Timestep should be less than CFL limit
        let dt = op.timestep();
        let dt_cfl = 1.0 / (C0 * (3.0f64).sqrt() / 0.001);
        assert!(dt < dt_cfl);
        assert!(dt > 0.0);
    }
}
