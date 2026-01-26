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

    fn create_test_grid() -> Grid {
        Grid::new(
            CoordinateSystem::Cartesian,
            vec![0.0, 0.001, 0.002, 0.003, 0.004, 0.005], // 1mm cells, 5 cells in x
            vec![0.0, 0.001, 0.002, 0.003, 0.004, 0.005], // 5 cells in y
            vec![0.0, 0.001, 0.002, 0.003, 0.004, 0.005], // 5 cells in z
        )
    }

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

    #[test]
    fn test_operator_creation() {
        let grid = create_test_grid();
        let op = Operator::new(grid, BoundaryConditions::default()).unwrap();

        // Verify basic properties
        assert!(op.timestep() > 0.0);
        let dims = op.dimensions();
        assert_eq!(dims.nx, 5);
        assert_eq!(dims.ny, 5);
        assert_eq!(dims.nz, 5);
    }

    #[test]
    fn test_e_field_coefficients_initialization() {
        let grid = create_test_grid();
        let op = Operator::new(grid, BoundaryConditions::default()).unwrap();
        let e_coeff = op.e_coefficients();

        // For vacuum, Ca should be 1.0
        assert!((e_coeff.ca[0].get(2, 2, 2) - 1.0).abs() < 1e-6);
        assert!((e_coeff.ca[1].get(2, 2, 2) - 1.0).abs() < 1e-6);
        assert!((e_coeff.ca[2].get(2, 2, 2) - 1.0).abs() < 1e-6);

        // Cb should be positive (dt / (eps0 * delta))
        assert!(e_coeff.cb[0].get(2, 2, 2) > 0.0);
        assert!(e_coeff.cb[1].get(2, 2, 2) > 0.0);
        assert!(e_coeff.cb[2].get(2, 2, 2) > 0.0);
    }

    #[test]
    fn test_h_field_coefficients_initialization() {
        let grid = create_test_grid();
        let op = Operator::new(grid, BoundaryConditions::default()).unwrap();
        let h_coeff = op.h_coefficients();

        // For vacuum, Da should be 1.0
        assert!((h_coeff.da[0].get(2, 2, 2) - 1.0).abs() < 1e-6);
        assert!((h_coeff.da[1].get(2, 2, 2) - 1.0).abs() < 1e-6);
        assert!((h_coeff.da[2].get(2, 2, 2) - 1.0).abs() < 1e-6);

        // Db should be negative (due to curl sign convention)
        assert!(h_coeff.db[0].get(2, 2, 2) < 0.0);
        assert!(h_coeff.db[1].get(2, 2, 2) < 0.0);
        assert!(h_coeff.db[2].get(2, 2, 2) < 0.0);
    }

    #[test]
    fn test_set_material() {
        let grid = create_test_grid();
        let mut op = Operator::new(grid, BoundaryConditions::default()).unwrap();

        // Set a dielectric region (eps_r = 4, lossless)
        op.set_material(4.0, 1.0, 0.0, 0.0, (1, 3, 1, 3, 1, 3));

        let e_coeff = op.e_coefficients();
        let h_coeff = op.h_coefficients();

        // Inside the region: Ca should still be 1.0 for lossless
        assert!((e_coeff.ca[0].get(2, 2, 2) - 1.0).abs() < 1e-6);

        // Cb should be reduced (1/eps_r factor)
        // Outside region (vacuum)
        let cb_vacuum = e_coeff.cb[0].get(0, 0, 0);
        // Inside region (dielectric)
        let cb_dielectric = e_coeff.cb[0].get(2, 2, 2);
        assert!(cb_dielectric < cb_vacuum);
        assert!((cb_dielectric / cb_vacuum - 0.25).abs() < 1e-5); // 1/4 for eps_r=4

        // H-coefficients should remain same for mu_r=1
        assert!((h_coeff.da[0].get(2, 2, 2) - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_set_material_with_losses() {
        let grid = create_test_grid();
        let mut op = Operator::new(grid, BoundaryConditions::default()).unwrap();

        // Set a lossy dielectric region
        op.set_material(2.0, 1.0, 0.01, 0.0, (1, 3, 1, 3, 1, 3));

        let e_coeff = op.e_coefficients();

        // With electric conductivity, Ca should be < 1.0
        let ca_lossy = e_coeff.ca[0].get(2, 2, 2);
        assert!(ca_lossy < 1.0);
        assert!(ca_lossy > 0.0);
    }

    #[test]
    fn test_set_material_magnetic() {
        let grid = create_test_grid();
        let mut op = Operator::new(grid, BoundaryConditions::default()).unwrap();

        // Set a magnetic material region with high magnetic conductivity
        // Need higher sigma_m for Da < 1.0 to be observable
        op.set_material(1.0, 4.0, 0.0, 1e6, (1, 3, 1, 3, 1, 3));

        let h_coeff = op.h_coefficients();

        // With high magnetic losses, Da should be < 1.0
        let da_lossy = h_coeff.da[0].get(2, 2, 2);
        assert!(da_lossy < 1.0);
        assert!(da_lossy > 0.0);

        // Db should be reduced due to mu_r = 4
        let db_outside = h_coeff.db[0].get(0, 0, 0);
        let db_inside = h_coeff.db[0].get(2, 2, 2);
        // Both negative, but inside should have smaller magnitude
        assert!(db_inside.abs() < db_outside.abs());
    }

    #[test]
    fn test_set_pec() {
        let grid = create_test_grid();
        let mut op = Operator::new(grid, BoundaryConditions::default()).unwrap();

        // Set a PEC region
        op.set_pec((1, 3, 1, 3, 1, 3));

        let e_coeff = op.e_coefficients();

        // Inside PEC: Ca and Cb should be 0 (E-field forced to zero)
        assert!((e_coeff.ca[0].get(2, 2, 2)).abs() < 1e-10);
        assert!((e_coeff.cb[0].get(2, 2, 2)).abs() < 1e-10);
        assert!((e_coeff.ca[1].get(2, 2, 2)).abs() < 1e-10);
        assert!((e_coeff.cb[1].get(2, 2, 2)).abs() < 1e-10);
        assert!((e_coeff.ca[2].get(2, 2, 2)).abs() < 1e-10);
        assert!((e_coeff.cb[2].get(2, 2, 2)).abs() < 1e-10);

        // Outside PEC: should still be vacuum
        assert!((e_coeff.ca[0].get(0, 0, 0) - 1.0).abs() < 1e-6);
        assert!(e_coeff.cb[0].get(0, 0, 0) > 0.0);
    }

    #[test]
    fn test_nyquist_timesteps() {
        let grid = create_test_grid();
        let op = Operator::new(grid, BoundaryConditions::default()).unwrap();

        // Test at 1 GHz
        let freq = 1e9;
        let steps = op.nyquist_timesteps(freq);

        // Period = 1ns, timestep is ~ps range, so should need many steps
        assert!(steps > 100);

        // Higher frequency should require fewer steps per period
        let steps_high = op.nyquist_timesteps(10e9);
        assert!(steps_high < steps);
    }

    #[test]
    fn test_grid_accessor() {
        let grid = create_test_grid();
        let op = Operator::new(grid, BoundaryConditions::default()).unwrap();

        let grid_ref = op.grid();
        let (dx, dy, dz) = grid_ref.cell_size();

        // All cell sizes should be 1mm
        assert!((dx - 0.001).abs() < 1e-10);
        assert!((dy - 0.001).abs() < 1e-10);
        assert!((dz - 0.001).abs() < 1e-10);
    }

    #[test]
    fn test_boundaries_accessor() {
        let grid = create_test_grid();
        let boundaries = BoundaryConditions::default();
        let op = Operator::new(grid, boundaries).unwrap();

        // Just verify we can access boundaries
        let _bc = op.boundaries();
    }

    #[test]
    fn test_timestep_cfl_condition() {
        // Test with anisotropic grid
        let grid = Grid::new(
            CoordinateSystem::Cartesian,
            vec![0.0, 0.001, 0.002, 0.003], // 1mm cells
            vec![0.0, 0.002, 0.004, 0.006], // 2mm cells
            vec![0.0, 0.0005, 0.001],       // 0.5mm cells
        );
        let op = Operator::new(grid, BoundaryConditions::default()).unwrap();

        // Timestep should be limited by smallest cell (0.5mm)
        let dt = op.timestep();
        let dt_max_smallest = 1.0 / (C0 * (1.0 / 0.001f64.powi(2) + 1.0 / 0.002f64.powi(2) + 1.0 / 0.0005f64.powi(2)).sqrt());

        // dt should be less than dt_max (with safety margin)
        assert!(dt < dt_max_smallest);
        assert!(dt > 0.9 * dt_max_smallest); // But not too small (we use 0.99 factor)
    }

    #[test]
    fn test_e_coefficient_clone() {
        let grid = create_test_grid();
        let op = Operator::new(grid, BoundaryConditions::default()).unwrap();

        let e_coeff = op.e_coefficients();
        let e_coeff_clone = e_coeff.clone();

        // Verify clone has same values
        assert!((e_coeff_clone.ca[0].get(2, 2, 2) - e_coeff.ca[0].get(2, 2, 2)).abs() < 1e-10);
        assert!((e_coeff_clone.cb[0].get(2, 2, 2) - e_coeff.cb[0].get(2, 2, 2)).abs() < 1e-10);
    }

    #[test]
    fn test_h_coefficient_clone() {
        let grid = create_test_grid();
        let op = Operator::new(grid, BoundaryConditions::default()).unwrap();

        let h_coeff = op.h_coefficients();
        let h_coeff_clone = h_coeff.clone();

        // Verify clone has same values
        assert!((h_coeff_clone.da[0].get(2, 2, 2) - h_coeff.da[0].get(2, 2, 2)).abs() < 1e-10);
        assert!((h_coeff_clone.db[0].get(2, 2, 2) - h_coeff.db[0].get(2, 2, 2)).abs() < 1e-10);
    }
}
