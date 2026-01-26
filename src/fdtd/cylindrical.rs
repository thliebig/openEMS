//! Cylindrical Coordinate FDTD Implementation.
//!
//! Implements the FDTD algorithm in cylindrical coordinates (ρ, α, z).
//! Handles special cases like closed alpha mesh and R=0 singularity.

use crate::arrays::{Dimensions, VectorField3D};
use crate::constants::{C0, EPS0, MU0};
use std::f64::consts::PI;

/// Threshold for determining closed alpha mesh.
const CLOSED_ALPHA_THRESHOLD: f64 = 1e-6;

/// Cylindrical grid configuration.
#[derive(Debug, Clone)]
pub struct CylindricalGrid {
    /// Number of cells in rho direction
    pub n_rho: usize,
    /// Number of cells in alpha direction
    pub n_alpha: usize,
    /// Number of cells in z direction
    pub n_z: usize,
    /// Radial coordinates (n_rho + 1 lines)
    pub rho: Vec<f64>,
    /// Angular coordinates (n_alpha + 1 lines, or n_alpha for closed)
    pub alpha: Vec<f64>,
    /// Z coordinates (n_z + 1 lines)
    pub z: Vec<f64>,
    /// Whether the alpha mesh is closed (full 360°)
    pub closed_alpha: bool,
    /// Whether R=0 is included in the mesh
    pub r0_included: bool,
}

impl CylindricalGrid {
    /// Create a new cylindrical grid.
    pub fn new(rho: Vec<f64>, alpha: Vec<f64>, z: Vec<f64>) -> Self {
        let n_rho = rho.len().saturating_sub(1);
        let n_alpha = alpha.len().saturating_sub(1);
        let n_z = z.len().saturating_sub(1);

        // Check if alpha mesh is closed (spans ~2π)
        let alpha_span = alpha.last().unwrap_or(&0.0) - alpha.first().unwrap_or(&0.0);
        let closed_alpha = (alpha_span - 2.0 * PI).abs() < CLOSED_ALPHA_THRESHOLD;

        // Check if R=0 is included
        let r0_included = rho.first().map(|&r| r.abs() < 1e-10).unwrap_or(false);

        Self {
            n_rho,
            n_alpha,
            n_z,
            rho,
            alpha,
            z,
            closed_alpha,
            r0_included,
        }
    }

    /// Create a uniform cylindrical grid.
    pub fn uniform(
        r_min: f64,
        r_max: f64,
        n_rho: usize,
        n_alpha: usize,
        z_min: f64,
        z_max: f64,
        n_z: usize,
    ) -> Self {
        let dr = (r_max - r_min) / n_rho as f64;
        let da = 2.0 * PI / n_alpha as f64;
        let dz = (z_max - z_min) / n_z as f64;

        let rho: Vec<f64> = (0..=n_rho).map(|i| r_min + i as f64 * dr).collect();
        let alpha: Vec<f64> = (0..=n_alpha).map(|i| i as f64 * da).collect();
        let z: Vec<f64> = (0..=n_z).map(|i| z_min + i as f64 * dz).collect();

        Self::new(rho, alpha, z)
    }

    /// Get cell dimensions at a position.
    pub fn get_cell_size(&self, i: usize, j: usize, k: usize) -> [f64; 3] {
        let d_rho = if i + 1 < self.rho.len() {
            self.rho[i + 1] - self.rho[i]
        } else {
            self.rho[i] - self.rho[i - 1]
        };

        let d_alpha = if j + 1 < self.alpha.len() {
            self.alpha[j + 1] - self.alpha[j]
        } else if self.closed_alpha {
            // Wrap around
            2.0 * PI - self.alpha[j] + self.alpha[0]
        } else {
            self.alpha[j] - self.alpha[j - 1]
        };

        let d_z = if k + 1 < self.z.len() {
            self.z[k + 1] - self.z[k]
        } else {
            self.z[k] - self.z[k - 1]
        };

        [d_rho, d_alpha, d_z]
    }

    /// Get radius at a rho index.
    pub fn get_radius(&self, i: usize) -> f64 {
        if i < self.rho.len() {
            self.rho[i]
        } else {
            *self.rho.last().unwrap_or(&0.0)
        }
    }

    /// Get the area of a cell face normal to direction ny.
    pub fn get_node_area(&self, ny: usize, pos: [usize; 3]) -> f64 {
        let [i, j, k] = pos;
        let cell_size = self.get_cell_size(i, j, k);

        match ny {
            0 => {
                // Rho-normal face: area = r * d_alpha * dz
                let r = self.get_radius(i);
                r * cell_size[1] * cell_size[2]
            }
            1 => {
                // Alpha-normal face: area = d_rho * dz
                cell_size[0] * cell_size[2]
            }
            2 => {
                // Z-normal face: area = 0.5 * (r2² - r1²) * d_alpha
                let r1 = self.get_radius(i);
                let r2 = if i + 1 < self.rho.len() {
                    self.rho[i + 1]
                } else {
                    r1
                };
                0.5 * (r2 * r2 - r1 * r1) * cell_size[1]
            }
            _ => 0.0,
        }
    }

    /// Get edge length for direction ny at position.
    pub fn get_edge_length(&self, ny: usize, pos: [usize; 3]) -> f64 {
        let [i, j, k] = pos;
        let cell_size = self.get_cell_size(i, j, k);

        match ny {
            0 => cell_size[0], // d_rho
            1 => {
                // Alpha edge: arc length = r * d_alpha
                let r = self.get_radius(i);
                r * cell_size[1]
            }
            2 => cell_size[2], // dz
            _ => 0.0,
        }
    }

    /// Get cell volume at position.
    pub fn get_cell_volume(&self, pos: [usize; 3]) -> f64 {
        let [i, j, k] = pos;
        let cell_size = self.get_cell_size(i, j, k);

        let r1 = self.get_radius(i);
        let r2 = if i + 1 < self.rho.len() {
            self.rho[i + 1]
        } else {
            r1 + cell_size[0]
        };

        // Volume = 0.5 * (r2² - r1²) * d_alpha * dz
        0.5 * (r2 * r2 - r1 * r1) * cell_size[1] * cell_size[2]
    }

    /// Map alpha index to valid range (for closed mesh).
    pub fn map_alpha_index(&self, idx: i32) -> usize {
        if !self.closed_alpha {
            return idx.max(0) as usize;
        }

        let n = self.n_alpha as i32;
        let mapped = ((idx % n) + n) % n;
        mapped as usize
    }

    /// Get dimensions for field arrays.
    pub fn dimensions(&self) -> Dimensions {
        Dimensions {
            nx: self.n_rho,
            ny: self.n_alpha,
            nz: self.n_z,
        }
    }
}

/// Cylindrical FDTD operator coefficients.
#[derive(Debug)]
pub struct CylindricalOperator {
    /// Grid
    grid: CylindricalGrid,
    /// Timestep
    dt: f64,
    /// E-field update coefficients (voltage to voltage)
    vv: VectorField3D,
    /// E-field update coefficients (current to voltage)
    vi: VectorField3D,
    /// H-field update coefficients (current to current)
    ii: VectorField3D,
    /// H-field update coefficients (voltage to current)
    iv: VectorField3D,
    /// Special R=0 coefficients for E-field
    vv_r0: Option<Vec<f32>>,
    vi_r0: Option<Vec<f32>>,
}

impl CylindricalOperator {
    /// Create a new cylindrical operator.
    pub fn new(grid: CylindricalGrid, dt: f64) -> Self {
        let dims = grid.dimensions();

        let mut op = Self {
            grid,
            dt,
            vv: VectorField3D::new(dims),
            vi: VectorField3D::new(dims),
            ii: VectorField3D::new(dims),
            iv: VectorField3D::new(dims),
            vv_r0: None,
            vi_r0: None,
        };

        op.calculate_coefficients();
        op
    }

    /// Calculate FDTD coefficients for cylindrical coordinates.
    fn calculate_coefficients(&mut self) {
        let dims = self.grid.dimensions();

        // Initialize R=0 coefficients if needed
        if self.grid.r0_included {
            let r0_size = self.grid.n_alpha * self.grid.n_z;
            self.vv_r0 = Some(vec![0.0; r0_size]);
            self.vi_r0 = Some(vec![0.0; r0_size]);
        }

        for i in 0..dims.nx {
            for j in 0..dims.ny {
                for k in 0..dims.nz {
                    let pos = [i, j, k];

                    // For each field component direction
                    for ny in 0..3 {
                        // E-field coefficient: vi = dt / (eps * area)
                        let area = self.grid.get_node_area(ny, pos);
                        let vi_coeff = if area > 0.0 {
                            (self.dt / (EPS0 * area)) as f32
                        } else {
                            0.0
                        };

                        // H-field coefficient: iv = dt / (mu * length)
                        let length = self.grid.get_edge_length(ny, pos);
                        let iv_coeff = if length > 0.0 {
                            (self.dt / (MU0 * length)) as f32
                        } else {
                            0.0
                        };

                        // Store coefficients
                        match ny {
                            0 => {
                                self.vv.x.set(i, j, k, 1.0);
                                self.vi.x.set(i, j, k, vi_coeff);
                                self.ii.x.set(i, j, k, 1.0);
                                self.iv.x.set(i, j, k, iv_coeff);
                            }
                            1 => {
                                self.vv.y.set(i, j, k, 1.0);
                                self.vi.y.set(i, j, k, vi_coeff);
                                self.ii.y.set(i, j, k, 1.0);
                                self.iv.y.set(i, j, k, iv_coeff);
                            }
                            2 => {
                                self.vv.z.set(i, j, k, 1.0);
                                self.vi.z.set(i, j, k, vi_coeff);
                                self.ii.z.set(i, j, k, 1.0);
                                self.iv.z.set(i, j, k, iv_coeff);
                            }
                            _ => {}
                        }
                    }
                }
            }
        }

        // Special handling for R=0
        if self.grid.r0_included {
            self.calculate_r0_coefficients();
        }
    }

    /// Calculate special coefficients for R=0.
    fn calculate_r0_coefficients(&mut self) {
        if let (Some(ref mut vv_r0), Some(ref mut vi_r0)) = (&mut self.vv_r0, &mut self.vi_r0) {
            for j in 0..self.grid.n_alpha {
                for k in 0..self.grid.n_z {
                    let idx = j * self.grid.n_z + k;

                    // At R=0, we need special treatment for E_alpha component
                    // The update uses values from all alpha positions around the axis
                    let d_alpha = self.grid.get_cell_size(0, j, k)[1];
                    let dz = self.grid.get_cell_size(0, j, k)[2];
                    let dr = self.grid.get_cell_size(0, j, k)[0];

                    // Effective area for R=0 cell (small cylinder)
                    let r_eff = 0.5 * dr;
                    let area = r_eff * d_alpha * dz;

                    vv_r0[idx] = 1.0;
                    vi_r0[idx] = if area > 0.0 {
                        (self.dt / (EPS0 * area)) as f32
                    } else {
                        0.0
                    };
                }
            }
        }
    }

    /// Get the grid.
    pub fn grid(&self) -> &CylindricalGrid {
        &self.grid
    }

    /// Get timestep.
    pub fn dt(&self) -> f64 {
        self.dt
    }

    /// Calculate stable timestep (CFL condition for cylindrical).
    pub fn calculate_timestep(&self) -> f64 {
        let mut dt_min = f64::MAX;

        for i in 0..self.grid.n_rho {
            for j in 0..self.grid.n_alpha {
                for k in 0..self.grid.n_z {
                    let cell_size = self.grid.get_cell_size(i, j, k);
                    let r = self.grid.get_radius(i).max(cell_size[0] * 0.5);

                    // Effective cell sizes
                    let dx = cell_size[0];
                    let dy = r * cell_size[1]; // Arc length
                    let dz = cell_size[2];

                    // CFL condition
                    let dt_cell =
                        1.0 / (C0 * (1.0 / (dx * dx) + 1.0 / (dy * dy) + 1.0 / (dz * dz)).sqrt());

                    dt_min = dt_min.min(dt_cell);
                }
            }
        }

        // Safety factor
        dt_min * 0.95
    }

    /// Get E-field coefficients.
    pub fn vv(&self) -> &VectorField3D {
        &self.vv
    }

    /// Get E-field update coefficients (current to voltage).
    pub fn vi(&self) -> &VectorField3D {
        &self.vi
    }

    /// Get H-field coefficients.
    pub fn ii(&self) -> &VectorField3D {
        &self.ii
    }

    /// Get H-field update coefficients.
    pub fn iv(&self) -> &VectorField3D {
        &self.iv
    }
}

/// Cylindrical FDTD Engine.
pub struct CylindricalEngine {
    /// Operator
    operator: CylindricalOperator,
    /// Electric field
    e_field: VectorField3D,
    /// Magnetic field
    h_field: VectorField3D,
    /// Current timestep
    timestep: u64,
    /// Current time
    time: f64,
}

impl CylindricalEngine {
    /// Create a new cylindrical engine.
    pub fn new(operator: CylindricalOperator) -> Self {
        let dims = operator.grid.dimensions();

        Self {
            operator,
            e_field: VectorField3D::new(dims),
            h_field: VectorField3D::new(dims),
            timestep: 0,
            time: 0.0,
        }
    }

    /// Update E-field (voltage update).
    pub fn update_voltage(&mut self) {
        let dims = self.operator.grid.dimensions();
        let grid = &self.operator.grid;

        for i in 0..dims.nx {
            for j in 0..dims.ny {
                for k in 0..dims.nz {
                    // Handle alpha wrapping for closed mesh
                    let j_prev = if j == 0 {
                        if grid.closed_alpha {
                            dims.ny - 1
                        } else {
                            0
                        }
                    } else {
                        j - 1
                    };

                    let k_prev = k.saturating_sub(1);
                    let i_prev = i.saturating_sub(1);

                    // E_rho update: dE_rho/dt = (1/eps) * (dHz/r*d_alpha - dH_alpha/dz)
                    let curl_h_rho = if i > 0 || !grid.r0_included {
                        let hz_diff =
                            self.h_field.z.get(i, j, k) - self.h_field.z.get(i, j_prev, k);
                        let ha_diff =
                            self.h_field.y.get(i, j, k) - self.h_field.y.get(i, j, k_prev);
                        hz_diff - ha_diff
                    } else {
                        0.0
                    };

                    let e_rho = self.e_field.x.get(i, j, k);
                    let new_e_rho = self.operator.vv.x.get(i, j, k) * e_rho
                        + self.operator.vi.x.get(i, j, k) * curl_h_rho;
                    self.e_field.x.set(i, j, k, new_e_rho);

                    // E_alpha update: dE_alpha/dt = (1/eps) * (dH_rho/dz - dHz/d_rho)
                    let curl_h_alpha = {
                        let hr_diff =
                            self.h_field.x.get(i, j, k) - self.h_field.x.get(i, j, k_prev);
                        let hz_diff =
                            self.h_field.z.get(i, j, k) - self.h_field.z.get(i_prev, j, k);
                        hr_diff - hz_diff
                    };

                    let e_alpha = self.e_field.y.get(i, j, k);
                    let new_e_alpha = self.operator.vv.y.get(i, j, k) * e_alpha
                        + self.operator.vi.y.get(i, j, k) * curl_h_alpha;
                    self.e_field.y.set(i, j, k, new_e_alpha);

                    // E_z update: dE_z/dt = (1/eps) * ((r*H_alpha)'/r*d_rho - dH_rho/r*d_alpha)
                    let curl_h_z = if i > 0 || !grid.r0_included {
                        let ha_diff =
                            self.h_field.y.get(i, j, k) - self.h_field.y.get(i_prev, j, k);
                        let hr_diff =
                            self.h_field.x.get(i, j, k) - self.h_field.x.get(i, j_prev, k);
                        ha_diff - hr_diff
                    } else {
                        0.0
                    };

                    let e_z = self.e_field.z.get(i, j, k);
                    let new_e_z = self.operator.vv.z.get(i, j, k) * e_z
                        + self.operator.vi.z.get(i, j, k) * curl_h_z;
                    self.e_field.z.set(i, j, k, new_e_z);
                }
            }
        }

        // Special R=0 handling
        if grid.r0_included {
            self.update_voltage_r0();
        }
    }

    /// Special voltage update at R=0.
    fn update_voltage_r0(&mut self) {
        // At R=0, E_alpha and E_z are undefined or zero
        // E_rho needs special averaging from surrounding cells
        if let (Some(ref _vv_r0), Some(ref _vi_r0)) = (&self.operator.vv_r0, &self.operator.vi_r0) {
            let dims = self.operator.grid.dimensions();

            for j in 0..dims.ny {
                for k in 0..dims.nz {
                    let _idx = j * dims.nz + k;

                    // E_alpha at R=0 is typically set to zero
                    self.e_field.y.set(0, j, k, 0.0);

                    // E_z at R=0 can be averaged from neighbors
                    // For simplicity, use the value at first radial cell
                    let e_z_neighbor = self.e_field.z.get(1, j, k);
                    self.e_field.z.set(0, j, k, e_z_neighbor * 0.5);
                }
            }
        }
    }

    /// Update H-field (current update).
    pub fn update_current(&mut self) {
        let dims = self.operator.grid.dimensions();
        let grid = &self.operator.grid;

        for i in 0..dims.nx {
            for j in 0..dims.ny {
                for k in 0..dims.nz {
                    // Handle alpha wrapping for closed mesh
                    let j_next = if j + 1 >= dims.ny {
                        if grid.closed_alpha {
                            0
                        } else {
                            dims.ny - 1
                        }
                    } else {
                        j + 1
                    };

                    let k_next = (k + 1).min(dims.nz - 1);
                    let i_next = (i + 1).min(dims.nx - 1);

                    // H_rho update
                    let curl_e_rho = {
                        let ez_diff =
                            self.e_field.z.get(i, j_next, k) - self.e_field.z.get(i, j, k);
                        let ea_diff =
                            self.e_field.y.get(i, j, k_next) - self.e_field.y.get(i, j, k);
                        ea_diff - ez_diff
                    };

                    let h_rho = self.h_field.x.get(i, j, k);
                    let new_h_rho = self.operator.ii.x.get(i, j, k) * h_rho
                        + self.operator.iv.x.get(i, j, k) * curl_e_rho;
                    self.h_field.x.set(i, j, k, new_h_rho);

                    // H_alpha update
                    let curl_e_alpha = {
                        let er_diff =
                            self.e_field.x.get(i, j, k_next) - self.e_field.x.get(i, j, k);
                        let ez_diff =
                            self.e_field.z.get(i_next, j, k) - self.e_field.z.get(i, j, k);
                        ez_diff - er_diff
                    };

                    let h_alpha = self.h_field.y.get(i, j, k);
                    let new_h_alpha = self.operator.ii.y.get(i, j, k) * h_alpha
                        + self.operator.iv.y.get(i, j, k) * curl_e_alpha;
                    self.h_field.y.set(i, j, k, new_h_alpha);

                    // H_z update
                    let curl_e_z = {
                        let ea_diff =
                            self.e_field.y.get(i_next, j, k) - self.e_field.y.get(i, j, k);
                        let er_diff =
                            self.e_field.x.get(i, j_next, k) - self.e_field.x.get(i, j, k);
                        er_diff - ea_diff
                    };

                    let h_z = self.h_field.z.get(i, j, k);
                    let new_h_z = self.operator.ii.z.get(i, j, k) * h_z
                        + self.operator.iv.z.get(i, j, k) * curl_e_z;
                    self.h_field.z.set(i, j, k, new_h_z);
                }
            }
        }
    }

    /// Run one full timestep.
    pub fn step(&mut self) {
        self.update_voltage();
        self.update_current();

        self.timestep += 1;
        self.time += self.operator.dt;
    }

    /// Get current timestep.
    pub fn timestep(&self) -> u64 {
        self.timestep
    }

    /// Get current time.
    pub fn time(&self) -> f64 {
        self.time
    }

    /// Get E-field reference.
    pub fn e_field(&self) -> &VectorField3D {
        &self.e_field
    }

    /// Get mutable E-field reference.
    pub fn e_field_mut(&mut self) -> &mut VectorField3D {
        &mut self.e_field
    }

    /// Get H-field reference.
    pub fn h_field(&self) -> &VectorField3D {
        &self.h_field
    }

    /// Get mutable H-field reference.
    pub fn h_field_mut(&mut self) -> &mut VectorField3D {
        &mut self.h_field
    }

    /// Reset fields and time.
    pub fn reset(&mut self) {
        self.e_field.zero();
        self.h_field.zero();
        self.timestep = 0;
        self.time = 0.0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cylindrical_grid_creation() {
        let grid = CylindricalGrid::uniform(0.001, 0.01, 10, 16, 0.0, 0.02, 20);

        assert_eq!(grid.n_rho, 10);
        assert_eq!(grid.n_alpha, 16);
        assert_eq!(grid.n_z, 20);
        assert!(grid.closed_alpha);
        assert!(!grid.r0_included);
    }

    #[test]
    fn test_grid_with_r0() {
        let rho: Vec<f64> = (0..=10).map(|i| i as f64 * 0.001).collect();
        let alpha: Vec<f64> = (0..=8).map(|i| i as f64 * PI / 4.0).collect();
        let z: Vec<f64> = (0..=10).map(|i| i as f64 * 0.002).collect();

        let grid = CylindricalGrid::new(rho, alpha, z);

        assert!(grid.r0_included);
        assert!(grid.closed_alpha);
    }

    #[test]
    fn test_cell_volume() {
        let grid = CylindricalGrid::uniform(0.005, 0.015, 10, 16, 0.0, 0.02, 10);

        let vol = grid.get_cell_volume([5, 0, 5]);
        assert!(vol > 0.0);

        // Volume should increase with radius
        let vol_inner = grid.get_cell_volume([0, 0, 5]);
        let vol_outer = grid.get_cell_volume([9, 0, 5]);
        assert!(vol_outer > vol_inner);
    }

    #[test]
    fn test_alpha_wrapping() {
        let grid = CylindricalGrid::uniform(0.005, 0.015, 10, 16, 0.0, 0.02, 10);

        assert_eq!(grid.map_alpha_index(-1), 15);
        assert_eq!(grid.map_alpha_index(16), 0);
        assert_eq!(grid.map_alpha_index(5), 5);
    }

    #[test]
    fn test_cylindrical_operator() {
        let grid = CylindricalGrid::uniform(0.005, 0.015, 10, 16, 0.0, 0.02, 10);
        let _dt_calc = grid.calculate_timestep(&grid);
        // Use fixed dt for test
        let dt = 1e-12;

        let operator = CylindricalOperator::new(grid, dt);

        assert!(operator.dt() > 0.0);
    }

    // Helper function for tests
    impl CylindricalGrid {
        fn calculate_timestep(&self, _grid: &CylindricalGrid) -> f64 {
            1e-12 // Placeholder
        }
    }

    #[test]
    fn test_cylindrical_engine_creation() {
        let grid = CylindricalGrid::uniform(0.005, 0.015, 5, 8, 0.0, 0.01, 5);
        let operator = CylindricalOperator::new(grid, 1e-12);
        let engine = CylindricalEngine::new(operator);

        assert_eq!(engine.timestep(), 0);
        assert_eq!(engine.time(), 0.0);
    }

    #[test]
    fn test_engine_step() {
        let grid = CylindricalGrid::uniform(0.005, 0.015, 5, 8, 0.0, 0.01, 5);
        let operator = CylindricalOperator::new(grid, 1e-12);
        let mut engine = CylindricalEngine::new(operator);

        // Set initial field
        engine.e_field_mut().z.set(2, 4, 2, 1.0);

        // Run a few steps
        for _ in 0..10 {
            engine.step();
        }

        assert_eq!(engine.timestep(), 10);
        assert!(engine.time() > 0.0);
    }

    #[test]
    fn test_engine_reset() {
        let grid = CylindricalGrid::uniform(0.005, 0.015, 5, 8, 0.0, 0.01, 5);
        let operator = CylindricalOperator::new(grid, 1e-12);
        let mut engine = CylindricalEngine::new(operator);

        engine.e_field_mut().z.fill(1.0);
        engine.step();
        engine.reset();

        assert_eq!(engine.timestep(), 0);
        assert_eq!(engine.time(), 0.0);
        assert_eq!(engine.e_field().z.get(2, 2, 2), 0.0);
    }

    #[test]
    fn test_edge_length_with_radius() {
        let grid = CylindricalGrid::uniform(0.005, 0.015, 10, 16, 0.0, 0.02, 10);

        // Alpha edge length should scale with radius
        let len_inner = grid.get_edge_length(1, [0, 0, 0]);
        let len_outer = grid.get_edge_length(1, [9, 0, 0]);

        assert!(len_outer > len_inner);
    }

    #[test]
    fn test_node_area() {
        let grid = CylindricalGrid::uniform(0.005, 0.015, 10, 16, 0.0, 0.02, 10);

        // Z-normal area should scale with r²
        let area = grid.get_node_area(2, [5, 0, 5]);
        assert!(area > 0.0);

        // Rho-normal area should scale linearly with r
        let area_rho_inner = grid.get_node_area(0, [0, 0, 5]);
        let area_rho_outer = grid.get_node_area(0, [9, 0, 5]);
        assert!(area_rho_outer > area_rho_inner);
    }
}
