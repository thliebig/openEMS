//! Cylindrical Multigrid FDTD Implementation.
//!
//! Implements hierarchical mesh refinement for cylindrical coordinates,
//! allowing higher accuracy near the axis (R=0) or other regions of interest.

use crate::arrays::{Dimensions, VectorField3D};
use crate::constants::{C0, EPS0, MU0};
use std::f64::consts::PI;

/// Multigrid level configuration.
#[derive(Debug, Clone)]
pub struct MultigridLevel {
    /// Level index (0 = coarsest)
    pub level: usize,
    /// Radial range [r_min, r_max]
    pub r_range: [f64; 2],
    /// Number of cells in rho direction
    pub n_rho: usize,
    /// Number of cells in alpha direction
    pub n_alpha: usize,
    /// Number of cells in z direction
    pub n_z: usize,
    /// Refinement factor relative to parent level
    pub refinement: usize,
}

/// Cylindrical multigrid configuration.
#[derive(Debug, Clone)]
pub struct MultigridConfig {
    /// Grid levels
    pub levels: Vec<MultigridLevel>,
    /// Total radial extent
    pub r_max: f64,
    /// Z range
    pub z_range: [f64; 2],
    /// Whether alpha mesh is closed
    pub closed_alpha: bool,
}

impl MultigridConfig {
    /// Create a new multigrid configuration.
    pub fn new(r_max: f64, z_range: [f64; 2]) -> Self {
        Self {
            levels: Vec::new(),
            r_max,
            z_range,
            closed_alpha: true,
        }
    }

    /// Add a refinement level.
    pub fn add_level(&mut self, r_range: [f64; 2], n_rho: usize, n_alpha: usize, n_z: usize, refinement: usize) {
        let level = self.levels.len();
        self.levels.push(MultigridLevel {
            level,
            r_range,
            n_rho,
            n_alpha,
            n_z,
            refinement,
        });
    }

    /// Create a simple two-level multigrid.
    pub fn two_level(
        r_max: f64,
        z_range: [f64; 2],
        coarse_cells: [usize; 3],
        fine_r_max: f64,
        fine_cells: [usize; 3],
    ) -> Self {
        let mut config = Self::new(r_max, z_range);

        // Coarse level
        config.add_level(
            [0.0, r_max],
            coarse_cells[0],
            coarse_cells[1],
            coarse_cells[2],
            1,
        );

        // Fine level near axis
        config.add_level(
            [0.0, fine_r_max],
            fine_cells[0],
            fine_cells[1],
            fine_cells[2],
            2,
        );

        config
    }
}

/// Multigrid operator for a single level.
#[derive(Debug)]
struct LevelOperator {
    /// Level configuration
    config: MultigridLevel,
    /// Radial coordinates
    rho: Vec<f64>,
    /// Angular coordinates
    alpha: Vec<f64>,
    /// Z coordinates
    z: Vec<f64>,
    /// E-field coefficients (vv, vi)
    vv: VectorField3D,
    vi: VectorField3D,
    /// H-field coefficients (ii, iv)
    ii: VectorField3D,
    iv: VectorField3D,
    /// Timestep
    dt: f64,
}

impl LevelOperator {
    fn new(config: MultigridLevel, z_range: [f64; 2], closed_alpha: bool, dt: f64) -> Self {
        let dims = Dimensions::new(config.n_rho, config.n_alpha, config.n_z);

        // Generate mesh lines
        let dr = (config.r_range[1] - config.r_range[0]) / config.n_rho as f64;
        let da = 2.0 * PI / config.n_alpha as f64;
        let dz = (z_range[1] - z_range[0]) / config.n_z as f64;

        let rho: Vec<f64> = (0..=config.n_rho)
            .map(|i| config.r_range[0] + i as f64 * dr)
            .collect();
        let alpha: Vec<f64> = (0..=config.n_alpha)
            .map(|i| i as f64 * da)
            .collect();
        let z: Vec<f64> = (0..=config.n_z)
            .map(|i| z_range[0] + i as f64 * dz)
            .collect();

        // Initialize coefficients
        let mut vv = VectorField3D::new(dims);
        let mut vi = VectorField3D::new(dims);
        let mut ii = VectorField3D::new(dims);
        let mut iv = VectorField3D::new(dims);

        // Calculate coefficients (simplified)
        for i in 0..config.n_rho {
            let r = rho[i] + dr * 0.5;
            for j in 0..config.n_alpha {
                for k in 0..config.n_z {
                    // VV coefficients (E-field update)
                    vv.x.set(i, j, k, 1.0);
                    vv.y.set(i, j, k, 1.0);
                    vv.z.set(i, j, k, 1.0);

                    // VI coefficients (dt / (eps * area))
                    let area_rho = r * da * dz;
                    let area_alpha = dr * dz;
                    let area_z = 0.5 * ((r + dr).powi(2) - r.powi(2)) * da;

                    vi.x.set(i, j, k, (dt / (EPS0 * area_rho)) as f32);
                    vi.y.set(i, j, k, (dt / (EPS0 * area_alpha)) as f32);
                    vi.z.set(i, j, k, (dt / (EPS0 * area_z)) as f32);

                    // II coefficients (H-field update)
                    ii.x.set(i, j, k, 1.0);
                    ii.y.set(i, j, k, 1.0);
                    ii.z.set(i, j, k, 1.0);

                    // IV coefficients (dt / (mu * length))
                    let len_rho = dr;
                    let len_alpha = r * da;
                    let len_z = dz;

                    iv.x.set(i, j, k, (dt / (MU0 * len_rho)) as f32);
                    iv.y.set(i, j, k, (dt / (MU0 * len_alpha.max(1e-10))) as f32);
                    iv.z.set(i, j, k, (dt / (MU0 * len_z)) as f32);
                }
            }
        }

        Self {
            config,
            rho,
            alpha,
            z,
            vv,
            vi,
            ii,
            iv,
            dt,
        }
    }

    fn dims(&self) -> Dimensions {
        self.vv.dims()
    }
}

/// Multigrid engine for a single level.
#[derive(Debug)]
struct LevelEngine {
    /// E-field
    e_field: VectorField3D,
    /// H-field
    h_field: VectorField3D,
}

impl LevelEngine {
    fn new(dims: Dimensions) -> Self {
        Self {
            e_field: VectorField3D::new(dims),
            h_field: VectorField3D::new(dims),
        }
    }
}

/// Cylindrical Multigrid FDTD.
pub struct CylindricalMultigrid {
    /// Configuration
    config: MultigridConfig,
    /// Level operators
    operators: Vec<LevelOperator>,
    /// Level engines
    engines: Vec<LevelEngine>,
    /// Timestep
    dt: f64,
    /// Current timestep number
    timestep: u64,
    /// Current time
    time: f64,
}

impl CylindricalMultigrid {
    /// Create a new cylindrical multigrid.
    pub fn new(config: MultigridConfig) -> Self {
        // Calculate stable timestep (CFL)
        let dt = Self::calculate_timestep(&config);

        // Create operators and engines for each level
        let mut operators = Vec::new();
        let mut engines = Vec::new();

        for level_config in &config.levels {
            let operator = LevelOperator::new(
                level_config.clone(),
                config.z_range,
                config.closed_alpha,
                dt,
            );
            let dims = operator.dims();
            engines.push(LevelEngine::new(dims));
            operators.push(operator);
        }

        Self {
            config,
            operators,
            engines,
            dt,
            timestep: 0,
            time: 0.0,
        }
    }

    /// Calculate stable timestep for the multigrid.
    fn calculate_timestep(config: &MultigridConfig) -> f64 {
        let mut dt_min = f64::MAX;

        for level in &config.levels {
            let dr = (level.r_range[1] - level.r_range[0]) / level.n_rho as f64;
            let da = 2.0 * PI / level.n_alpha as f64;
            let dz = (config.z_range[1] - config.z_range[0]) / level.n_z as f64;

            // Use minimum radius for arc length (excluding R=0)
            let r_min = level.r_range[0].max(dr * 0.5);
            let d_alpha = r_min * da;

            // CFL condition
            let dt_level = 1.0 / (C0 * (1.0 / (dr * dr) + 1.0 / (d_alpha * d_alpha) + 1.0 / (dz * dz)).sqrt());
            dt_min = dt_min.min(dt_level);
        }

        dt_min * 0.95 // Safety factor
    }

    /// Get number of levels.
    pub fn num_levels(&self) -> usize {
        self.operators.len()
    }

    /// Get timestep.
    pub fn dt(&self) -> f64 {
        self.dt
    }

    /// Get current timestep number.
    pub fn timestep(&self) -> u64 {
        self.timestep
    }

    /// Get current time.
    pub fn time(&self) -> f64 {
        self.time
    }

    /// Update E-field on a single level.
    fn update_e_level(&mut self, level: usize) {
        let operator = &self.operators[level];
        let engine = &mut self.engines[level];
        let dims = operator.dims();

        for i in 0..dims.nx {
            for j in 0..dims.ny {
                for k in 0..dims.nz {
                    // Get neighboring indices with wrapping for closed alpha
                    let j_prev = if j == 0 {
                        if self.config.closed_alpha {
                            dims.ny - 1
                        } else {
                            0
                        }
                    } else {
                        j - 1
                    };
                    let k_prev = k.saturating_sub(1);
                    let i_prev = i.saturating_sub(1);

                    // E_rho update
                    let curl_h_rho = engine.h_field.z.get(i, j, k)
                        - engine.h_field.z.get(i, j_prev, k)
                        - (engine.h_field.y.get(i, j, k) - engine.h_field.y.get(i, j, k_prev));
                    let new_e_rho = operator.vv.x.get(i, j, k) * engine.e_field.x.get(i, j, k)
                        + operator.vi.x.get(i, j, k) * curl_h_rho;
                    engine.e_field.x.set(i, j, k, new_e_rho);

                    // E_alpha update
                    let curl_h_alpha = engine.h_field.x.get(i, j, k)
                        - engine.h_field.x.get(i, j, k_prev)
                        - (engine.h_field.z.get(i, j, k) - engine.h_field.z.get(i_prev, j, k));
                    let new_e_alpha = operator.vv.y.get(i, j, k) * engine.e_field.y.get(i, j, k)
                        + operator.vi.y.get(i, j, k) * curl_h_alpha;
                    engine.e_field.y.set(i, j, k, new_e_alpha);

                    // E_z update
                    let curl_h_z = engine.h_field.y.get(i, j, k)
                        - engine.h_field.y.get(i_prev, j, k)
                        - (engine.h_field.x.get(i, j, k) - engine.h_field.x.get(i, j_prev, k));
                    let new_e_z = operator.vv.z.get(i, j, k) * engine.e_field.z.get(i, j, k)
                        + operator.vi.z.get(i, j, k) * curl_h_z;
                    engine.e_field.z.set(i, j, k, new_e_z);
                }
            }
        }
    }

    /// Update H-field on a single level.
    fn update_h_level(&mut self, level: usize) {
        let operator = &self.operators[level];
        let engine = &mut self.engines[level];
        let dims = operator.dims();

        for i in 0..dims.nx {
            for j in 0..dims.ny {
                for k in 0..dims.nz {
                    let j_next = if j + 1 >= dims.ny {
                        if self.config.closed_alpha {
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
                    let curl_e_rho = engine.e_field.y.get(i, j, k_next)
                        - engine.e_field.y.get(i, j, k)
                        - (engine.e_field.z.get(i, j_next, k) - engine.e_field.z.get(i, j, k));
                    let new_h_rho = operator.ii.x.get(i, j, k) * engine.h_field.x.get(i, j, k)
                        + operator.iv.x.get(i, j, k) * curl_e_rho;
                    engine.h_field.x.set(i, j, k, new_h_rho);

                    // H_alpha update
                    let curl_e_alpha = engine.e_field.z.get(i_next, j, k)
                        - engine.e_field.z.get(i, j, k)
                        - (engine.e_field.x.get(i, j, k_next) - engine.e_field.x.get(i, j, k));
                    let new_h_alpha = operator.ii.y.get(i, j, k) * engine.h_field.y.get(i, j, k)
                        + operator.iv.y.get(i, j, k) * curl_e_alpha;
                    engine.h_field.y.set(i, j, k, new_h_alpha);

                    // H_z update
                    let curl_e_z = engine.e_field.x.get(i, j_next, k)
                        - engine.e_field.x.get(i, j, k)
                        - (engine.e_field.y.get(i_next, j, k) - engine.e_field.y.get(i, j, k));
                    let new_h_z = operator.ii.z.get(i, j, k) * engine.h_field.z.get(i, j, k)
                        + operator.iv.z.get(i, j, k) * curl_e_z;
                    engine.h_field.z.set(i, j, k, new_h_z);
                }
            }
        }
    }

    /// Interpolate fields from coarse to fine level.
    fn interpolate_coarse_to_fine(&mut self, coarse_level: usize, fine_level: usize) {
        // Simple injection interpolation at boundaries
        // In practice, would use higher-order interpolation
        let fine_dims = self.operators[fine_level].dims();
        let coarse_dims = self.operators[coarse_level].dims();

        // For boundary cells of fine level, interpolate from coarse level
        // This is simplified - full implementation would handle proper spatial mapping
        let refinement = self.config.levels[fine_level].refinement;

        // Interpolate E-field at outer boundary of fine region
        let i_boundary = fine_dims.nx - 1;
        let i_coarse = i_boundary / refinement;

        if i_coarse < coarse_dims.nx {
            for j in 0..fine_dims.ny.min(coarse_dims.ny) {
                for k in 0..fine_dims.nz.min(coarse_dims.nz) {
                    // Simple copy (injection)
                    let e_coarse = self.engines[coarse_level].e_field.x.get(i_coarse, j, k);
                    self.engines[fine_level].e_field.x.set(i_boundary, j, k, e_coarse);
                }
            }
        }
    }

    /// Restrict fields from fine to coarse level.
    fn restrict_fine_to_coarse(&mut self, fine_level: usize, coarse_level: usize) {
        // Simple averaging restriction
        let fine_dims = self.operators[fine_level].dims();
        let coarse_dims = self.operators[coarse_level].dims();
        let refinement = self.config.levels[fine_level].refinement;

        // Average fine level values to coarse level
        for i in 0..coarse_dims.nx.min(fine_dims.nx / refinement) {
            for j in 0..coarse_dims.ny.min(fine_dims.ny) {
                for k in 0..coarse_dims.nz.min(fine_dims.nz) {
                    let i_fine = i * refinement;

                    // Average E-field
                    let mut sum = 0.0f32;
                    let mut count = 0;
                    for di in 0..refinement {
                        if i_fine + di < fine_dims.nx {
                            sum += self.engines[fine_level].e_field.x.get(i_fine + di, j, k);
                            count += 1;
                        }
                    }
                    if count > 0 {
                        self.engines[coarse_level].e_field.x.set(i, j, k, sum / count as f32);
                    }
                }
            }
        }
    }

    /// Run one full timestep on all levels.
    pub fn step(&mut self) {
        let num_levels = self.num_levels();

        // Update from coarsest to finest
        for level in 0..num_levels {
            // Multiple sub-steps for finer levels
            let sub_steps = if level == 0 {
                1
            } else {
                self.config.levels[level].refinement
            };

            for _ in 0..sub_steps {
                self.update_e_level(level);
                self.update_h_level(level);
            }

            // Interpolate to finer level
            if level + 1 < num_levels {
                self.interpolate_coarse_to_fine(level, level + 1);
            }
        }

        // Restrict from finest to coarsest
        for level in (1..num_levels).rev() {
            self.restrict_fine_to_coarse(level, level - 1);
        }

        self.timestep += 1;
        self.time += self.dt;
    }

    /// Get E-field for a level.
    pub fn e_field(&self, level: usize) -> Option<&VectorField3D> {
        self.engines.get(level).map(|e| &e.e_field)
    }

    /// Get H-field for a level.
    pub fn h_field(&self, level: usize) -> Option<&VectorField3D> {
        self.engines.get(level).map(|e| &e.h_field)
    }

    /// Get mutable E-field for a level.
    pub fn e_field_mut(&mut self, level: usize) -> Option<&mut VectorField3D> {
        self.engines.get_mut(level).map(|e| &mut e.e_field)
    }

    /// Reset all fields.
    pub fn reset(&mut self) {
        for engine in &mut self.engines {
            engine.e_field.zero();
            engine.h_field.zero();
        }
        self.timestep = 0;
        self.time = 0.0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_multigrid_config() {
        let config = MultigridConfig::two_level(
            0.01,
            [0.0, 0.02],
            [10, 16, 20],
            0.005,
            [10, 16, 20],
        );

        assert_eq!(config.levels.len(), 2);
        assert_eq!(config.levels[0].level, 0);
        assert_eq!(config.levels[1].level, 1);
    }

    #[test]
    fn test_multigrid_creation() {
        let config = MultigridConfig::two_level(
            0.01,
            [0.0, 0.02],
            [5, 8, 10],
            0.005,
            [5, 8, 10],
        );

        let multigrid = CylindricalMultigrid::new(config);

        assert_eq!(multigrid.num_levels(), 2);
        assert!(multigrid.dt() > 0.0);
    }

    #[test]
    fn test_multigrid_step() {
        let config = MultigridConfig::two_level(
            0.01,
            [0.0, 0.02],
            [5, 8, 10],
            0.005,
            [5, 8, 10],
        );

        let mut multigrid = CylindricalMultigrid::new(config);

        // Set initial field
        if let Some(e) = multigrid.e_field_mut(0) {
            e.z.set(2, 4, 5, 1.0);
        }

        // Run a few steps
        for _ in 0..5 {
            multigrid.step();
        }

        assert_eq!(multigrid.timestep(), 5);
        assert!(multigrid.time() > 0.0);
    }

    #[test]
    fn test_multigrid_reset() {
        let config = MultigridConfig::two_level(
            0.01,
            [0.0, 0.02],
            [5, 8, 10],
            0.005,
            [5, 8, 10],
        );

        let mut multigrid = CylindricalMultigrid::new(config);

        multigrid.step();
        multigrid.reset();

        assert_eq!(multigrid.timestep(), 0);
        assert_eq!(multigrid.time(), 0.0);
    }
}
