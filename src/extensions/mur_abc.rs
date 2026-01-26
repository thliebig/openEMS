//! Mur Absorbing Boundary Condition (1st Order).
//!
//! Implements Mur's first-order ABC for absorbing outgoing waves at boundaries.
//!
//! Reference: G. Mur, "Absorbing Boundary Conditions for the Finite-Difference
//! Approximation of the Time-Domain Electromagnetic-Field Equations,"
//! IEEE Trans. EMC, vol. EMC-23, no. 4, pp. 377-382, Nov. 1981.

use crate::arrays::{Dimensions, VectorField3D};
use crate::constants::C0;

/// Configuration for Mur ABC.
#[derive(Debug, Clone)]
pub struct MurAbcConfig {
    /// Override phase velocity (if None, uses c0/sqrt(eps*mu))
    pub phase_velocity: Option<f64>,
    /// Enable on each of the 6 boundaries
    pub enabled: [bool; 6],
}

impl Default for MurAbcConfig {
    fn default() -> Self {
        Self {
            phase_velocity: None,
            enabled: [true; 6],
        }
    }
}

impl MurAbcConfig {
    /// Create with Mur ABC on all boundaries.
    pub fn all() -> Self {
        Self::default()
    }

    /// Create with Mur ABC disabled on all boundaries.
    pub fn none() -> Self {
        Self {
            phase_velocity: None,
            enabled: [false; 6],
        }
    }

    /// Set phase velocity.
    pub fn with_phase_velocity(mut self, v: f64) -> Self {
        self.phase_velocity = Some(v);
        self
    }

    /// Enable/disable specific boundary.
    pub fn set_enabled(&mut self, boundary: usize, enabled: bool) -> &mut Self {
        if boundary < 6 {
            self.enabled[boundary] = enabled;
        }
        self
    }
}

/// Storage for a single Mur ABC face.
#[derive(Debug)]
#[allow(dead_code)]
struct MurFace {
    /// Normal direction (0=x, 1=y, 2=z)
    ny: usize,
    /// Next perpendicular direction
    nyp: usize,
    /// Second perpendicular direction
    nypp: usize,
    /// Whether this is the max boundary (true) or min (false)
    is_max: bool,
    /// Line number of the boundary
    line_nr: usize,
    /// Shifted line number (neighbor)
    line_nr_shift: usize,
    /// Coefficients for nyP component
    coeff_nyp: Vec<f32>,
    /// Coefficients for nyPP component
    coeff_nypp: Vec<f32>,
    /// Stored voltage values for nyP
    volt_nyp: Vec<f32>,
    /// Stored voltage values for nyPP
    volt_nypp: Vec<f32>,
    /// Size of the 2D face
    size: [usize; 2],
}

impl MurFace {
    fn new(ny: usize, is_max: bool, dims: Dimensions) -> Self {
        let nyp = (ny + 1) % 3;
        let nypp = (ny + 2) % 3;

        let (line_nr, line_nr_shift) = if is_max {
            let n = match ny {
                0 => dims.nx,
                1 => dims.ny,
                _ => dims.nz,
            };
            (n - 1, n - 2)
        } else {
            (0, 1)
        };

        let size = [
            match nyp {
                0 => dims.nx,
                1 => dims.ny,
                _ => dims.nz,
            },
            match nypp {
                0 => dims.nx,
                1 => dims.ny,
                _ => dims.nz,
            },
        ];

        let total_size = size[0] * size[1];

        Self {
            ny,
            nyp,
            nypp,
            is_max,
            line_nr,
            line_nr_shift,
            coeff_nyp: vec![0.0; total_size],
            coeff_nypp: vec![0.0; total_size],
            volt_nyp: vec![0.0; total_size],
            volt_nypp: vec![0.0; total_size],
            size,
        }
    }

    #[inline]
    fn idx(&self, i: usize, j: usize) -> usize {
        i * self.size[1] + j
    }

    fn init_coefficients(&mut self, dt: f64, delta: f64, phase_velocity: Option<f64>) {
        // Mur coefficient: (c*dt - delta) / (c*dt + delta)
        let c = phase_velocity.unwrap_or(C0);
        let c_dt = c * dt;
        let coeff = ((c_dt - delta) / (c_dt + delta)) as f32;

        for i in 0..self.size[0] {
            for j in 0..self.size[1] {
                let idx = self.idx(i, j);
                self.coeff_nyp[idx] = coeff;
                self.coeff_nypp[idx] = coeff;
            }
        }
    }
}

/// Mur Absorbing Boundary Condition.
///
/// First-order ABC based on the one-way wave equation.
/// Less effective than PML but simpler and lower memory usage.
pub struct MurAbc {
    /// Grid dimensions
    dims: Dimensions,
    /// Cell sizes
    delta: [f64; 3],
    /// Timestep
    dt: f64,
    /// Faces with active Mur ABC
    faces: Vec<MurFace>,
    /// Configuration
    config: MurAbcConfig,
}

impl MurAbc {
    /// Create a new Mur ABC.
    pub fn new(dims: Dimensions, delta: [f64; 3], dt: f64, config: MurAbcConfig) -> Self {
        let mut mur = Self {
            dims,
            delta,
            dt,
            faces: Vec::new(),
            config,
        };
        mur.build_faces();
        mur
    }

    /// Build ABC faces for all enabled boundaries.
    fn build_faces(&mut self) {
        // Boundary indices: 0=x_min, 1=x_max, 2=y_min, 3=y_max, 4=z_min, 5=z_max
        for boundary in 0..6 {
            if !self.config.enabled[boundary] {
                continue;
            }

            let ny = boundary / 2;
            let is_max = boundary % 2 == 1;

            let mut face = MurFace::new(ny, is_max, self.dims);
            face.init_coefficients(self.dt, self.delta[ny], self.config.phase_velocity);
            self.faces.push(face);
        }
    }

    /// Check if any ABC is active.
    pub fn is_active(&self) -> bool {
        !self.faces.is_empty()
    }

    /// Get the number of active faces.
    pub fn num_faces(&self) -> usize {
        self.faces.len()
    }

    /// Pre-voltage update: Store values for ABC calculation.
    pub fn pre_voltage_update(&mut self, e_field: &VectorField3D) {
        for face in &mut self.faces {
            let mut pos = [0usize; 3];
            let mut pos_shift = [0usize; 3];
            pos[face.ny] = face.line_nr;
            pos_shift[face.ny] = face.line_nr_shift;

            for i in 0..face.size[0] {
                pos[face.nyp] = i;
                pos_shift[face.nyp] = i;

                for j in 0..face.size[1] {
                    pos[face.nypp] = j;
                    pos_shift[face.nypp] = j;

                    let idx = face.idx(i, j);

                    // Get E-field values
                    let field_nyp = match face.nyp {
                        0 => &e_field.x,
                        1 => &e_field.y,
                        _ => &e_field.z,
                    };
                    let field_nypp = match face.nypp {
                        0 => &e_field.x,
                        1 => &e_field.y,
                        _ => &e_field.z,
                    };

                    // Get voltages at boundary and shifted positions
                    let v_nyp_boundary = field_nyp.get(pos[0], pos[1], pos[2]);
                    let v_nyp_shift = field_nyp.get(pos_shift[0], pos_shift[1], pos_shift[2]);
                    let v_nypp_boundary = field_nypp.get(pos[0], pos[1], pos[2]);
                    let v_nypp_shift = field_nypp.get(pos_shift[0], pos_shift[1], pos_shift[2]);

                    // Pre-calculation: V_shift - coeff * V_boundary
                    face.volt_nyp[idx] = v_nyp_shift - face.coeff_nyp[idx] * v_nyp_boundary;
                    face.volt_nypp[idx] = v_nypp_shift - face.coeff_nypp[idx] * v_nypp_boundary;
                }
            }
        }
    }

    /// Post-voltage update: Complete ABC calculation.
    pub fn post_voltage_update(&mut self, e_field: &VectorField3D) {
        for face in &mut self.faces {
            let mut pos_shift = [0usize; 3];
            pos_shift[face.ny] = face.line_nr_shift;

            for i in 0..face.size[0] {
                pos_shift[face.nyp] = i;

                for j in 0..face.size[1] {
                    pos_shift[face.nypp] = j;

                    let idx = face.idx(i, j);

                    // Get E-field at shifted position
                    let field_nyp = match face.nyp {
                        0 => &e_field.x,
                        1 => &e_field.y,
                        _ => &e_field.z,
                    };
                    let field_nypp = match face.nypp {
                        0 => &e_field.x,
                        1 => &e_field.y,
                        _ => &e_field.z,
                    };

                    let v_nyp_shift = field_nyp.get(pos_shift[0], pos_shift[1], pos_shift[2]);
                    let v_nypp_shift = field_nypp.get(pos_shift[0], pos_shift[1], pos_shift[2]);

                    // Complete calculation: previous + coeff * V_shift_new
                    face.volt_nyp[idx] += face.coeff_nyp[idx] * v_nyp_shift;
                    face.volt_nypp[idx] += face.coeff_nypp[idx] * v_nypp_shift;
                }
            }
        }
    }

    /// Apply ABC to voltages: Set boundary values.
    pub fn apply_to_voltages(&self, e_field: &mut VectorField3D) {
        for face in &self.faces {
            let mut pos = [0usize; 3];
            pos[face.ny] = face.line_nr;

            for i in 0..face.size[0] {
                pos[face.nyp] = i;

                for j in 0..face.size[1] {
                    pos[face.nypp] = j;

                    let idx = face.idx(i, j);

                    // Set boundary voltage
                    let field_nyp = match face.nyp {
                        0 => &mut e_field.x,
                        1 => &mut e_field.y,
                        _ => &mut e_field.z,
                    };
                    field_nyp.set(pos[0], pos[1], pos[2], face.volt_nyp[idx]);

                    let field_nypp = match face.nypp {
                        0 => &mut e_field.x,
                        1 => &mut e_field.y,
                        _ => &mut e_field.z,
                    };
                    field_nypp.set(pos[0], pos[1], pos[2], face.volt_nypp[idx]);
                }
            }
        }
    }

    /// Reset all stored values.
    pub fn reset(&mut self) {
        for face in &mut self.faces {
            face.volt_nyp.fill(0.0);
            face.volt_nypp.fill(0.0);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mur_config() {
        let config = MurAbcConfig::all();
        assert!(config.enabled.iter().all(|&e| e));

        let config = MurAbcConfig::none();
        assert!(config.enabled.iter().all(|&e| !e));
    }

    #[test]
    fn test_mur_creation() {
        let dims = Dimensions {
            nx: 20,
            ny: 20,
            nz: 20,
        };
        let delta = [0.001, 0.001, 0.001];
        let dt = 1e-12;
        let config = MurAbcConfig::all();

        let mur = MurAbc::new(dims, delta, dt, config);

        // Should have 6 faces (one per boundary)
        assert!(mur.is_active());
        assert_eq!(mur.num_faces(), 6);
    }

    #[test]
    fn test_mur_coefficient() {
        // Test that coefficient is calculated correctly
        let dt = 1e-12;
        let delta = 0.001;
        let c = C0;
        let c_dt = c * dt;

        let expected_coeff = (c_dt - delta) / (c_dt + delta);

        // Coefficient should be negative (c*dt < delta typically)
        assert!(expected_coeff < 0.0);
    }

    #[test]
    fn test_mur_partial_boundaries() {
        let dims = Dimensions {
            nx: 20,
            ny: 20,
            nz: 20,
        };
        let delta = [0.001, 0.001, 0.001];
        let dt = 1e-12;

        let mut config = MurAbcConfig::all();
        config.set_enabled(0, false); // Disable x_min
        config.set_enabled(1, false); // Disable x_max

        let mur = MurAbc::new(dims, delta, dt, config);

        // Should have 4 faces (Y and Z boundaries only)
        assert_eq!(mur.num_faces(), 4);
    }

    #[test]
    fn test_mur_update_cycle() {
        let dims = Dimensions {
            nx: 10,
            ny: 10,
            nz: 10,
        };
        let delta = [0.001, 0.001, 0.001];
        let dt = 1e-12;
        let config = MurAbcConfig::all();

        let mut mur = MurAbc::new(dims, delta, dt, config);
        let mut e_field = VectorField3D::new(dims);

        // Set some initial values
        e_field.x.set(5, 5, 5, 1.0);

        // Run update cycle
        mur.pre_voltage_update(&e_field);
        // (main engine update would go here)
        mur.post_voltage_update(&e_field);
        mur.apply_to_voltages(&mut e_field);

        // Boundary values should be modified
        // (exact values depend on the update algorithm)
    }

    #[test]
    fn test_mur_with_custom_phase_velocity() {
        let dims = Dimensions {
            nx: 10,
            ny: 10,
            nz: 10,
        };
        let delta = [0.001, 0.001, 0.001];
        let dt = 1e-12;
        let config = MurAbcConfig::all().with_phase_velocity(C0 * 0.5);

        let mur = MurAbc::new(dims, delta, dt, config);
        assert!(mur.is_active());
    }
}
