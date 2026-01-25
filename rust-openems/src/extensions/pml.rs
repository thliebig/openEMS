//! Perfectly Matched Layer (PML) absorbing boundary condition.
//!
//! Implements the Uniaxial PML (UPML) for absorbing outgoing waves
//! at the boundaries of the simulation domain.

use crate::arrays::{Dimensions, Field3D};
use crate::constants::{EPS0, MU0};

/// PML configuration parameters.
#[derive(Debug, Clone)]
pub struct PmlConfig {
    /// Number of PML layers
    pub layers: usize,
    /// Polynomial grading order
    pub grading_order: f64,
    /// Maximum conductivity (automatically calculated if None)
    pub sigma_max: Option<f64>,
    /// Target reflection coefficient (used to calculate sigma_max)
    pub reflection_coeff: f64,
}

impl Default for PmlConfig {
    fn default() -> Self {
        Self {
            layers: 8,
            grading_order: 3.0,
            sigma_max: None,
            reflection_coeff: 1e-6, // -120 dB reflection
        }
    }
}

impl PmlConfig {
    /// Create PML with specified number of layers.
    pub fn new(layers: usize) -> Self {
        Self {
            layers,
            ..Default::default()
        }
    }

    /// Calculate optimal sigma_max for given parameters.
    pub fn calculate_sigma_max(&self, delta: f64) -> f64 {
        if let Some(sigma) = self.sigma_max {
            return sigma;
        }

        // sigma_max = -(m+1) * ln(R) / (2 * eta * d)
        // where m = grading_order, R = reflection_coeff, d = PML thickness
        let eta = (MU0 / EPS0).sqrt(); // Free space impedance
        let d = self.layers as f64 * delta;
        -(self.grading_order + 1.0) * self.reflection_coeff.ln() / (2.0 * eta * d)
    }
}

/// PML absorbing boundary condition.
pub struct Pml {
    /// Configuration
    config: PmlConfig,
    /// Conductivity profile in each direction (graded)
    sigma_e: [Vec<f64>; 3],
    sigma_m: [Vec<f64>; 3],
    /// Split-field components for PML (Psi fields)
    psi_ex_y: Option<Field3D>,
    psi_ex_z: Option<Field3D>,
    psi_ey_x: Option<Field3D>,
    psi_ey_z: Option<Field3D>,
    psi_ez_x: Option<Field3D>,
    psi_ez_y: Option<Field3D>,
    psi_hx_y: Option<Field3D>,
    psi_hx_z: Option<Field3D>,
    psi_hy_x: Option<Field3D>,
    psi_hy_z: Option<Field3D>,
    psi_hz_x: Option<Field3D>,
    psi_hz_y: Option<Field3D>,
}

impl Pml {
    /// Create a new PML with the given configuration.
    pub fn new(config: PmlConfig, dims: Dimensions, cell_sizes: (f64, f64, f64)) -> Self {
        let (dx, dy, dz) = cell_sizes;
        let n = config.layers;

        // Calculate sigma profiles (polynomial grading)
        let sigma_x = Self::calculate_sigma_profile(&config, n, dx);
        let sigma_y = Self::calculate_sigma_profile(&config, n, dy);
        let sigma_z = Self::calculate_sigma_profile(&config, n, dz);

        // Create split-field storage for boundary regions
        // For now, create full-size fields (optimization: only allocate boundary regions)
        let create_field = || Some(Field3D::new(dims));

        Self {
            config,
            sigma_e: [sigma_x.clone(), sigma_y.clone(), sigma_z.clone()],
            sigma_m: [sigma_x, sigma_y, sigma_z],
            psi_ex_y: create_field(),
            psi_ex_z: create_field(),
            psi_ey_x: create_field(),
            psi_ey_z: create_field(),
            psi_ez_x: create_field(),
            psi_ez_y: create_field(),
            psi_hx_y: create_field(),
            psi_hx_z: create_field(),
            psi_hy_x: create_field(),
            psi_hy_z: create_field(),
            psi_hz_x: create_field(),
            psi_hz_y: create_field(),
        }
    }

    /// Calculate the conductivity profile for PML.
    fn calculate_sigma_profile(config: &PmlConfig, n: usize, delta: f64) -> Vec<f64> {
        let sigma_max = config.calculate_sigma_max(delta);
        let m = config.grading_order;

        (0..n)
            .map(|i| {
                // Polynomial grading: sigma(x) = sigma_max * (x/d)^m
                let x = (i as f64 + 0.5) / n as f64;
                sigma_max * x.powf(m)
            })
            .collect()
    }

    /// Get the number of PML layers.
    pub fn layers(&self) -> usize {
        self.config.layers
    }

    /// Update PML auxiliary fields for E-field update.
    ///
    /// This should be called after the H-field update and before the E-field update.
    pub fn update_e_fields(&mut self, _dims: Dimensions, _dt: f64) {
        // CPML update equations for Psi fields
        // Psi_ex_y[n+1] = b_y * Psi_ex_y[n] + a_y * (Hz[i,j,k] - Hz[i,j-1,k])
        // etc.

        // TODO: Implement full CPML update
        // For now, this is a placeholder
    }

    /// Update PML auxiliary fields for H-field update.
    pub fn update_h_fields(&mut self, _dims: Dimensions, _dt: f64) {
        // Similar CPML update for magnetic field Psi components
        // TODO: Implement full CPML update
    }

    /// Apply PML corrections to E-field update.
    pub fn apply_e_corrections(
        &self,
        _e_field: &mut crate::arrays::VectorField3D,
        _dims: Dimensions,
    ) {
        // E = E + Cb * (Psi_contributions)
        // TODO: Implement full CPML corrections
    }

    /// Apply PML corrections to H-field update.
    pub fn apply_h_corrections(
        &self,
        _h_field: &mut crate::arrays::VectorField3D,
        _dims: Dimensions,
    ) {
        // H = H + Db * (Psi_contributions)
        // TODO: Implement full CPML corrections
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pml_config() {
        let config = PmlConfig::new(8);
        assert_eq!(config.layers, 8);

        let sigma = config.calculate_sigma_max(0.001);
        assert!(sigma > 0.0);
    }

    #[test]
    fn test_sigma_profile() {
        let config = PmlConfig::new(10);
        let profile = Pml::calculate_sigma_profile(&config, 10, 0.001);

        // Profile should be monotonically increasing
        for i in 1..profile.len() {
            assert!(profile[i] >= profile[i - 1]);
        }
    }
}
