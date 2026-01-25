//! Near-Field to Far-Field (NF2FF) transformation.
//!
//! Computes far-field radiation patterns from near-field data
//! using the equivalence principle.

use num_complex::Complex64;
use std::f64::consts::PI;

use crate::constants::{C0, Z0};

/// Far-field result data.
#[derive(Debug, Clone)]
pub struct FarFieldResult {
    /// Theta angles (radians)
    pub theta: Vec<f64>,
    /// Phi angles (radians)
    pub phi: Vec<f64>,
    /// E-theta component (complex)
    pub e_theta: Vec<Vec<Complex64>>,
    /// E-phi component (complex)
    pub e_phi: Vec<Vec<Complex64>>,
    /// Frequency (Hz)
    pub frequency: f64,
}

impl FarFieldResult {
    /// Compute directivity pattern (linear scale).
    pub fn directivity(&self) -> Vec<Vec<f64>> {
        let mut dir = vec![vec![0.0; self.phi.len()]; self.theta.len()];
        let mut max_power = 0.0f64;

        // Find maximum radiation intensity
        for (i, row) in self.e_theta.iter().enumerate() {
            for (j, et) in row.iter().enumerate() {
                let ep = self.e_phi[i][j];
                let power = et.norm_sqr() + ep.norm_sqr();
                if power > max_power {
                    max_power = power;
                }
            }
        }

        // Normalize
        if max_power > 0.0 {
            for (i, row) in self.e_theta.iter().enumerate() {
                for (j, et) in row.iter().enumerate() {
                    let ep = self.e_phi[i][j];
                    let power = et.norm_sqr() + ep.norm_sqr();
                    dir[i][j] = power / max_power;
                }
            }
        }

        dir
    }

    /// Compute directivity in dB.
    pub fn directivity_db(&self) -> Vec<Vec<f64>> {
        self.directivity()
            .iter()
            .map(|row| row.iter().map(|&d| 10.0 * d.log10()).collect())
            .collect()
    }

    /// Get the maximum directivity value.
    pub fn max_directivity(&self) -> f64 {
        // TODO: Proper directivity calculation with integration
        1.0
    }
}

/// NF2FF calculation engine.
pub struct Nf2ff {
    /// Frequency for calculation
    frequency: f64,
    /// Wavenumber k = 2*pi*f/c
    k: f64,
    /// Near-field surface data (to be filled from simulation)
    surfaces: Vec<NearFieldSurface>,
}

/// Near-field surface data.
pub struct NearFieldSurface {
    /// Surface normal direction (0=x, 1=y, 2=z)
    pub normal_dir: usize,
    /// Position along normal
    pub position: f64,
    /// Grid of (E, H) field values
    pub e_tangential: Vec<Vec<Complex64>>,
    pub h_tangential: Vec<Vec<Complex64>>,
    /// Surface area element
    pub ds: f64,
}

impl Nf2ff {
    /// Create a new NF2FF calculator.
    pub fn new(frequency: f64) -> Self {
        let k = 2.0 * PI * frequency / C0;
        Self {
            frequency,
            k,
            surfaces: Vec::new(),
        }
    }

    /// Add a near-field surface.
    pub fn add_surface(&mut self, surface: NearFieldSurface) {
        self.surfaces.push(surface);
    }

    /// Calculate far-field at given angles.
    ///
    /// # Arguments
    /// * `theta` - Elevation angles (radians, 0 to pi)
    /// * `phi` - Azimuth angles (radians, 0 to 2*pi)
    pub fn calculate(&self, theta: &[f64], phi: &[f64]) -> FarFieldResult {
        let n_theta = theta.len();
        let n_phi = phi.len();

        let mut e_theta = vec![vec![Complex64::new(0.0, 0.0); n_phi]; n_theta];
        let mut e_phi = vec![vec![Complex64::new(0.0, 0.0); n_phi]; n_theta];

        // For each observation direction
        for (it, &th) in theta.iter().enumerate() {
            for (ip, &ph) in phi.iter().enumerate() {
                // Unit vector in observation direction
                let r_hat = (th.sin() * ph.cos(), th.sin() * ph.sin(), th.cos());

                // Theta unit vector
                let th_hat = (th.cos() * ph.cos(), th.cos() * ph.sin(), -th.sin());

                // Phi unit vector
                let ph_hat = (-ph.sin(), ph.cos(), 0.0);

                // Integrate over all surfaces
                let (n_th, n_ph) = self.integrate_surfaces(r_hat, th_hat, ph_hat);

                // Far-field: E = -j*k*eta/(4*pi) * N
                let factor =
                    Complex64::new(0.0, -self.k * Z0 / (4.0 * PI));

                e_theta[it][ip] = factor * n_th;
                e_phi[it][ip] = factor * n_ph;
            }
        }

        FarFieldResult {
            theta: theta.to_vec(),
            phi: phi.to_vec(),
            e_theta,
            e_phi,
            frequency: self.frequency,
        }
    }

    /// Integrate equivalent currents over all surfaces.
    fn integrate_surfaces(
        &self,
        _r_hat: (f64, f64, f64),
        _th_hat: (f64, f64, f64),
        _ph_hat: (f64, f64, f64),
    ) -> (Complex64, Complex64) {
        let n_theta = Complex64::new(0.0, 0.0);
        let n_phi = Complex64::new(0.0, 0.0);

        // TODO: Implement proper surface integration
        // This is a placeholder for the full NF2FF algorithm

        (n_theta, n_phi)
    }

    /// Calculate radiation pattern in principal planes.
    pub fn calculate_principal_planes(&self, n_points: usize) -> (FarFieldResult, FarFieldResult) {
        let theta: Vec<f64> = (0..n_points)
            .map(|i| PI * i as f64 / (n_points - 1) as f64)
            .collect();

        // E-plane (phi = 0)
        let e_plane = self.calculate(&theta, &[0.0]);

        // H-plane (phi = pi/2)
        let h_plane = self.calculate(&theta, &[PI / 2.0]);

        (e_plane, h_plane)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nf2ff_creation() {
        let nf2ff = Nf2ff::new(1e9);
        assert!((nf2ff.k - 2.0 * PI * 1e9 / C0).abs() < 1e-6);
    }
}
