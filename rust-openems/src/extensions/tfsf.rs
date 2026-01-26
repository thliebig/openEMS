//! Total-Field/Scattered-Field (TF/SF) Boundary.
//!
//! Implements the TF/SF technique for injecting plane waves into the
//! FDTD domain while separating total and scattered field regions.

use crate::arrays::{Dimensions, VectorField3D};
use crate::constants::C0;
use std::f64::consts::PI;

/// Direction of wave propagation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PropagationDirection {
    /// Propagating in +x direction
    XPlus,
    /// Propagating in -x direction
    XMinus,
    /// Propagating in +y direction
    YPlus,
    /// Propagating in -y direction
    YMinus,
    /// Propagating in +z direction
    ZPlus,
    /// Propagating in -z direction
    ZMinus,
}

impl PropagationDirection {
    /// Get the axis index (0, 1, 2 for x, y, z).
    pub fn axis(&self) -> usize {
        match self {
            Self::XPlus | Self::XMinus => 0,
            Self::YPlus | Self::YMinus => 1,
            Self::ZPlus | Self::ZMinus => 2,
        }
    }

    /// Check if propagating in positive direction.
    pub fn is_positive(&self) -> bool {
        matches!(self, Self::XPlus | Self::YPlus | Self::ZPlus)
    }
}

/// Polarization of the incident wave.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Polarization {
    /// E-field along first perpendicular axis
    Transverse1,
    /// E-field along second perpendicular axis
    Transverse2,
}

/// Configuration for TF/SF boundary.
#[derive(Debug, Clone)]
pub struct TfsfConfig {
    /// Propagation direction
    pub direction: PropagationDirection,
    /// Polarization
    pub polarization: Polarization,
    /// Center frequency (Hz)
    pub frequency: f64,
    /// Amplitude (V/m)
    pub amplitude: f64,
    /// TF/SF boundary position (layer number from boundary)
    pub boundary_layer: usize,
    /// Phase velocity (defaults to c0)
    pub phase_velocity: Option<f64>,
}

impl TfsfConfig {
    /// Create a TF/SF configuration.
    pub fn new(direction: PropagationDirection, frequency: f64, amplitude: f64) -> Self {
        Self {
            direction,
            polarization: Polarization::Transverse1,
            frequency,
            amplitude,
            boundary_layer: 10,
            phase_velocity: None,
        }
    }

    /// Set polarization.
    pub fn with_polarization(mut self, pol: Polarization) -> Self {
        self.polarization = pol;
        self
    }

    /// Set boundary layer position.
    pub fn with_boundary_layer(mut self, layer: usize) -> Self {
        self.boundary_layer = layer;
        self
    }

    /// Set phase velocity.
    pub fn with_phase_velocity(mut self, v: f64) -> Self {
        self.phase_velocity = Some(v);
        self
    }
}

/// Internal storage for TF/SF field values.
#[derive(Debug)]
struct TfsfStorage {
    /// E-field component on source plane
    e_src: Vec<f32>,
    /// H-field component on source plane
    h_src: Vec<f32>,
    /// 1D auxiliary E-field for plane wave
    e_1d: Vec<f32>,
    /// 1D auxiliary H-field for plane wave
    h_1d: Vec<f32>,
}

/// Total-Field/Scattered-Field boundary.
pub struct TfsfBoundary {
    /// Grid dimensions
    dims: Dimensions,
    /// Cell sizes
    delta: [f64; 3],
    /// Timestep
    dt: f64,
    /// Configuration
    config: TfsfConfig,
    /// Storage for field values
    storage: TfsfStorage,
    /// Propagation axis
    prop_axis: usize,
    /// E-field axis
    e_axis: usize,
    /// H-field axis
    h_axis: usize,
    /// Boundary position
    boundary_pos: usize,
    /// Phase velocity
    phase_velocity: f64,
    /// Current timestep
    timestep: u64,
}

impl TfsfBoundary {
    /// Create a new TF/SF boundary.
    pub fn new(dims: Dimensions, delta: [f64; 3], dt: f64, config: TfsfConfig) -> Self {
        let prop_axis = config.direction.axis();
        let v = config.phase_velocity.unwrap_or(C0);

        // Determine E and H field axes based on polarization
        let (e_axis, h_axis) = match config.polarization {
            Polarization::Transverse1 => ((prop_axis + 1) % 3, (prop_axis + 2) % 3),
            Polarization::Transverse2 => ((prop_axis + 2) % 3, (prop_axis + 1) % 3),
        };

        // Boundary position
        let boundary_pos = config.boundary_layer;

        // Calculate plane size
        let plane_size = match prop_axis {
            0 => dims.ny * dims.nz,
            1 => dims.nx * dims.nz,
            _ => dims.nx * dims.ny,
        };

        // 1D grid size for auxiliary propagation
        let aux_size = match prop_axis {
            0 => dims.nx,
            1 => dims.ny,
            _ => dims.nz,
        };

        let storage = TfsfStorage {
            e_src: vec![0.0; plane_size],
            h_src: vec![0.0; plane_size],
            e_1d: vec![0.0; aux_size + 1],
            h_1d: vec![0.0; aux_size + 1],
        };

        Self {
            dims,
            delta,
            dt,
            config,
            storage,
            prop_axis,
            e_axis,
            h_axis,
            boundary_pos,
            phase_velocity: v,
            timestep: 0,
        }
    }

    /// Check if TF/SF is active.
    pub fn is_active(&self) -> bool {
        true
    }

    /// Calculate incident E-field at given time.
    #[allow(dead_code)]
    fn incident_e(&self, t: f64) -> f64 {
        let omega = 2.0 * PI * self.config.frequency;
        self.config.amplitude * (omega * t).sin()
    }

    /// Calculate incident H-field at given time.
    fn incident_h(&self, t: f64) -> f64 {
        let omega = 2.0 * PI * self.config.frequency;
        // H = E / eta, with 90-degree phase shift for plane wave
        let eta = (crate::constants::MU0 / crate::constants::EPS0).sqrt();
        (self.config.amplitude / eta) * (omega * t).sin()
    }

    /// Update 1D auxiliary grid.
    fn update_1d_grid(&mut self) {
        let dt = self.dt;
        let dx = self.delta[self.prop_axis];
        let c = self.phase_velocity;

        // 1D FDTD update
        // H update
        for i in 0..self.storage.h_1d.len() - 1 {
            let db = -(dt / (crate::constants::MU0 * dx)) as f32;
            self.storage.h_1d[i] += db * (self.storage.e_1d[i + 1] - self.storage.e_1d[i]);
        }

        // Inject source at left boundary
        let t = self.timestep as f64 * dt;
        self.storage.h_1d[0] = self.incident_h(t) as f32;

        // E update
        for i in 1..self.storage.e_1d.len() {
            let cb = (dt / (crate::constants::EPS0 * dx)) as f32;
            self.storage.e_1d[i] += cb * (self.storage.h_1d[i] - self.storage.h_1d[i - 1]);
        }

        // ABC at right boundary (simple first-order Mur)
        let coeff = ((c * dt - dx) / (c * dt + dx)) as f32;
        let n = self.storage.e_1d.len() - 1;
        self.storage.e_1d[n] = coeff * self.storage.e_1d[n - 1];
    }

    /// Pre-update for H-field: Add incident field.
    pub fn pre_update_h(&mut self, h_field: &mut VectorField3D) {
        // Update 1D auxiliary grid
        self.update_1d_grid();

        // Get incident field at boundary
        let e_inc = if self.boundary_pos < self.storage.e_1d.len() {
            self.storage.e_1d[self.boundary_pos]
        } else {
            0.0
        };

        // Add correction to H-field at TF/SF boundary
        // This subtracts the incident field contribution from the scattered field region
        let db = -(self.dt / (crate::constants::MU0 * self.delta[self.prop_axis])) as f32;

        let h_field_comp = match self.h_axis {
            0 => &mut h_field.x,
            1 => &mut h_field.y,
            _ => &mut h_field.z,
        };

        // Apply correction on the boundary plane
        self.apply_h_correction(h_field_comp, e_inc, db);
    }

    /// Apply H-field correction at boundary.
    fn apply_h_correction(&self, h_field: &mut crate::arrays::Field3D, e_inc: f32, db: f32) {
        let pos = self.boundary_pos;

        match self.prop_axis {
            0 => {
                // X propagation: correct H on YZ plane at x=pos
                for j in 0..self.dims.ny {
                    for k in 0..self.dims.nz {
                        let h_old = h_field.get(pos, j, k);
                        h_field.set(pos, j, k, h_old - db * e_inc);
                    }
                }
            }
            1 => {
                // Y propagation: correct H on XZ plane
                for i in 0..self.dims.nx {
                    for k in 0..self.dims.nz {
                        let h_old = h_field.get(i, pos, k);
                        h_field.set(i, pos, k, h_old - db * e_inc);
                    }
                }
            }
            _ => {
                // Z propagation: correct H on XY plane
                for i in 0..self.dims.nx {
                    for j in 0..self.dims.ny {
                        let h_old = h_field.get(i, j, pos);
                        h_field.set(i, j, pos, h_old - db * e_inc);
                    }
                }
            }
        }
    }

    /// Pre-update for E-field: Add incident field.
    pub fn pre_update_e(&mut self, e_field: &mut VectorField3D) {
        // Get incident H-field at boundary
        let h_inc = if self.boundary_pos < self.storage.h_1d.len() {
            self.storage.h_1d[self.boundary_pos]
        } else {
            0.0
        };

        // Add correction to E-field
        let cb = (self.dt / (crate::constants::EPS0 * self.delta[self.prop_axis])) as f32;

        let e_field_comp = match self.e_axis {
            0 => &mut e_field.x,
            1 => &mut e_field.y,
            _ => &mut e_field.z,
        };

        self.apply_e_correction(e_field_comp, h_inc, cb);

        self.timestep += 1;
    }

    /// Apply E-field correction at boundary.
    fn apply_e_correction(&self, e_field: &mut crate::arrays::Field3D, h_inc: f32, cb: f32) {
        let pos = self.boundary_pos;

        match self.prop_axis {
            0 => {
                for j in 0..self.dims.ny {
                    for k in 0..self.dims.nz {
                        let e_old = e_field.get(pos, j, k);
                        e_field.set(pos, j, k, e_old + cb * h_inc);
                    }
                }
            }
            1 => {
                for i in 0..self.dims.nx {
                    for k in 0..self.dims.nz {
                        let e_old = e_field.get(i, pos, k);
                        e_field.set(i, pos, k, e_old + cb * h_inc);
                    }
                }
            }
            _ => {
                for i in 0..self.dims.nx {
                    for j in 0..self.dims.ny {
                        let e_old = e_field.get(i, j, pos);
                        e_field.set(i, j, pos, e_old + cb * h_inc);
                    }
                }
            }
        }
    }

    /// Reset the TF/SF boundary.
    pub fn reset(&mut self) {
        self.storage.e_src.fill(0.0);
        self.storage.h_src.fill(0.0);
        self.storage.e_1d.fill(0.0);
        self.storage.h_1d.fill(0.0);
        self.timestep = 0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tfsf_config() {
        let config = TfsfConfig::new(PropagationDirection::ZPlus, 1e9, 1.0);
        assert_eq!(config.direction, PropagationDirection::ZPlus);
        assert_eq!(config.frequency, 1e9);
    }

    #[test]
    fn test_propagation_direction() {
        assert_eq!(PropagationDirection::XPlus.axis(), 0);
        assert_eq!(PropagationDirection::YMinus.axis(), 1);
        assert_eq!(PropagationDirection::ZPlus.axis(), 2);

        assert!(PropagationDirection::XPlus.is_positive());
        assert!(!PropagationDirection::XMinus.is_positive());
    }

    #[test]
    fn test_tfsf_creation() {
        let dims = Dimensions {
            nx: 50,
            ny: 50,
            nz: 50,
        };
        let delta = [0.001, 0.001, 0.001];
        let dt = 1e-12;

        let config = TfsfConfig::new(PropagationDirection::ZPlus, 1e9, 1.0);
        let tfsf = TfsfBoundary::new(dims, delta, dt, config);

        assert!(tfsf.is_active());
    }

    #[test]
    fn test_incident_field() {
        let dims = Dimensions {
            nx: 20,
            ny: 20,
            nz: 20,
        };
        let delta = [0.001, 0.001, 0.001];
        let dt = 1e-12;

        let config = TfsfConfig::new(PropagationDirection::ZPlus, 1e9, 1.0);
        let tfsf = TfsfBoundary::new(dims, delta, dt, config);

        let e = tfsf.incident_e(0.0);
        assert_eq!(e, 0.0); // sin(0) = 0

        let t = 0.25e-9; // Quarter period for 1 GHz
        let e = tfsf.incident_e(t);
        assert!((e - 1.0).abs() < 0.01); // Should be near maximum
    }

    #[test]
    fn test_tfsf_update() {
        let dims = Dimensions {
            nx: 20,
            ny: 20,
            nz: 20,
        };
        let delta = [0.001, 0.001, 0.001];
        let dt = 1e-12;

        let config = TfsfConfig::new(PropagationDirection::ZPlus, 1e9, 1.0).with_boundary_layer(5);
        let mut tfsf = TfsfBoundary::new(dims, delta, dt, config);

        let mut e_field = VectorField3D::new(dims);
        let mut h_field = VectorField3D::new(dims);

        // Run a few updates
        for _ in 0..10 {
            tfsf.pre_update_h(&mut h_field);
            tfsf.pre_update_e(&mut e_field);
        }

        // Fields should have been modified
    }

    #[test]
    fn test_tfsf_reset() {
        let dims = Dimensions {
            nx: 20,
            ny: 20,
            nz: 20,
        };
        let delta = [0.001, 0.001, 0.001];
        let dt = 1e-12;

        let config = TfsfConfig::new(PropagationDirection::XPlus, 1e9, 1.0);
        let mut tfsf = TfsfBoundary::new(dims, delta, dt, config);

        // Modify state
        tfsf.storage.e_1d[5] = 1.0;
        tfsf.timestep = 100;

        tfsf.reset();
        assert_eq!(tfsf.storage.e_1d[5], 0.0);
        assert_eq!(tfsf.timestep, 0);
    }
}
