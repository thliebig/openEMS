//! Lumped RLC Circuit Elements.
//!
//! Implements localized R, L, C elements within the FDTD grid.
//! These can be used to model discrete circuit components, ports,
//! and terminations.

use crate::arrays::{Dimensions, VectorField3D};
use crate::constants::EPS0;

/// Type of RLC element connection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RlcType {
    /// Parallel RLC (elements in parallel)
    Parallel,
    /// Series RLC (elements in series)
    Series,
}

/// Configuration for a single lumped RLC element.
#[derive(Debug, Clone)]
pub struct RlcElement {
    /// Position in grid
    pub position: [usize; 3],
    /// Direction (0=x, 1=y, 2=z)
    pub direction: usize,
    /// Resistance (Ohms), None = infinite (open)
    pub resistance: Option<f64>,
    /// Inductance (Henries), None = 0 (short)
    pub inductance: Option<f64>,
    /// Capacitance (Farads), None = 0 (open)
    pub capacitance: Option<f64>,
    /// Connection type
    pub rlc_type: RlcType,
}

impl RlcElement {
    /// Create a pure resistor.
    pub fn resistor(position: [usize; 3], direction: usize, resistance: f64) -> Self {
        Self {
            position,
            direction,
            resistance: Some(resistance),
            inductance: None,
            capacitance: None,
            rlc_type: RlcType::Series,
        }
    }

    /// Create a pure inductor.
    pub fn inductor(position: [usize; 3], direction: usize, inductance: f64) -> Self {
        Self {
            position,
            direction,
            resistance: None,
            inductance: Some(inductance),
            capacitance: None,
            rlc_type: RlcType::Series,
        }
    }

    /// Create a pure capacitor.
    pub fn capacitor(position: [usize; 3], direction: usize, capacitance: f64) -> Self {
        Self {
            position,
            direction,
            resistance: None,
            inductance: None,
            capacitance: Some(capacitance),
            rlc_type: RlcType::Series,
        }
    }

    /// Create a parallel RLC circuit.
    pub fn parallel_rlc(
        position: [usize; 3],
        direction: usize,
        r: Option<f64>,
        l: Option<f64>,
        c: Option<f64>,
    ) -> Self {
        Self {
            position,
            direction,
            resistance: r,
            inductance: l,
            capacitance: c,
            rlc_type: RlcType::Parallel,
        }
    }

    /// Create a series RLC circuit.
    pub fn series_rlc(
        position: [usize; 3],
        direction: usize,
        r: Option<f64>,
        l: Option<f64>,
        c: Option<f64>,
    ) -> Self {
        Self {
            position,
            direction,
            resistance: r,
            inductance: l,
            capacitance: c,
            rlc_type: RlcType::Series,
        }
    }
}

/// Internal storage for RLC update coefficients.
#[derive(Debug)]
#[allow(dead_code)]
struct RlcCoeffs {
    /// Position
    position: [usize; 3],
    /// Direction
    direction: usize,
    /// Voltage update coefficient (v_new = v_vv * v_old + ...)
    v_vv: f32,
    /// Voltage from current coefficient
    v_vi: f32,
    /// Voltage from integral coefficient (for inductance)
    v_il: f32,
    /// Current update coefficient
    i_ii: f32,
    /// Current from voltage coefficient
    i_iv: f32,
    /// Accumulated integral for inductance
    i_int: f32,
    /// Previous current for capacitance
    i_prev: f32,
}

/// Lumped RLC element extension.
#[allow(dead_code)]
pub struct LumpedRlc {
    /// Grid dimensions
    dims: Dimensions,
    /// Timestep
    dt: f64,
    /// Cell sizes
    delta: [f64; 3],
    /// RLC element coefficients
    elements: Vec<RlcCoeffs>,
}

impl LumpedRlc {
    /// Create a new lumped RLC extension.
    pub fn new(dims: Dimensions, delta: [f64; 3], dt: f64, elements: Vec<RlcElement>) -> Self {
        let mut rlc = Self {
            dims,
            dt,
            delta,
            elements: Vec::new(),
        };

        for elem in elements {
            rlc.add_element(elem);
        }

        rlc
    }

    /// Add an RLC element.
    fn add_element(&mut self, elem: RlcElement) {
        let dt = self.dt;
        let dl = self.delta[elem.direction];
        let da = dl * dl; // Edge area (simplified as square cell)

        // Calculate update coefficients based on RLC type
        let (v_vv, v_vi, v_il, i_ii, i_iv) = match elem.rlc_type {
            RlcType::Parallel => {
                self.calc_parallel_coeffs(&elem, dt, dl, da)
            }
            RlcType::Series => {
                self.calc_series_coeffs(&elem, dt, dl, da)
            }
        };

        self.elements.push(RlcCoeffs {
            position: elem.position,
            direction: elem.direction,
            v_vv,
            v_vi,
            v_il,
            i_ii,
            i_iv,
            i_int: 0.0,
            i_prev: 0.0,
        });
    }

    /// Calculate coefficients for parallel RLC.
    fn calc_parallel_coeffs(
        &self,
        elem: &RlcElement,
        dt: f64,
        dl: f64,
        da: f64,
    ) -> (f32, f32, f32, f32, f32) {
        // Parallel RLC: admittances add
        // Y = 1/R + 1/(j*omega*L) + j*omega*C
        // Discretized: I = (1/R)*V + (dt/L)*integral(V) + C*dV/dt

        let g = elem.resistance.map(|r| 1.0 / r).unwrap_or(0.0); // Conductance
        let inv_l = elem.inductance.map(|l| dt / l).unwrap_or(0.0);
        let c = elem.capacitance.unwrap_or(0.0);

        // Effective parameters normalized by cell geometry
        let g_eff = g * dl / da;
        let c_eff = c * dl / da + EPS0 * dl / da;

        // Update coefficients
        let denom = 2.0 * c_eff + g_eff * dt;
        let v_vv = ((2.0 * c_eff - g_eff * dt) / denom) as f32;
        let v_vi = (2.0 * dt / denom) as f32;
        let v_il = (inv_l * dl / da) as f32;

        (v_vv, v_vi, v_il, 1.0, 0.0)
    }

    /// Calculate coefficients for series RLC.
    fn calc_series_coeffs(
        &self,
        elem: &RlcElement,
        dt: f64,
        dl: f64,
        da: f64,
    ) -> (f32, f32, f32, f32, f32) {
        // Series RLC: impedances add
        // Z = R + j*omega*L + 1/(j*omega*C)
        // Discretized: V = R*I + L*dI/dt + (1/C)*integral(I)

        let r = elem.resistance.unwrap_or(0.0);
        let l = elem.inductance.unwrap_or(0.0);
        let inv_c = elem.capacitance.map(|c| dt / c).unwrap_or(0.0);

        // Effective parameters
        let r_eff = r * da / dl;
        let l_eff = l * da / dl;

        // Current update coefficients
        let denom = 2.0 * l_eff + r_eff * dt;
        if denom.abs() < 1e-20 {
            // No inductance or resistance, pure capacitor
            return (1.0, 0.0, 0.0, 1.0, (inv_c * da / dl) as f32);
        }

        let i_ii = ((2.0 * l_eff - r_eff * dt) / denom) as f32;
        let i_iv = (2.0 * dt / denom) as f32;
        let v_il = (inv_c * da / dl) as f32;

        (1.0, 0.0, v_il, i_ii, i_iv)
    }

    /// Check if any RLC elements are defined.
    pub fn is_active(&self) -> bool {
        !self.elements.is_empty()
    }

    /// Get number of RLC elements.
    pub fn num_elements(&self) -> usize {
        self.elements.len()
    }

    /// Apply RLC updates after E-field update.
    pub fn post_update_e(&mut self, e_field: &mut VectorField3D, _h_field: &VectorField3D) {
        for elem in &mut self.elements {
            let pos = elem.position;
            let dir = elem.direction;

            let e = match dir {
                0 => e_field.x.get(pos[0], pos[1], pos[2]),
                1 => e_field.y.get(pos[0], pos[1], pos[2]),
                _ => e_field.z.get(pos[0], pos[1], pos[2]),
            };

            // Apply RLC modification
            let e_new = elem.v_vv * e + elem.v_il * elem.i_int;

            // Update integral for inductance
            elem.i_int += e * self.dt as f32;

            match dir {
                0 => e_field.x.set(pos[0], pos[1], pos[2], e_new),
                1 => e_field.y.set(pos[0], pos[1], pos[2], e_new),
                _ => e_field.z.set(pos[0], pos[1], pos[2], e_new),
            }
        }
    }

    /// Apply RLC updates after H-field update.
    pub fn post_update_h(&mut self, _h_field: &mut VectorField3D, e_field: &VectorField3D) {
        for elem in &mut self.elements {
            let pos = elem.position;
            let dir = elem.direction;

            // For series RLC, update current through voltage
            let e = match dir {
                0 => e_field.x.get(pos[0], pos[1], pos[2]),
                1 => e_field.y.get(pos[0], pos[1], pos[2]),
                _ => e_field.z.get(pos[0], pos[1], pos[2]),
            };

            // Store for capacitor update
            elem.i_prev = e;
        }
    }

    /// Reset all internal state.
    pub fn reset(&mut self) {
        for elem in &mut self.elements {
            elem.i_int = 0.0;
            elem.i_prev = 0.0;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rlc_element_creation() {
        let r = RlcElement::resistor([5, 5, 5], 2, 50.0);
        assert_eq!(r.resistance, Some(50.0));
        assert_eq!(r.inductance, None);

        let l = RlcElement::inductor([5, 5, 5], 2, 1e-9);
        assert_eq!(l.inductance, Some(1e-9));

        let c = RlcElement::capacitor([5, 5, 5], 2, 1e-12);
        assert_eq!(c.capacitance, Some(1e-12));
    }

    #[test]
    fn test_lumped_rlc_creation() {
        let dims = Dimensions { nx: 20, ny: 20, nz: 20 };
        let delta = [0.001, 0.001, 0.001];
        let dt = 1e-12;

        let elements = vec![
            RlcElement::resistor([10, 10, 10], 2, 50.0),
            RlcElement::capacitor([10, 10, 11], 2, 1e-12),
        ];

        let rlc = LumpedRlc::new(dims, delta, dt, elements);

        assert!(rlc.is_active());
        assert_eq!(rlc.num_elements(), 2);
    }

    #[test]
    fn test_rlc_update() {
        let dims = Dimensions { nx: 10, ny: 10, nz: 10 };
        let delta = [0.001, 0.001, 0.001];
        let dt = 1e-12;

        let elements = vec![RlcElement::resistor([5, 5, 5], 0, 50.0)];
        let mut rlc = LumpedRlc::new(dims, delta, dt, elements);

        let mut e_field = VectorField3D::new(dims);
        let h_field = VectorField3D::new(dims);

        e_field.x.set(5, 5, 5, 1.0);

        rlc.post_update_e(&mut e_field, &h_field);

        // Resistor should modify the field
        let _e = e_field.x.get(5, 5, 5);
        // Value depends on coefficients
    }

    #[test]
    fn test_parallel_rlc() {
        let elem = RlcElement::parallel_rlc(
            [5, 5, 5],
            2,
            Some(50.0),  // R
            Some(1e-9),  // L
            Some(1e-12), // C
        );

        assert_eq!(elem.rlc_type, RlcType::Parallel);
    }

    #[test]
    fn test_series_rlc() {
        let elem = RlcElement::series_rlc(
            [5, 5, 5],
            2,
            Some(50.0),
            Some(1e-9),
            Some(1e-12),
        );

        assert_eq!(elem.rlc_type, RlcType::Series);
    }

    #[test]
    fn test_rlc_reset() {
        let dims = Dimensions { nx: 10, ny: 10, nz: 10 };
        let delta = [0.001, 0.001, 0.001];
        let dt = 1e-12;

        let elements = vec![RlcElement::inductor([5, 5, 5], 0, 1e-9)];
        let mut rlc = LumpedRlc::new(dims, delta, dt, elements);

        // Modify state
        rlc.elements[0].i_int = 1.0;

        rlc.reset();
        assert_eq!(rlc.elements[0].i_int, 0.0);
    }
}
