//! Field probes for extracting data from simulations.

use crate::arrays::VectorField3D;
use crate::fdtd::Engine;
use num_complex::Complex64;

/// Collection of processing operations.
pub struct Processing {
    /// Voltage probes
    pub voltage_probes: Vec<VoltageProbe>,
    /// Field probes
    pub field_probes: Vec<FieldProbe>,
}

impl Processing {
    /// Create a new processing collection.
    pub fn new() -> Self {
        Self {
            voltage_probes: Vec::new(),
            field_probes: Vec::new(),
        }
    }

    /// Add a voltage probe.
    pub fn add_voltage_probe(&mut self, probe: VoltageProbe) {
        self.voltage_probes.push(probe);
    }

    /// Add a field probe.
    pub fn add_field_probe(&mut self, probe: FieldProbe) {
        self.field_probes.push(probe);
    }

    /// Process all probes for current timestep.
    pub fn process(&mut self, engine: &Engine, time: f64) {
        for probe in &mut self.voltage_probes {
            probe.sample(engine.e_field(), time);
        }
        for probe in &mut self.field_probes {
            probe.sample(engine.e_field(), engine.h_field(), time);
        }
    }
}

impl Default for Processing {
    fn default() -> Self {
        Self::new()
    }
}

/// Voltage probe (line integral of E-field).
pub struct VoltageProbe {
    /// Probe name
    pub name: String,
    /// Start position
    pub start: (usize, usize, usize),
    /// End position
    pub end: (usize, usize, usize),
    /// Direction (0=x, 1=y, 2=z)
    pub direction: usize,
    /// Time samples
    pub times: Vec<f64>,
    /// Voltage samples
    pub values: Vec<f64>,
}

impl VoltageProbe {
    /// Create a new voltage probe.
    pub fn new(
        name: &str,
        start: (usize, usize, usize),
        end: (usize, usize, usize),
        direction: usize,
    ) -> Self {
        Self {
            name: name.to_string(),
            start,
            end,
            direction,
            times: Vec::new(),
            values: Vec::new(),
        }
    }

    /// Sample the voltage at current timestep.
    pub fn sample(&mut self, e_field: &VectorField3D, time: f64) {
        let field = e_field.component(self.direction);

        // Simple implementation: sum E-field values along the line
        // TODO: Proper line integral with weighting
        let mut voltage = 0.0f32;
        let (i0, j0, k0) = self.start;
        let (i1, j1, k1) = self.end;

        match self.direction {
            0 => {
                // X-direction
                for i in i0..=i1 {
                    voltage += field.get(i, j0, k0);
                }
            }
            1 => {
                // Y-direction
                for j in j0..=j1 {
                    voltage += field.get(i0, j, k0);
                }
            }
            2 => {
                // Z-direction
                for k in k0..=k1 {
                    voltage += field.get(i0, j0, k);
                }
            }
            _ => {}
        }

        self.times.push(time);
        self.values.push(voltage as f64);
    }

    /// Get the frequency-domain response using FFT.
    pub fn frequency_response(&self, frequencies: &[f64], dt: f64) -> Vec<Complex64> {
        use crate::tools::signal;

        let spectrum = signal::fft(&self.values);
        let n = self.values.len();

        // Interpolate to requested frequencies
        frequencies
            .iter()
            .map(|&f| {
                let bin = (f * n as f64 * dt).round() as usize;
                if bin < spectrum.len() {
                    spectrum[bin]
                } else {
                    Complex64::new(0.0, 0.0)
                }
            })
            .collect()
    }
}

/// Field probe for point measurements.
pub struct FieldProbe {
    /// Probe name
    pub name: String,
    /// Position
    pub position: (usize, usize, usize),
    /// Time samples
    pub times: Vec<f64>,
    /// E-field samples (Ex, Ey, Ez)
    pub e_values: Vec<(f64, f64, f64)>,
    /// H-field samples (Hx, Hy, Hz)
    pub h_values: Vec<(f64, f64, f64)>,
}

impl FieldProbe {
    /// Create a new field probe.
    pub fn new(name: &str, position: (usize, usize, usize)) -> Self {
        Self {
            name: name.to_string(),
            position,
            times: Vec::new(),
            e_values: Vec::new(),
            h_values: Vec::new(),
        }
    }

    /// Sample fields at current timestep.
    pub fn sample(&mut self, e_field: &VectorField3D, h_field: &VectorField3D, time: f64) {
        let (i, j, k) = self.position;

        let ex = e_field.x.get(i, j, k) as f64;
        let ey = e_field.y.get(i, j, k) as f64;
        let ez = e_field.z.get(i, j, k) as f64;

        let hx = h_field.x.get(i, j, k) as f64;
        let hy = h_field.y.get(i, j, k) as f64;
        let hz = h_field.z.get(i, j, k) as f64;

        self.times.push(time);
        self.e_values.push((ex, ey, ez));
        self.h_values.push((hx, hy, hz));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arrays::Dimensions;

    #[test]
    fn test_voltage_probe() {
        let mut probe = VoltageProbe::new("test", (5, 5, 0), (5, 5, 10), 2);

        // Create test field
        let dims = Dimensions::new(10, 10, 20);
        let mut e_field = VectorField3D::new(dims);
        e_field.z.fill(1.0);

        probe.sample(&e_field, 0.0);

        assert_eq!(probe.times.len(), 1);
        assert!(probe.values[0] > 0.0);
    }
}
