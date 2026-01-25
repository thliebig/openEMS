//! Field probes for extracting data from simulations.

use crate::arrays::VectorField3D;
use crate::fdtd::Engine;
use num_complex::Complex64;

/// Collection of processing operations.
pub struct Processing {
    /// Voltage probes
    pub voltage_probes: Vec<VoltageProbe>,
    /// Current probes
    pub current_probes: Vec<CurrentProbe>,
    /// Field probes
    pub field_probes: Vec<FieldProbe>,
}

impl Processing {
    /// Create a new processing collection.
    pub fn new() -> Self {
        Self {
            voltage_probes: Vec::new(),
            current_probes: Vec::new(),
            field_probes: Vec::new(),
        }
    }

    /// Add a voltage probe.
    pub fn add_voltage_probe(&mut self, probe: VoltageProbe) {
        self.voltage_probes.push(probe);
    }

    /// Add a current probe.
    pub fn add_current_probe(&mut self, probe: CurrentProbe) {
        self.current_probes.push(probe);
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
        for probe in &mut self.current_probes {
            probe.sample(engine.h_field(), time);
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

/// Current probe (area integral of H-field / surface current).
pub struct CurrentProbe {
    /// Probe name
    pub name: String,
    /// Start position of measurement surface
    pub start: (usize, usize, usize),
    /// End position of measurement surface
    pub end: (usize, usize, usize),
    /// Normal direction (0=x, 1=y, 2=z)
    pub normal_direction: usize,
    /// Time samples
    pub times: Vec<f64>,
    /// Current samples
    pub values: Vec<f64>,
}

impl CurrentProbe {
    /// Create a new current probe.
    pub fn new(
        name: &str,
        start: (usize, usize, usize),
        end: (usize, usize, usize),
        normal_direction: usize,
    ) -> Self {
        Self {
            name: name.to_string(),
            start,
            end,
            normal_direction,
            times: Vec::new(),
            values: Vec::new(),
        }
    }

    /// Sample the current at current timestep.
    ///
    /// Current is computed as the line integral of H around the measurement surface.
    /// For a surface with normal in the x-direction, we integrate Hy and Hz.
    pub fn sample(&mut self, h_field: &VectorField3D, time: f64) {
        let (i0, j0, k0) = self.start;
        let (i1, j1, k1) = self.end;

        let mut current = 0.0f32;

        // The current through a surface is the integral of (curl H) · dA
        // By Stokes' theorem, this equals the line integral of H around the boundary
        // For simplicity, we sum the H-field components tangent to the surface

        match self.normal_direction {
            0 => {
                // X-normal surface: integrate Hy and Hz around boundary
                // This is a simplified calculation - proper implementation needs
                // edge contributions
                for j in j0..=j1 {
                    for k in k0..=k1 {
                        // Use H at the surface position
                        current += h_field.y.get(i0, j, k);
                        current -= h_field.z.get(i0, j, k);
                    }
                }
            }
            1 => {
                // Y-normal surface: integrate Hx and Hz
                for i in i0..=i1 {
                    for k in k0..=k1 {
                        current += h_field.z.get(i, j0, k);
                        current -= h_field.x.get(i, j0, k);
                    }
                }
            }
            2 => {
                // Z-normal surface: integrate Hx and Hy
                for i in i0..=i1 {
                    for j in j0..=j1 {
                        current += h_field.x.get(i, j, k0);
                        current -= h_field.y.get(i, j, k0);
                    }
                }
            }
            _ => {}
        }

        self.times.push(time);
        self.values.push(current as f64);
    }

    /// Get the frequency-domain response using FFT.
    pub fn frequency_response(&self, frequencies: &[f64], dt: f64) -> Vec<Complex64> {
        use crate::tools::signal;

        let spectrum = signal::fft(&self.values);
        let n = self.values.len();

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

    #[test]
    fn test_current_probe() {
        let mut probe = CurrentProbe::new("test_current", (5, 0, 0), (5, 10, 10), 0);

        // Create test field
        let dims = Dimensions::new(10, 20, 20);
        let mut h_field = VectorField3D::new(dims);
        h_field.y.fill(1.0);

        probe.sample(&h_field, 0.0);

        assert_eq!(probe.times.len(), 1);
        // Current should be non-zero due to Hy
        assert!(probe.values[0].abs() > 0.0);
    }

    #[test]
    fn test_field_probe() {
        let mut probe = FieldProbe::new("test_field", (5, 5, 5));

        let dims = Dimensions::new(10, 10, 10);
        let mut e_field = VectorField3D::new(dims);
        let mut h_field = VectorField3D::new(dims);

        e_field.x.set(5, 5, 5, 1.0);
        e_field.y.set(5, 5, 5, 2.0);
        e_field.z.set(5, 5, 5, 3.0);
        h_field.x.set(5, 5, 5, 0.1);

        probe.sample(&e_field, &h_field, 0.0);

        assert_eq!(probe.times.len(), 1);
        // Use approximate comparison due to f32/f64 conversion
        let (ex, ey, ez) = probe.e_values[0];
        assert!((ex - 1.0).abs() < 1e-6);
        assert!((ey - 2.0).abs() < 1e-6);
        assert!((ez - 3.0).abs() < 1e-6);
        assert!((probe.h_values[0].0 - 0.1).abs() < 1e-6);
    }

    #[test]
    fn test_processing_collection() {
        let mut processing = Processing::new();

        processing.add_voltage_probe(VoltageProbe::new("v1", (0, 0, 0), (0, 0, 5), 2));
        processing.add_current_probe(CurrentProbe::new("i1", (5, 0, 0), (5, 5, 5), 0));
        processing.add_field_probe(FieldProbe::new("f1", (2, 2, 2)));

        assert_eq!(processing.voltage_probes.len(), 1);
        assert_eq!(processing.current_probes.len(), 1);
        assert_eq!(processing.field_probes.len(), 1);
    }
}
