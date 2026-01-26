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

    #[test]
    fn test_processing_default() {
        let processing = Processing::default();
        assert!(processing.voltage_probes.is_empty());
        assert!(processing.current_probes.is_empty());
        assert!(processing.field_probes.is_empty());
    }

    #[test]
    fn test_voltage_probe_x_direction() {
        let mut probe = VoltageProbe::new("test_x", (0, 5, 5), (9, 5, 5), 0);

        let dims = Dimensions::new(10, 10, 10);
        let mut e_field = VectorField3D::new(dims);
        e_field.x.fill(2.0);

        probe.sample(&e_field, 1.0e-9);

        assert_eq!(probe.times.len(), 1);
        assert!((probe.times[0] - 1.0e-9).abs() < 1e-15);
        // Should sum 10 cells with value 2.0
        assert!((probe.values[0] - 20.0).abs() < 1e-6);
    }

    #[test]
    fn test_voltage_probe_y_direction() {
        let mut probe = VoltageProbe::new("test_y", (5, 0, 5), (5, 9, 5), 1);

        let dims = Dimensions::new(10, 10, 10);
        let mut e_field = VectorField3D::new(dims);
        e_field.y.fill(3.0);

        probe.sample(&e_field, 2.0e-9);

        assert_eq!(probe.times.len(), 1);
        assert!((probe.times[0] - 2.0e-9).abs() < 1e-15);
        // Should sum 10 cells with value 3.0
        assert!((probe.values[0] - 30.0).abs() < 1e-6);
    }

    #[test]
    fn test_voltage_probe_multiple_samples() {
        let mut probe = VoltageProbe::new("test_multi", (5, 5, 0), (5, 5, 4), 2);

        let dims = Dimensions::new(10, 10, 10);
        let mut e_field = VectorField3D::new(dims);

        // Sample multiple times with different field values
        e_field.z.fill(1.0);
        probe.sample(&e_field, 0.0);

        e_field.z.fill(2.0);
        probe.sample(&e_field, 1e-9);

        e_field.z.fill(0.5);
        probe.sample(&e_field, 2e-9);

        assert_eq!(probe.times.len(), 3);
        assert_eq!(probe.values.len(), 3);
        assert!((probe.values[0] - 5.0).abs() < 1e-6); // 5 cells * 1.0
        assert!((probe.values[1] - 10.0).abs() < 1e-6); // 5 cells * 2.0
        assert!((probe.values[2] - 2.5).abs() < 1e-6); // 5 cells * 0.5
    }

    #[test]
    fn test_voltage_probe_z_direction_single_point() {
        // Test a single point probe (start == end)
        let mut probe = VoltageProbe::new("test_single", (5, 5, 5), (5, 5, 5), 2);

        let dims = Dimensions::new(10, 10, 10);
        let mut e_field = VectorField3D::new(dims);
        e_field.z.set(5, 5, 5, 2.5);

        probe.sample(&e_field, 0.0);

        assert_eq!(probe.times.len(), 1);
        // Single point should return that point's value
        assert!((probe.values[0] - 2.5).abs() < 1e-5);
    }

    #[test]
    fn test_current_probe_y_normal() {
        let mut probe = CurrentProbe::new("test_y_normal", (0, 5, 0), (9, 5, 9), 1);

        let dims = Dimensions::new(10, 10, 10);
        let mut h_field = VectorField3D::new(dims);
        h_field.z.fill(1.0);
        h_field.x.fill(0.5);

        probe.sample(&h_field, 0.0);

        assert_eq!(probe.times.len(), 1);
        // Current = sum(Hz - Hx) = 10*10*(1.0 - 0.5) = 50
        assert!((probe.values[0] - 50.0).abs() < 1e-5);
    }

    #[test]
    fn test_current_probe_z_normal() {
        let mut probe = CurrentProbe::new("test_z_normal", (0, 0, 5), (9, 9, 5), 2);

        let dims = Dimensions::new(10, 10, 10);
        let mut h_field = VectorField3D::new(dims);
        h_field.x.fill(2.0);
        h_field.y.fill(1.0);

        probe.sample(&h_field, 0.0);

        assert_eq!(probe.times.len(), 1);
        // Current = sum(Hx - Hy) = 10*10*(2.0 - 1.0) = 100
        assert!((probe.values[0] - 100.0).abs() < 1e-5);
    }

    #[test]
    fn test_current_probe_invalid_direction() {
        let mut probe = CurrentProbe::new("test_invalid", (5, 0, 0), (5, 4, 4), 10); // Invalid

        let dims = Dimensions::new(10, 10, 10);
        let mut h_field = VectorField3D::new(dims);
        h_field.y.fill(1.0);

        probe.sample(&h_field, 0.0);

        // Should sample but current should be 0
        assert_eq!(probe.times.len(), 1);
        assert!((probe.values[0]).abs() < 1e-10);
    }

    #[test]
    fn test_field_probe_all_components() {
        let mut probe = FieldProbe::new("test_all", (3, 4, 5));

        let dims = Dimensions::new(10, 10, 10);
        let mut e_field = VectorField3D::new(dims);
        let mut h_field = VectorField3D::new(dims);

        // Set specific values at probe location
        e_field.x.set(3, 4, 5, 1.0);
        e_field.y.set(3, 4, 5, 2.0);
        e_field.z.set(3, 4, 5, 3.0);
        h_field.x.set(3, 4, 5, 0.1);
        h_field.y.set(3, 4, 5, 0.2);
        h_field.z.set(3, 4, 5, 0.3);

        probe.sample(&e_field, &h_field, 1e-9);

        let (ex, ey, ez) = probe.e_values[0];
        let (hx, hy, hz) = probe.h_values[0];

        assert!((ex - 1.0).abs() < 1e-6);
        assert!((ey - 2.0).abs() < 1e-6);
        assert!((ez - 3.0).abs() < 1e-6);
        assert!((hx - 0.1).abs() < 1e-6);
        assert!((hy - 0.2).abs() < 1e-6);
        assert!((hz - 0.3).abs() < 1e-6);
    }

    #[test]
    fn test_field_probe_multiple_samples() {
        let mut probe = FieldProbe::new("test_multi", (2, 2, 2));

        let dims = Dimensions::new(5, 5, 5);
        let mut e_field = VectorField3D::new(dims);
        let mut h_field = VectorField3D::new(dims);

        // First sample
        e_field.x.set(2, 2, 2, 1.0);
        h_field.x.set(2, 2, 2, 0.1);
        probe.sample(&e_field, &h_field, 0.0);

        // Second sample
        e_field.x.set(2, 2, 2, 2.0);
        h_field.x.set(2, 2, 2, 0.2);
        probe.sample(&e_field, &h_field, 1e-9);

        assert_eq!(probe.times.len(), 2);
        assert_eq!(probe.e_values.len(), 2);
        assert_eq!(probe.h_values.len(), 2);

        assert!((probe.e_values[0].0 - 1.0).abs() < 1e-6);
        assert!((probe.e_values[1].0 - 2.0).abs() < 1e-6);
        assert!((probe.h_values[0].0 - 0.1).abs() < 1e-6);
        assert!((probe.h_values[1].0 - 0.2).abs() < 1e-6);
    }

    #[test]
    fn test_voltage_probe_frequency_response() {
        let mut probe = VoltageProbe::new("test_freq", (5, 5, 0), (5, 5, 4), 2);

        // Generate a constant signal (should have strong DC component)
        let n_samples = 16;
        let dt = 1e-9;

        for i in 0..n_samples {
            let t = i as f64 * dt;
            probe.times.push(t);
            probe.values.push(1.0); // Constant value
        }

        // Get frequency response at a few frequencies
        let test_frequencies = vec![0.0, 1e8, 2e8];
        let response = probe.frequency_response(&test_frequencies, dt);

        assert_eq!(response.len(), 3);
        // DC component should be strong for constant signal
        assert!(response[0].norm() > 0.5);
    }

    #[test]
    fn test_current_probe_frequency_response() {
        let mut probe = CurrentProbe::new("test_freq", (5, 0, 0), (5, 4, 4), 0);

        // Generate a constant signal (should have DC component)
        let n_samples = 8;
        let dt = 1e-9;

        for i in 0..n_samples {
            let t = i as f64 * dt;
            probe.times.push(t);
            probe.values.push(1.0); // Constant value
        }

        let test_frequencies = vec![0.0, 1e8];
        let response = probe.frequency_response(&test_frequencies, dt);

        assert_eq!(response.len(), 2);
        // DC component should be strong for constant signal
        assert!(response[0].norm() > 0.5);
    }

    #[test]
    fn test_voltage_probe_frequency_response_out_of_range() {
        let mut probe = VoltageProbe::new("test_oor", (5, 5, 0), (5, 5, 4), 2);

        // Add some samples
        for i in 0..4 {
            probe.times.push(i as f64 * 1e-9);
            probe.values.push(1.0);
        }

        // Request frequency that would be out of spectrum range
        let dt = 1e-9;
        let very_high_freq = vec![1e12]; // Much higher than Nyquist
        let response = probe.frequency_response(&very_high_freq, dt);

        // Should return zero for out-of-range frequency
        assert_eq!(response.len(), 1);
        assert!((response[0].re).abs() < 1e-10);
        assert!((response[0].im).abs() < 1e-10);
    }

    #[test]
    fn test_probe_name_storage() {
        let v_probe = VoltageProbe::new("voltage_1", (0, 0, 0), (0, 0, 5), 2);
        let c_probe = CurrentProbe::new("current_1", (5, 0, 0), (5, 5, 5), 0);
        let f_probe = FieldProbe::new("field_1", (2, 2, 2));

        assert_eq!(v_probe.name, "voltage_1");
        assert_eq!(c_probe.name, "current_1");
        assert_eq!(f_probe.name, "field_1");
    }

    #[test]
    fn test_probe_positions() {
        let v_probe = VoltageProbe::new("v", (1, 2, 3), (4, 5, 6), 0);
        assert_eq!(v_probe.start, (1, 2, 3));
        assert_eq!(v_probe.end, (4, 5, 6));
        assert_eq!(v_probe.direction, 0);

        let c_probe = CurrentProbe::new("c", (10, 20, 30), (40, 50, 60), 2);
        assert_eq!(c_probe.start, (10, 20, 30));
        assert_eq!(c_probe.end, (40, 50, 60));
        assert_eq!(c_probe.normal_direction, 2);

        let f_probe = FieldProbe::new("f", (7, 8, 9));
        assert_eq!(f_probe.position, (7, 8, 9));
    }
}
