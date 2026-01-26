//! Frequency-domain field processing.
//!
//! Accumulates field data in the frequency domain using DFT.
//! This allows extraction of field phasors at specific frequencies.

use crate::arrays::{Dimensions, Field3D, VectorField3D};
use num_complex::Complex64;
use std::f64::consts::PI;

/// Accumulate a single field component into a complex accumulator.
fn accumulate_field_component(
    field: &Field3D,
    acc: &mut [Complex64],
    factor: Complex64,
    dims: Dimensions,
    sub_sample: [usize; 3],
) {
    for i in (0..dims.nx).step_by(sub_sample[0]) {
        for j in (0..dims.ny).step_by(sub_sample[1]) {
            for k in (0..dims.nz).step_by(sub_sample[2]) {
                let idx = dims.to_linear(i, j, k);
                let value = field.get(i, j, k) as f64;
                acc[idx] += Complex64::new(value, 0.0) * factor;
            }
        }
    }
}

/// Frequency domain field accumulator.
///
/// Accumulates time-domain field samples into frequency-domain phasors
/// using the DFT formula: X(f) = sum(x(t) * exp(-j*2*pi*f*t) * dt)
#[derive(Debug)]
pub struct FrequencyDomainFieldDump {
    /// Field dimensions
    dims: Dimensions,
    /// Frequencies to sample (Hz)
    frequencies: Vec<f64>,
    /// Accumulated E-field (complex, per frequency)
    e_field_fd: Vec<[Vec<Complex64>; 3]>,
    /// Accumulated H-field (complex, per frequency)
    h_field_fd: Vec<[Vec<Complex64>; 3]>,
    /// Number of samples accumulated
    num_samples: u64,
    /// Dump type
    dump_type: FdDumpType,
    /// Sub-sampling rate
    sub_sample: [usize; 3],
    /// Whether to dump E-field
    dump_e: bool,
    /// Whether to dump H-field
    dump_h: bool,
}

/// Type of frequency domain dump.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FdDumpType {
    /// E-field only
    EField,
    /// H-field only
    HField,
    /// Both E and H fields
    Both,
    /// Current density (J-field)
    JField,
    /// Total current density (rot H)
    RotH,
}

impl FrequencyDomainFieldDump {
    /// Create a new frequency domain field dump.
    pub fn new(dims: Dimensions, frequencies: Vec<f64>, dump_type: FdDumpType) -> Self {
        let size = dims.total();
        let num_freqs = frequencies.len();

        let (dump_e, dump_h) = match dump_type {
            FdDumpType::EField | FdDumpType::JField => (true, false),
            FdDumpType::HField => (false, true),
            FdDumpType::Both | FdDumpType::RotH => (true, true),
        };

        // Initialize complex field storage
        let e_field_fd = if dump_e {
            (0..num_freqs)
                .map(|_| {
                    [
                        vec![Complex64::new(0.0, 0.0); size],
                        vec![Complex64::new(0.0, 0.0); size],
                        vec![Complex64::new(0.0, 0.0); size],
                    ]
                })
                .collect()
        } else {
            Vec::new()
        };

        let h_field_fd = if dump_h {
            (0..num_freqs)
                .map(|_| {
                    [
                        vec![Complex64::new(0.0, 0.0); size],
                        vec![Complex64::new(0.0, 0.0); size],
                        vec![Complex64::new(0.0, 0.0); size],
                    ]
                })
                .collect()
        } else {
            Vec::new()
        };

        Self {
            dims,
            frequencies,
            e_field_fd,
            h_field_fd,
            num_samples: 0,
            dump_type,
            sub_sample: [1, 1, 1],
            dump_e,
            dump_h,
        }
    }

    /// Set sub-sampling rate.
    pub fn set_sub_sample(&mut self, rate: [usize; 3]) {
        self.sub_sample = [rate[0].max(1), rate[1].max(1), rate[2].max(1)];
    }

    /// Get frequencies.
    pub fn frequencies(&self) -> &[f64] {
        &self.frequencies
    }

    /// Get number of frequencies.
    pub fn num_frequencies(&self) -> usize {
        self.frequencies.len()
    }

    /// Get number of accumulated samples.
    pub fn num_samples(&self) -> u64 {
        self.num_samples
    }

    /// Accumulate a time-domain sample.
    ///
    /// This performs one step of the DFT accumulation:
    /// X(f) += x(t) * exp(-j*2*pi*f*t) * 2 * dt
    pub fn accumulate(
        &mut self,
        e_field: &VectorField3D,
        h_field: &VectorField3D,
        time: f64,
        dt: f64,
    ) {
        let scale = 2.0 * dt;
        let dims = self.dims;
        let sub_sample = self.sub_sample;

        for (freq_idx, &freq) in self.frequencies.iter().enumerate() {
            // exp(-j*2*pi*f*t)
            let phase = -2.0 * PI * freq * time;
            let exp_factor = Complex64::new(phase.cos(), phase.sin()) * scale;

            // Accumulate E-field
            if self.dump_e {
                accumulate_field_component(&e_field.x, &mut self.e_field_fd[freq_idx][0], exp_factor, dims, sub_sample);
                accumulate_field_component(&e_field.y, &mut self.e_field_fd[freq_idx][1], exp_factor, dims, sub_sample);
                accumulate_field_component(&e_field.z, &mut self.e_field_fd[freq_idx][2], exp_factor, dims, sub_sample);
            }

            // Accumulate H-field
            if self.dump_h {
                accumulate_field_component(&h_field.x, &mut self.h_field_fd[freq_idx][0], exp_factor, dims, sub_sample);
                accumulate_field_component(&h_field.y, &mut self.h_field_fd[freq_idx][1], exp_factor, dims, sub_sample);
                accumulate_field_component(&h_field.z, &mut self.h_field_fd[freq_idx][2], exp_factor, dims, sub_sample);
            }
        }

        self.num_samples += 1;
    }

    /// Get accumulated E-field for a specific frequency.
    pub fn get_e_field(&self, freq_idx: usize) -> Option<&[Vec<Complex64>; 3]> {
        if self.dump_e && freq_idx < self.e_field_fd.len() {
            Some(&self.e_field_fd[freq_idx])
        } else {
            None
        }
    }

    /// Get accumulated H-field for a specific frequency.
    pub fn get_h_field(&self, freq_idx: usize) -> Option<&[Vec<Complex64>; 3]> {
        if self.dump_h && freq_idx < self.h_field_fd.len() {
            Some(&self.h_field_fd[freq_idx])
        } else {
            None
        }
    }

    /// Get E-field magnitude at a specific position and frequency.
    pub fn get_e_magnitude(&self, freq_idx: usize, i: usize, j: usize, k: usize) -> f64 {
        if let Some(e) = self.get_e_field(freq_idx) {
            let idx = self.dims.to_linear(i, j, k);
            (e[0][idx].norm_sqr() + e[1][idx].norm_sqr() + e[2][idx].norm_sqr()).sqrt()
        } else {
            0.0
        }
    }

    /// Get H-field magnitude at a specific position and frequency.
    pub fn get_h_magnitude(&self, freq_idx: usize, i: usize, j: usize, k: usize) -> f64 {
        if let Some(h) = self.get_h_field(freq_idx) {
            let idx = self.dims.to_linear(i, j, k);
            (h[0][idx].norm_sqr() + h[1][idx].norm_sqr() + h[2][idx].norm_sqr()).sqrt()
        } else {
            0.0
        }
    }

    /// Extract E-field magnitude as a 3D field for a specific frequency.
    pub fn extract_e_magnitude(&self, freq_idx: usize) -> Option<Field3D> {
        if let Some(e) = self.get_e_field(freq_idx) {
            let mut mag = Field3D::new(self.dims);

            for i in 0..self.dims.nx {
                for j in 0..self.dims.ny {
                    for k in 0..self.dims.nz {
                        let idx = self.dims.to_linear(i, j, k);
                        let m = (e[0][idx].norm_sqr()
                            + e[1][idx].norm_sqr()
                            + e[2][idx].norm_sqr())
                        .sqrt();
                        mag.set(i, j, k, m as f32);
                    }
                }
            }

            Some(mag)
        } else {
            None
        }
    }

    /// Extract E-field phase (of Ex component) as a 3D field for a specific frequency.
    pub fn extract_e_phase(&self, freq_idx: usize, component: usize) -> Option<Field3D> {
        if let Some(e) = self.get_e_field(freq_idx) {
            let comp = component.min(2);
            let mut phase = Field3D::new(self.dims);

            for i in 0..self.dims.nx {
                for j in 0..self.dims.ny {
                    for k in 0..self.dims.nz {
                        let idx = self.dims.to_linear(i, j, k);
                        let p = e[comp][idx].arg();
                        phase.set(i, j, k, p as f32);
                    }
                }
            }

            Some(phase)
        } else {
            None
        }
    }

    /// Reset accumulation.
    pub fn reset(&mut self) {
        for e in &mut self.e_field_fd {
            for comp in e.iter_mut() {
                for c in comp.iter_mut() {
                    *c = Complex64::new(0.0, 0.0);
                }
            }
        }

        for h in &mut self.h_field_fd {
            for comp in h.iter_mut() {
                for c in comp.iter_mut() {
                    *c = Complex64::new(0.0, 0.0);
                }
            }
        }

        self.num_samples = 0;
    }
}

/// Time-domain field dump configuration.
#[derive(Debug, Clone)]
pub struct TimeDomainFieldDump {
    /// Dump interval (timesteps)
    pub interval: u64,
    /// Sub-sampling rate
    pub sub_sample: [usize; 3],
    /// Output path prefix
    pub path_prefix: String,
    /// Padding for timestep numbers
    pub padding: usize,
    /// Current file index
    file_index: u64,
}

impl TimeDomainFieldDump {
    /// Create a new time-domain field dump.
    pub fn new(interval: u64, path_prefix: &str) -> Self {
        Self {
            interval,
            sub_sample: [1, 1, 1],
            path_prefix: path_prefix.to_string(),
            padding: 8,
            file_index: 0,
        }
    }

    /// Set sub-sampling rate.
    pub fn set_sub_sample(&mut self, rate: [usize; 3]) {
        self.sub_sample = rate;
    }

    /// Check if dump is due at this timestep.
    pub fn should_dump(&self, timestep: u64) -> bool {
        self.interval > 0 && timestep % self.interval == 0
    }

    /// Get output filename for current dump.
    pub fn get_filename(&self) -> String {
        format!(
            "{}_{:0width$}",
            self.path_prefix,
            self.file_index,
            width = self.padding
        )
    }

    /// Increment file index.
    pub fn increment(&mut self) {
        self.file_index += 1;
    }

    /// Get current file index.
    pub fn file_index(&self) -> u64 {
        self.file_index
    }

    /// Reset file index.
    pub fn reset(&mut self) {
        self.file_index = 0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fd_dump_creation() {
        let dims = Dimensions::new(10, 10, 10);
        let freqs = vec![1e9, 2e9, 3e9];

        let dump = FrequencyDomainFieldDump::new(dims, freqs.clone(), FdDumpType::Both);

        assert_eq!(dump.num_frequencies(), 3);
        assert_eq!(dump.num_samples(), 0);
    }

    #[test]
    fn test_fd_accumulation() {
        let dims = Dimensions::new(5, 5, 5);
        let freqs = vec![1e9];

        let mut dump = FrequencyDomainFieldDump::new(dims, freqs, FdDumpType::EField);

        let mut e_field = VectorField3D::new(dims);
        let h_field = VectorField3D::new(dims);

        // Set uniform E-field
        e_field.x.fill(1.0);

        // Accumulate several samples
        let dt = 1e-12;
        for i in 0..10 {
            dump.accumulate(&e_field, &h_field, i as f64 * dt, dt);
        }

        assert_eq!(dump.num_samples(), 10);

        // Check that something was accumulated
        let e_mag = dump.get_e_magnitude(0, 2, 2, 2);
        assert!(e_mag > 0.0);
    }

    #[test]
    fn test_fd_reset() {
        let dims = Dimensions::new(5, 5, 5);
        let freqs = vec![1e9];

        let mut dump = FrequencyDomainFieldDump::new(dims, freqs, FdDumpType::EField);

        let mut e_field = VectorField3D::new(dims);
        let h_field = VectorField3D::new(dims);
        e_field.x.fill(1.0);

        dump.accumulate(&e_field, &h_field, 0.0, 1e-12);
        assert_eq!(dump.num_samples(), 1);

        dump.reset();
        assert_eq!(dump.num_samples(), 0);
    }

    #[test]
    fn test_extract_magnitude() {
        let dims = Dimensions::new(5, 5, 5);
        let freqs = vec![1e9];

        let mut dump = FrequencyDomainFieldDump::new(dims, freqs, FdDumpType::EField);

        let mut e_field = VectorField3D::new(dims);
        let h_field = VectorField3D::new(dims);
        e_field.x.fill(1.0);

        dump.accumulate(&e_field, &h_field, 0.0, 1e-12);

        let mag = dump.extract_e_magnitude(0);
        assert!(mag.is_some());

        let mag_field = mag.unwrap();
        assert!(mag_field.get(2, 2, 2) > 0.0);
    }

    #[test]
    fn test_td_dump() {
        let mut dump = TimeDomainFieldDump::new(100, "field_dump");

        assert!(dump.should_dump(0));
        assert!(!dump.should_dump(50));
        assert!(dump.should_dump(100));

        assert_eq!(dump.get_filename(), "field_dump_00000000");

        dump.increment();
        assert_eq!(dump.get_filename(), "field_dump_00000001");
    }

    #[test]
    fn test_sub_sampling() {
        let dims = Dimensions::new(10, 10, 10);
        let freqs = vec![1e9];

        let mut dump = FrequencyDomainFieldDump::new(dims, freqs, FdDumpType::EField);
        dump.set_sub_sample([2, 2, 2]);

        // Sub-sampling should still work
        let mut e_field = VectorField3D::new(dims);
        let h_field = VectorField3D::new(dims);
        e_field.x.fill(1.0);

        dump.accumulate(&e_field, &h_field, 0.0, 1e-12);
        assert_eq!(dump.num_samples(), 1);
    }
}
