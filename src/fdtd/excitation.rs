//! Excitation sources for FDTD simulations.
//!
//! Provides various time-domain excitation waveforms for FDTD simulations.

use std::f64::consts::PI;

/// Excitation waveform type
#[derive(Debug, Clone)]
pub enum ExcitationType {
    /// Gaussian pulse: exp(-((t - t0) / tau)^2)
    Gaussian {
        /// Center time
        t0: f64,
        /// Pulse width (time constant)
        tau: f64,
    },
    /// Sinusoidal: sin(2*pi*f*t)
    Sinusoidal {
        /// Frequency (Hz)
        frequency: f64,
    },
    /// Modulated Gaussian: sin(2*pi*f*t) * exp(-((t - t0) / tau)^2)
    GaussianModulated {
        /// Center frequency (Hz)
        frequency: f64,
        /// Center time
        t0: f64,
        /// Pulse width
        tau: f64,
    },
    /// Dirac delta (impulse)
    Dirac,
    /// Step function (Heaviside)
    Step,
    /// Custom waveform from sampled data
    Custom {
        /// Time samples
        times: Vec<f64>,
        /// Amplitude samples
        values: Vec<f64>,
    },
}

/// Excitation source configuration
#[derive(Debug, Clone)]
pub struct Excitation {
    /// Excitation type/waveform
    pub excitation_type: ExcitationType,
    /// Polarization direction (0=x, 1=y, 2=z)
    pub direction: usize,
    /// Source position (i, j, k)
    pub position: (usize, usize, usize),
    /// Amplitude
    pub amplitude: f64,
    /// Is this a soft source? (adds to field instead of replacing)
    pub soft_source: bool,
}

impl Excitation {
    /// Create a Gaussian pulse excitation.
    ///
    /// # Arguments
    /// * `center_freq` - Center frequency of the excitation
    /// * `bandwidth` - Fractional bandwidth (0.0 to 1.0)
    /// * `direction` - Polarization direction (0=x, 1=y, 2=z)
    /// * `position` - Source position (i, j, k)
    pub fn gaussian(
        center_freq: f64,
        bandwidth: f64,
        direction: usize,
        position: (usize, usize, usize),
    ) -> Self {
        // Calculate Gaussian parameters from frequency specs
        let f_max = center_freq * (1.0 + 0.5 * bandwidth);
        let tau = 1.0 / (PI * f_max);
        let t0 = 4.0 * tau; // Start 4 time constants after t=0

        Self {
            excitation_type: ExcitationType::Gaussian { t0, tau },
            direction,
            position,
            amplitude: 1.0,
            soft_source: true,
        }
    }

    /// Create a modulated Gaussian pulse (common for narrowband excitation).
    pub fn gaussian_modulated(
        center_freq: f64,
        bandwidth: f64,
        direction: usize,
        position: (usize, usize, usize),
    ) -> Self {
        let tau = 1.0 / (PI * center_freq * bandwidth);
        let t0 = 4.0 * tau;

        Self {
            excitation_type: ExcitationType::GaussianModulated {
                frequency: center_freq,
                t0,
                tau,
            },
            direction,
            position,
            amplitude: 1.0,
            soft_source: true,
        }
    }

    /// Create a sinusoidal (CW) excitation.
    pub fn sinusoidal(frequency: f64, direction: usize, position: (usize, usize, usize)) -> Self {
        Self {
            excitation_type: ExcitationType::Sinusoidal { frequency },
            direction,
            position,
            amplitude: 1.0,
            soft_source: true,
        }
    }

    /// Set the amplitude.
    pub fn with_amplitude(mut self, amplitude: f64) -> Self {
        self.amplitude = amplitude;
        self
    }

    /// Set as hard source (replaces field value).
    pub fn hard_source(mut self) -> Self {
        self.soft_source = false;
        self
    }

    /// Evaluate the excitation at a given time.
    pub fn evaluate(&self, t: f64) -> f64 {
        let value = match &self.excitation_type {
            ExcitationType::Gaussian { t0, tau } => {
                let arg = (t - t0) / tau;
                (-arg * arg).exp()
            }
            ExcitationType::Sinusoidal { frequency } => (2.0 * PI * frequency * t).sin(),
            ExcitationType::GaussianModulated { frequency, t0, tau } => {
                let arg = (t - t0) / tau;
                let envelope = (-arg * arg).exp();
                let carrier = (2.0 * PI * frequency * t).sin();
                envelope * carrier
            }
            ExcitationType::Dirac => {
                if t.abs() < 1e-15 {
                    1.0
                } else {
                    0.0
                }
            }
            ExcitationType::Step => {
                if t >= 0.0 {
                    1.0
                } else {
                    0.0
                }
            }
            ExcitationType::Custom { times, values } => {
                // Linear interpolation
                if t <= times[0] {
                    values[0]
                } else if t >= *times.last().unwrap() {
                    *values.last().unwrap()
                } else {
                    // Find bracketing indices
                    let mut i = 0;
                    while i < times.len() - 1 && times[i + 1] < t {
                        i += 1;
                    }
                    let t0 = times[i];
                    let t1 = times[i + 1];
                    let v0 = values[i];
                    let v1 = values[i + 1];
                    let alpha = (t - t0) / (t1 - t0);
                    v0 + alpha * (v1 - v0)
                }
            }
        };

        self.amplitude * value
    }

    /// Get the signal duration (time until signal is effectively zero).
    pub fn duration(&self) -> f64 {
        match &self.excitation_type {
            ExcitationType::Gaussian { t0, tau }
            | ExcitationType::GaussianModulated { t0, tau, .. } => {
                t0 + 6.0 * tau // 6 sigma covers >99.9% of energy
            }
            ExcitationType::Sinusoidal { .. } => f64::INFINITY,
            ExcitationType::Dirac | ExcitationType::Step => f64::INFINITY,
            ExcitationType::Custom { times, .. } => *times.last().unwrap_or(&0.0),
        }
    }

    /// Check if the excitation is still active at time t.
    pub fn is_active(&self, t: f64) -> bool {
        match &self.excitation_type {
            ExcitationType::Gaussian { t0, tau }
            | ExcitationType::GaussianModulated { t0, tau, .. } => t <= t0 + 6.0 * tau,
            ExcitationType::Sinusoidal { .. } => true,
            ExcitationType::Dirac => t.abs() < 1e-15,
            ExcitationType::Step => true,
            ExcitationType::Custom { times, .. } => t <= *times.last().unwrap_or(&0.0),
        }
    }
}

/// Create excitation for a rectangular waveguide TE10 mode.
#[allow(dead_code)]
pub fn waveguide_te10(
    _a: f64,         // waveguide width
    freq: f64,       // frequency
    position: usize, // z-position
    ny: usize,       // grid size in y
) -> Vec<Excitation> {
    // TE10 mode: Ey varies as sin(pi*x/a), Ex = 0
    // For simplicity, create multiple point sources that approximate this
    let mut excitations = Vec::new();

    for j in 1..ny - 1 {
        let y_norm = j as f64 / (ny - 1) as f64;
        let amplitude = (PI * y_norm).sin(); // sin(pi*y/a) profile

        if amplitude > 0.01 {
            // Only include significant amplitudes
            let exc = Excitation::sinusoidal(freq, 1, (0, j, position)).with_amplitude(amplitude);
            excitations.push(exc);
        }
    }

    excitations
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gaussian_excitation() {
        let exc = Excitation::gaussian(1e9, 0.5, 2, (10, 10, 10));

        // At t=t0, should be maximum
        if let ExcitationType::Gaussian { t0, .. } = exc.excitation_type {
            let val = exc.evaluate(t0);
            assert!((val - 1.0).abs() < 1e-10);
        }

        // At t=0, should be very small
        let val_t0 = exc.evaluate(0.0);
        assert!(val_t0 < 0.01);
    }

    #[test]
    fn test_sinusoidal_excitation() {
        let freq = 1e9;
        let exc = Excitation::sinusoidal(freq, 0, (5, 5, 5));

        // Check periodicity
        let period = 1.0 / freq;
        let v1 = exc.evaluate(0.0);
        let v2 = exc.evaluate(period);
        assert!((v1 - v2).abs() < 1e-10);

        // Check value at quarter period
        let v_quarter = exc.evaluate(period / 4.0);
        assert!((v_quarter - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_gaussian_modulated_excitation() {
        let exc = Excitation::gaussian_modulated(1e9, 0.5, 1, (5, 5, 5));

        // Verify it's a GaussianModulated type
        assert!(matches!(
            exc.excitation_type,
            ExcitationType::GaussianModulated { .. }
        ));
        assert_eq!(exc.direction, 1);
        assert_eq!(exc.position, (5, 5, 5));
        assert_eq!(exc.amplitude, 1.0);
        assert!(exc.soft_source);

        // At t=0, envelope should be very small (t0 = 4*tau)
        let val_t0 = exc.evaluate(0.0);
        assert!(val_t0.abs() < 0.01);

        // At t=t0, envelope should be at maximum, but carrier modulates it
        if let ExcitationType::GaussianModulated { t0, frequency, .. } = exc.excitation_type {
            // Check that the signal has oscillations
            let v1 = exc.evaluate(t0);
            let v2 = exc.evaluate(t0 + 0.5 / frequency); // Half period later
                                                         // Signal should change sign due to carrier
            assert!(v1 * v2 <= 0.0 || v1.abs() < 0.1 || v2.abs() < 0.1);
        }
    }

    #[test]
    fn test_with_amplitude() {
        let exc = Excitation::gaussian(1e9, 0.5, 2, (10, 10, 10)).with_amplitude(2.5);

        assert!((exc.amplitude - 2.5).abs() < 1e-10);

        // At t=t0, value should be amplitude
        if let ExcitationType::Gaussian { t0, .. } = exc.excitation_type {
            let val = exc.evaluate(t0);
            assert!((val - 2.5).abs() < 1e-10);
        }
    }

    #[test]
    fn test_hard_source() {
        let exc = Excitation::gaussian(1e9, 0.5, 2, (10, 10, 10)).hard_source();

        assert!(!exc.soft_source);
    }

    #[test]
    fn test_dirac_excitation() {
        let exc = Excitation {
            excitation_type: ExcitationType::Dirac,
            direction: 0,
            position: (5, 5, 5),
            amplitude: 1.0,
            soft_source: true,
        };

        // At t=0, should be 1.0
        let val_at_zero = exc.evaluate(0.0);
        assert!((val_at_zero - 1.0).abs() < 1e-10);

        // At any other time, should be 0
        let val_at_1ns = exc.evaluate(1e-9);
        assert!(val_at_1ns.abs() < 1e-10);

        let val_at_neg = exc.evaluate(-1e-9);
        assert!(val_at_neg.abs() < 1e-10);
    }

    #[test]
    fn test_step_excitation() {
        let exc = Excitation {
            excitation_type: ExcitationType::Step,
            direction: 1,
            position: (5, 5, 5),
            amplitude: 3.0,
            soft_source: true,
        };

        // At t >= 0, should be amplitude
        let val_at_zero = exc.evaluate(0.0);
        assert!((val_at_zero - 3.0).abs() < 1e-10);

        let val_at_pos = exc.evaluate(1e-9);
        assert!((val_at_pos - 3.0).abs() < 1e-10);

        // At t < 0, should be 0
        let val_at_neg = exc.evaluate(-1e-9);
        assert!(val_at_neg.abs() < 1e-10);
    }

    #[test]
    fn test_custom_excitation() {
        let times = vec![0.0, 1.0, 2.0, 3.0];
        let values = vec![0.0, 1.0, 0.5, 0.0];

        let exc = Excitation {
            excitation_type: ExcitationType::Custom {
                times: times.clone(),
                values: values.clone(),
            },
            direction: 2,
            position: (5, 5, 5),
            amplitude: 1.0,
            soft_source: true,
        };

        // At exact sample points
        assert!((exc.evaluate(0.0) - 0.0).abs() < 1e-10);
        assert!((exc.evaluate(1.0) - 1.0).abs() < 1e-10);
        assert!((exc.evaluate(2.0) - 0.5).abs() < 1e-10);
        assert!((exc.evaluate(3.0) - 0.0).abs() < 1e-10);

        // Linear interpolation between points
        let val_half = exc.evaluate(0.5); // Between 0.0 and 1.0
        assert!((val_half - 0.5).abs() < 1e-10);

        let val_1_5 = exc.evaluate(1.5); // Between 1.0 and 0.5
        assert!((val_1_5 - 0.75).abs() < 1e-10);

        // Before first sample
        let val_neg = exc.evaluate(-1.0);
        assert!((val_neg - 0.0).abs() < 1e-10);

        // After last sample
        let val_after = exc.evaluate(5.0);
        assert!((val_after - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_duration_gaussian() {
        let exc = Excitation::gaussian(1e9, 0.5, 2, (10, 10, 10));

        if let ExcitationType::Gaussian { t0, tau } = exc.excitation_type {
            let expected_duration = t0 + 6.0 * tau;
            assert!((exc.duration() - expected_duration).abs() < 1e-20);
        }
    }

    #[test]
    fn test_duration_gaussian_modulated() {
        let exc = Excitation::gaussian_modulated(1e9, 0.5, 1, (5, 5, 5));

        if let ExcitationType::GaussianModulated { t0, tau, .. } = exc.excitation_type {
            let expected_duration = t0 + 6.0 * tau;
            assert!((exc.duration() - expected_duration).abs() < 1e-20);
        }
    }

    #[test]
    fn test_duration_sinusoidal() {
        let exc = Excitation::sinusoidal(1e9, 0, (5, 5, 5));
        assert!(exc.duration().is_infinite());
    }

    #[test]
    fn test_duration_dirac() {
        let exc = Excitation {
            excitation_type: ExcitationType::Dirac,
            direction: 0,
            position: (5, 5, 5),
            amplitude: 1.0,
            soft_source: true,
        };
        assert!(exc.duration().is_infinite());
    }

    #[test]
    fn test_duration_step() {
        let exc = Excitation {
            excitation_type: ExcitationType::Step,
            direction: 0,
            position: (5, 5, 5),
            amplitude: 1.0,
            soft_source: true,
        };
        assert!(exc.duration().is_infinite());
    }

    #[test]
    fn test_duration_custom() {
        let times = vec![0.0, 1.0, 2.0, 3.0];
        let values = vec![0.0, 1.0, 0.5, 0.0];

        let exc = Excitation {
            excitation_type: ExcitationType::Custom { times, values },
            direction: 2,
            position: (5, 5, 5),
            amplitude: 1.0,
            soft_source: true,
        };

        assert!((exc.duration() - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_duration_custom_empty() {
        let exc = Excitation {
            excitation_type: ExcitationType::Custom {
                times: vec![],
                values: vec![],
            },
            direction: 2,
            position: (5, 5, 5),
            amplitude: 1.0,
            soft_source: true,
        };

        assert!((exc.duration() - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_is_active_gaussian() {
        let exc = Excitation::gaussian(1e9, 0.5, 2, (10, 10, 10));

        if let ExcitationType::Gaussian { t0, tau } = exc.excitation_type {
            let end_time = t0 + 6.0 * tau;

            assert!(exc.is_active(0.0));
            assert!(exc.is_active(t0));
            assert!(exc.is_active(end_time));
            assert!(!exc.is_active(end_time + 1.0));
        }
    }

    #[test]
    fn test_is_active_gaussian_modulated() {
        let exc = Excitation::gaussian_modulated(1e9, 0.5, 1, (5, 5, 5));

        if let ExcitationType::GaussianModulated { t0, tau, .. } = exc.excitation_type {
            let end_time = t0 + 6.0 * tau;

            assert!(exc.is_active(0.0));
            assert!(exc.is_active(t0));
            assert!(exc.is_active(end_time));
            assert!(!exc.is_active(end_time + 1.0));
        }
    }

    #[test]
    fn test_is_active_sinusoidal() {
        let exc = Excitation::sinusoidal(1e9, 0, (5, 5, 5));

        assert!(exc.is_active(0.0));
        assert!(exc.is_active(1e10)); // Always active
    }

    #[test]
    fn test_is_active_dirac() {
        let exc = Excitation {
            excitation_type: ExcitationType::Dirac,
            direction: 0,
            position: (5, 5, 5),
            amplitude: 1.0,
            soft_source: true,
        };

        assert!(exc.is_active(0.0));
        assert!(!exc.is_active(1e-9));
        assert!(!exc.is_active(-1e-9));
    }

    #[test]
    fn test_is_active_step() {
        let exc = Excitation {
            excitation_type: ExcitationType::Step,
            direction: 0,
            position: (5, 5, 5),
            amplitude: 1.0,
            soft_source: true,
        };

        // Step is always active
        assert!(exc.is_active(0.0));
        assert!(exc.is_active(1e10));
    }

    #[test]
    fn test_is_active_custom() {
        let times = vec![0.0, 1.0, 2.0, 3.0];
        let values = vec![0.0, 1.0, 0.5, 0.0];

        let exc = Excitation {
            excitation_type: ExcitationType::Custom { times, values },
            direction: 2,
            position: (5, 5, 5),
            amplitude: 1.0,
            soft_source: true,
        };

        assert!(exc.is_active(0.0));
        assert!(exc.is_active(1.5));
        assert!(exc.is_active(3.0));
        assert!(!exc.is_active(4.0));
    }

    #[test]
    fn test_waveguide_te10_creation() {
        let a = 0.01; // 10mm waveguide width
        let freq = 10e9; // 10 GHz
        let position = 5;
        let ny = 20;

        let excitations = waveguide_te10(a, freq, position, ny);

        // Should create multiple excitations
        assert!(!excitations.is_empty());

        // All excitations should be sinusoidal
        for exc in &excitations {
            assert!(matches!(
                exc.excitation_type,
                ExcitationType::Sinusoidal { .. }
            ));
            // Direction should be y (1) for TE10 mode
            assert_eq!(exc.direction, 1);
            // z-position should match
            assert_eq!(exc.position.2, position);
        }

        // Center excitations should have higher amplitude (sin profile)
        let mid_idx = excitations.len() / 2;
        let center_amp = excitations[mid_idx].amplitude;
        let edge_amp = excitations[0].amplitude;
        assert!(center_amp > edge_amp);
    }

    #[test]
    fn test_custom_excitation_with_amplitude() {
        let times = vec![0.0, 1.0, 2.0];
        let values = vec![0.0, 1.0, 0.0];

        let exc = Excitation {
            excitation_type: ExcitationType::Custom { times, values },
            direction: 2,
            position: (5, 5, 5),
            amplitude: 5.0,
            soft_source: true,
        };

        // Value at t=1.0 should be amplitude * 1.0 = 5.0
        assert!((exc.evaluate(1.0) - 5.0).abs() < 1e-10);

        // Value at t=0.5 should be amplitude * 0.5 = 2.5
        assert!((exc.evaluate(0.5) - 2.5).abs() < 1e-10);
    }
}
