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
}
