//! Steady-State Detection Extension.
//!
//! Monitors field values over time to detect when a simulation has reached
//! steady state (periodic behavior), allowing early termination.

use crate::arrays::{Dimensions, VectorField3D};
use std::collections::VecDeque;

/// Configuration for steady-state detection.
#[derive(Debug, Clone)]
pub struct SteadyStateConfig {
    /// Number of periods to check for convergence
    pub num_periods: usize,
    /// Convergence threshold (relative difference)
    pub threshold: f64,
    /// Minimum number of timesteps before checking
    pub min_timesteps: u64,
    /// Check interval (timesteps between checks)
    pub check_interval: u64,
}

impl Default for SteadyStateConfig {
    fn default() -> Self {
        Self {
            num_periods: 3,
            threshold: 1e-6,
            min_timesteps: 1000,
            check_interval: 100,
        }
    }
}

impl SteadyStateConfig {
    /// Create with custom threshold.
    pub fn with_threshold(mut self, threshold: f64) -> Self {
        self.threshold = threshold;
        self
    }

    /// Create with custom minimum timesteps.
    pub fn with_min_timesteps(mut self, min: u64) -> Self {
        self.min_timesteps = min;
        self
    }
}

/// Probe point for monitoring.
#[derive(Debug, Clone)]
struct ProbePoint {
    /// Position in grid
    position: [usize; 3],
    /// Field component (0=Ex, 1=Ey, 2=Ez, 3=Hx, 4=Hy, 5=Hz)
    component: usize,
    /// History of values
    history: VecDeque<f32>,
    /// Maximum history size
    max_history: usize,
}

impl ProbePoint {
    fn new(position: [usize; 3], component: usize, max_history: usize) -> Self {
        Self {
            position,
            component,
            history: VecDeque::with_capacity(max_history),
            max_history,
        }
    }

    fn record(&mut self, value: f32) {
        if self.history.len() >= self.max_history {
            self.history.pop_front();
        }
        self.history.push_back(value);
    }

    fn get_value(&self, e_field: &VectorField3D, h_field: &VectorField3D) -> f32 {
        let pos = self.position;
        match self.component {
            0 => e_field.x.get(pos[0], pos[1], pos[2]),
            1 => e_field.y.get(pos[0], pos[1], pos[2]),
            2 => e_field.z.get(pos[0], pos[1], pos[2]),
            3 => h_field.x.get(pos[0], pos[1], pos[2]),
            4 => h_field.y.get(pos[0], pos[1], pos[2]),
            _ => h_field.z.get(pos[0], pos[1], pos[2]),
        }
    }
}

/// Steady-state detection result.
#[derive(Debug, Clone)]
pub struct SteadyStateResult {
    /// Whether steady state was detected
    pub converged: bool,
    /// Maximum relative difference found
    pub max_difference: f64,
    /// Timestep at which check was performed
    pub timestep: u64,
}

/// Steady-State Detection Extension.
///
/// Monitors E and H field probes to detect when the simulation
/// has reached a periodic steady state.
pub struct SteadyStateDetector {
    /// Grid dimensions
    dims: Dimensions,
    /// Configuration
    config: SteadyStateConfig,
    /// Probe points
    probes: Vec<ProbePoint>,
    /// Current timestep
    timestep: u64,
    /// Last check result
    last_result: Option<SteadyStateResult>,
    /// Period estimate (in timesteps)
    period_estimate: Option<u64>,
}

impl SteadyStateDetector {
    /// Create a new steady-state detector.
    pub fn new(dims: Dimensions, config: SteadyStateConfig) -> Self {
        Self {
            dims,
            config,
            probes: Vec::new(),
            timestep: 0,
            last_result: None,
            period_estimate: None,
        }
    }

    /// Add a probe point for monitoring.
    pub fn add_probe(&mut self, position: [usize; 3], component: usize) {
        // Keep enough history for period detection
        let max_history = self.config.num_periods * 1000; // Assume max 1000 samples per period
        self.probes.push(ProbePoint::new(position, component, max_history));
    }

    /// Add default probes at strategic locations.
    pub fn add_default_probes(&mut self) {
        // Add probes at center and near boundaries
        let cx = self.dims.nx / 2;
        let cy = self.dims.ny / 2;
        let cz = self.dims.nz / 2;

        // Center Ez
        self.add_probe([cx, cy, cz], 2);

        // Near X boundaries
        self.add_probe([self.dims.nx / 4, cy, cz], 2);
        self.add_probe([3 * self.dims.nx / 4, cy, cz], 2);

        // Add Hx probe
        self.add_probe([cx, cy, cz], 3);
    }

    /// Set the expected period in timesteps.
    pub fn set_period(&mut self, period: u64) {
        self.period_estimate = Some(period);
    }

    /// Record current field values.
    pub fn record(&mut self, e_field: &VectorField3D, h_field: &VectorField3D) {
        for probe in &mut self.probes {
            let value = probe.get_value(e_field, h_field);
            probe.record(value);
        }
        self.timestep += 1;
    }

    /// Check if steady state has been reached.
    pub fn check(&mut self) -> SteadyStateResult {
        // Don't check if not enough timesteps
        if self.timestep < self.config.min_timesteps {
            return SteadyStateResult {
                converged: false,
                max_difference: f64::INFINITY,
                timestep: self.timestep,
            };
        }

        // Need period estimate for comparison
        let period = match self.period_estimate {
            Some(p) => p as usize,
            None => {
                // Try to estimate period from autocorrelation
                self.estimate_period().unwrap_or(100)
            }
        };

        if period == 0 {
            return SteadyStateResult {
                converged: false,
                max_difference: f64::INFINITY,
                timestep: self.timestep,
            };
        }

        // Compare values across multiple periods
        let mut max_diff = 0.0f64;

        for probe in &self.probes {
            let n = probe.history.len();
            if n < period * self.config.num_periods {
                return SteadyStateResult {
                    converged: false,
                    max_difference: f64::INFINITY,
                    timestep: self.timestep,
                };
            }

            // Compare most recent period with previous periods
            for i in 0..period {
                let current = probe.history[n - period + i] as f64;
                let previous = probe.history[n - 2 * period + i] as f64;

                let diff = if current.abs() > 1e-20 {
                    ((current - previous) / current).abs()
                } else if previous.abs() > 1e-20 {
                    ((current - previous) / previous).abs()
                } else {
                    0.0
                };

                max_diff = max_diff.max(diff);
            }
        }

        let converged = max_diff < self.config.threshold;

        let result = SteadyStateResult {
            converged,
            max_difference: max_diff,
            timestep: self.timestep,
        };

        self.last_result = Some(result.clone());
        result
    }

    /// Estimate the period from probe history using autocorrelation.
    fn estimate_period(&self) -> Option<usize> {
        if self.probes.is_empty() {
            return None;
        }

        let probe = &self.probes[0];
        let n = probe.history.len();

        if n < 100 {
            return None;
        }

        // Simple zero-crossing period detection
        let mut zero_crossings = Vec::new();
        let data: Vec<f64> = probe.history.iter().map(|&x| x as f64).collect();

        for i in 1..n {
            if (data[i - 1] < 0.0 && data[i] >= 0.0) || (data[i - 1] >= 0.0 && data[i] < 0.0) {
                zero_crossings.push(i);
            }
        }

        if zero_crossings.len() < 4 {
            return None;
        }

        // Calculate average period between zero crossings (half period)
        let mut sum = 0;
        for i in 1..zero_crossings.len() {
            sum += zero_crossings[i] - zero_crossings[i - 1];
        }
        let avg_half_period = sum / (zero_crossings.len() - 1);

        Some(avg_half_period * 2)
    }

    /// Get the last check result.
    pub fn last_result(&self) -> Option<&SteadyStateResult> {
        self.last_result.as_ref()
    }

    /// Check if steady state has been achieved.
    pub fn is_converged(&self) -> bool {
        self.last_result.as_ref().map(|r| r.converged).unwrap_or(false)
    }

    /// Reset the detector.
    pub fn reset(&mut self) {
        for probe in &mut self.probes {
            probe.history.clear();
        }
        self.timestep = 0;
        self.last_result = None;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_config() {
        let config = SteadyStateConfig::default()
            .with_threshold(1e-5)
            .with_min_timesteps(500);

        assert_eq!(config.threshold, 1e-5);
        assert_eq!(config.min_timesteps, 500);
    }

    #[test]
    fn test_detector_creation() {
        let dims = Dimensions { nx: 20, ny: 20, nz: 20 };
        let config = SteadyStateConfig::default();
        let detector = SteadyStateDetector::new(dims, config);

        assert!(!detector.is_converged());
    }

    #[test]
    fn test_add_probes() {
        let dims = Dimensions { nx: 20, ny: 20, nz: 20 };
        let config = SteadyStateConfig::default();
        let mut detector = SteadyStateDetector::new(dims, config);

        detector.add_probe([10, 10, 10], 2); // Ez
        detector.add_probe([5, 5, 5], 3);    // Hx

        assert_eq!(detector.probes.len(), 2);
    }

    #[test]
    fn test_default_probes() {
        let dims = Dimensions { nx: 40, ny: 40, nz: 40 };
        let config = SteadyStateConfig::default();
        let mut detector = SteadyStateDetector::new(dims, config);

        detector.add_default_probes();

        assert!(detector.probes.len() >= 4);
    }

    #[test]
    fn test_record() {
        let dims = Dimensions { nx: 20, ny: 20, nz: 20 };
        let config = SteadyStateConfig::default();
        let mut detector = SteadyStateDetector::new(dims, config);

        detector.add_probe([10, 10, 10], 2);

        let e_field = VectorField3D::new(dims);
        let h_field = VectorField3D::new(dims);

        detector.record(&e_field, &h_field);
        assert_eq!(detector.timestep, 1);

        detector.record(&e_field, &h_field);
        assert_eq!(detector.timestep, 2);
    }

    #[test]
    fn test_not_converged_early() {
        let dims = Dimensions { nx: 20, ny: 20, nz: 20 };
        let config = SteadyStateConfig::default().with_min_timesteps(100);
        let mut detector = SteadyStateDetector::new(dims, config);

        detector.add_probe([10, 10, 10], 2);

        let e_field = VectorField3D::new(dims);
        let h_field = VectorField3D::new(dims);

        // Record less than min_timesteps
        for _ in 0..50 {
            detector.record(&e_field, &h_field);
        }

        let result = detector.check();
        assert!(!result.converged);
    }

    #[test]
    fn test_reset() {
        let dims = Dimensions { nx: 20, ny: 20, nz: 20 };
        let config = SteadyStateConfig::default();
        let mut detector = SteadyStateDetector::new(dims, config);

        detector.add_probe([10, 10, 10], 2);
        detector.timestep = 100;

        detector.reset();

        assert_eq!(detector.timestep, 0);
        assert!(detector.last_result.is_none());
    }
}
