//! Timestep calculation and information.

use crate::constants::C0;
use crate::geometry::Grid;

/// Information about the simulation timestep.
#[derive(Debug, Clone)]
pub struct TimestepInfo {
    /// Timestep value in seconds
    pub dt: f64,
    /// Maximum stable timestep (CFL limit)
    pub dt_max: f64,
    /// Nyquist frequency for this timestep
    pub nyquist_freq: f64,
    /// Number of cells in the simulation
    pub num_cells: usize,
    /// Estimated memory usage in bytes
    pub memory_bytes: usize,
}

impl TimestepInfo {
    /// Calculate timestep information for a grid.
    pub fn calculate(grid: &Grid) -> Self {
        let (dx, dy, dz) = grid.cell_size();
        let dims = grid.dimensions();

        // CFL condition for 3D: dt <= 1 / (c * sqrt(1/dx^2 + 1/dy^2 + 1/dz^2))
        let inv_sum = 1.0 / (dx * dx) + 1.0 / (dy * dy) + 1.0 / (dz * dz);
        let dt_max = 1.0 / (C0 * inv_sum.sqrt());

        // Use 99% of maximum for stability margin
        let dt = 0.99 * dt_max;

        // Nyquist frequency
        let nyquist_freq = 1.0 / (2.0 * dt);

        let num_cells = dims.total();

        // Memory estimate: 6 field components (E and H) + 12 coefficient arrays
        // Each is f32 = 4 bytes
        let fields_per_cell = 6 + 12; // E(3) + H(3) + Ca(3) + Cb(3) + Da(3) + Db(3)
        let memory_bytes = num_cells * fields_per_cell * 4;

        Self {
            dt,
            dt_max,
            nyquist_freq,
            num_cells,
            memory_bytes,
        }
    }

    /// Calculate number of timesteps needed for given signal duration.
    pub fn timesteps_for_duration(&self, duration: f64) -> usize {
        (duration / self.dt).ceil() as usize
    }

    /// Calculate number of timesteps per period at given frequency.
    pub fn timesteps_per_period(&self, frequency: f64) -> f64 {
        1.0 / (frequency * self.dt)
    }

    /// Format memory size for display.
    pub fn memory_display(&self) -> String {
        let bytes = self.memory_bytes as f64;
        if bytes < 1024.0 {
            format!("{:.0} B", bytes)
        } else if bytes < 1024.0 * 1024.0 {
            format!("{:.2} KB", bytes / 1024.0)
        } else if bytes < 1024.0 * 1024.0 * 1024.0 {
            format!("{:.2} MB", bytes / (1024.0 * 1024.0))
        } else {
            format!("{:.2} GB", bytes / (1024.0 * 1024.0 * 1024.0))
        }
    }
}
