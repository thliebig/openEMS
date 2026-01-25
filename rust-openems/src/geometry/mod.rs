//! Geometry and grid definitions for FDTD simulations.
//!
//! This module provides grid definitions and coordinate system support
//! for both Cartesian and cylindrical FDTD simulations.

use crate::arrays::Dimensions;

/// Supported coordinate systems
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CoordinateSystem {
    /// Standard Cartesian (x, y, z)
    Cartesian,
    /// Cylindrical (r, phi, z)
    Cylindrical,
}

impl Default for CoordinateSystem {
    fn default() -> Self {
        Self::Cartesian
    }
}

/// Grid definition for FDTD simulation.
///
/// Stores the mesh lines in each direction and provides methods
/// for coordinate transformations and cell size calculations.
#[derive(Clone)]
pub struct Grid {
    /// Coordinate system
    coord_system: CoordinateSystem,
    /// X (or R) coordinate lines
    x_lines: Vec<f64>,
    /// Y (or Phi) coordinate lines
    y_lines: Vec<f64>,
    /// Z coordinate lines
    z_lines: Vec<f64>,
    /// Unit scaling factor (default 1.0 = meters)
    unit: f64,
}

impl Grid {
    /// Create a new grid with explicit mesh lines.
    pub fn new(
        coord_system: CoordinateSystem,
        x_lines: Vec<f64>,
        y_lines: Vec<f64>,
        z_lines: Vec<f64>,
    ) -> Self {
        Self {
            coord_system,
            x_lines,
            y_lines,
            z_lines,
            unit: 1.0,
        }
    }

    /// Create a uniform Cartesian grid.
    ///
    /// # Arguments
    /// * `nx` - Number of cells in x direction
    /// * `ny` - Number of cells in y direction
    /// * `nz` - Number of cells in z direction
    /// * `delta` - Cell size in all directions (meters)
    pub fn uniform(nx: usize, ny: usize, nz: usize, delta: f64) -> Self {
        let x_lines: Vec<f64> = (0..=nx).map(|i| i as f64 * delta).collect();
        let y_lines: Vec<f64> = (0..=ny).map(|i| i as f64 * delta).collect();
        let z_lines: Vec<f64> = (0..=nz).map(|i| i as f64 * delta).collect();

        Self::new(CoordinateSystem::Cartesian, x_lines, y_lines, z_lines)
    }

    /// Create a Cartesian grid from mesh line vectors.
    pub fn cartesian(x_lines: Vec<f64>, y_lines: Vec<f64>, z_lines: Vec<f64>) -> Self {
        Self::new(CoordinateSystem::Cartesian, x_lines, y_lines, z_lines)
    }

    /// Set the unit scaling factor (e.g., 1e-3 for millimeters).
    pub fn set_unit(&mut self, unit: f64) {
        self.unit = unit;
    }

    /// Get the coordinate system.
    #[inline]
    pub fn coord_system(&self) -> CoordinateSystem {
        self.coord_system
    }

    /// Get grid dimensions (number of cells in each direction).
    pub fn dimensions(&self) -> Dimensions {
        // Number of cells = number of lines - 1
        Dimensions::new(
            self.x_lines.len().saturating_sub(1).max(1),
            self.y_lines.len().saturating_sub(1).max(1),
            self.z_lines.len().saturating_sub(1).max(1),
        )
    }

    /// Get X/R mesh lines.
    #[inline]
    pub fn x_lines(&self) -> &[f64] {
        &self.x_lines
    }

    /// Get Y/Phi mesh lines.
    #[inline]
    pub fn y_lines(&self) -> &[f64] {
        &self.y_lines
    }

    /// Get Z mesh lines.
    #[inline]
    pub fn z_lines(&self) -> &[f64] {
        &self.z_lines
    }

    /// Get the minimum cell size in each direction.
    pub fn cell_size(&self) -> (f64, f64, f64) {
        let dx = Self::min_delta(&self.x_lines) * self.unit;
        let dy = Self::min_delta(&self.y_lines) * self.unit;
        let dz = Self::min_delta(&self.z_lines) * self.unit;
        (dx, dy, dz)
    }

    /// Get the cell size at a specific location.
    pub fn cell_size_at(&self, i: usize, j: usize, k: usize) -> (f64, f64, f64) {
        let dx = (self.x_lines[i + 1] - self.x_lines[i]) * self.unit;
        let dy = (self.y_lines[j + 1] - self.y_lines[j]) * self.unit;
        let dz = (self.z_lines[k + 1] - self.z_lines[k]) * self.unit;
        (dx, dy, dz)
    }

    /// Get the cell center coordinates.
    pub fn cell_center(&self, i: usize, j: usize, k: usize) -> (f64, f64, f64) {
        let x = 0.5 * (self.x_lines[i] + self.x_lines[i + 1]) * self.unit;
        let y = 0.5 * (self.y_lines[j] + self.y_lines[j + 1]) * self.unit;
        let z = 0.5 * (self.z_lines[k] + self.z_lines[k + 1]) * self.unit;
        (x, y, z)
    }

    /// Get the total simulation volume.
    pub fn volume(&self) -> f64 {
        let (x_min, x_max) = (*self.x_lines.first().unwrap(), *self.x_lines.last().unwrap());
        let (y_min, y_max) = (*self.y_lines.first().unwrap(), *self.y_lines.last().unwrap());
        let (z_min, z_max) = (*self.z_lines.first().unwrap(), *self.z_lines.last().unwrap());

        let unit3 = self.unit * self.unit * self.unit;
        (x_max - x_min) * (y_max - y_min) * (z_max - z_min) * unit3
    }

    /// Get the simulation domain bounds.
    pub fn bounds(&self) -> ((f64, f64), (f64, f64), (f64, f64)) {
        (
            (
                *self.x_lines.first().unwrap() * self.unit,
                *self.x_lines.last().unwrap() * self.unit,
            ),
            (
                *self.y_lines.first().unwrap() * self.unit,
                *self.y_lines.last().unwrap() * self.unit,
            ),
            (
                *self.z_lines.first().unwrap() * self.unit,
                *self.z_lines.last().unwrap() * self.unit,
            ),
        )
    }

    /// Add mesh lines in X direction.
    pub fn add_x_lines(&mut self, lines: &[f64]) {
        self.x_lines.extend_from_slice(lines);
        self.x_lines.sort_by(|a, b| a.partial_cmp(b).unwrap());
        self.x_lines.dedup();
    }

    /// Add mesh lines in Y direction.
    pub fn add_y_lines(&mut self, lines: &[f64]) {
        self.y_lines.extend_from_slice(lines);
        self.y_lines.sort_by(|a, b| a.partial_cmp(b).unwrap());
        self.y_lines.dedup();
    }

    /// Add mesh lines in Z direction.
    pub fn add_z_lines(&mut self, lines: &[f64]) {
        self.z_lines.extend_from_slice(lines);
        self.z_lines.sort_by(|a, b| a.partial_cmp(b).unwrap());
        self.z_lines.dedup();
    }

    /// Smooth mesh lines with maximum ratio between adjacent cells.
    pub fn smooth_mesh(&mut self, max_ratio: f64) {
        self.x_lines = Self::smooth_lines(&self.x_lines, max_ratio);
        self.y_lines = Self::smooth_lines(&self.y_lines, max_ratio);
        self.z_lines = Self::smooth_lines(&self.z_lines, max_ratio);
    }

    /// Smooth a single set of mesh lines.
    fn smooth_lines(lines: &[f64], max_ratio: f64) -> Vec<f64> {
        if lines.len() < 3 {
            return lines.to_vec();
        }

        let mut result = lines.to_vec();
        let mut changed = true;

        while changed {
            changed = false;
            let mut new_lines = Vec::new();
            new_lines.push(result[0]);

            for i in 0..result.len() - 1 {
                let delta_curr = result[i + 1] - result[i];

                if i > 0 {
                    let delta_prev = result[i] - result[i - 1];
                    let ratio = delta_curr / delta_prev;

                    if ratio > max_ratio {
                        // Insert intermediate point
                        let mid = result[i] + delta_prev * max_ratio;
                        if mid < result[i + 1] {
                            new_lines.push(mid);
                            changed = true;
                        }
                    }
                }

                new_lines.push(result[i + 1]);
            }

            result = new_lines;
        }

        result
    }

    /// Find minimum delta in a set of lines.
    fn min_delta(lines: &[f64]) -> f64 {
        if lines.len() < 2 {
            return 1.0; // Default
        }

        let mut min = f64::MAX;
        for i in 0..lines.len() - 1 {
            let delta = lines[i + 1] - lines[i];
            if delta > 0.0 && delta < min {
                min = delta;
            }
        }
        min
    }

    /// Find the cell index containing a given coordinate.
    pub fn find_cell(&self, x: f64, y: f64, z: f64) -> Option<(usize, usize, usize)> {
        let x_scaled = x / self.unit;
        let y_scaled = y / self.unit;
        let z_scaled = z / self.unit;

        let i = Self::find_index(&self.x_lines, x_scaled)?;
        let j = Self::find_index(&self.y_lines, y_scaled)?;
        let k = Self::find_index(&self.z_lines, z_scaled)?;

        Some((i, j, k))
    }

    /// Find index in sorted lines array.
    fn find_index(lines: &[f64], value: f64) -> Option<usize> {
        if value < lines[0] || value > *lines.last()? {
            return None;
        }

        for i in 0..lines.len() - 1 {
            if value >= lines[i] && value <= lines[i + 1] {
                return Some(i);
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_uniform_grid() {
        let grid = Grid::uniform(10, 20, 30, 0.001);
        let dims = grid.dimensions();

        assert_eq!(dims.nx, 10);
        assert_eq!(dims.ny, 20);
        assert_eq!(dims.nz, 30);

        let (dx, dy, dz) = grid.cell_size();
        assert!((dx - 0.001).abs() < 1e-10);
        assert!((dy - 0.001).abs() < 1e-10);
        assert!((dz - 0.001).abs() < 1e-10);
    }

    #[test]
    fn test_find_cell() {
        let grid = Grid::uniform(10, 10, 10, 0.001);

        let cell = grid.find_cell(0.0055, 0.0055, 0.0055);
        assert_eq!(cell, Some((5, 5, 5)));

        let outside = grid.find_cell(0.02, 0.0, 0.0);
        assert_eq!(outside, None);
    }
}
