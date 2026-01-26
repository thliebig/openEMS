//! Geometry and grid definitions for FDTD simulations.
//!
//! This module provides grid definitions and coordinate system support
//! for both Cartesian and cylindrical FDTD simulations.

use crate::arrays::Dimensions;

/// Supported coordinate systems
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum CoordinateSystem {
    /// Standard Cartesian (x, y, z)
    #[default]
    Cartesian,
    /// Cylindrical (r, phi, z)
    Cylindrical,
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
        let (x_min, x_max) = (
            *self.x_lines.first().unwrap(),
            *self.x_lines.last().unwrap(),
        );
        let (y_min, y_max) = (
            *self.y_lines.first().unwrap(),
            *self.y_lines.last().unwrap(),
        );
        let (z_min, z_max) = (
            *self.z_lines.first().unwrap(),
            *self.z_lines.last().unwrap(),
        );

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

        (0..lines.len() - 1).find(|&i| value >= lines[i] && value <= lines[i + 1])
    }

    /// Get delta (cell size) in X direction at index.
    #[inline]
    pub fn delta_x(&self, i: usize) -> f64 {
        if i + 1 < self.x_lines.len() {
            (self.x_lines[i + 1] - self.x_lines[i]) * self.unit
        } else if i > 0 {
            (self.x_lines[i] - self.x_lines[i - 1]) * self.unit
        } else {
            self.unit
        }
    }

    /// Get delta (cell size) in Y direction at index.
    #[inline]
    pub fn delta_y(&self, j: usize) -> f64 {
        if j + 1 < self.y_lines.len() {
            (self.y_lines[j + 1] - self.y_lines[j]) * self.unit
        } else if j > 0 {
            (self.y_lines[j] - self.y_lines[j - 1]) * self.unit
        } else {
            self.unit
        }
    }

    /// Get delta (cell size) in Z direction at index.
    #[inline]
    pub fn delta_z(&self, k: usize) -> f64 {
        if k + 1 < self.z_lines.len() {
            (self.z_lines[k + 1] - self.z_lines[k]) * self.unit
        } else if k > 0 {
            (self.z_lines[k] - self.z_lines[k - 1]) * self.unit
        } else {
            self.unit
        }
    }

    /// Get X coordinate at index.
    #[inline]
    pub fn x_line(&self, i: usize) -> f64 {
        self.x_lines.get(i).copied().unwrap_or(0.0) * self.unit
    }

    /// Get Y coordinate at index.
    #[inline]
    pub fn y_line(&self, j: usize) -> f64 {
        self.y_lines.get(j).copied().unwrap_or(0.0) * self.unit
    }

    /// Get Z coordinate at index.
    #[inline]
    pub fn z_line(&self, k: usize) -> f64 {
        self.z_lines.get(k).copied().unwrap_or(0.0) * self.unit
    }

    /// Find cell index in X direction for a coordinate.
    pub fn find_cell_x(&self, x: f64) -> usize {
        let x_scaled = x / self.unit;
        for i in 0..self.x_lines.len().saturating_sub(1) {
            if x_scaled >= self.x_lines[i] && x_scaled <= self.x_lines[i + 1] {
                return i;
            }
        }
        self.x_lines.len().saturating_sub(2)
    }

    /// Find cell index in Y direction for a coordinate.
    pub fn find_cell_y(&self, y: f64) -> usize {
        let y_scaled = y / self.unit;
        for j in 0..self.y_lines.len().saturating_sub(1) {
            if y_scaled >= self.y_lines[j] && y_scaled <= self.y_lines[j + 1] {
                return j;
            }
        }
        self.y_lines.len().saturating_sub(2)
    }

    /// Find cell index in Z direction for a coordinate.
    pub fn find_cell_z(&self, z: f64) -> usize {
        let z_scaled = z / self.unit;
        for k in 0..self.z_lines.len().saturating_sub(1) {
            if z_scaled >= self.z_lines[k] && z_scaled <= self.z_lines[k + 1] {
                return k;
            }
        }
        self.z_lines.len().saturating_sub(2)
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

    #[test]
    fn test_cartesian_constructor() {
        let x_lines = vec![0.0, 0.001, 0.003, 0.006];
        let y_lines = vec![0.0, 0.002, 0.004];
        let z_lines = vec![0.0, 0.001, 0.002, 0.003, 0.004];

        let grid = Grid::cartesian(x_lines.clone(), y_lines.clone(), z_lines.clone());

        assert_eq!(grid.coord_system(), CoordinateSystem::Cartesian);

        let dims = grid.dimensions();
        assert_eq!(dims.nx, 3); // 4 lines = 3 cells
        assert_eq!(dims.ny, 2); // 3 lines = 2 cells
        assert_eq!(dims.nz, 4); // 5 lines = 4 cells
    }

    #[test]
    fn test_set_unit_scaling() {
        let mut grid = Grid::uniform(10, 10, 10, 1.0); // 1 meter cells

        // Set unit to millimeters (1mm = 0.001m)
        grid.set_unit(0.001);

        // Cell size should now be in millimeters
        let (dx, dy, dz) = grid.cell_size();
        assert!((dx - 0.001).abs() < 1e-10);
        assert!((dy - 0.001).abs() < 1e-10);
        assert!((dz - 0.001).abs() < 1e-10);
    }

    #[test]
    fn test_coord_system_getter() {
        let grid_cartesian = Grid::uniform(10, 10, 10, 0.001);
        assert_eq!(grid_cartesian.coord_system(), CoordinateSystem::Cartesian);

        let grid_cylindrical = Grid::new(
            CoordinateSystem::Cylindrical,
            vec![0.0, 0.001, 0.002],
            vec![0.0, std::f64::consts::FRAC_PI_2, std::f64::consts::PI],
            vec![0.0, 0.001, 0.002],
        );
        assert_eq!(
            grid_cylindrical.coord_system(),
            CoordinateSystem::Cylindrical
        );
    }

    #[test]
    fn test_line_getters() {
        let x_lines = vec![0.0, 0.001, 0.003, 0.006];
        let y_lines = vec![0.0, 0.002, 0.004];
        let z_lines = vec![0.0, 0.001, 0.002, 0.003, 0.004];

        let grid = Grid::cartesian(x_lines.clone(), y_lines.clone(), z_lines.clone());

        assert_eq!(grid.x_lines(), &x_lines[..]);
        assert_eq!(grid.y_lines(), &y_lines[..]);
        assert_eq!(grid.z_lines(), &z_lines[..]);
    }

    #[test]
    fn test_cell_size_at_position() {
        // Non-uniform grid
        let x_lines = vec![0.0, 0.001, 0.003, 0.007]; // dx = [0.001, 0.002, 0.004]
        let y_lines = vec![0.0, 0.002, 0.005]; // dy = [0.002, 0.003]
        let z_lines = vec![0.0, 0.001, 0.002, 0.004]; // dz = [0.001, 0.001, 0.002]

        let grid = Grid::cartesian(x_lines, y_lines, z_lines);

        // Cell at (0, 0, 0)
        let (dx, dy, dz) = grid.cell_size_at(0, 0, 0);
        assert!((dx - 0.001).abs() < 1e-10);
        assert!((dy - 0.002).abs() < 1e-10);
        assert!((dz - 0.001).abs() < 1e-10);

        // Cell at (1, 1, 2)
        let (dx, dy, dz) = grid.cell_size_at(1, 1, 2);
        assert!((dx - 0.002).abs() < 1e-10);
        assert!((dy - 0.003).abs() < 1e-10);
        assert!((dz - 0.002).abs() < 1e-10);
    }

    #[test]
    fn test_cell_center() {
        let grid = Grid::uniform(10, 10, 10, 0.001);

        // Cell (0, 0, 0) center should be at (0.0005, 0.0005, 0.0005)
        let (cx, cy, cz) = grid.cell_center(0, 0, 0);
        assert!((cx - 0.0005).abs() < 1e-10);
        assert!((cy - 0.0005).abs() < 1e-10);
        assert!((cz - 0.0005).abs() < 1e-10);

        // Cell (5, 5, 5) center should be at (0.0055, 0.0055, 0.0055)
        let (cx, cy, cz) = grid.cell_center(5, 5, 5);
        assert!((cx - 0.0055).abs() < 1e-10);
        assert!((cy - 0.0055).abs() < 1e-10);
        assert!((cz - 0.0055).abs() < 1e-10);
    }

    #[test]
    fn test_grid_volume() {
        let grid = Grid::uniform(10, 10, 10, 0.001);

        // Volume = 0.01 * 0.01 * 0.01 = 1e-6 m^3
        let volume = grid.volume();
        assert!((volume - 1e-6).abs() < 1e-15);
    }

    #[test]
    fn test_grid_volume_with_unit() {
        let mut grid = Grid::uniform(10, 10, 10, 10.0); // 10mm cells (in mm)

        // Without scaling: 100mm x 100mm x 100mm = 1e6 mm^3
        let volume_mm = grid.volume();
        assert!((volume_mm - 1e6).abs() < 1e-6);

        // Set unit to meters (1mm = 0.001m)
        grid.set_unit(0.001);

        // With scaling: 0.1m x 0.1m x 0.1m = 0.001 m^3
        let volume_m = grid.volume();
        assert!((volume_m - 0.001).abs() < 1e-10);
    }

    #[test]
    fn test_grid_bounds() {
        let grid = Grid::uniform(10, 20, 30, 0.001);

        let ((x_min, x_max), (y_min, y_max), (z_min, z_max)) = grid.bounds();

        assert!((x_min - 0.0).abs() < 1e-10);
        assert!((x_max - 0.01).abs() < 1e-10);

        assert!((y_min - 0.0).abs() < 1e-10);
        assert!((y_max - 0.02).abs() < 1e-10);

        assert!((z_min - 0.0).abs() < 1e-10);
        assert!((z_max - 0.03).abs() < 1e-10);
    }

    #[test]
    fn test_add_x_lines() {
        let mut grid = Grid::uniform(10, 10, 10, 0.001);

        // Add some new X lines
        grid.add_x_lines(&[0.0005, 0.0015, 0.0025]);

        // Check that lines were added and sorted
        let x_lines = grid.x_lines();
        assert!(x_lines.contains(&0.0005));
        assert!(x_lines.contains(&0.0015));
        assert!(x_lines.contains(&0.0025));

        // Check sorting
        for i in 0..x_lines.len() - 1 {
            assert!(x_lines[i] < x_lines[i + 1]);
        }
    }

    #[test]
    fn test_add_y_lines() {
        let mut grid = Grid::uniform(10, 10, 10, 0.001);

        grid.add_y_lines(&[0.0005, 0.0035]);

        let y_lines = grid.y_lines();
        assert!(y_lines.contains(&0.0005));
        assert!(y_lines.contains(&0.0035));

        // Verify sorted
        for i in 0..y_lines.len() - 1 {
            assert!(y_lines[i] < y_lines[i + 1]);
        }
    }

    #[test]
    fn test_add_z_lines() {
        let mut grid = Grid::uniform(10, 10, 10, 0.001);

        grid.add_z_lines(&[0.0005, 0.0045, 0.0085]);

        let z_lines = grid.z_lines();
        assert!(z_lines.contains(&0.0005));
        assert!(z_lines.contains(&0.0045));
        assert!(z_lines.contains(&0.0085));

        // Verify sorted
        for i in 0..z_lines.len() - 1 {
            assert!(z_lines[i] < z_lines[i + 1]);
        }
    }

    #[test]
    fn test_add_lines_dedup() {
        let mut grid = Grid::uniform(10, 10, 10, 0.001);

        let original_count = grid.x_lines().len();

        // Add duplicate lines (0.001 already exists)
        grid.add_x_lines(&[0.001, 0.002, 0.002]);

        // Should not add duplicates
        assert_eq!(grid.x_lines().len(), original_count);
    }

    #[test]
    fn test_smooth_mesh() {
        // Create a grid with large ratio between adjacent cells
        let x_lines = vec![0.0, 0.001, 0.01, 0.02]; // Large jump from 0.001 to 0.01
        let y_lines = vec![0.0, 0.001, 0.002];
        let z_lines = vec![0.0, 0.001, 0.002];

        let mut grid = Grid::cartesian(x_lines, y_lines, z_lines);

        // Smooth with max ratio of 2.0
        grid.smooth_mesh(2.0);

        // Check that new lines were added to reduce ratio
        let x_lines = grid.x_lines();
        assert!(x_lines.len() > 4); // More lines than original

        // Verify max ratio constraint
        for i in 1..x_lines.len() - 1 {
            let delta_prev = x_lines[i] - x_lines[i - 1];
            let delta_curr = x_lines[i + 1] - x_lines[i];

            if delta_prev > 0.0 && delta_curr > 0.0 {
                let ratio = delta_curr / delta_prev;
                // Allow small tolerance
                assert!(ratio <= 2.1, "Ratio {} exceeds max at index {}", ratio, i);
            }
        }
    }

    #[test]
    fn test_smooth_mesh_small_grid() {
        // Grid with less than 3 lines should not change
        let x_lines = vec![0.0, 1.0];
        let y_lines = vec![0.0, 1.0];
        let z_lines = vec![0.0, 1.0];

        let mut grid = Grid::cartesian(x_lines.clone(), y_lines.clone(), z_lines.clone());
        grid.smooth_mesh(2.0);

        assert_eq!(grid.x_lines().len(), 2);
        assert_eq!(grid.y_lines().len(), 2);
        assert_eq!(grid.z_lines().len(), 2);
    }

    #[test]
    fn test_delta_x_edge_cases() {
        let grid = Grid::uniform(10, 10, 10, 0.001);

        // Normal case
        let dx = grid.delta_x(5);
        assert!((dx - 0.001).abs() < 1e-10);

        // Last index - should use previous delta
        let dx_last = grid.delta_x(10);
        assert!((dx_last - 0.001).abs() < 1e-10);

        // Index 0
        let dx_first = grid.delta_x(0);
        assert!((dx_first - 0.001).abs() < 1e-10);
    }

    #[test]
    fn test_delta_y_edge_cases() {
        let grid = Grid::uniform(10, 10, 10, 0.001);

        // Normal case
        let dy = grid.delta_y(5);
        assert!((dy - 0.001).abs() < 1e-10);

        // Last index
        let dy_last = grid.delta_y(10);
        assert!((dy_last - 0.001).abs() < 1e-10);
    }

    #[test]
    fn test_delta_z_edge_cases() {
        let grid = Grid::uniform(10, 10, 10, 0.001);

        // Normal case
        let dz = grid.delta_z(5);
        assert!((dz - 0.001).abs() < 1e-10);

        // Last index
        let dz_last = grid.delta_z(10);
        assert!((dz_last - 0.001).abs() < 1e-10);
    }

    #[test]
    fn test_delta_single_line_grid() {
        // Edge case: grid with only one line
        let grid = Grid::cartesian(vec![0.0], vec![0.0], vec![0.0]);

        // Should return unit (default 1.0)
        let dx = grid.delta_x(0);
        assert!((dx - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_x_line_out_of_bounds() {
        let grid = Grid::uniform(10, 10, 10, 0.001);

        // Valid index
        let x = grid.x_line(5);
        assert!((x - 0.005).abs() < 1e-10);

        // Out of bounds - should return 0.0
        let x_oob = grid.x_line(100);
        assert!((x_oob - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_y_line_out_of_bounds() {
        let grid = Grid::uniform(10, 10, 10, 0.001);

        // Valid index
        let y = grid.y_line(5);
        assert!((y - 0.005).abs() < 1e-10);

        // Out of bounds
        let y_oob = grid.y_line(100);
        assert!((y_oob - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_z_line_out_of_bounds() {
        let grid = Grid::uniform(10, 10, 10, 0.001);

        // Valid index
        let z = grid.z_line(5);
        assert!((z - 0.005).abs() < 1e-10);

        // Out of bounds
        let z_oob = grid.z_line(100);
        assert!((z_oob - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_find_cell_out_of_range() {
        let grid = Grid::uniform(10, 10, 10, 0.001);

        // Before grid start
        let before = grid.find_cell(-0.001, 0.005, 0.005);
        assert_eq!(before, None);

        // After grid end
        let after_x = grid.find_cell(0.015, 0.005, 0.005);
        assert_eq!(after_x, None);

        let after_y = grid.find_cell(0.005, 0.015, 0.005);
        assert_eq!(after_y, None);

        let after_z = grid.find_cell(0.005, 0.005, 0.015);
        assert_eq!(after_z, None);
    }

    #[test]
    fn test_find_cell_at_boundary() {
        let grid = Grid::uniform(10, 10, 10, 0.001);

        // At grid origin
        let origin = grid.find_cell(0.0, 0.0, 0.0);
        assert_eq!(origin, Some((0, 0, 0)));

        // At grid end (exactly)
        let end = grid.find_cell(0.01, 0.01, 0.01);
        assert_eq!(end, Some((9, 9, 9)));
    }

    #[test]
    fn test_find_cell_x() {
        let grid = Grid::uniform(10, 10, 10, 0.001);

        assert_eq!(grid.find_cell_x(0.0), 0);
        assert_eq!(grid.find_cell_x(0.0055), 5);
        assert_eq!(grid.find_cell_x(0.0095), 9);

        // Beyond grid should clamp to last cell
        let clamped = grid.find_cell_x(0.02);
        assert_eq!(clamped, 9);
    }

    #[test]
    fn test_find_cell_y() {
        let grid = Grid::uniform(10, 10, 10, 0.001);

        assert_eq!(grid.find_cell_y(0.0), 0);
        assert_eq!(grid.find_cell_y(0.0055), 5);

        // Beyond grid
        let clamped = grid.find_cell_y(0.02);
        assert_eq!(clamped, 9);
    }

    #[test]
    fn test_find_cell_z() {
        let grid = Grid::uniform(10, 10, 10, 0.001);

        assert_eq!(grid.find_cell_z(0.0), 0);
        assert_eq!(grid.find_cell_z(0.0055), 5);

        // Beyond grid
        let clamped = grid.find_cell_z(0.02);
        assert_eq!(clamped, 9);
    }

    #[test]
    fn test_dimensions_with_single_lines() {
        // Edge case: minimal grid
        let grid = Grid::cartesian(vec![0.0, 1.0], vec![0.0, 1.0], vec![0.0, 1.0]);

        let dims = grid.dimensions();
        assert_eq!(dims.nx, 1);
        assert_eq!(dims.ny, 1);
        assert_eq!(dims.nz, 1);
    }

    #[test]
    fn test_new_constructor() {
        let grid = Grid::new(
            CoordinateSystem::Cartesian,
            vec![0.0, 0.1, 0.2],
            vec![0.0, 0.1, 0.2, 0.3],
            vec![0.0, 0.1],
        );

        assert_eq!(grid.coord_system(), CoordinateSystem::Cartesian);
        assert_eq!(grid.dimensions().nx, 2);
        assert_eq!(grid.dimensions().ny, 3);
        assert_eq!(grid.dimensions().nz, 1);
    }

    #[test]
    fn test_min_delta_function() {
        // Via cell_size which uses min_delta
        let x_lines = vec![0.0, 0.001, 0.003, 0.01]; // min delta = 0.001
        let y_lines = vec![0.0, 0.002, 0.003]; // min delta = 0.001
        let z_lines = vec![0.0, 0.005, 0.006]; // min delta = 0.001

        let grid = Grid::cartesian(x_lines, y_lines, z_lines);
        let (dx, dy, dz) = grid.cell_size();

        assert!((dx - 0.001).abs() < 1e-10);
        assert!((dy - 0.001).abs() < 1e-10);
        assert!((dz - 0.001).abs() < 1e-10);
    }

    #[test]
    fn test_coordinate_system_default() {
        let cs = CoordinateSystem::default();
        assert_eq!(cs, CoordinateSystem::Cartesian);
    }

    #[test]
    fn test_bounds_with_unit() {
        let mut grid = Grid::uniform(10, 10, 10, 10.0); // 10 unit cells

        // Set unit scaling
        grid.set_unit(0.001); // Convert to mm

        let ((x_min, x_max), (y_min, y_max), (z_min, z_max)) = grid.bounds();

        // Bounds should be scaled
        assert!((x_min - 0.0).abs() < 1e-10);
        assert!((x_max - 0.1).abs() < 1e-10); // 100 * 0.001 = 0.1
        assert!((y_min - 0.0).abs() < 1e-10);
        assert!((y_max - 0.1).abs() < 1e-10);
        assert!((z_min - 0.0).abs() < 1e-10);
        assert!((z_max - 0.1).abs() < 1e-10);
    }

    #[test]
    fn test_cell_center_with_unit() {
        let mut grid = Grid::uniform(10, 10, 10, 10.0); // 10 unit cells
        grid.set_unit(0.001); // mm to m

        let (cx, cy, cz) = grid.cell_center(0, 0, 0);

        // Center of first cell: 0.5 * 10 * 0.001 = 0.005
        assert!((cx - 0.005).abs() < 1e-10);
        assert!((cy - 0.005).abs() < 1e-10);
        assert!((cz - 0.005).abs() < 1e-10);
    }
}
