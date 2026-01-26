//! Engine Interface abstraction for field access.
//!
//! Provides a high-level interface to access simulation fields with
//! interpolation support.

use crate::arrays::{Dimensions, VectorField3D};
use crate::constants::{EPS0, MU0};
use crate::geometry::Grid;

/// Interpolation type for field access.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum InterpolationType {
    /// No interpolation (nearest neighbor)
    #[default]
    NoInterpolation,
    /// Interpolate to node position
    NodeInterpolate,
    /// Interpolate to cell center
    CellInterpolate,
}

/// Engine interface for accessing simulation fields.
///
/// This provides a unified interface to access E and H fields,
/// calculate derived quantities like J-field and energy,
/// and perform interpolation for field values at arbitrary positions.
pub struct EngineInterface<'a> {
    /// Reference to E-field
    e_field: &'a VectorField3D,
    /// Reference to H-field
    h_field: &'a VectorField3D,
    /// Grid information
    grid: &'a Grid,
    /// Operator coefficients (optional, for derived fields)
    op_eps: Option<&'a VectorField3D>,
    op_mu: Option<&'a VectorField3D>,
    /// Current timestep
    timestep: u64,
    /// Current time
    time: f64,
    /// Timestep size
    dt: f64,
    /// Interpolation type
    interp_type: InterpolationType,
}

impl<'a> EngineInterface<'a> {
    /// Create a new engine interface.
    pub fn new(
        e_field: &'a VectorField3D,
        h_field: &'a VectorField3D,
        grid: &'a Grid,
        dt: f64,
    ) -> Self {
        Self {
            e_field,
            h_field,
            grid,
            op_eps: None,
            op_mu: None,
            timestep: 0,
            time: 0.0,
            dt,
            interp_type: InterpolationType::default(),
        }
    }

    /// Set operator epsilon coefficients for D-field calculation.
    pub fn with_epsilon(mut self, eps: &'a VectorField3D) -> Self {
        self.op_eps = Some(eps);
        self
    }

    /// Set operator mu coefficients for B-field calculation.
    pub fn with_mu(mut self, mu: &'a VectorField3D) -> Self {
        self.op_mu = Some(mu);
        self
    }

    /// Set interpolation type.
    pub fn with_interpolation(mut self, interp: InterpolationType) -> Self {
        self.interp_type = interp;
        self
    }

    /// Update time information.
    pub fn set_time(&mut self, timestep: u64, time: f64) {
        self.timestep = timestep;
        self.time = time;
    }

    /// Get current timestep.
    pub fn timestep(&self) -> u64 {
        self.timestep
    }

    /// Get current simulation time.
    pub fn time(&self) -> f64 {
        self.time
    }

    /// Get timestep size.
    pub fn dt(&self) -> f64 {
        self.dt
    }

    /// Get grid dimensions.
    pub fn dims(&self) -> Dimensions {
        self.e_field.dims()
    }

    /// Get E-field at grid position (no interpolation).
    pub fn get_e_field_raw(&self, i: usize, j: usize, k: usize) -> [f64; 3] {
        [
            self.e_field.x.get(i, j, k) as f64,
            self.e_field.y.get(i, j, k) as f64,
            self.e_field.z.get(i, j, k) as f64,
        ]
    }

    /// Get H-field at grid position (no interpolation).
    pub fn get_h_field_raw(&self, i: usize, j: usize, k: usize) -> [f64; 3] {
        [
            self.h_field.x.get(i, j, k) as f64,
            self.h_field.y.get(i, j, k) as f64,
            self.h_field.z.get(i, j, k) as f64,
        ]
    }

    /// Get E-field at physical position with interpolation.
    pub fn get_e_field(&self, pos: [f64; 3]) -> [f64; 3] {
        match self.interp_type {
            InterpolationType::NoInterpolation => {
                let (i, j, k) = self.snap_to_grid(pos);
                self.get_e_field_raw(i, j, k)
            }
            InterpolationType::NodeInterpolate => {
                self.interpolate_e_field_to_node(pos)
            }
            InterpolationType::CellInterpolate => {
                self.interpolate_e_field_to_cell(pos)
            }
        }
    }

    /// Get H-field at physical position with interpolation.
    pub fn get_h_field(&self, pos: [f64; 3]) -> [f64; 3] {
        match self.interp_type {
            InterpolationType::NoInterpolation => {
                let (i, j, k) = self.snap_to_grid(pos);
                self.get_h_field_raw(i, j, k)
            }
            InterpolationType::NodeInterpolate => {
                self.interpolate_h_field_to_node(pos)
            }
            InterpolationType::CellInterpolate => {
                self.interpolate_h_field_to_cell(pos)
            }
        }
    }

    /// Calculate J-field (current density) from E-field.
    /// J = sigma * E (for conductive materials)
    /// Returns curl(H) for general case.
    pub fn get_j_field(&self, pos: [f64; 3]) -> [f64; 3] {
        // For simplicity, return curl(H) which equals J + dD/dt
        let (i, j, k) = self.snap_to_grid(pos);
        self.calc_curl_h(i, j, k)
    }

    /// Calculate D-field (electric flux density) from E-field.
    /// D = epsilon * E
    pub fn get_d_field(&self, pos: [f64; 3]) -> [f64; 3] {
        let e = self.get_e_field(pos);

        if let Some(eps) = self.op_eps {
            let (i, j, k) = self.snap_to_grid(pos);
            [
                e[0] * eps.x.get(i, j, k) as f64,
                e[1] * eps.y.get(i, j, k) as f64,
                e[2] * eps.z.get(i, j, k) as f64,
            ]
        } else {
            // Free space
            [e[0] * EPS0, e[1] * EPS0, e[2] * EPS0]
        }
    }

    /// Calculate B-field (magnetic flux density) from H-field.
    /// B = mu * H
    pub fn get_b_field(&self, pos: [f64; 3]) -> [f64; 3] {
        let h = self.get_h_field(pos);

        if let Some(mu) = self.op_mu {
            let (i, j, k) = self.snap_to_grid(pos);
            [
                h[0] * mu.x.get(i, j, k) as f64,
                h[1] * mu.y.get(i, j, k) as f64,
                h[2] * mu.z.get(i, j, k) as f64,
            ]
        } else {
            // Free space
            [h[0] * MU0, h[1] * MU0, h[2] * MU0]
        }
    }

    /// Calculate total current density from curl(H).
    /// rot(H) = J + dD/dt
    pub fn get_rot_h_field(&self, pos: [f64; 3]) -> [f64; 3] {
        let (i, j, k) = self.snap_to_grid(pos);
        self.calc_curl_h(i, j, k)
    }

    /// Calculate voltage integral along a line.
    /// V = -integral(E . dl)
    pub fn calc_voltage_integral(&self, start: [f64; 3], stop: [f64; 3]) -> f64 {
        let dims = self.dims();
        let mut voltage = 0.0;

        // Integrate along x, then y, then z (axis-aligned path)
        let (i_start, j_start, k_start) = self.snap_to_grid(start);
        let (i_stop, j_stop, k_stop) = self.snap_to_grid(stop);

        // X-direction
        let i_min = i_start.min(i_stop);
        let i_max = i_start.max(i_stop);
        for i in i_min..i_max {
            let i_clamped = i.min(dims.nx - 1);
            let e_x = self.e_field.x.get(i_clamped, j_start, k_start) as f64;
            let dx = self.grid.delta_x(i_clamped);
            voltage -= e_x * dx;
        }

        // Y-direction
        let j_min = j_start.min(j_stop);
        let j_max = j_start.max(j_stop);
        for j in j_min..j_max {
            let j_clamped = j.min(dims.ny - 1);
            let e_y = self.e_field.y.get(i_stop, j_clamped, k_start) as f64;
            let dy = self.grid.delta_y(j_clamped);
            voltage -= e_y * dy;
        }

        // Z-direction
        let k_min = k_start.min(k_stop);
        let k_max = k_start.max(k_stop);
        for k in k_min..k_max {
            let k_clamped = k.min(dims.nz - 1);
            let e_z = self.e_field.z.get(i_stop, j_stop, k_clamped) as f64;
            let dz = self.grid.delta_z(k_clamped);
            voltage -= e_z * dz;
        }

        voltage
    }

    /// Calculate current integral through a surface.
    /// I = integral(H . dl) around the surface
    pub fn calc_current_integral(
        &self,
        start: [f64; 3],
        stop: [f64; 3],
        normal_dir: usize,
    ) -> f64 {
        let dims = self.dims();
        let mut current = 0.0;

        let (i_start, j_start, k_start) = self.snap_to_grid(start);
        let (i_stop, j_stop, k_stop) = self.snap_to_grid(stop);

        // Integrate H around the contour based on normal direction
        match normal_dir {
            0 => {
                // X-normal surface: integrate Hy and Hz
                for j in j_start.min(j_stop)..j_start.max(j_stop) {
                    for k in k_start.min(k_stop)..k_start.max(k_stop) {
                        let j_c = j.min(dims.ny - 1);
                        let k_c = k.min(dims.nz - 1);
                        let h_y = self.h_field.y.get(i_start, j_c, k_c) as f64;
                        let h_z = self.h_field.z.get(i_start, j_c, k_c) as f64;
                        let dy = self.grid.delta_y(j_c);
                        let dz = self.grid.delta_z(k_c);
                        current += h_y * dz - h_z * dy;
                    }
                }
            }
            1 => {
                // Y-normal surface: integrate Hx and Hz
                for i in i_start.min(i_stop)..i_start.max(i_stop) {
                    for k in k_start.min(k_stop)..k_start.max(k_stop) {
                        let i_c = i.min(dims.nx - 1);
                        let k_c = k.min(dims.nz - 1);
                        let h_x = self.h_field.x.get(i_c, j_start, k_c) as f64;
                        let h_z = self.h_field.z.get(i_c, j_start, k_c) as f64;
                        let dx = self.grid.delta_x(i_c);
                        let dz = self.grid.delta_z(k_c);
                        current += h_z * dx - h_x * dz;
                    }
                }
            }
            _ => {
                // Z-normal surface: integrate Hx and Hy
                for i in i_start.min(i_stop)..i_start.max(i_stop) {
                    for j in j_start.min(j_stop)..j_start.max(j_stop) {
                        let i_c = i.min(dims.nx - 1);
                        let j_c = j.min(dims.ny - 1);
                        let h_x = self.h_field.x.get(i_c, j_c, k_start) as f64;
                        let h_y = self.h_field.y.get(i_c, j_c, k_start) as f64;
                        let dx = self.grid.delta_x(i_c);
                        let dy = self.grid.delta_y(j_c);
                        current += h_x * dy - h_y * dx;
                    }
                }
            }
        }

        current
    }

    /// Calculate total electromagnetic energy (rough estimate).
    pub fn calc_fast_energy(&self) -> f64 {
        let e_energy = self.e_field.energy();
        let h_energy = self.h_field.energy();

        // W = 0.5 * (eps * E^2 + mu * H^2)
        0.5 * (EPS0 * e_energy + MU0 * h_energy)
    }

    /// Calculate total E-field energy.
    pub fn calc_e_energy(&self) -> f64 {
        0.5 * EPS0 * self.e_field.energy()
    }

    /// Calculate total H-field energy.
    pub fn calc_h_energy(&self) -> f64 {
        0.5 * MU0 * self.h_field.energy()
    }

    /// Snap physical position to grid indices.
    fn snap_to_grid(&self, pos: [f64; 3]) -> (usize, usize, usize) {
        let dims = self.dims();
        let i = self.grid.find_cell_x(pos[0]).min(dims.nx - 1);
        let j = self.grid.find_cell_y(pos[1]).min(dims.ny - 1);
        let k = self.grid.find_cell_z(pos[2]).min(dims.nz - 1);
        (i, j, k)
    }

    /// Calculate curl of H at grid position.
    fn calc_curl_h(&self, i: usize, j: usize, k: usize) -> [f64; 3] {
        let dims = self.dims();

        // dHz/dy - dHy/dz
        let dhz_dy = if j + 1 < dims.ny {
            (self.h_field.z.get(i, j + 1, k) - self.h_field.z.get(i, j, k)) as f64
                / self.grid.delta_y(j)
        } else {
            0.0
        };
        let dhy_dz = if k + 1 < dims.nz {
            (self.h_field.y.get(i, j, k + 1) - self.h_field.y.get(i, j, k)) as f64
                / self.grid.delta_z(k)
        } else {
            0.0
        };
        let curl_x = dhz_dy - dhy_dz;

        // dHx/dz - dHz/dx
        let dhx_dz = if k + 1 < dims.nz {
            (self.h_field.x.get(i, j, k + 1) - self.h_field.x.get(i, j, k)) as f64
                / self.grid.delta_z(k)
        } else {
            0.0
        };
        let dhz_dx = if i + 1 < dims.nx {
            (self.h_field.z.get(i + 1, j, k) - self.h_field.z.get(i, j, k)) as f64
                / self.grid.delta_x(i)
        } else {
            0.0
        };
        let curl_y = dhx_dz - dhz_dx;

        // dHy/dx - dHx/dy
        let dhy_dx = if i + 1 < dims.nx {
            (self.h_field.y.get(i + 1, j, k) - self.h_field.y.get(i, j, k)) as f64
                / self.grid.delta_x(i)
        } else {
            0.0
        };
        let dhx_dy = if j + 1 < dims.ny {
            (self.h_field.x.get(i, j + 1, k) - self.h_field.x.get(i, j, k)) as f64
                / self.grid.delta_y(j)
        } else {
            0.0
        };
        let curl_z = dhy_dx - dhx_dy;

        [curl_x, curl_y, curl_z]
    }

    /// Trilinear interpolation for E-field to node position.
    fn interpolate_e_field_to_node(&self, pos: [f64; 3]) -> [f64; 3] {
        let dims = self.dims();
        let (i, j, k) = self.snap_to_grid(pos);

        // Get fractional position within cell
        let fx = (pos[0] - self.grid.x_line(i)) / self.grid.delta_x(i);
        let fy = (pos[1] - self.grid.y_line(j)) / self.grid.delta_y(j);
        let fz = (pos[2] - self.grid.z_line(k)) / self.grid.delta_z(k);

        let fx = fx.clamp(0.0, 1.0);
        let fy = fy.clamp(0.0, 1.0);
        let fz = fz.clamp(0.0, 1.0);

        // Trilinear interpolation for each component
        let mut result = [0.0; 3];
        let fields = [&self.e_field.x, &self.e_field.y, &self.e_field.z];

        for (res, field) in result.iter_mut().zip(fields.iter()) {
            let i1 = (i + 1).min(dims.nx - 1);
            let j1 = (j + 1).min(dims.ny - 1);
            let k1 = (k + 1).min(dims.nz - 1);

            let c000 = field.get(i, j, k) as f64;
            let c100 = field.get(i1, j, k) as f64;
            let c010 = field.get(i, j1, k) as f64;
            let c110 = field.get(i1, j1, k) as f64;
            let c001 = field.get(i, j, k1) as f64;
            let c101 = field.get(i1, j, k1) as f64;
            let c011 = field.get(i, j1, k1) as f64;
            let c111 = field.get(i1, j1, k1) as f64;

            // Interpolate along x
            let c00 = c000 * (1.0 - fx) + c100 * fx;
            let c01 = c001 * (1.0 - fx) + c101 * fx;
            let c10 = c010 * (1.0 - fx) + c110 * fx;
            let c11 = c011 * (1.0 - fx) + c111 * fx;

            // Interpolate along y
            let c0 = c00 * (1.0 - fy) + c10 * fy;
            let c1 = c01 * (1.0 - fy) + c11 * fy;

            // Interpolate along z
            *res = c0 * (1.0 - fz) + c1 * fz;
        }

        result
    }

    /// Trilinear interpolation for H-field to node position.
    fn interpolate_h_field_to_node(&self, pos: [f64; 3]) -> [f64; 3] {
        let dims = self.dims();
        let (i, j, k) = self.snap_to_grid(pos);

        let fx = (pos[0] - self.grid.x_line(i)) / self.grid.delta_x(i);
        let fy = (pos[1] - self.grid.y_line(j)) / self.grid.delta_y(j);
        let fz = (pos[2] - self.grid.z_line(k)) / self.grid.delta_z(k);

        let fx = fx.clamp(0.0, 1.0);
        let fy = fy.clamp(0.0, 1.0);
        let fz = fz.clamp(0.0, 1.0);

        let mut result = [0.0; 3];
        let fields = [&self.h_field.x, &self.h_field.y, &self.h_field.z];

        for (res, field) in result.iter_mut().zip(fields.iter()) {
            let i1 = (i + 1).min(dims.nx - 1);
            let j1 = (j + 1).min(dims.ny - 1);
            let k1 = (k + 1).min(dims.nz - 1);

            let c000 = field.get(i, j, k) as f64;
            let c100 = field.get(i1, j, k) as f64;
            let c010 = field.get(i, j1, k) as f64;
            let c110 = field.get(i1, j1, k) as f64;
            let c001 = field.get(i, j, k1) as f64;
            let c101 = field.get(i1, j, k1) as f64;
            let c011 = field.get(i, j1, k1) as f64;
            let c111 = field.get(i1, j1, k1) as f64;

            let c00 = c000 * (1.0 - fx) + c100 * fx;
            let c01 = c001 * (1.0 - fx) + c101 * fx;
            let c10 = c010 * (1.0 - fx) + c110 * fx;
            let c11 = c011 * (1.0 - fx) + c111 * fx;

            let c0 = c00 * (1.0 - fy) + c10 * fy;
            let c1 = c01 * (1.0 - fy) + c11 * fy;

            *res = c0 * (1.0 - fz) + c1 * fz;
        }

        result
    }

    /// Interpolation to cell center (average of 8 corners).
    fn interpolate_e_field_to_cell(&self, pos: [f64; 3]) -> [f64; 3] {
        // For cell-centered interpolation, shift position by half a cell
        let (i, j, k) = self.snap_to_grid(pos);

        let dx = self.grid.delta_x(i) * 0.5;
        let dy = self.grid.delta_y(j) * 0.5;
        let dz = self.grid.delta_z(k) * 0.5;

        let shifted_pos = [pos[0] + dx, pos[1] + dy, pos[2] + dz];
        self.interpolate_e_field_to_node(shifted_pos)
    }

    /// Interpolation to cell center for H-field.
    fn interpolate_h_field_to_cell(&self, pos: [f64; 3]) -> [f64; 3] {
        let (i, j, k) = self.snap_to_grid(pos);

        let dx = self.grid.delta_x(i) * 0.5;
        let dy = self.grid.delta_y(j) * 0.5;
        let dz = self.grid.delta_z(k) * 0.5;

        let shifted_pos = [pos[0] + dx, pos[1] + dy, pos[2] + dz];
        self.interpolate_h_field_to_node(shifted_pos)
    }
}

/// Mutable engine interface for field modifications.
#[allow(dead_code)]
pub struct EngineInterfaceMut<'a> {
    /// Mutable reference to E-field
    e_field: &'a mut VectorField3D,
    /// Mutable reference to H-field
    h_field: &'a mut VectorField3D,
    /// Grid information
    grid: &'a Grid,
    /// Current timestep
    timestep: u64,
    /// Current time
    time: f64,
    /// Timestep size
    dt: f64,
}

impl<'a> EngineInterfaceMut<'a> {
    /// Create a new mutable engine interface.
    pub fn new(
        e_field: &'a mut VectorField3D,
        h_field: &'a mut VectorField3D,
        grid: &'a Grid,
        dt: f64,
    ) -> Self {
        Self {
            e_field,
            h_field,
            grid,
            timestep: 0,
            time: 0.0,
            dt,
        }
    }

    /// Set E-field value at position.
    pub fn set_e_field(&mut self, i: usize, j: usize, k: usize, value: [f32; 3]) {
        self.e_field.x.set(i, j, k, value[0]);
        self.e_field.y.set(i, j, k, value[1]);
        self.e_field.z.set(i, j, k, value[2]);
    }

    /// Set H-field value at position.
    pub fn set_h_field(&mut self, i: usize, j: usize, k: usize, value: [f32; 3]) {
        self.h_field.x.set(i, j, k, value[0]);
        self.h_field.y.set(i, j, k, value[1]);
        self.h_field.z.set(i, j, k, value[2]);
    }

    /// Add to E-field at position.
    pub fn add_e_field(&mut self, i: usize, j: usize, k: usize, value: [f32; 3]) {
        self.e_field.x.add(i, j, k, value[0]);
        self.e_field.y.add(i, j, k, value[1]);
        self.e_field.z.add(i, j, k, value[2]);
    }

    /// Get E-field reference.
    pub fn e_field(&self) -> &VectorField3D {
        self.e_field
    }

    /// Get mutable E-field reference.
    pub fn e_field_mut(&mut self) -> &mut VectorField3D {
        self.e_field
    }

    /// Get H-field reference.
    pub fn h_field(&self) -> &VectorField3D {
        self.h_field
    }

    /// Get mutable H-field reference.
    pub fn h_field_mut(&mut self) -> &mut VectorField3D {
        self.h_field
    }

    /// Get grid reference.
    pub fn grid(&self) -> &Grid {
        self.grid
    }

    /// Update time information.
    pub fn set_time(&mut self, timestep: u64, time: f64) {
        self.timestep = timestep;
        self.time = time;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_grid() -> Grid {
        Grid::uniform(10, 10, 10, 0.001)
    }

    #[test]
    fn test_engine_interface_creation() {
        let dims = Dimensions::new(10, 10, 10);
        let e_field = VectorField3D::new(dims);
        let h_field = VectorField3D::new(dims);
        let grid = create_test_grid();

        let interface = EngineInterface::new(&e_field, &h_field, &grid, 1e-12);
        assert_eq!(interface.timestep(), 0);
        assert_eq!(interface.time(), 0.0);
    }

    #[test]
    fn test_field_access() {
        let dims = Dimensions::new(10, 10, 10);
        let mut e_field = VectorField3D::new(dims);
        let h_field = VectorField3D::new(dims);
        let grid = create_test_grid();

        e_field.x.set(5, 5, 5, 1.0);
        e_field.y.set(5, 5, 5, 2.0);
        e_field.z.set(5, 5, 5, 3.0);

        let interface = EngineInterface::new(&e_field, &h_field, &grid, 1e-12);
        let e = interface.get_e_field_raw(5, 5, 5);

        assert!((e[0] - 1.0).abs() < 1e-6);
        assert!((e[1] - 2.0).abs() < 1e-6);
        assert!((e[2] - 3.0).abs() < 1e-6);
    }

    #[test]
    fn test_energy_calculation() {
        let dims = Dimensions::new(10, 10, 10);
        let mut e_field = VectorField3D::new(dims);
        let h_field = VectorField3D::new(dims);
        let grid = create_test_grid();

        // Set some field values
        e_field.x.fill(1.0);

        let interface = EngineInterface::new(&e_field, &h_field, &grid, 1e-12);
        let energy = interface.calc_fast_energy();

        assert!(energy > 0.0);
    }

    #[test]
    fn test_interpolation_types() {
        let dims = Dimensions::new(10, 10, 10);
        let mut e_field = VectorField3D::new(dims);
        let h_field = VectorField3D::new(dims);
        let grid = create_test_grid();

        e_field.x.fill(1.0);

        let interface = EngineInterface::new(&e_field, &h_field, &grid, 1e-12)
            .with_interpolation(InterpolationType::NodeInterpolate);

        let e = interface.get_e_field([0.005, 0.005, 0.005]);
        assert!(e[0] > 0.0);
    }

    #[test]
    fn test_voltage_integral() {
        let dims = Dimensions::new(10, 10, 10);
        let mut e_field = VectorField3D::new(dims);
        let h_field = VectorField3D::new(dims);
        let grid = create_test_grid();

        // Set uniform E-field in x-direction
        e_field.x.fill(100.0); // 100 V/m

        let interface = EngineInterface::new(&e_field, &h_field, &grid, 1e-12);

        // Integrate over 5 cells in x-direction (5 * 0.001 = 0.005 m)
        // Expected voltage = -100 * 0.005 = -0.5 V
        let voltage = interface.calc_voltage_integral([0.0, 0.005, 0.005], [0.005, 0.005, 0.005]);

        // Approximate check (depends on grid discretization)
        assert!(voltage < 0.0);
    }

    #[test]
    fn test_mutable_interface() {
        let dims = Dimensions::new(10, 10, 10);
        let mut e_field = VectorField3D::new(dims);
        let mut h_field = VectorField3D::new(dims);
        let grid = create_test_grid();

        {
            let mut interface = EngineInterfaceMut::new(&mut e_field, &mut h_field, &grid, 1e-12);
            interface.set_e_field(5, 5, 5, [1.0, 2.0, 3.0]);
        }

        assert!((e_field.x.get(5, 5, 5) - 1.0).abs() < 1e-6);
        assert!((e_field.y.get(5, 5, 5) - 2.0).abs() < 1e-6);
        assert!((e_field.z.get(5, 5, 5) - 3.0).abs() < 1e-6);
    }
}
