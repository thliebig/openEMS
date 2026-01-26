//! Local Absorbing Boundary Condition Extension.
//!
//! Implements first-order Mur absorbing boundary conditions on arbitrary
//! 2D sheet primitives within the simulation domain.
//!
//! Two modes are supported:
//! - MUR_1ST: Standard first-order Mur BC
//! - MUR_1ST_SA: First-order Mur BC with Super Absorption (SIBC)
//!
//! References:
//! - G. Mur, "Absorbing Boundary Conditions for the Finite-Difference
//!   Approximation of the Time-Domain Electromagnetic-Field Equations,"
//!   IEEE Trans. EMC, vol. EMC-23, no. 4, pp. 377-382, Nov. 1981.
//! - Betz, V. T. and Mittra, R. "Absorbing boundary conditions for the
//!   finite-difference time-domain analysis of guided-wave structures."
//!   Coordinated Science Laboratory Report no. UILU-ENG-93-2243 (1993).

use crate::arrays::{Field3D, VectorField3D};
use crate::constants::C0;

/// Type of absorbing boundary condition.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum AbcType {
    /// Undefined/disabled
    #[default]
    Undefined,
    /// Mur's BC, 1st order
    Mur1st,
    /// Mur's BC, 1st order with Super Absorption
    Mur1stSa,
}

impl AbcType {
    /// Check if this type applies to currents (H-field).
    pub fn applies_to_currents(&self) -> bool {
        matches!(self, AbcType::Mur1stSa)
    }
}

/// Configuration for a local absorbing BC sheet.
#[derive(Debug, Clone)]
pub struct LocalAbsorbingBcConfig {
    /// Type of ABC
    pub abc_type: AbcType,
    /// Phase velocity (if None, uses C0)
    pub phase_velocity: Option<f64>,
    /// Sheet bounding box start indices [x, y, z]
    pub start: [usize; 3],
    /// Sheet bounding box stop indices [x, y, z]
    pub stop: [usize; 3],
    /// Whether normal direction points positive
    pub normal_sign_positive: bool,
}

impl Default for LocalAbsorbingBcConfig {
    fn default() -> Self {
        Self {
            abc_type: AbcType::Undefined,
            phase_velocity: None,
            start: [0; 3],
            stop: [0; 3],
            normal_sign_positive: true,
        }
    }
}

impl LocalAbsorbingBcConfig {
    /// Create a new configuration.
    pub fn new(abc_type: AbcType, start: [usize; 3], stop: [usize; 3]) -> Self {
        Self {
            abc_type,
            phase_velocity: None,
            start,
            stop,
            normal_sign_positive: true,
        }
    }

    /// Set phase velocity.
    pub fn with_phase_velocity(mut self, v: f64) -> Self {
        self.phase_velocity = Some(v);
        self
    }

    /// Set normal sign direction.
    pub fn with_normal_sign(mut self, positive: bool) -> Self {
        self.normal_sign_positive = positive;
        self
    }

    /// Determine the normal direction from the sheet geometry.
    /// Returns None if the geometry is not a valid sheet.
    pub fn detect_normal_direction(&self) -> Option<usize> {
        let mut normal_dir = None;
        let mut sheet_count = 0;

        for dim in 0..3 {
            if self.start[dim] == self.stop[dim] {
                normal_dir = Some(dim);
                sheet_count += 1;
            }
        }

        // Valid sheet has exactly one collapsed dimension
        if sheet_count == 1 {
            normal_dir
        } else {
            None
        }
    }
}

/// 2D array storage for coefficient and field values.
#[derive(Debug, Clone)]
struct Array2D {
    data: Vec<f32>,
    size: [usize; 2],
}

impl Array2D {
    fn new(size: [usize; 2]) -> Self {
        Self {
            data: vec![0.0; size[0] * size[1]],
            size,
        }
    }

    #[inline]
    fn idx(&self, i: usize, j: usize) -> usize {
        i * self.size[1] + j
    }

    #[inline]
    fn get(&self, i: usize, j: usize) -> f32 {
        self.data[self.idx(i, j)]
    }

    #[inline]
    fn set(&mut self, i: usize, j: usize, value: f32) {
        let idx = self.idx(i, j);
        self.data[idx] = value;
    }

    fn fill(&mut self, value: f32) {
        self.data.fill(value);
    }
}

/// A single local absorbing BC sheet.
#[derive(Debug)]
#[allow(dead_code)]
pub struct LocalAbsorbingSheet {
    /// ABC type
    abc_type: AbcType,
    /// Normal direction (0=x, 1=y, 2=z)
    ny: usize,
    /// First tangential direction (ny + 1) % 3
    ny_p: usize,
    /// Second tangential direction (ny + 2) % 3
    ny_pp: usize,
    /// Sheet position start indices
    pos_start: [usize; 3],
    /// Sheet position stop indices
    pos_stop: [usize; 3],
    /// Number of lines in each tangential direction
    num_lines: [usize; 2],
    /// Whether normal sign is positive
    normal_sign_positive: bool,
    /// Position for voltage boundary
    pos_ny0_shift_v: usize,
    /// Position for current boundary
    pos_ny0_i: usize,
    /// Shifted position for current
    pos_ny0_shift_i: usize,
    /// K1 coefficient for ny_p direction
    k1_ny_p: Array2D,
    /// K1 coefficient for ny_pp direction
    k1_ny_pp: Array2D,
    /// K2 coefficient for ny_p direction (super absorption only)
    k2_ny_p: Array2D,
    /// K2 coefficient for ny_pp direction (super absorption only)
    k2_ny_pp: Array2D,
    /// Stored voltage values for ny_p
    v_ny_p: Array2D,
    /// Stored voltage values for ny_pp
    v_ny_pp: Array2D,
    /// Stored current values for ny_p (super absorption only)
    i_ny_p: Array2D,
    /// Stored current values for ny_pp (super absorption only)
    i_ny_pp: Array2D,
    /// Starting timestep (activation delay)
    start_ts: u64,
}

impl LocalAbsorbingSheet {
    /// Create a new local absorbing sheet.
    pub fn new(
        config: &LocalAbsorbingBcConfig,
        dt: f64,
        get_delta: impl Fn(usize, [usize; 3]) -> f64,
    ) -> Option<Self> {
        let ny = config.detect_normal_direction()?;
        let ny_p = (ny + 1) % 3;
        let ny_pp = (ny + 2) % 3;

        let num_lines = [
            config.stop[ny_p] - config.start[ny_p] + 1,
            config.stop[ny_pp] - config.start[ny_pp] + 1,
        ];

        let normal_sign_positive = config.normal_sign_positive;
        let pos_start = config.start;
        let pos_stop = config.stop;

        // Calculate shifted positions
        let pos_ny0_shift_v = if normal_sign_positive {
            pos_start[ny] + 1
        } else {
            pos_start[ny].saturating_sub(1)
        };

        let pos_ny0_i = if normal_sign_positive {
            pos_start[ny]
        } else {
            pos_start[ny].saturating_sub(1)
        };

        let pos_ny0_shift_i = if normal_sign_positive {
            pos_start[ny] + 1
        } else {
            pos_start[ny].saturating_sub(2)
        };

        // Initialize coefficient arrays
        let mut k1_ny_p = Array2D::new(num_lines);
        let mut k1_ny_pp = Array2D::new(num_lines);
        let mut k2_ny_p = Array2D::new(num_lines);
        let mut k2_ny_pp = Array2D::new(num_lines);

        // Calculate phase velocity * dt
        let phase_velocity = config.phase_velocity.unwrap_or(C0);
        let vt = phase_velocity * dt;

        // Calculate coefficients for each cell
        let mut pos = [0usize; 3];
        pos[ny] = pos_start[ny];

        for i in 0..num_lines[0] {
            pos[ny_p] = pos_start[ny_p] + i;
            for j in 0..num_lines[1] {
                pos[ny_pp] = pos_start[ny_pp] + j;

                // Get cell size in normal direction
                let delta = get_delta(ny, pos);

                // K1 = (v*dt - delta) / (v*dt + delta)
                let k1 = ((vt - delta) / (vt + delta)) as f32;
                k1_ny_p.set(i, j, k1);
                k1_ny_pp.set(i, j, k1);

                // K2 = v*dt / delta (for super absorption)
                if config.abc_type == AbcType::Mur1stSa {
                    let k2 = (vt / delta) as f32;
                    k2_ny_p.set(i, j, k2);
                    k2_ny_pp.set(i, j, k2);
                }
            }
        }

        // Initialize storage arrays
        let v_ny_p = Array2D::new(num_lines);
        let v_ny_pp = Array2D::new(num_lines);
        let i_ny_p = Array2D::new(num_lines);
        let i_ny_pp = Array2D::new(num_lines);

        Some(Self {
            abc_type: config.abc_type,
            ny,
            ny_p,
            ny_pp,
            pos_start,
            pos_stop,
            num_lines,
            normal_sign_positive,
            pos_ny0_shift_v,
            pos_ny0_i,
            pos_ny0_shift_i,
            k1_ny_p,
            k1_ny_pp,
            k2_ny_p,
            k2_ny_pp,
            v_ny_p,
            v_ny_pp,
            i_ny_p,
            i_ny_pp,
            start_ts: 0,
        })
    }

    /// Set starting timestep for activation delay.
    pub fn set_start_timestep(&mut self, ts: u64) {
        self.start_ts = ts;
    }

    /// Check if this sheet is active at the given timestep.
    #[inline]
    pub fn is_active(&self, timestep: u64) -> bool {
        timestep >= self.start_ts
    }

    /// Get the number of cells in this sheet.
    pub fn num_cells(&self) -> usize {
        self.num_lines[0] * self.num_lines[1]
    }

    /// Get field component accessor.
    #[inline]
    fn get_field_component<'a>(&self, field: &'a VectorField3D, direction: usize) -> &'a Field3D {
        match direction {
            0 => &field.x,
            1 => &field.y,
            _ => &field.z,
        }
    }

    /// Get mutable field component accessor.
    #[inline]
    fn get_field_component_mut<'a>(
        &self,
        field: &'a mut VectorField3D,
        direction: usize,
    ) -> &'a mut Field3D {
        match direction {
            0 => &mut field.x,
            1 => &mut field.y,
            _ => &mut field.z,
        }
    }

    /// Pre-voltage update: Store E(i+s,n) - K1*E(i,n).
    pub fn pre_voltage_update(&mut self, e_field: &VectorField3D, timestep: u64) {
        if !self.is_active(timestep) {
            return;
        }

        let mut pos = [0usize; 3];
        let mut pos_shift = [0usize; 3];

        pos[self.ny] = self.pos_start[self.ny];
        pos_shift[self.ny] = self.pos_ny0_shift_v;

        for i in 0..self.num_lines[0] {
            pos[self.ny_p] = self.pos_start[self.ny_p] + i;
            pos_shift[self.ny_p] = pos[self.ny_p];

            for j in 0..self.num_lines[1] {
                pos[self.ny_pp] = self.pos_start[self.ny_pp] + j;
                pos_shift[self.ny_pp] = pos[self.ny_pp];

                // Get E-field values
                let field_ny_p = self.get_field_component(e_field, self.ny_p);
                let field_ny_pp = self.get_field_component(e_field, self.ny_pp);

                let e_shift_ny_p = field_ny_p.get(pos_shift[0], pos_shift[1], pos_shift[2]);
                let e_boundary_ny_p = field_ny_p.get(pos[0], pos[1], pos[2]);
                let e_shift_ny_pp = field_ny_pp.get(pos_shift[0], pos_shift[1], pos_shift[2]);
                let e_boundary_ny_pp = field_ny_pp.get(pos[0], pos[1], pos[2]);

                // E(i+s,n) - K1*E(i,n)
                self.v_ny_p.set(
                    i,
                    j,
                    e_shift_ny_p - self.k1_ny_p.get(i, j) * e_boundary_ny_p,
                );
                self.v_ny_pp.set(
                    i,
                    j,
                    e_shift_ny_pp - self.k1_ny_pp.get(i, j) * e_boundary_ny_pp,
                );
            }
        }
    }

    /// Post-voltage update: Add K1*E(i+s,n+1).
    pub fn post_voltage_update(&mut self, e_field: &VectorField3D, timestep: u64) {
        if !self.is_active(timestep) {
            return;
        }

        let mut pos_shift = [0usize; 3];
        pos_shift[self.ny] = self.pos_ny0_shift_v;

        for i in 0..self.num_lines[0] {
            pos_shift[self.ny_p] = self.pos_start[self.ny_p] + i;

            for j in 0..self.num_lines[1] {
                pos_shift[self.ny_pp] = self.pos_start[self.ny_pp] + j;

                // Get E-field at shifted position (after FDTD update)
                let field_ny_p = self.get_field_component(e_field, self.ny_p);
                let field_ny_pp = self.get_field_component(e_field, self.ny_pp);

                let e_shift_ny_p = field_ny_p.get(pos_shift[0], pos_shift[1], pos_shift[2]);
                let e_shift_ny_pp = field_ny_pp.get(pos_shift[0], pos_shift[1], pos_shift[2]);

                // E(i+s,n) - K1*E(i,n) + K1*E(i+s,n+1)
                let prev_ny_p = self.v_ny_p.get(i, j);
                let prev_ny_pp = self.v_ny_pp.get(i, j);
                self.v_ny_p
                    .set(i, j, prev_ny_p + self.k1_ny_p.get(i, j) * e_shift_ny_p);
                self.v_ny_pp
                    .set(i, j, prev_ny_pp + self.k1_ny_pp.get(i, j) * e_shift_ny_pp);
            }
        }
    }

    /// Apply ABC to voltages: Set boundary E-field.
    pub fn apply_to_voltages(&self, e_field: &mut VectorField3D, timestep: u64) {
        if !self.is_active(timestep) {
            return;
        }

        let mut pos = [0usize; 3];
        pos[self.ny] = self.pos_start[self.ny];

        for i in 0..self.num_lines[0] {
            pos[self.ny_p] = self.pos_start[self.ny_p] + i;

            for j in 0..self.num_lines[1] {
                pos[self.ny_pp] = self.pos_start[self.ny_pp] + j;

                // Set E-field at boundary
                let field_ny_p = self.get_field_component_mut(e_field, self.ny_p);
                field_ny_p.set(pos[0], pos[1], pos[2], self.v_ny_p.get(i, j));

                let field_ny_pp = self.get_field_component_mut(e_field, self.ny_pp);
                field_ny_pp.set(pos[0], pos[1], pos[2], self.v_ny_pp.get(i, j));
            }
        }
    }

    /// Pre-current update (super absorption only): Store H(i+s,n) - K1*H(i,n).
    pub fn pre_current_update(&mut self, h_field: &VectorField3D, timestep: u64) {
        if !self.is_active(timestep) || self.abc_type != AbcType::Mur1stSa {
            return;
        }

        let mut pos = [0usize; 3];
        let mut pos_shift = [0usize; 3];

        pos[self.ny] = self.pos_ny0_i;
        pos_shift[self.ny] = self.pos_ny0_shift_i;

        // For H-field, we iterate one less than E-field due to dual grid
        let num_lines_0 = self.num_lines[0].saturating_sub(1);
        let num_lines_1 = self.num_lines[1].saturating_sub(1);

        for i in 0..num_lines_0 {
            pos[self.ny_p] = self.pos_start[self.ny_p] + i;
            pos_shift[self.ny_p] = pos[self.ny_p];

            for j in 0..num_lines_1 {
                pos[self.ny_pp] = self.pos_start[self.ny_pp] + j;
                pos_shift[self.ny_pp] = pos[self.ny_pp];

                // Get H-field values
                let field_ny_p = self.get_field_component(h_field, self.ny_p);
                let field_ny_pp = self.get_field_component(h_field, self.ny_pp);

                let h_shift_ny_p = field_ny_p.get(pos_shift[0], pos_shift[1], pos_shift[2]);
                let h_boundary_ny_p = field_ny_p.get(pos[0], pos[1], pos[2]);
                let h_shift_ny_pp = field_ny_pp.get(pos_shift[0], pos_shift[1], pos_shift[2]);
                let h_boundary_ny_pp = field_ny_pp.get(pos[0], pos[1], pos[2]);

                // H(i+s,n) - K1*H(i,n)
                self.i_ny_p.set(
                    i,
                    j,
                    h_shift_ny_p - self.k1_ny_p.get(i, j) * h_boundary_ny_p,
                );
                self.i_ny_pp.set(
                    i,
                    j,
                    h_shift_ny_pp - self.k1_ny_pp.get(i, j) * h_boundary_ny_pp,
                );
            }
        }
    }

    /// Post-current update (super absorption only): Add K1*H(i+s,n+1).
    pub fn post_current_update(&mut self, h_field: &VectorField3D, timestep: u64) {
        if !self.is_active(timestep) || self.abc_type != AbcType::Mur1stSa {
            return;
        }

        let mut pos_shift = [0usize; 3];
        pos_shift[self.ny] = self.pos_ny0_shift_i;

        let num_lines_0 = self.num_lines[0].saturating_sub(1);
        let num_lines_1 = self.num_lines[1].saturating_sub(1);

        for i in 0..num_lines_0 {
            pos_shift[self.ny_p] = self.pos_start[self.ny_p] + i;

            for j in 0..num_lines_1 {
                pos_shift[self.ny_pp] = self.pos_start[self.ny_pp] + j;

                // Get H-field at shifted position (after FDTD update)
                let field_ny_p = self.get_field_component(h_field, self.ny_p);
                let field_ny_pp = self.get_field_component(h_field, self.ny_pp);

                let h_shift_ny_p = field_ny_p.get(pos_shift[0], pos_shift[1], pos_shift[2]);
                let h_shift_ny_pp = field_ny_pp.get(pos_shift[0], pos_shift[1], pos_shift[2]);

                // H(i+s,n) - K1*H(i,n) + K1*H(i+s,n+1)
                let prev_ny_p = self.i_ny_p.get(i, j);
                let prev_ny_pp = self.i_ny_pp.get(i, j);
                self.i_ny_p
                    .set(i, j, prev_ny_p + self.k1_ny_p.get(i, j) * h_shift_ny_p);
                self.i_ny_pp
                    .set(i, j, prev_ny_pp + self.k1_ny_pp.get(i, j) * h_shift_ny_pp);
            }
        }
    }

    /// Apply super absorption to currents: H(n+1) = (K2*Hsa + Hc)/(K2+1).
    pub fn apply_to_currents(&self, h_field: &mut VectorField3D, timestep: u64) {
        if !self.is_active(timestep) || self.abc_type != AbcType::Mur1stSa {
            return;
        }

        let mut pos = [0usize; 3];
        pos[self.ny] = self.pos_ny0_i;

        let num_lines_0 = self.num_lines[0].saturating_sub(1);
        let num_lines_1 = self.num_lines[1].saturating_sub(1);

        for i in 0..num_lines_0 {
            pos[self.ny_p] = self.pos_start[self.ny_p] + i;

            for j in 0..num_lines_1 {
                pos[self.ny_pp] = self.pos_start[self.ny_pp] + j;

                // Get current H-field values (after FDTD update)
                let k2_ny_p = self.k2_ny_p.get(i, j);
                let k2_ny_pp = self.k2_ny_pp.get(i, j);
                let h_sa_ny_p = self.i_ny_p.get(i, j);
                let h_sa_ny_pp = self.i_ny_pp.get(i, j);

                // Apply super absorption: H = (K2*Hsa + Hc)/(K2+1)
                let field_ny_p = self.get_field_component_mut(h_field, self.ny_p);
                let h_c_ny_p = field_ny_p.get(pos[0], pos[1], pos[2]);
                let h_new_ny_p = (k2_ny_p * h_sa_ny_p + h_c_ny_p) / (k2_ny_p + 1.0);
                field_ny_p.set(pos[0], pos[1], pos[2], h_new_ny_p);

                let field_ny_pp = self.get_field_component_mut(h_field, self.ny_pp);
                let h_c_ny_pp = field_ny_pp.get(pos[0], pos[1], pos[2]);
                let h_new_ny_pp = (k2_ny_pp * h_sa_ny_pp + h_c_ny_pp) / (k2_ny_pp + 1.0);
                field_ny_pp.set(pos[0], pos[1], pos[2], h_new_ny_pp);
            }
        }
    }

    /// Reset all stored values.
    pub fn reset(&mut self) {
        self.v_ny_p.fill(0.0);
        self.v_ny_pp.fill(0.0);
        self.i_ny_p.fill(0.0);
        self.i_ny_pp.fill(0.0);
    }
}

/// Manager for multiple local absorbing BC sheets.
#[derive(Debug)]
pub struct LocalAbsorbingBc {
    /// Collection of sheets
    sheets: Vec<LocalAbsorbingSheet>,
    /// Total cell count
    total_cells: usize,
}

impl LocalAbsorbingBc {
    /// Create a new empty local absorbing BC manager.
    pub fn new() -> Self {
        Self {
            sheets: Vec::new(),
            total_cells: 0,
        }
    }

    /// Add a sheet to the manager.
    pub fn add_sheet(&mut self, sheet: LocalAbsorbingSheet) {
        self.total_cells += sheet.num_cells();
        self.sheets.push(sheet);
    }

    /// Get the number of sheets.
    pub fn num_sheets(&self) -> usize {
        self.sheets.len()
    }

    /// Get the total number of cells across all sheets.
    pub fn total_cells(&self) -> usize {
        self.total_cells
    }

    /// Check if any sheets are registered.
    pub fn is_active(&self) -> bool {
        !self.sheets.is_empty()
    }

    /// Pre-voltage update for all sheets.
    pub fn pre_voltage_update(&mut self, e_field: &VectorField3D, timestep: u64) {
        for sheet in &mut self.sheets {
            sheet.pre_voltage_update(e_field, timestep);
        }
    }

    /// Post-voltage update for all sheets.
    pub fn post_voltage_update(&mut self, e_field: &VectorField3D, timestep: u64) {
        for sheet in &mut self.sheets {
            sheet.post_voltage_update(e_field, timestep);
        }
    }

    /// Apply voltage ABC to all sheets.
    pub fn apply_to_voltages(&self, e_field: &mut VectorField3D, timestep: u64) {
        for sheet in &self.sheets {
            sheet.apply_to_voltages(e_field, timestep);
        }
    }

    /// Pre-current update for all sheets (super absorption only).
    pub fn pre_current_update(&mut self, h_field: &VectorField3D, timestep: u64) {
        for sheet in &mut self.sheets {
            sheet.pre_current_update(h_field, timestep);
        }
    }

    /// Post-current update for all sheets (super absorption only).
    pub fn post_current_update(&mut self, h_field: &VectorField3D, timestep: u64) {
        for sheet in &mut self.sheets {
            sheet.post_current_update(h_field, timestep);
        }
    }

    /// Apply current ABC to all sheets (super absorption only).
    pub fn apply_to_currents(&self, h_field: &mut VectorField3D, timestep: u64) {
        for sheet in &self.sheets {
            sheet.apply_to_currents(h_field, timestep);
        }
    }

    /// Reset all sheets.
    pub fn reset(&mut self) {
        for sheet in &mut self.sheets {
            sheet.reset();
        }
    }

    /// Get sheet statistics.
    pub fn show_stats(&self) -> String {
        format!(
            "Local Absorbing BC: {} sheets, {} total cells",
            self.sheets.len(),
            self.total_cells
        )
    }
}

impl Default for LocalAbsorbingBc {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arrays::Dimensions;

    #[test]
    fn test_abc_type() {
        assert!(!AbcType::Mur1st.applies_to_currents());
        assert!(AbcType::Mur1stSa.applies_to_currents());
        assert!(!AbcType::Undefined.applies_to_currents());
    }

    #[test]
    fn test_abc_type_default() {
        let abc_type = AbcType::default();
        assert_eq!(abc_type, AbcType::Undefined);
    }

    #[test]
    fn test_abc_type_debug() {
        let abc_type = AbcType::Mur1st;
        let debug_str = format!("{:?}", abc_type);
        assert!(debug_str.contains("Mur1st"));
    }

    #[test]
    fn test_abc_type_clone() {
        let abc_type = AbcType::Mur1stSa;
        let cloned = abc_type;
        assert_eq!(cloned, AbcType::Mur1stSa);
    }

    #[test]
    fn test_config_normal_detection() {
        // X-normal sheet (y-z plane)
        let config = LocalAbsorbingBcConfig::new(AbcType::Mur1st, [5, 0, 0], [5, 10, 10]);
        assert_eq!(config.detect_normal_direction(), Some(0));

        // Y-normal sheet (x-z plane)
        let config = LocalAbsorbingBcConfig::new(AbcType::Mur1st, [0, 5, 0], [10, 5, 10]);
        assert_eq!(config.detect_normal_direction(), Some(1));

        // Z-normal sheet (x-y plane)
        let config = LocalAbsorbingBcConfig::new(AbcType::Mur1st, [0, 0, 5], [10, 10, 5]);
        assert_eq!(config.detect_normal_direction(), Some(2));

        // Not a sheet (3D box)
        let config = LocalAbsorbingBcConfig::new(AbcType::Mur1st, [0, 0, 0], [10, 10, 10]);
        assert_eq!(config.detect_normal_direction(), None);

        // Line (two collapsed dimensions)
        let config = LocalAbsorbingBcConfig::new(AbcType::Mur1st, [5, 5, 0], [5, 5, 10]);
        assert_eq!(config.detect_normal_direction(), None);
    }

    #[test]
    fn test_config_default() {
        let config = LocalAbsorbingBcConfig::default();
        assert_eq!(config.abc_type, AbcType::Undefined);
        assert_eq!(config.phase_velocity, None);
        assert_eq!(config.start, [0, 0, 0]);
        assert_eq!(config.stop, [0, 0, 0]);
        assert!(config.normal_sign_positive);
    }

    #[test]
    fn test_config_with_phase_velocity() {
        let config = LocalAbsorbingBcConfig::new(AbcType::Mur1st, [5, 0, 0], [5, 10, 10])
            .with_phase_velocity(2e8);

        assert_eq!(config.phase_velocity, Some(2e8));
    }

    #[test]
    fn test_config_with_normal_sign() {
        let config = LocalAbsorbingBcConfig::new(AbcType::Mur1st, [5, 0, 0], [5, 10, 10])
            .with_normal_sign(false);

        assert!(!config.normal_sign_positive);
    }

    #[test]
    fn test_config_builder_chain() {
        let config = LocalAbsorbingBcConfig::new(AbcType::Mur1stSa, [0, 5, 0], [10, 5, 10])
            .with_phase_velocity(1.5e8)
            .with_normal_sign(false);

        assert_eq!(config.abc_type, AbcType::Mur1stSa);
        assert_eq!(config.phase_velocity, Some(1.5e8));
        assert!(!config.normal_sign_positive);
        assert_eq!(config.detect_normal_direction(), Some(1));
    }

    #[test]
    fn test_config_clone() {
        let config = LocalAbsorbingBcConfig::new(AbcType::Mur1st, [5, 0, 0], [5, 10, 10])
            .with_phase_velocity(2e8);

        let cloned = config.clone();
        assert_eq!(cloned.abc_type, config.abc_type);
        assert_eq!(cloned.phase_velocity, config.phase_velocity);
        assert_eq!(cloned.start, config.start);
        assert_eq!(cloned.stop, config.stop);
    }

    #[test]
    fn test_sheet_creation() {
        let config = LocalAbsorbingBcConfig::new(AbcType::Mur1st, [5, 0, 0], [5, 10, 10]);

        let dt = 1e-12;
        let get_delta = |_ny: usize, _pos: [usize; 3]| 0.001;

        let sheet = LocalAbsorbingSheet::new(&config, dt, get_delta);
        assert!(sheet.is_some());

        let sheet = sheet.unwrap();
        assert_eq!(sheet.ny, 0);
        assert_eq!(sheet.ny_p, 1);
        assert_eq!(sheet.ny_pp, 2);
        assert_eq!(sheet.num_lines, [11, 11]);
    }

    #[test]
    fn test_sheet_creation_y_normal() {
        let config = LocalAbsorbingBcConfig::new(AbcType::Mur1st, [0, 5, 0], [10, 5, 10]);

        let dt = 1e-12;
        let get_delta = |_ny: usize, _pos: [usize; 3]| 0.001;

        let sheet = LocalAbsorbingSheet::new(&config, dt, get_delta).unwrap();
        assert_eq!(sheet.ny, 1);
        assert_eq!(sheet.ny_p, 2);
        assert_eq!(sheet.ny_pp, 0);
    }

    #[test]
    fn test_sheet_creation_z_normal() {
        let config = LocalAbsorbingBcConfig::new(AbcType::Mur1st, [0, 0, 5], [10, 10, 5]);

        let dt = 1e-12;
        let get_delta = |_ny: usize, _pos: [usize; 3]| 0.001;

        let sheet = LocalAbsorbingSheet::new(&config, dt, get_delta).unwrap();
        assert_eq!(sheet.ny, 2);
        assert_eq!(sheet.ny_p, 0);
        assert_eq!(sheet.ny_pp, 1);
    }

    #[test]
    fn test_sheet_creation_invalid_geometry() {
        // 3D box - not a valid sheet
        let config = LocalAbsorbingBcConfig::new(AbcType::Mur1st, [0, 0, 0], [10, 10, 10]);

        let dt = 1e-12;
        let get_delta = |_ny: usize, _pos: [usize; 3]| 0.001;

        let sheet = LocalAbsorbingSheet::new(&config, dt, get_delta);
        assert!(sheet.is_none());
    }

    #[test]
    fn test_sheet_creation_negative_normal() {
        let config = LocalAbsorbingBcConfig::new(AbcType::Mur1st, [5, 0, 0], [5, 10, 10])
            .with_normal_sign(false);

        let dt = 1e-12;
        let get_delta = |_ny: usize, _pos: [usize; 3]| 0.001;

        let sheet = LocalAbsorbingSheet::new(&config, dt, get_delta).unwrap();
        assert!(!sheet.normal_sign_positive);
    }

    #[test]
    fn test_sheet_super_absorption() {
        let config = LocalAbsorbingBcConfig::new(AbcType::Mur1stSa, [5, 0, 0], [5, 10, 10]);

        let dt = 1e-12;
        let get_delta = |_ny: usize, _pos: [usize; 3]| 0.001;

        let sheet = LocalAbsorbingSheet::new(&config, dt, get_delta);
        assert!(sheet.is_some());

        let sheet = sheet.unwrap();
        assert_eq!(sheet.abc_type, AbcType::Mur1stSa);
    }

    #[test]
    fn test_sheet_num_cells() {
        let config = LocalAbsorbingBcConfig::new(AbcType::Mur1st, [5, 0, 0], [5, 9, 14]);

        let dt = 1e-12;
        let get_delta = |_ny: usize, _pos: [usize; 3]| 0.001;

        let sheet = LocalAbsorbingSheet::new(&config, dt, get_delta).unwrap();
        // 10 x 15 = 150 cells
        assert_eq!(sheet.num_cells(), 150);
    }

    #[test]
    fn test_manager() {
        let mut manager = LocalAbsorbingBc::new();
        assert!(!manager.is_active());
        assert_eq!(manager.num_sheets(), 0);

        let config = LocalAbsorbingBcConfig::new(AbcType::Mur1st, [5, 0, 0], [5, 9, 9]);

        let dt = 1e-12;
        let get_delta = |_ny: usize, _pos: [usize; 3]| 0.001;

        if let Some(sheet) = LocalAbsorbingSheet::new(&config, dt, get_delta) {
            manager.add_sheet(sheet);
        }

        assert!(manager.is_active());
        assert_eq!(manager.num_sheets(), 1);
        assert_eq!(manager.total_cells(), 100); // 10 x 10
    }

    #[test]
    fn test_manager_default() {
        let manager = LocalAbsorbingBc::default();
        assert!(!manager.is_active());
        assert_eq!(manager.num_sheets(), 0);
        assert_eq!(manager.total_cells(), 0);
    }

    #[test]
    fn test_manager_multiple_sheets() {
        let mut manager = LocalAbsorbingBc::new();

        let dt = 1e-12;
        let get_delta = |_ny: usize, _pos: [usize; 3]| 0.001;

        // Add first sheet (X-normal)
        let config1 = LocalAbsorbingBcConfig::new(AbcType::Mur1st, [5, 0, 0], [5, 9, 9]);
        if let Some(sheet) = LocalAbsorbingSheet::new(&config1, dt, get_delta) {
            manager.add_sheet(sheet);
        }

        // Add second sheet (Y-normal)
        let config2 = LocalAbsorbingBcConfig::new(AbcType::Mur1st, [0, 5, 0], [9, 5, 9]);
        if let Some(sheet) = LocalAbsorbingSheet::new(&config2, dt, get_delta) {
            manager.add_sheet(sheet);
        }

        assert_eq!(manager.num_sheets(), 2);
        assert_eq!(manager.total_cells(), 200); // 100 + 100
    }

    #[test]
    fn test_manager_show_stats() {
        let mut manager = LocalAbsorbingBc::new();

        let dt = 1e-12;
        let get_delta = |_ny: usize, _pos: [usize; 3]| 0.001;

        let config = LocalAbsorbingBcConfig::new(AbcType::Mur1st, [5, 0, 0], [5, 9, 9]);
        if let Some(sheet) = LocalAbsorbingSheet::new(&config, dt, get_delta) {
            manager.add_sheet(sheet);
        }

        let stats = manager.show_stats();
        assert!(stats.contains("1 sheets"));
        assert!(stats.contains("100 total cells"));
    }

    #[test]
    fn test_voltage_update_cycle() {
        let dims = Dimensions {
            nx: 20,
            ny: 20,
            nz: 20,
        };
        let mut e_field = VectorField3D::new(dims);

        // Set some initial values
        e_field.y.set(5, 5, 5, 1.0);
        e_field.z.set(5, 5, 5, 1.0);

        let config = LocalAbsorbingBcConfig::new(AbcType::Mur1st, [5, 0, 0], [5, 19, 19]);

        let dt = 1e-12;
        let get_delta = |_ny: usize, _pos: [usize; 3]| 0.001;

        let mut sheet = LocalAbsorbingSheet::new(&config, dt, get_delta).unwrap();

        // Run voltage update cycle
        sheet.pre_voltage_update(&e_field, 0);
        // (FDTD update would go here)
        sheet.post_voltage_update(&e_field, 0);
        sheet.apply_to_voltages(&mut e_field, 0);

        // Boundary values should be modified
        // The exact values depend on the ABC algorithm
    }

    #[test]
    fn test_current_update_cycle_super_absorption() {
        let dims = Dimensions {
            nx: 20,
            ny: 20,
            nz: 20,
        };
        let mut h_field = VectorField3D::new(dims);

        // Set some initial values
        h_field.y.set(5, 5, 5, 0.1);
        h_field.z.set(5, 5, 5, 0.1);

        let config = LocalAbsorbingBcConfig::new(
            AbcType::Mur1stSa, // Super absorption for current updates
            [5, 0, 0],
            [5, 19, 19],
        );

        let dt = 1e-12;
        let get_delta = |_ny: usize, _pos: [usize; 3]| 0.001;

        let mut sheet = LocalAbsorbingSheet::new(&config, dt, get_delta).unwrap();

        // Run current update cycle (only applies to Mur1stSa)
        sheet.pre_current_update(&h_field, 0);
        // (FDTD update would go here)
        sheet.post_current_update(&h_field, 0);
        sheet.apply_to_currents(&mut h_field, 0);

        // For super absorption, H-field should be modified
    }

    #[test]
    fn test_current_update_skipped_for_mur1st() {
        let dims = Dimensions {
            nx: 20,
            ny: 20,
            nz: 20,
        };
        let mut h_field = VectorField3D::new(dims);

        h_field.y.fill(1.0);
        let original_value = h_field.y.get(5, 5, 5);

        let config = LocalAbsorbingBcConfig::new(
            AbcType::Mur1st, // Not super absorption
            [5, 0, 0],
            [5, 19, 19],
        );

        let dt = 1e-12;
        let get_delta = |_ny: usize, _pos: [usize; 3]| 0.001;

        let mut sheet = LocalAbsorbingSheet::new(&config, dt, get_delta).unwrap();

        // Current updates should be skipped for Mur1st
        sheet.pre_current_update(&h_field, 0);
        sheet.post_current_update(&h_field, 0);
        sheet.apply_to_currents(&mut h_field, 0);

        // H-field should be unchanged
        assert!((h_field.y.get(5, 5, 5) - original_value).abs() < 1e-10);
    }

    #[test]
    fn test_manager_voltage_update_cycle() {
        let dims = Dimensions {
            nx: 20,
            ny: 20,
            nz: 20,
        };
        let mut e_field = VectorField3D::new(dims);
        e_field.y.fill(1.0);

        let mut manager = LocalAbsorbingBc::new();

        let dt = 1e-12;
        let get_delta = |_ny: usize, _pos: [usize; 3]| 0.001;

        let config = LocalAbsorbingBcConfig::new(AbcType::Mur1st, [5, 0, 0], [5, 19, 19]);
        if let Some(sheet) = LocalAbsorbingSheet::new(&config, dt, get_delta) {
            manager.add_sheet(sheet);
        }

        // Run full voltage update cycle through manager
        manager.pre_voltage_update(&e_field, 0);
        manager.post_voltage_update(&e_field, 0);
        manager.apply_to_voltages(&mut e_field, 0);
    }

    #[test]
    fn test_manager_current_update_cycle() {
        let dims = Dimensions {
            nx: 20,
            ny: 20,
            nz: 20,
        };
        let mut h_field = VectorField3D::new(dims);
        h_field.y.fill(0.1);

        let mut manager = LocalAbsorbingBc::new();

        let dt = 1e-12;
        let get_delta = |_ny: usize, _pos: [usize; 3]| 0.001;

        let config = LocalAbsorbingBcConfig::new(AbcType::Mur1stSa, [5, 0, 0], [5, 19, 19]);
        if let Some(sheet) = LocalAbsorbingSheet::new(&config, dt, get_delta) {
            manager.add_sheet(sheet);
        }

        // Run full current update cycle through manager
        manager.pre_current_update(&h_field, 0);
        manager.post_current_update(&h_field, 0);
        manager.apply_to_currents(&mut h_field, 0);
    }

    #[test]
    fn test_activation_delay() {
        let config = LocalAbsorbingBcConfig::new(AbcType::Mur1st, [5, 0, 0], [5, 9, 9]);

        let dt = 1e-12;
        let get_delta = |_ny: usize, _pos: [usize; 3]| 0.001;

        let mut sheet = LocalAbsorbingSheet::new(&config, dt, get_delta).unwrap();
        sheet.set_start_timestep(100);

        assert!(!sheet.is_active(0));
        assert!(!sheet.is_active(99));
        assert!(sheet.is_active(100));
        assert!(sheet.is_active(1000));
    }

    #[test]
    fn test_inactive_sheet_skips_updates() {
        let dims = Dimensions {
            nx: 20,
            ny: 20,
            nz: 20,
        };
        let mut e_field = VectorField3D::new(dims);
        e_field.y.fill(1.0);
        let original_value = e_field.y.get(5, 5, 5);

        let config = LocalAbsorbingBcConfig::new(AbcType::Mur1st, [5, 0, 0], [5, 19, 19]);

        let dt = 1e-12;
        let get_delta = |_ny: usize, _pos: [usize; 3]| 0.001;

        let mut sheet = LocalAbsorbingSheet::new(&config, dt, get_delta).unwrap();
        sheet.set_start_timestep(100); // Activate at timestep 100

        // Run updates at timestep 0 (inactive)
        sheet.pre_voltage_update(&e_field, 0);
        sheet.post_voltage_update(&e_field, 0);
        sheet.apply_to_voltages(&mut e_field, 0);

        // E-field should be unchanged since sheet is inactive
        assert!((e_field.y.get(5, 5, 5) - original_value).abs() < 1e-10);
    }

    #[test]
    fn test_coefficient_calculation() {
        // Test that K1 is calculated correctly
        let c = C0;
        let dt = 1e-12;
        let delta = 0.001;
        let vt = c * dt;

        let expected_k1 = (vt - delta) / (vt + delta);
        let expected_k2 = vt / delta;

        // K1 should be negative (vt < delta typically at these parameters)
        assert!(expected_k1 < 0.0);

        // K2 should be positive
        assert!(expected_k2 > 0.0);
    }

    #[test]
    fn test_coefficient_with_custom_phase_velocity() {
        // With slower phase velocity
        let v = C0 / 2.0;
        let dt = 1e-12;
        let delta = 0.001;
        let vt = v * dt;

        let k1_slow = (vt - delta) / (vt + delta);

        // K1 should be more negative with slower phase velocity
        let vt_c0 = C0 * dt;
        let k1_c0 = (vt_c0 - delta) / (vt_c0 + delta);

        assert!(k1_slow < k1_c0);
    }

    #[test]
    fn test_reset() {
        let config = LocalAbsorbingBcConfig::new(AbcType::Mur1stSa, [5, 0, 0], [5, 9, 9]);

        let dt = 1e-12;
        let get_delta = |_ny: usize, _pos: [usize; 3]| 0.001;

        let mut manager = LocalAbsorbingBc::new();
        if let Some(sheet) = LocalAbsorbingSheet::new(&config, dt, get_delta) {
            manager.add_sheet(sheet);
        }

        // Reset should complete without panic
        manager.reset();
    }

    #[test]
    fn test_array2d_operations() {
        let mut arr = Array2D::new([5, 10]);

        // Test initial values (should be zero)
        assert!((arr.get(0, 0)).abs() < 1e-10);
        assert!((arr.get(4, 9)).abs() < 1e-10);

        // Test set and get
        arr.set(2, 3, 1.5);
        assert!((arr.get(2, 3) - 1.5).abs() < 1e-6);

        // Test fill
        arr.fill(2.0);
        assert!((arr.get(0, 0) - 2.0).abs() < 1e-6);
        assert!((arr.get(4, 9) - 2.0).abs() < 1e-6);
    }

    #[test]
    fn test_array2d_index_calculation() {
        let arr = Array2D::new([3, 4]);

        // idx(i, j) = i * size[1] + j
        assert_eq!(arr.idx(0, 0), 0);
        assert_eq!(arr.idx(0, 1), 1);
        assert_eq!(arr.idx(1, 0), 4);
        assert_eq!(arr.idx(2, 3), 11);
    }

    #[test]
    fn test_sheet_with_variable_delta() {
        let config = LocalAbsorbingBcConfig::new(AbcType::Mur1st, [5, 0, 0], [5, 4, 4]);

        let dt = 1e-12;
        // Variable cell size based on position
        let get_delta = |_ny: usize, pos: [usize; 3]| 0.001 * (1.0 + 0.1 * pos[1] as f64);

        let sheet = LocalAbsorbingSheet::new(&config, dt, get_delta);
        assert!(sheet.is_some());
    }

    #[test]
    fn test_point_sheet() {
        // Single point (degenerate case)
        let config = LocalAbsorbingBcConfig::new(AbcType::Mur1st, [5, 5, 5], [5, 5, 5]);

        let dt = 1e-12;
        let get_delta = |_ny: usize, _pos: [usize; 3]| 0.001;

        // This should fail - all three dimensions collapsed (point, not sheet)
        let sheet = LocalAbsorbingSheet::new(&config, dt, get_delta);
        assert!(sheet.is_none());
    }
}
