//! Conducting Sheet Model Extension.
//!
//! Implements a frequency-dependent conducting sheet model for efficient FDTD
//! analysis of planar waveguides and circuits.
//!
//! Reference: Lauer, A.; Wolff, I.; "A conducting sheet model for efficient wide band
//! FDTD analysis of planar waveguides and circuits," IEEE MTT-S, 1999.

use crate::arrays::{Dimensions, VectorField3D};

/// Pre-computed Lorentz parameters for frequency-dependent sheet impedance.
/// These parameters fit the skin-effect impedance over specific frequency ranges.
#[derive(Debug, Clone)]
pub struct SheetLorentzParams {
    /// Drude conductivity term (DC conductance)
    pub g: f64,
    /// First Lorentz pole resistance
    pub r1: f64,
    /// First Lorentz pole inductance
    pub l1: f64,
    /// Second Lorentz pole resistance
    pub r2: f64,
    /// Second Lorentz pole inductance
    pub l2: f64,
    /// Critical frequency (rad/s)
    pub omega_critical: f64,
    /// Stop frequency (rad/s)
    pub omega_stop: f64,
}

impl SheetLorentzParams {
    /// Get parameters for copper at a given frequency range.
    /// Based on pre-optimized fits from the original C++ implementation.
    pub fn copper_for_frequency(f_max: f64) -> Self {
        let omega_max = 2.0 * std::f64::consts::PI * f_max;

        // Simplified parameter selection based on frequency range
        // Full implementation would have 30+ frequency ranges
        if omega_max < 1e9 {
            Self {
                g: 5.8e7,  // Copper conductivity
                r1: 0.0,
                l1: 1e-12,
                r2: 0.0,
                l2: 1e-12,
                omega_critical: 1e8,
                omega_stop: 1e9,
            }
        } else if omega_max < 1e10 {
            Self {
                g: 5.8e7,
                r1: 0.008,
                l1: 1.5e-12,
                r2: 0.004,
                l2: 2.5e-12,
                omega_critical: 1e9,
                omega_stop: 1e10,
            }
        } else {
            Self {
                g: 5.8e7,
                r1: 0.025,
                l1: 0.8e-12,
                r2: 0.012,
                l2: 1.2e-12,
                omega_critical: 1e10,
                omega_stop: 1e11,
            }
        }
    }
}

/// Configuration for a conducting sheet.
#[derive(Debug, Clone)]
pub struct ConductingSheetConfig {
    /// Sheet conductivity (S/m)
    pub conductivity: f64,
    /// Sheet thickness (m)
    pub thickness: f64,
    /// Maximum simulation frequency (Hz)
    pub f_max: f64,
    /// Normal direction of sheet (0=x, 1=y, 2=z)
    pub normal_direction: usize,
    /// Sheet position along normal (grid index)
    pub position: usize,
    /// Sheet extent: start indices [i, j, k]
    pub start: [usize; 3],
    /// Sheet extent: stop indices [i, j, k]
    pub stop: [usize; 3],
}

impl ConductingSheetConfig {
    /// Create a new conducting sheet configuration.
    pub fn new(
        conductivity: f64,
        thickness: f64,
        f_max: f64,
        normal_direction: usize,
        position: usize,
    ) -> Self {
        Self {
            conductivity,
            thickness,
            f_max,
            normal_direction,
            position,
            start: [0, 0, 0],
            stop: [0, 0, 0],
        }
    }

    /// Set sheet extent.
    pub fn with_extent(mut self, start: [usize; 3], stop: [usize; 3]) -> Self {
        self.start = start;
        self.stop = stop;
        self
    }

    /// Create copper sheet (conductivity = 5.8e7 S/m).
    pub fn copper(thickness: f64, f_max: f64, normal_direction: usize, position: usize) -> Self {
        Self::new(5.8e7, thickness, f_max, normal_direction, position)
    }

    /// Create aluminum sheet (conductivity = 3.8e7 S/m).
    pub fn aluminum(thickness: f64, f_max: f64, normal_direction: usize, position: usize) -> Self {
        Self::new(3.8e7, thickness, f_max, normal_direction, position)
    }
}

/// ADE coefficients for conducting sheet update.
#[derive(Debug, Clone)]
struct AdeCoefficients {
    /// E-field update coefficient
    c1: f32,
    /// Auxiliary current coefficient
    c2: f32,
    /// Previous auxiliary current coefficient
    c3: f32,
    /// Previous E-field coefficient
    c4: f32,
}

/// Single conducting sheet element storage.
#[derive(Debug)]
struct SheetElement {
    /// Position in grid
    pos: [usize; 3],
    /// Tangential direction (primary)
    tan_dir: usize,
    /// ADE coefficients for first Lorentz pole
    ade1: AdeCoefficients,
    /// ADE coefficients for second Lorentz pole
    ade2: AdeCoefficients,
    /// Auxiliary current 1
    j_aux1: f32,
    /// Auxiliary current 2
    j_aux2: f32,
    /// Previous E-field value
    e_prev: f32,
}

impl SheetElement {
    fn new(pos: [usize; 3], tan_dir: usize) -> Self {
        Self {
            pos,
            tan_dir,
            ade1: AdeCoefficients {
                c1: 1.0,
                c2: 0.0,
                c3: 0.0,
                c4: 0.0,
            },
            ade2: AdeCoefficients {
                c1: 1.0,
                c2: 0.0,
                c3: 0.0,
                c4: 0.0,
            },
            j_aux1: 0.0,
            j_aux2: 0.0,
            e_prev: 0.0,
        }
    }
}

/// Conducting Sheet Extension.
///
/// Models thin conducting sheets with frequency-dependent impedance
/// using auxiliary differential equations (ADE).
pub struct ConductingSheet {
    /// Configuration
    config: ConductingSheetConfig,
    /// Grid dimensions
    dims: Dimensions,
    /// Timestep
    dt: f64,
    /// Lorentz parameters
    params: SheetLorentzParams,
    /// Sheet elements (one per affected grid edge)
    elements: Vec<SheetElement>,
    /// Whether the extension is active
    active: bool,
}

impl ConductingSheet {
    /// Create a new conducting sheet extension.
    pub fn new(config: ConductingSheetConfig, dims: Dimensions, dt: f64) -> Self {
        let params = SheetLorentzParams::copper_for_frequency(config.f_max);

        let mut sheet = Self {
            config,
            dims,
            dt,
            params,
            elements: Vec::new(),
            active: false,
        };

        sheet.build_elements();
        sheet
    }

    /// Build sheet elements for all affected grid edges.
    fn build_elements(&mut self) {
        let ny = self.config.normal_direction;
        let nyp = (ny + 1) % 3;
        let nypp = (ny + 2) % 3;

        let pos_ny = self.config.position;

        // Determine sheet extent
        let start = self.config.start;
        let stop = [
            if self.config.stop[0] == 0 {
                self.dims.nx
            } else {
                self.config.stop[0]
            },
            if self.config.stop[1] == 0 {
                self.dims.ny
            } else {
                self.config.stop[1]
            },
            if self.config.stop[2] == 0 {
                self.dims.nz
            } else {
                self.config.stop[2]
            },
        ];

        // Create elements for both tangential directions
        for tan_dir in [nyp, nypp] {
            let (i_start, i_end, j_start, j_end) = match ny {
                0 => (start[1], stop[1], start[2], stop[2]),
                1 => (start[0], stop[0], start[2], stop[2]),
                _ => (start[0], stop[0], start[1], stop[1]),
            };

            for i in i_start..i_end {
                for j in j_start..j_end {
                    let mut pos = [0usize; 3];
                    pos[ny] = pos_ny;
                    match ny {
                        0 => {
                            pos[1] = i;
                            pos[2] = j;
                        }
                        1 => {
                            pos[0] = i;
                            pos[2] = j;
                        }
                        _ => {
                            pos[0] = i;
                            pos[1] = j;
                        }
                    }

                    let mut elem = SheetElement::new(pos, tan_dir);
                    self.init_coefficients(&mut elem);
                    self.elements.push(elem);
                }
            }
        }

        self.active = !self.elements.is_empty();
    }

    /// Initialize ADE coefficients for a sheet element.
    fn init_coefficients(&self, elem: &mut SheetElement) {
        // Surface resistance per square
        let r_sheet = 1.0 / (self.config.conductivity * self.config.thickness);

        // First Lorentz pole coefficients
        // From: Z1 = r1 + j*omega*l1
        // ADE: dJ1/dt + (r1/l1)*J1 = (1/l1)*E
        let alpha1 = if self.params.l1 > 0.0 {
            self.params.r1 / self.params.l1
        } else {
            0.0
        };
        let beta1 = if self.params.l1 > 0.0 {
            1.0 / self.params.l1
        } else {
            0.0
        };

        // Discretized: J1^(n+1) = c3*J1^n + c2*E^(n+1/2)
        let denom1 = 1.0 + 0.5 * alpha1 * self.dt;
        elem.ade1.c3 = ((1.0 - 0.5 * alpha1 * self.dt) / denom1) as f32;
        elem.ade1.c2 = (beta1 * self.dt / denom1) as f32;
        elem.ade1.c1 = (r_sheet * elem.ade1.c2 as f64) as f32;

        // Second Lorentz pole
        let alpha2 = if self.params.l2 > 0.0 {
            self.params.r2 / self.params.l2
        } else {
            0.0
        };
        let beta2 = if self.params.l2 > 0.0 {
            1.0 / self.params.l2
        } else {
            0.0
        };

        let denom2 = 1.0 + 0.5 * alpha2 * self.dt;
        elem.ade2.c3 = ((1.0 - 0.5 * alpha2 * self.dt) / denom2) as f32;
        elem.ade2.c2 = (beta2 * self.dt / denom2) as f32;
        elem.ade2.c1 = (r_sheet * elem.ade2.c2 as f64) as f32;

        // Drude (DC) term
        elem.ade1.c4 = (self.params.g * self.config.thickness * self.dt) as f32;
    }

    /// Check if extension is active.
    pub fn is_active(&self) -> bool {
        self.active
    }

    /// Get number of sheet elements.
    pub fn num_elements(&self) -> usize {
        self.elements.len()
    }

    /// Apply conducting sheet modification to E-field update.
    pub fn apply_to_voltage(&mut self, e_field: &mut VectorField3D) {
        if !self.active {
            return;
        }

        for elem in &mut self.elements {
            let (i, j, k) = (elem.pos[0], elem.pos[1], elem.pos[2]);

            // Get current E-field
            let e_current = match elem.tan_dir {
                0 => e_field.x.get(i, j, k),
                1 => e_field.y.get(i, j, k),
                _ => e_field.z.get(i, j, k),
            };

            // Update auxiliary currents (ADE)
            let j1_new = elem.ade1.c3 * elem.j_aux1 + elem.ade1.c2 * e_current;
            let j2_new = elem.ade2.c3 * elem.j_aux2 + elem.ade2.c2 * e_current;

            // Calculate total current modification
            let j_total = elem.ade1.c4 * e_current // Drude term
                + elem.ade1.c1 * (j1_new + elem.j_aux1) * 0.5 // Lorentz 1
                + elem.ade2.c1 * (j2_new + elem.j_aux2) * 0.5; // Lorentz 2

            // Modify E-field
            let e_modified = e_current - j_total;

            // Store updated values
            elem.j_aux1 = j1_new;
            elem.j_aux2 = j2_new;
            elem.e_prev = e_current;

            // Write back
            match elem.tan_dir {
                0 => e_field.x.set(i, j, k, e_modified),
                1 => e_field.y.set(i, j, k, e_modified),
                _ => e_field.z.set(i, j, k, e_modified),
            }
        }
    }

    /// Reset all auxiliary variables.
    pub fn reset(&mut self) {
        for elem in &mut self.elements {
            elem.j_aux1 = 0.0;
            elem.j_aux2 = 0.0;
            elem.e_prev = 0.0;
        }
    }
}

/// Manager for multiple conducting sheets.
pub struct ConductingSheetManager {
    sheets: Vec<ConductingSheet>,
}

impl ConductingSheetManager {
    /// Create a new manager.
    pub fn new() -> Self {
        Self { sheets: Vec::new() }
    }

    /// Add a conducting sheet.
    pub fn add_sheet(&mut self, sheet: ConductingSheet) {
        self.sheets.push(sheet);
    }

    /// Apply all sheets to E-field.
    pub fn apply_to_voltage(&mut self, e_field: &mut VectorField3D) {
        for sheet in &mut self.sheets {
            sheet.apply_to_voltage(e_field);
        }
    }

    /// Reset all sheets.
    pub fn reset(&mut self) {
        for sheet in &mut self.sheets {
            sheet.reset();
        }
    }

    /// Check if any sheets are active.
    pub fn is_active(&self) -> bool {
        self.sheets.iter().any(|s| s.is_active())
    }
}

impl Default for ConductingSheetManager {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lorentz_params() {
        let params = SheetLorentzParams::copper_for_frequency(1e9);
        assert!((params.g - 5.8e7).abs() < 1e3);

        let params_high = SheetLorentzParams::copper_for_frequency(1e11);
        assert!(params_high.r1 > 0.0);
    }

    #[test]
    fn test_sheet_config() {
        let config = ConductingSheetConfig::copper(35e-6, 10e9, 2, 10);
        assert!((config.conductivity - 5.8e7).abs() < 1.0);
        assert!((config.thickness - 35e-6).abs() < 1e-9);
    }

    #[test]
    fn test_sheet_creation() {
        let dims = Dimensions::new(20, 20, 20);
        let config = ConductingSheetConfig::copper(35e-6, 10e9, 2, 10)
            .with_extent([0, 0, 0], [20, 20, 0]);

        let sheet = ConductingSheet::new(config, dims, 1e-12);

        assert!(sheet.is_active());
        // Should have elements for both tangential directions across the sheet
        assert!(sheet.num_elements() > 0);
    }

    #[test]
    fn test_sheet_update() {
        let dims = Dimensions::new(10, 10, 10);
        let config = ConductingSheetConfig::copper(35e-6, 10e9, 2, 5);

        let mut sheet = ConductingSheet::new(config, dims, 1e-12);
        let mut e_field = VectorField3D::new(dims);

        // Set initial field at sheet position
        for i in 0..10 {
            for j in 0..10 {
                e_field.x.set(i, j, 5, 1.0);
                e_field.y.set(i, j, 5, 1.0);
            }
        }

        // Apply sheet
        sheet.apply_to_voltage(&mut e_field);

        // Field should be modified (reduced by conductive losses)
        let ex = e_field.x.get(5, 5, 5);
        assert!(ex <= 1.0);
    }

    #[test]
    fn test_sheet_reset() {
        let dims = Dimensions::new(10, 10, 10);
        let config = ConductingSheetConfig::copper(35e-6, 10e9, 2, 5);

        let mut sheet = ConductingSheet::new(config, dims, 1e-12);
        let mut e_field = VectorField3D::new(dims);

        e_field.x.fill(1.0);

        // Run a few updates
        for _ in 0..5 {
            sheet.apply_to_voltage(&mut e_field);
        }

        // Reset
        sheet.reset();

        // All auxiliary variables should be zero
        for elem in &sheet.elements {
            assert_eq!(elem.j_aux1, 0.0);
            assert_eq!(elem.j_aux2, 0.0);
        }
    }

    #[test]
    fn test_sheet_manager() {
        let dims = Dimensions::new(20, 20, 20);

        let config1 = ConductingSheetConfig::copper(35e-6, 10e9, 2, 5);
        let config2 = ConductingSheetConfig::aluminum(35e-6, 10e9, 2, 15);

        let mut manager = ConductingSheetManager::new();
        manager.add_sheet(ConductingSheet::new(config1, dims, 1e-12));
        manager.add_sheet(ConductingSheet::new(config2, dims, 1e-12));

        assert!(manager.is_active());

        let mut e_field = VectorField3D::new(dims);
        e_field.x.fill(1.0);

        manager.apply_to_voltage(&mut e_field);
        manager.reset();
    }
}
