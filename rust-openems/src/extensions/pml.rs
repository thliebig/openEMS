//! Perfectly Matched Layer (PML) absorbing boundary condition.
//!
//! Implements the Uniaxial PML (UPML) for absorbing outgoing waves
//! at the boundaries of the simulation domain.
//!
//! Based on: Allen Taflove, "Computational Electrodynamics - The FDTD Method",
//! Third Edition, Chapter 7.8, pages 297-300.

use crate::arrays::{Dimensions, Field3D, VectorField3D};
use crate::constants::{EPS0, Z0};

/// PML configuration parameters.
#[derive(Debug, Clone)]
pub struct PmlConfig {
    /// Number of PML layers
    pub layers: usize,
    /// Polynomial grading order (typically 3-4)
    pub grading_order: f64,
    /// Maximum conductivity (automatically calculated if None)
    pub sigma_max: Option<f64>,
    /// Target reflection coefficient (used to calculate sigma_max)
    pub reflection_coeff: f64,
}

impl Default for PmlConfig {
    fn default() -> Self {
        Self {
            layers: 8,
            grading_order: 3.0,
            sigma_max: None,
            reflection_coeff: 1e-6, // -120 dB reflection
        }
    }
}

impl PmlConfig {
    /// Create PML with specified number of layers.
    pub fn new(layers: usize) -> Self {
        Self {
            layers,
            ..Default::default()
        }
    }

    /// Create PML with custom parameters.
    pub fn with_params(layers: usize, grading_order: f64, reflection_coeff: f64) -> Self {
        Self {
            layers,
            grading_order,
            sigma_max: None,
            reflection_coeff,
        }
    }

    /// Calculate optimal sigma_max for given parameters.
    pub fn calculate_sigma_max(&self, delta: f64) -> f64 {
        if let Some(sigma) = self.sigma_max {
            return sigma;
        }

        // sigma_max = -(m+1) * ln(R) / (2 * eta * d)
        // where m = grading_order, R = reflection_coeff, d = PML thickness
        let d = self.layers as f64 * delta;
        -(self.grading_order + 1.0) * self.reflection_coeff.ln() / (2.0 * Z0 * d)
    }
}

/// Boundary identifier for PML regions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PmlBoundary {
    /// X minimum boundary (-x direction)
    XMin = 0,
    /// X maximum boundary (+x direction)
    XMax = 1,
    /// Y minimum boundary (-y direction)
    YMin = 2,
    /// Y maximum boundary (+y direction)
    YMax = 3,
    /// Z minimum boundary (-z direction)
    ZMin = 4,
    /// Z maximum boundary (+z direction)
    ZMax = 5,
}

/// UPML boundary configuration for all 6 faces.
#[derive(Debug, Clone)]
pub struct PmlBoundaries {
    /// PML layers for each boundary (0 = disabled)
    pub layers: [usize; 6],
    /// Grading order for each boundary
    pub grading_order: f64,
    /// Target reflection coefficient
    pub reflection_coeff: f64,
}

impl Default for PmlBoundaries {
    fn default() -> Self {
        Self {
            layers: [8; 6],
            grading_order: 3.0,
            reflection_coeff: 1e-6,
        }
    }
}

impl PmlBoundaries {
    /// Create with uniform PML on all boundaries.
    pub fn uniform(layers: usize) -> Self {
        Self {
            layers: [layers; 6],
            ..Default::default()
        }
    }

    /// Disable PML on specific boundary.
    pub fn disable(&mut self, boundary: PmlBoundary) -> &mut Self {
        self.layers[boundary as usize] = 0;
        self
    }

    /// Set layers for specific boundary.
    pub fn set_layers(&mut self, boundary: PmlBoundary, layers: usize) -> &mut Self {
        self.layers[boundary as usize] = layers;
        self
    }

    /// Check if any PML is enabled.
    pub fn any_enabled(&self) -> bool {
        self.layers.iter().any(|&l| l > 0)
    }
}

/// UPML operator coefficients for a single PML region.
#[derive(Debug)]
struct PmlRegion {
    /// Start position in global grid
    start: [usize; 3],
    /// Number of cells in each direction
    size: [usize; 3],

    /// Coefficients for voltage update (vv * V_old - vvfo * flux_old)
    vv: [Field3D; 3],
    vvfo: [Field3D; 3],
    vvfn: [Field3D; 3],

    /// Coefficients for current update (ii * I_old - iifo * flux_old)
    ii: [Field3D; 3],
    iifo: [Field3D; 3],
    iifn: [Field3D; 3],

    /// Flux storage
    volt_flux: [Field3D; 3],
    curr_flux: [Field3D; 3],
}

impl PmlRegion {
    fn new(start: [usize; 3], size: [usize; 3]) -> Self {
        let dims = Dimensions {
            nx: size[0],
            ny: size[1],
            nz: size[2],
        };

        let create_field = || Field3D::new(dims);
        let create_fields = || [create_field(), create_field(), create_field()];

        Self {
            start,
            size,
            vv: create_fields(),
            vvfo: create_fields(),
            vvfn: create_fields(),
            ii: create_fields(),
            iifo: create_fields(),
            iifn: create_fields(),
            volt_flux: create_fields(),
            curr_flux: create_fields(),
        }
    }

    fn dims(&self) -> Dimensions {
        Dimensions {
            nx: self.size[0],
            ny: self.size[1],
            nz: self.size[2],
        }
    }
}

/// UPML absorbing boundary condition.
///
/// Implements the Uniaxial PML using flux-based auxiliary differential
/// equations as described in Taflove's FDTD textbook.
pub struct Upml {
    /// Global grid dimensions
    dims: Dimensions,
    /// Cell sizes
    delta: [f64; 3],
    /// Timestep
    dt: f64,
    /// PML regions (one per enabled boundary)
    regions: Vec<PmlRegion>,
    /// Boundary configuration
    boundaries: PmlBoundaries,
}

impl Upml {
    /// Create a new UPML with the given configuration.
    pub fn new(
        dims: Dimensions,
        delta: [f64; 3],
        dt: f64,
        boundaries: PmlBoundaries,
    ) -> Self {
        let mut upml = Self {
            dims,
            delta,
            dt,
            regions: Vec::new(),
            boundaries,
        };
        upml.build_regions();
        upml
    }

    /// Build PML regions for all enabled boundaries.
    fn build_regions(&mut self) {
        let layers = &self.boundaries.layers;
        let nx = self.dims.nx;
        let ny = self.dims.ny;
        let nz = self.dims.nz;

        // X-min boundary (full Y-Z plane)
        if layers[0] > 0 && layers[0] < nx {
            let start = [0, 0, 0];
            let size = [layers[0] + 1, ny, nz];
            let mut region = PmlRegion::new(start, size);
            self.init_region_coefficients(&mut region, 0, true);
            self.regions.push(region);
        }

        // X-max boundary (full Y-Z plane)
        if layers[1] > 0 && layers[1] < nx {
            let start = [nx - layers[1] - 1, 0, 0];
            let size = [layers[1] + 1, ny, nz];
            let mut region = PmlRegion::new(start, size);
            self.init_region_coefficients(&mut region, 0, false);
            self.regions.push(region);
        }

        // Y-min boundary (excluding X-PML corners)
        let x_start = if layers[0] > 0 { layers[0] + 1 } else { 0 };
        let x_end = if layers[1] > 0 { nx - layers[1] - 1 } else { nx };

        if layers[2] > 0 && layers[2] < ny && x_start < x_end {
            let start = [x_start, 0, 0];
            let size = [x_end - x_start, layers[2] + 1, nz];
            let mut region = PmlRegion::new(start, size);
            self.init_region_coefficients(&mut region, 1, true);
            self.regions.push(region);
        }

        // Y-max boundary (excluding X-PML corners)
        if layers[3] > 0 && layers[3] < ny && x_start < x_end {
            let start = [x_start, ny - layers[3] - 1, 0];
            let size = [x_end - x_start, layers[3] + 1, nz];
            let mut region = PmlRegion::new(start, size);
            self.init_region_coefficients(&mut region, 1, false);
            self.regions.push(region);
        }

        // Y range for Z-PML (excluding Y-PML edges)
        let y_start = if layers[2] > 0 { layers[2] + 1 } else { 0 };
        let y_end = if layers[3] > 0 { ny - layers[3] - 1 } else { ny };

        // Z-min boundary (excluding X and Y PML corners/edges)
        if layers[4] > 0 && layers[4] < nz && x_start < x_end && y_start < y_end {
            let start = [x_start, y_start, 0];
            let size = [x_end - x_start, y_end - y_start, layers[4] + 1];
            let mut region = PmlRegion::new(start, size);
            self.init_region_coefficients(&mut region, 2, true);
            self.regions.push(region);
        }

        // Z-max boundary (excluding X and Y PML corners/edges)
        if layers[5] > 0 && layers[5] < nz && x_start < x_end && y_start < y_end {
            let start = [x_start, y_start, nz - layers[5] - 1];
            let size = [x_end - x_start, y_end - y_start, layers[5] + 1];
            let mut region = PmlRegion::new(start, size);
            self.init_region_coefficients(&mut region, 2, false);
            self.regions.push(region);
        }
    }

    /// Initialize coefficients for a PML region.
    fn init_region_coefficients(&self, region: &mut PmlRegion, direction: usize, is_min: bool) {
        let dt = self.dt;
        let m = self.boundaries.grading_order;
        let layers = if is_min {
            self.boundaries.layers[direction * 2]
        } else {
            self.boundaries.layers[direction * 2 + 1]
        };

        let delta = self.delta[direction];
        let width = layers as f64 * delta;

        // Calculate sigma_max using polynomial grading
        let sigma_max = -(m + 1.0) * self.boundaries.reflection_coeff.ln() / (2.0 * Z0 * width);

        let _dims = region.dims();

        for li in 0..region.size[0] {
            for lj in 0..region.size[1] {
                for lk in 0..region.size[2] {
                    let gi = li + region.start[0];
                    let gj = lj + region.start[1];
                    let gk = lk + region.start[2];

                    let pos = [gi, gj, gk];

                    for n in 0..3 {
                        // Calculate kappa (conductivity) for voltage and current
                        let (_kappa_v, _kappa_i) = self.calculate_kappa(
                            n,
                            direction,
                            pos,
                            is_min,
                            layers,
                            width,
                            sigma_max,
                            m,
                        );

                        let _np = (n + 1) % 3;
                        let npp = (n + 2) % 3;

                        // Calculate kappa for perpendicular directions
                        let (kappa_v_arr, kappa_i_arr) = self.calculate_all_kappa(
                            pos, is_min, direction, layers, width, sigma_max, m,
                        );

                        // Check if we need PML here
                        let kappa_sum_v: f64 = kappa_v_arr.iter().sum();
                        let kappa_sum_i: f64 = kappa_i_arr.iter().sum();

                        if kappa_sum_v > 0.0 {
                            // Voltage coefficients (eq. 7.88 from Taflove)
                            // vv:   (2*eps0 - kappa_npp*dt) / (2*eps0 + kappa_npp*dt)
                            // vvfn: (2*eps0 + kappa_n*dt) / (2*eps0 + kappa_npp*dt) / eps_r
                            // vvfo: (2*eps0 - kappa_n*dt) / (2*eps0 + kappa_npp*dt) / eps_r
                            let denom = 2.0 * EPS0 + kappa_v_arr[npp] * dt;
                            region.vv[n].set(li, lj, lk,
                                ((2.0 * EPS0 - kappa_v_arr[npp] * dt) / denom) as f32);
                            region.vvfn[n].set(li, lj, lk,
                                ((2.0 * EPS0 + kappa_v_arr[n] * dt) / denom) as f32);
                            region.vvfo[n].set(li, lj, lk,
                                ((2.0 * EPS0 - kappa_v_arr[n] * dt) / denom) as f32);
                        } else {
                            // No PML absorption, passthrough
                            region.vv[n].set(li, lj, lk, 1.0);
                            region.vvfn[n].set(li, lj, lk, 1.0);
                            region.vvfo[n].set(li, lj, lk, 0.0);
                        }

                        if kappa_sum_i > 0.0 {
                            // Current coefficients (eq. 7.90 from Taflove)
                            let denom = 2.0 * EPS0 + kappa_i_arr[npp] * dt;
                            region.ii[n].set(li, lj, lk,
                                ((2.0 * EPS0 - kappa_i_arr[npp] * dt) / denom) as f32);
                            region.iifn[n].set(li, lj, lk,
                                ((2.0 * EPS0 + kappa_i_arr[n] * dt) / denom) as f32);
                            region.iifo[n].set(li, lj, lk,
                                ((2.0 * EPS0 - kappa_i_arr[n] * dt) / denom) as f32);
                        } else {
                            region.ii[n].set(li, lj, lk, 1.0);
                            region.iifn[n].set(li, lj, lk, 1.0);
                            region.iifo[n].set(li, lj, lk, 0.0);
                        }
                    }
                }
            }
        }
    }

    /// Calculate conductivity (kappa) at a given position.
    #[allow(clippy::too_many_arguments)]
    fn calculate_kappa(
        &self,
        component: usize,
        direction: usize,
        pos: [usize; 3],
        is_min: bool,
        layers: usize,
        width: f64,
        sigma_max: f64,
        m: f64,
    ) -> (f64, f64) {
        let delta = self.delta[direction];

        // Calculate depth into PML
        let depth = if is_min {
            // Distance from outer boundary
            let boundary_pos = layers;
            if pos[direction] <= boundary_pos {
                (boundary_pos - pos[direction]) as f64 * delta
            } else {
                0.0
            }
        } else {
            // Distance from outer boundary (high side)
            let boundary_start = if direction == 0 {
                self.dims.nx - layers - 1
            } else if direction == 1 {
                self.dims.ny - layers - 1
            } else {
                self.dims.nz - layers - 1
            };
            if pos[direction] >= boundary_start {
                (pos[direction] - boundary_start) as f64 * delta
            } else {
                0.0
            }
        };

        if depth <= 0.0 {
            return (0.0, 0.0);
        }

        // Polynomial grading: sigma(x) = sigma_max * (x/d)^m
        let normalized_depth = depth / width;
        let kappa_v = sigma_max * normalized_depth.powf(m);

        // Current grading is offset by half a cell
        let depth_i = depth + if component == direction { -0.5 } else { 0.5 } * delta;
        let normalized_depth_i = (depth_i / width).clamp(0.0, 1.0);
        let kappa_i = sigma_max * normalized_depth_i.powf(m);

        (kappa_v, kappa_i)
    }

    /// Calculate all three kappa components for a position.
    #[allow(clippy::too_many_arguments)]
    fn calculate_all_kappa(
        &self,
        pos: [usize; 3],
        is_min: bool,
        pml_direction: usize,
        layers: usize,
        width: f64,
        sigma_max: f64,
        m: f64,
    ) -> ([f64; 3], [f64; 3]) {
        let mut kappa_v = [0.0; 3];
        let mut kappa_i = [0.0; 3];

        for n in 0..3 {
            let (kv, ki) = self.calculate_kappa(
                n, pml_direction, pos, is_min, layers, width, sigma_max, m,
            );
            kappa_v[n] = kv;
            kappa_i[n] = ki;
        }

        (kappa_v, kappa_i)
    }

    /// Get the number of PML regions.
    pub fn num_regions(&self) -> usize {
        self.regions.len()
    }

    /// Check if any PML is active.
    pub fn is_active(&self) -> bool {
        !self.regions.is_empty()
    }

    /// Pre-voltage update: Swap voltage with flux and compute intermediate.
    ///
    /// This implements the first part of eq. 7.88 from Taflove.
    pub fn pre_voltage_update(&mut self, e_field: &mut VectorField3D) {
        for region in &mut self.regions {
            for li in 0..region.size[0] {
                for lj in 0..region.size[1] {
                    for lk in 0..region.size[2] {
                        let gi = li + region.start[0];
                        let gj = lj + region.start[1];
                        let gk = lk + region.start[2];

                        // Check bounds
                        if gi >= self.dims.nx || gj >= self.dims.ny || gk >= self.dims.nz {
                            continue;
                        }

                        for n in 0..3 {
                            let field = match n {
                                0 => &mut e_field.x,
                                1 => &mut e_field.y,
                                _ => &mut e_field.z,
                            };

                            let v_old = field.get(gi, gj, gk);
                            let flux_old = region.volt_flux[n].get(li, lj, lk);

                            let vv = region.vv[n].get(li, lj, lk);
                            let vvfo = region.vvfo[n].get(li, lj, lk);

                            // f_help = vv * V_old - vvfo * flux_old
                            let f_help = vv * v_old - vvfo * flux_old;

                            // Store flux_old in field temporarily
                            field.set(gi, gj, gk, flux_old);

                            // Store f_help for post-update
                            region.volt_flux[n].set(li, lj, lk, f_help);
                        }
                    }
                }
            }
        }
    }

    /// Post-voltage update: Complete the voltage calculation.
    ///
    /// This implements the second part of eq. 7.88 from Taflove.
    pub fn post_voltage_update(&mut self, e_field: &mut VectorField3D) {
        for region in &mut self.regions {
            for li in 0..region.size[0] {
                for lj in 0..region.size[1] {
                    for lk in 0..region.size[2] {
                        let gi = li + region.start[0];
                        let gj = lj + region.start[1];
                        let gk = lk + region.start[2];

                        if gi >= self.dims.nx || gj >= self.dims.ny || gk >= self.dims.nz {
                            continue;
                        }

                        for n in 0..3 {
                            let field = match n {
                                0 => &mut e_field.x,
                                1 => &mut e_field.y,
                                _ => &mut e_field.z,
                            };

                            // f_help was stored in volt_flux
                            let f_help = region.volt_flux[n].get(li, lj, lk);

                            // New flux is in the field (after main engine update)
                            let flux_new = field.get(gi, gj, gk);

                            // Store new flux
                            region.volt_flux[n].set(li, lj, lk, flux_new);

                            // Calculate new voltage: V_new = f_help + vvfn * flux_new
                            let vvfn = region.vvfn[n].get(li, lj, lk);
                            field.set(gi, gj, gk, f_help + vvfn * flux_new);
                        }
                    }
                }
            }
        }
    }

    /// Pre-current update: Swap current with flux and compute intermediate.
    ///
    /// This implements the first part of eq. 7.90 from Taflove.
    pub fn pre_current_update(&mut self, h_field: &mut VectorField3D) {
        for region in &mut self.regions {
            for li in 0..region.size[0] {
                for lj in 0..region.size[1] {
                    for lk in 0..region.size[2] {
                        let gi = li + region.start[0];
                        let gj = lj + region.start[1];
                        let gk = lk + region.start[2];

                        if gi >= self.dims.nx || gj >= self.dims.ny || gk >= self.dims.nz {
                            continue;
                        }

                        for n in 0..3 {
                            let field = match n {
                                0 => &mut h_field.x,
                                1 => &mut h_field.y,
                                _ => &mut h_field.z,
                            };

                            let i_old = field.get(gi, gj, gk);
                            let flux_old = region.curr_flux[n].get(li, lj, lk);

                            let ii = region.ii[n].get(li, lj, lk);
                            let iifo = region.iifo[n].get(li, lj, lk);

                            let f_help = ii * i_old - iifo * flux_old;

                            field.set(gi, gj, gk, flux_old);
                            region.curr_flux[n].set(li, lj, lk, f_help);
                        }
                    }
                }
            }
        }
    }

    /// Post-current update: Complete the current calculation.
    ///
    /// This implements the second part of eq. 7.90 from Taflove.
    pub fn post_current_update(&mut self, h_field: &mut VectorField3D) {
        for region in &mut self.regions {
            for li in 0..region.size[0] {
                for lj in 0..region.size[1] {
                    for lk in 0..region.size[2] {
                        let gi = li + region.start[0];
                        let gj = lj + region.start[1];
                        let gk = lk + region.start[2];

                        if gi >= self.dims.nx || gj >= self.dims.ny || gk >= self.dims.nz {
                            continue;
                        }

                        for n in 0..3 {
                            let field = match n {
                                0 => &mut h_field.x,
                                1 => &mut h_field.y,
                                _ => &mut h_field.z,
                            };

                            let f_help = region.curr_flux[n].get(li, lj, lk);
                            let flux_new = field.get(gi, gj, gk);

                            region.curr_flux[n].set(li, lj, lk, flux_new);

                            let iifn = region.iifn[n].get(li, lj, lk);
                            field.set(gi, gj, gk, f_help + iifn * flux_new);
                        }
                    }
                }
            }
        }
    }

    /// Reset all flux storage to zero.
    pub fn reset(&mut self) {
        for region in &mut self.regions {
            for n in 0..3 {
                region.volt_flux[n].clear();
                region.curr_flux[n].clear();
            }
        }
    }
}

// Keep the old Pml struct for backwards compatibility
pub use Upml as Pml;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pml_config() {
        let config = PmlConfig::new(8);
        assert_eq!(config.layers, 8);

        let sigma = config.calculate_sigma_max(0.001);
        assert!(sigma > 0.0);
    }

    #[test]
    fn test_pml_boundaries() {
        let mut boundaries = PmlBoundaries::uniform(10);
        assert!(boundaries.any_enabled());

        boundaries.disable(PmlBoundary::XMin);
        assert_eq!(boundaries.layers[0], 0);
        assert!(boundaries.any_enabled());
    }

    #[test]
    fn test_upml_creation() {
        let dims = Dimensions { nx: 50, ny: 50, nz: 50 };
        let delta = [0.001, 0.001, 0.001];
        let dt = 1e-12;
        let boundaries = PmlBoundaries::uniform(8);

        let upml = Upml::new(dims, delta, dt, boundaries);

        // Should have 6 regions (one per boundary)
        assert!(upml.is_active());
        assert_eq!(upml.num_regions(), 6);
    }

    #[test]
    fn test_upml_disabled_boundary() {
        let dims = Dimensions { nx: 50, ny: 50, nz: 50 };
        let delta = [0.001, 0.001, 0.001];
        let dt = 1e-12;

        let mut boundaries = PmlBoundaries::uniform(8);
        boundaries.disable(PmlBoundary::XMin);
        boundaries.disable(PmlBoundary::XMax);

        let upml = Upml::new(dims, delta, dt, boundaries);

        // Should have 4 regions (Y and Z boundaries only)
        assert!(upml.is_active());
        assert_eq!(upml.num_regions(), 4);
    }

    #[test]
    fn test_kappa_grading() {
        let dims = Dimensions { nx: 50, ny: 50, nz: 50 };
        let delta = [0.001, 0.001, 0.001];
        let dt = 1e-12;
        let boundaries = PmlBoundaries::uniform(8);

        let upml = Upml::new(dims, delta, dt, boundaries);

        // Kappa should be zero outside PML
        let (kv, ki) = upml.calculate_kappa(0, 0, [25, 25, 25], true, 8, 0.008, 1000.0, 3.0);
        assert_eq!(kv, 0.0);
        assert_eq!(ki, 0.0);

        // Kappa should be non-zero inside PML
        let (kv, ki) = upml.calculate_kappa(0, 0, [4, 25, 25], true, 8, 0.008, 1000.0, 3.0);
        assert!(kv > 0.0);
    }

    #[test]
    fn test_upml_field_update() {
        let dims = Dimensions { nx: 20, ny: 20, nz: 20 };
        let delta = [0.001, 0.001, 0.001];
        let dt = 1e-12;
        let boundaries = PmlBoundaries::uniform(4);

        let mut upml = Upml::new(dims, delta, dt, boundaries);

        let mut e_field = VectorField3D::new(dims);
        let mut h_field = VectorField3D::new(dims);

        // Set some initial field values
        e_field.z.set(10, 10, 10, 1.0);

        // Run through the PML update cycle
        upml.pre_current_update(&mut h_field);
        // (main engine H update would go here)
        upml.post_current_update(&mut h_field);

        upml.pre_voltage_update(&mut e_field);
        // (main engine E update would go here)
        upml.post_voltage_update(&mut e_field);

        // Field in center should be unchanged (not in PML)
        let e_center = e_field.z.get(10, 10, 10);
        assert!((e_center - 1.0).abs() < 1e-6);
    }
}
