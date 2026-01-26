//! SAR (Specific Absorption Rate) calculation.
//!
//! Implements SAR calculation for electromagnetic exposure assessment.
//!
//! Reference: IEEE Std C95.3, IEEE Std 62704-1

use crate::arrays::{Dimensions, Field3D};
use num_complex::Complex32;
use std::f64::consts::PI;
use std::str::FromStr;

/// SAR averaging method.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SarAveragingMethod {
    /// IEEE C95.3 method (ICNIRP standard)
    IeeeC95_3,
    /// IEEE 62704 method (IEC/IEEE standard, strictest)
    #[default]
    Ieee62704,
    /// Simple averaging (basic cubic averaging)
    Simple,
    /// Local SAR (no averaging, point values)
    Local,
}

impl FromStr for SarAveragingMethod {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "ieee_c95_3" | "c95.3" | "icnirp" => Ok(Self::IeeeC95_3),
            "ieee_62704" | "62704" | "iec" => Ok(Self::Ieee62704),
            "simple" => Ok(Self::Simple),
            "local" => Ok(Self::Local),
            _ => Err(format!("Unknown SAR averaging method: {}", s)),
        }
    }
}

/// Configuration for SAR calculation.
#[derive(Debug, Clone)]
pub struct SarConfig {
    /// Averaging method
    pub method: SarAveragingMethod,
    /// Averaging mass in kg (typically 0.001 for 1g or 0.01 for 10g)
    pub averaging_mass: f64,
    /// Mass tolerance for cubical fitting (relative)
    pub mass_tolerance: f64,
    /// Maximum iterations for mass fitting
    pub max_mass_iterations: u32,
    /// Maximum background ratio for valid averaging
    pub max_bg_ratio: f64,
    /// Use cell conductivity (true) or E*J (false)
    pub use_cell_conductivity: bool,
}

impl Default for SarConfig {
    fn default() -> Self {
        Self {
            method: SarAveragingMethod::Ieee62704,
            averaging_mass: 0.01, // 10g averaging
            mass_tolerance: 0.05, // 5% tolerance
            max_mass_iterations: 100,
            max_bg_ratio: 1.0, // No background allowed
            use_cell_conductivity: true,
        }
    }
}

impl SarConfig {
    /// Create config for 1g averaging.
    pub fn one_gram() -> Self {
        Self {
            averaging_mass: 0.001,
            ..Default::default()
        }
    }

    /// Create config for 10g averaging.
    pub fn ten_gram() -> Self {
        Self {
            averaging_mass: 0.01,
            ..Default::default()
        }
    }

    /// Create config for local (no averaging) SAR.
    pub fn local() -> Self {
        Self {
            method: SarAveragingMethod::Local,
            averaging_mass: 0.0,
            ..Default::default()
        }
    }
}

/// Frequency-domain field storage for SAR calculation.
#[derive(Debug)]
pub struct FrequencyDomainFields {
    /// Complex E-field (Ex, Ey, Ez)
    pub e_field: [Vec<Complex32>; 3],
    /// Complex J-field (current density) (Jx, Jy, Jz)
    pub j_field: Option<[Vec<Complex32>; 3]>,
    /// Dimensions
    dims: Dimensions,
}

impl FrequencyDomainFields {
    /// Create new frequency domain field storage.
    pub fn new(dims: Dimensions) -> Self {
        let size = dims.total();
        Self {
            e_field: [
                vec![Complex32::new(0.0, 0.0); size],
                vec![Complex32::new(0.0, 0.0); size],
                vec![Complex32::new(0.0, 0.0); size],
            ],
            j_field: None,
            dims,
        }
    }

    /// Create with J-field storage.
    pub fn with_j_field(dims: Dimensions) -> Self {
        let size = dims.total();
        Self {
            e_field: [
                vec![Complex32::new(0.0, 0.0); size],
                vec![Complex32::new(0.0, 0.0); size],
                vec![Complex32::new(0.0, 0.0); size],
            ],
            j_field: Some([
                vec![Complex32::new(0.0, 0.0); size],
                vec![Complex32::new(0.0, 0.0); size],
                vec![Complex32::new(0.0, 0.0); size],
            ]),
            dims,
        }
    }

    #[inline]
    fn idx(&self, i: usize, j: usize, k: usize) -> usize {
        i * self.dims.ny * self.dims.nz + j * self.dims.nz + k
    }

    /// Get E-field component at position.
    pub fn get_e(&self, ny: usize, i: usize, j: usize, k: usize) -> Complex32 {
        self.e_field[ny][self.idx(i, j, k)]
    }

    /// Set E-field component at position.
    pub fn set_e(&mut self, ny: usize, i: usize, j: usize, k: usize, val: Complex32) {
        let idx = self.idx(i, j, k);
        self.e_field[ny][idx] = val;
    }

    /// Accumulate time-domain sample to frequency domain.
    ///
    /// Uses DFT: field_fd += field_td * exp(-j*2*pi*f*t) * 2 * dt
    pub fn accumulate_sample(
        &mut self,
        e_field: &crate::arrays::VectorField3D,
        time: f64,
        frequency: f64,
        dt: f64,
    ) {
        let exp_jwt = Complex32::new(
            (2.0 * PI * frequency * time).cos() as f32,
            -(2.0 * PI * frequency * time).sin() as f32,
        );
        let scale = Complex32::new(2.0 * dt as f32, 0.0);
        let factor = exp_jwt * scale;

        for i in 0..self.dims.nx {
            for j in 0..self.dims.ny {
                for k in 0..self.dims.nz {
                    let idx = self.idx(i, j, k);

                    // E-field
                    let ex = Complex32::new(e_field.x.get(i, j, k), 0.0);
                    let ey = Complex32::new(e_field.y.get(i, j, k), 0.0);
                    let ez = Complex32::new(e_field.z.get(i, j, k), 0.0);

                    self.e_field[0][idx] += ex * factor;
                    self.e_field[1][idx] += ey * factor;
                    self.e_field[2][idx] += ez * factor;
                }
            }
        }
    }
}

/// SAR calculator.
pub struct SarCalculation {
    /// Configuration
    config: SarConfig,
    /// Grid dimensions
    dims: Dimensions,
    /// Cell widths in each direction
    cell_width: [Vec<f64>; 3],
    /// Cell volumes (optional, computed if not provided)
    cell_volume: Option<Field3D>,
    /// Cell densities (kg/m³)
    cell_density: Field3D,
    /// Cell conductivities (S/m)
    cell_conductivity: Option<Field3D>,
    /// Voxel used flags for averaging
    vx_used: Vec<bool>,
    /// Voxel valid flags for averaging
    vx_valid: Vec<bool>,
    /// Statistics
    stats: SarStats,
}

/// SAR calculation statistics.
#[derive(Debug, Default, Clone)]
pub struct SarStats {
    /// Number of valid voxels
    pub valid: usize,
    /// Number of used voxels
    pub used: usize,
    /// Number of unused voxels
    pub unused: usize,
    /// Number of air voxels (zero density)
    pub air_voxels: usize,
    /// Peak local SAR (W/kg)
    pub peak_local_sar: f64,
    /// Peak averaged SAR (W/kg)
    pub peak_averaged_sar: f64,
    /// Total absorbed power (W)
    pub total_power: f64,
}

impl SarCalculation {
    /// Create a new SAR calculator.
    pub fn new(
        dims: Dimensions,
        cell_width: [Vec<f64>; 3],
        cell_density: Field3D,
        config: SarConfig,
    ) -> Self {
        let size = dims.total();
        Self {
            config,
            dims,
            cell_width,
            cell_volume: None,
            cell_density,
            cell_conductivity: None,
            vx_used: vec![false; size],
            vx_valid: vec![true; size],
            stats: SarStats::default(),
        }
    }

    /// Set cell conductivity distribution.
    pub fn set_conductivity(&mut self, conductivity: Field3D) {
        self.cell_conductivity = Some(conductivity);
    }

    /// Set cell volumes (optional, for speedup).
    pub fn set_cell_volumes(&mut self, volumes: Field3D) {
        self.cell_volume = Some(volumes);
    }

    /// Reset calculation state.
    pub fn reset(&mut self) {
        self.vx_used.fill(false);
        self.vx_valid.fill(true);
        self.stats = SarStats::default();
    }

    /// Get cell volume at position.
    fn cell_volume(&self, pos: [usize; 3]) -> f64 {
        if let Some(ref vol) = self.cell_volume {
            vol.get(pos[0], pos[1], pos[2]) as f64
        } else {
            // Compute from cell widths
            let dx = if pos[0] < self.cell_width[0].len() {
                self.cell_width[0][pos[0]]
            } else {
                self.cell_width[0].last().copied().unwrap_or(1.0)
            };
            let dy = if pos[1] < self.cell_width[1].len() {
                self.cell_width[1][pos[1]]
            } else {
                self.cell_width[1].last().copied().unwrap_or(1.0)
            };
            let dz = if pos[2] < self.cell_width[2].len() {
                self.cell_width[2][pos[2]]
            } else {
                self.cell_width[2].last().copied().unwrap_or(1.0)
            };
            dx * dy * dz
        }
    }

    /// Get cell mass at position.
    fn cell_mass(&self, pos: [usize; 3]) -> f64 {
        let density = self.cell_density.get(pos[0], pos[1], pos[2]) as f64;
        density * self.cell_volume(pos)
    }

    /// Calculate local power density at position.
    ///
    /// P = 0.5 * σ * |E|² (using conductivity)
    /// or P = 0.5 * |E| * |J| (using current density)
    fn calc_local_power_density(&self, pos: [usize; 3], fd_fields: &FrequencyDomainFields) -> f64 {
        let ex = fd_fields.get_e(0, pos[0], pos[1], pos[2]);
        let ey = fd_fields.get_e(1, pos[0], pos[1], pos[2]);
        let ez = fd_fields.get_e(2, pos[0], pos[1], pos[2]);

        let e_mag_sq = ex.norm_sqr() + ey.norm_sqr() + ez.norm_sqr();

        if self.config.use_cell_conductivity {
            if let Some(ref cond) = self.cell_conductivity {
                let sigma = cond.get(pos[0], pos[1], pos[2]) as f64;
                0.5 * sigma * e_mag_sq as f64
            } else {
                0.0
            }
        } else if let Some(ref j_field) = fd_fields.j_field {
            // Use E*J formulation
            let idx = pos[0] * self.dims.ny * self.dims.nz + pos[1] * self.dims.nz + pos[2];
            let jx = j_field[0][idx];
            let jy = j_field[1][idx];
            let jz = j_field[2][idx];

            // P = 0.5 * Re(E · J*)
            let power = 0.5
                * ((ex * jx.conj()).re + (ey * jy.conj()).re + (ez * jz.conj()).re) as f64;
            power.max(0.0)
        } else {
            0.0
        }
    }

    /// Calculate local SAR (no averaging).
    pub fn calc_local_sar(&mut self, fd_fields: &FrequencyDomainFields) -> Field3D {
        let mut sar = Field3D::new(self.dims);

        self.stats.peak_local_sar = 0.0;
        self.stats.total_power = 0.0;
        self.stats.air_voxels = 0;

        for i in 0..self.dims.nx {
            for j in 0..self.dims.ny {
                for k in 0..self.dims.nz {
                    let pos = [i, j, k];
                    let density = self.cell_density.get(i, j, k) as f64;

                    if density <= 0.0 {
                        self.stats.air_voxels += 1;
                        continue;
                    }

                    let power_density = self.calc_local_power_density(pos, fd_fields);
                    let local_sar = power_density / density;

                    sar.set(i, j, k, local_sar as f32);

                    if local_sar > self.stats.peak_local_sar {
                        self.stats.peak_local_sar = local_sar;
                    }

                    self.stats.total_power += power_density * self.cell_volume(pos);
                }
            }
        }

        sar
    }

    /// Calculate averaged SAR using cubical mass fitting.
    pub fn calc_averaged_sar(&mut self, fd_fields: &FrequencyDomainFields) -> Field3D {
        // First calculate local SAR
        let local_sar = self.calc_local_sar(fd_fields);

        if self.config.method == SarAveragingMethod::Local {
            return local_sar;
        }

        let mut avg_sar = Field3D::new(self.dims);
        self.stats.peak_averaged_sar = 0.0;

        // Simple cubic averaging implementation
        let target_mass = self.config.averaging_mass;

        for i in 0..self.dims.nx {
            for j in 0..self.dims.ny {
                for k in 0..self.dims.nz {
                    let density = self.cell_density.get(i, j, k) as f64;
                    if density <= 0.0 {
                        continue;
                    }

                    // Find cubic region with target mass
                    let (avg_value, _mass) =
                        self.calc_cubic_average(&local_sar, [i, j, k], target_mass);

                    avg_sar.set(i, j, k, avg_value as f32);

                    if avg_value > self.stats.peak_averaged_sar {
                        self.stats.peak_averaged_sar = avg_value;
                    }
                }
            }
        }

        avg_sar
    }

    /// Calculate cubic average around a position.
    fn calc_cubic_average(
        &self,
        local_sar: &Field3D,
        center: [usize; 3],
        target_mass: f64,
    ) -> (f64, f64) {
        // Iteratively expand cubic region until target mass is reached
        let mut half_size = 1usize;
        let max_half_size = self
            .dims
            .nx
            .min(self.dims.ny)
            .min(self.dims.nz)
            .min(20);

        let mut best_sar = local_sar.get(center[0], center[1], center[2]) as f64;
        let mut best_mass = self.cell_mass(center);

        for _ in 0..self.config.max_mass_iterations {
            let (sum_sar, total_mass) = self.sum_cubic_region(local_sar, center, half_size);

            if total_mass >= target_mass * (1.0 - self.config.mass_tolerance) {
                // Close enough to target mass
                if total_mass > 0.0 {
                    best_sar = sum_sar / total_mass;
                    best_mass = total_mass;
                }
                break;
            }

            best_sar = if total_mass > 0.0 {
                sum_sar / total_mass
            } else {
                0.0
            };
            best_mass = total_mass;

            half_size += 1;
            if half_size > max_half_size {
                break;
            }
        }

        (best_sar, best_mass)
    }

    /// Sum SAR values in a cubic region.
    fn sum_cubic_region(
        &self,
        local_sar: &Field3D,
        center: [usize; 3],
        half_size: usize,
    ) -> (f64, f64) {
        let i_start = center[0].saturating_sub(half_size);
        let i_end = (center[0] + half_size + 1).min(self.dims.nx);
        let j_start = center[1].saturating_sub(half_size);
        let j_end = (center[1] + half_size + 1).min(self.dims.ny);
        let k_start = center[2].saturating_sub(half_size);
        let k_end = (center[2] + half_size + 1).min(self.dims.nz);

        let mut sum_sar = 0.0;
        let mut total_mass = 0.0;

        for i in i_start..i_end {
            for j in j_start..j_end {
                for k in k_start..k_end {
                    let pos = [i, j, k];
                    let mass = self.cell_mass(pos);
                    let sar = local_sar.get(i, j, k) as f64;

                    sum_sar += sar * mass;
                    total_mass += mass;
                }
            }
        }

        (sum_sar, total_mass)
    }

    /// Calculate total absorbed power.
    pub fn calc_total_power(&self, fd_fields: &FrequencyDomainFields) -> f64 {
        let mut total_power = 0.0;

        for i in 0..self.dims.nx {
            for j in 0..self.dims.ny {
                for k in 0..self.dims.nz {
                    let pos = [i, j, k];
                    let power_density = self.calc_local_power_density(pos, fd_fields);
                    total_power += power_density * self.cell_volume(pos);
                }
            }
        }

        total_power
    }

    /// Get calculation statistics.
    pub fn stats(&self) -> &SarStats {
        &self.stats
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sar_config() {
        let config = SarConfig::default();
        assert_eq!(config.method, SarAveragingMethod::Ieee62704);
        assert!((config.averaging_mass - 0.01).abs() < 1e-6);

        let config = SarConfig::one_gram();
        assert!((config.averaging_mass - 0.001).abs() < 1e-6);
    }

    #[test]
    fn test_averaging_method_parse() {
        assert_eq!(
            SarAveragingMethod::from_str("ieee_c95_3"),
            Ok(SarAveragingMethod::IeeeC95_3)
        );
        assert_eq!(
            SarAveragingMethod::from_str("62704"),
            Ok(SarAveragingMethod::Ieee62704)
        );
        assert_eq!(
            SarAveragingMethod::from_str("local"),
            Ok(SarAveragingMethod::Local)
        );
    }

    #[test]
    fn test_fd_fields_creation() {
        let dims = Dimensions::new(10, 10, 10);
        let fields = FrequencyDomainFields::new(dims);

        assert_eq!(fields.e_field[0].len(), 1000);
        assert!(fields.j_field.is_none());

        let fields = FrequencyDomainFields::with_j_field(dims);
        assert!(fields.j_field.is_some());
    }

    #[test]
    fn test_sar_calculation_creation() {
        let dims = Dimensions::new(10, 10, 10);
        let cell_width = [
            vec![0.001; 10],
            vec![0.001; 10],
            vec![0.001; 10],
        ];

        let mut density = Field3D::new(dims);
        density.fill(1000.0); // 1000 kg/m³ (water-like)

        let sar_calc = SarCalculation::new(dims, cell_width, density, SarConfig::default());
        assert_eq!(sar_calc.dims.nx, 10);
    }

    #[test]
    fn test_local_sar_calculation() {
        let dims = Dimensions::new(5, 5, 5);
        let cell_width = [
            vec![0.001; 5],
            vec![0.001; 5],
            vec![0.001; 5],
        ];

        let mut density = Field3D::new(dims);
        density.fill(1000.0);

        let mut conductivity = Field3D::new(dims);
        conductivity.fill(0.5); // 0.5 S/m

        let mut sar_calc = SarCalculation::new(dims, cell_width, density, SarConfig::local());
        sar_calc.set_conductivity(conductivity);

        // Create frequency domain fields with uniform E-field
        let mut fd_fields = FrequencyDomainFields::new(dims);
        for i in 0..dims.nx {
            for j in 0..dims.ny {
                for k in 0..dims.nz {
                    // E = 100 V/m
                    fd_fields.set_e(0, i, j, k, Complex32::new(100.0, 0.0));
                }
            }
        }

        let sar = sar_calc.calc_local_sar(&fd_fields);

        // SAR = 0.5 * σ * |E|² / ρ = 0.5 * 0.5 * 10000 / 1000 = 2.5 W/kg
        let expected_sar = 2.5;
        let calculated_sar = sar.get(2, 2, 2) as f64;
        assert!(
            (calculated_sar - expected_sar).abs() < 0.1,
            "Expected SAR ~{}, got {}",
            expected_sar,
            calculated_sar
        );
    }

    #[test]
    fn test_averaged_sar() {
        let dims = Dimensions::new(10, 10, 10);
        let cell_width = [
            vec![0.001; 10],
            vec![0.001; 10],
            vec![0.001; 10],
        ];

        let mut density = Field3D::new(dims);
        density.fill(1000.0);

        let mut conductivity = Field3D::new(dims);
        conductivity.fill(0.5);

        let mut config = SarConfig::default();
        config.averaging_mass = 1e-6; // Very small for test grid

        let mut sar_calc = SarCalculation::new(dims, cell_width, density, config);
        sar_calc.set_conductivity(conductivity);

        let mut fd_fields = FrequencyDomainFields::new(dims);
        // Non-uniform E-field - higher in center
        for i in 0..dims.nx {
            for j in 0..dims.ny {
                for k in 0..dims.nz {
                    let dist = ((i as f32 - 5.0).powi(2)
                        + (j as f32 - 5.0).powi(2)
                        + (k as f32 - 5.0).powi(2))
                    .sqrt();
                    let e_val = 100.0 * (-dist / 3.0).exp();
                    fd_fields.set_e(0, i, j, k, Complex32::new(e_val, 0.0));
                }
            }
        }

        let avg_sar = sar_calc.calc_averaged_sar(&fd_fields);

        // Averaged SAR should be smoothed compared to local
        // Just verify it computed without panicking
        let center_sar = avg_sar.get(5, 5, 5);
        assert!(center_sar > 0.0);
    }

    #[test]
    fn test_total_power() {
        let dims = Dimensions::new(5, 5, 5);
        let cell_width = [
            vec![0.001; 5], // 1mm cells
            vec![0.001; 5],
            vec![0.001; 5],
        ];

        let mut density = Field3D::new(dims);
        density.fill(1000.0);

        let mut conductivity = Field3D::new(dims);
        conductivity.fill(0.5);

        let mut sar_calc = SarCalculation::new(dims, cell_width, density, SarConfig::default());
        sar_calc.set_conductivity(conductivity);

        let mut fd_fields = FrequencyDomainFields::new(dims);
        for i in 0..dims.nx {
            for j in 0..dims.ny {
                for k in 0..dims.nz {
                    fd_fields.set_e(0, i, j, k, Complex32::new(100.0, 0.0));
                }
            }
        }

        let total_power = sar_calc.calc_total_power(&fd_fields);

        // Power density = 0.5 * σ * |E|² = 0.5 * 0.5 * 10000 = 2500 W/m³
        // Total volume = 125 * (0.001)³ = 1.25e-7 m³
        // Total power = 2500 * 1.25e-7 = 3.125e-4 W
        let expected_power = 2500.0 * 125.0 * 1e-9;
        assert!(
            (total_power - expected_power).abs() / expected_power < 0.01,
            "Expected power ~{}, got {}",
            expected_power,
            total_power
        );
    }

    #[test]
    fn test_sar_stats() {
        let dims = Dimensions::new(5, 5, 5);
        let cell_width = [
            vec![0.001; 5],
            vec![0.001; 5],
            vec![0.001; 5],
        ];

        let mut density = Field3D::new(dims);
        // Half filled with tissue, half air
        for i in 0..dims.nx {
            for j in 0..dims.ny {
                for k in 0..dims.nz {
                    if i < 3 {
                        density.set(i, j, k, 1000.0);
                    } else {
                        density.set(i, j, k, 0.0);
                    }
                }
            }
        }

        let mut conductivity = Field3D::new(dims);
        conductivity.fill(0.5);

        let mut sar_calc = SarCalculation::new(dims, cell_width, density, SarConfig::local());
        sar_calc.set_conductivity(conductivity);

        let fd_fields = FrequencyDomainFields::new(dims);
        let _sar = sar_calc.calc_local_sar(&fd_fields);

        let stats = sar_calc.stats();
        assert_eq!(stats.air_voxels, 50); // 2 * 5 * 5 = 50 air voxels
    }
}
