//! Physical and mathematical constants for electromagnetic simulations.
//!
//! All values are in SI units unless otherwise noted.

use std::f64::consts::PI;

/// Speed of light in vacuum (m/s)
pub const C0: f64 = 299_792_458.0;

/// Permittivity of free space (F/m)
pub const EPS0: f64 = 8.854_187_817e-12;

/// Permeability of free space (H/m)
pub const MU0: f64 = 1.256_637_062e-6;

/// Impedance of free space (Ohm)
/// Z0 = sqrt(MU0/EPS0) = MU0 * C0
pub const Z0: f64 = 376.730_313_668;

/// Pi constant
pub const PI_CONST: f64 = PI;

/// 2 * Pi
pub const TWO_PI: f64 = 2.0 * PI;

/// Elementary charge (C)
pub const ELEMENTARY_CHARGE: f64 = 1.602_176_634e-19;

/// Electron mass (kg)
pub const ELECTRON_MASS: f64 = 9.109_383_7e-31;

/// Boltzmann constant (J/K)
pub const BOLTZMANN: f64 = 1.380_649e-23;

/// Planck constant (J*s)
pub const PLANCK: f64 = 6.626_070_15e-34;

/// Reduced Planck constant (J*s)
pub const HBAR: f64 = 1.054_571_817e-34;

/// Convert frequency to angular frequency
#[inline]
pub fn freq_to_omega(freq: f64) -> f64 {
    TWO_PI * freq
}

/// Convert angular frequency to frequency
#[inline]
pub fn omega_to_freq(omega: f64) -> f64 {
    omega / TWO_PI
}

/// Convert frequency to wavelength in free space
#[inline]
pub fn freq_to_wavelength(freq: f64) -> f64 {
    C0 / freq
}

/// Convert wavelength to frequency
#[inline]
pub fn wavelength_to_freq(wavelength: f64) -> f64 {
    C0 / wavelength
}

/// Calculate wave impedance for given permittivity and permeability
#[inline]
pub fn wave_impedance(eps_r: f64, mu_r: f64) -> f64 {
    Z0 * (mu_r / eps_r).sqrt()
}

/// Calculate phase velocity for given permittivity and permeability
#[inline]
pub fn phase_velocity(eps_r: f64, mu_r: f64) -> f64 {
    C0 / (eps_r * mu_r).sqrt()
}

/// Calculate refractive index
#[inline]
pub fn refractive_index(eps_r: f64, mu_r: f64) -> f64 {
    (eps_r * mu_r).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_z0_calculation() {
        // Z0 should equal MU0 * C0
        let z0_calc = MU0 * C0;
        assert!((Z0 - z0_calc).abs() < 1e-6);
    }

    #[test]
    fn test_eps0_mu0_c0_relation() {
        // c0 = 1/sqrt(eps0 * mu0)
        let c0_calc = 1.0 / (EPS0 * MU0).sqrt();
        assert!((C0 - c0_calc).abs() < 1.0);
    }

    #[test]
    fn test_frequency_wavelength() {
        let freq = 1e9; // 1 GHz
        let wavelength = freq_to_wavelength(freq);
        assert!((wavelength - 0.299792458).abs() < 1e-6);
    }
}
