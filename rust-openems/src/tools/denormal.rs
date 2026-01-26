//! Denormal number handling for high-performance FDTD.
//!
//! Denormal (subnormal) floating-point numbers can cause significant
//! performance degradation in FDTD simulations. This module provides
//! utilities to detect and handle denormals.

use std::sync::atomic::{AtomicBool, Ordering};

/// Global flag indicating if flush-to-zero mode is enabled.
static FTZ_ENABLED: AtomicBool = AtomicBool::new(false);

/// Denormal handling modes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DenormalMode {
    /// Flush denormals to zero (FTZ)
    FlushToZero,
    /// Keep denormals (IEEE compliant, slower)
    KeepDenormals,
    /// Platform default
    Default,
}

/// RAII guard for denormal mode.
pub struct DenormalGuard {
    #[cfg(target_arch = "x86_64")]
    original_mxcsr: u32,
    #[cfg(not(target_arch = "x86_64"))]
    _marker: std::marker::PhantomData<()>,
}

impl DenormalGuard {
    /// Enable flush-to-zero mode and return a guard that restores the original state.
    pub fn flush_to_zero() -> Self {
        #[cfg(target_arch = "x86_64")]
        {
            let original = unsafe { get_mxcsr() };
            unsafe { enable_ftz_daz() };
            FTZ_ENABLED.store(true, Ordering::SeqCst);
            Self {
                original_mxcsr: original,
            }
        }

        #[cfg(not(target_arch = "x86_64"))]
        {
            // On non-x86 platforms, we can't easily control FTZ mode
            // Just return a no-op guard
            Self {
                _marker: std::marker::PhantomData,
            }
        }
    }

    /// Check if FTZ mode is currently enabled globally.
    pub fn is_ftz_enabled() -> bool {
        FTZ_ENABLED.load(Ordering::SeqCst)
    }
}

impl Drop for DenormalGuard {
    fn drop(&mut self) {
        #[cfg(target_arch = "x86_64")]
        {
            unsafe { set_mxcsr(self.original_mxcsr) };
            FTZ_ENABLED.store(false, Ordering::SeqCst);
        }
    }
}

// x86_64 specific implementations
#[cfg(target_arch = "x86_64")]
mod x86 {
    /// MXCSR bit for Flush-To-Zero
    pub const FTZ_BIT: u32 = 1 << 15;
    /// MXCSR bit for Denormals-Are-Zero
    pub const DAZ_BIT: u32 = 1 << 6;

    /// Get current MXCSR register value.
    #[inline]
    pub unsafe fn get_mxcsr() -> u32 {
        let mut mxcsr: u32 = 0;
        std::arch::asm!(
            "stmxcsr [{}]",
            in(reg) &mut mxcsr,
            options(nostack)
        );
        mxcsr
    }

    /// Set MXCSR register value.
    #[inline]
    pub unsafe fn set_mxcsr(value: u32) {
        std::arch::asm!(
            "ldmxcsr [{}]",
            in(reg) &value,
            options(nostack)
        );
    }

    /// Enable Flush-To-Zero and Denormals-Are-Zero modes.
    #[inline]
    pub unsafe fn enable_ftz_daz() {
        let mxcsr = get_mxcsr();
        set_mxcsr(mxcsr | FTZ_BIT | DAZ_BIT);
    }

    /// Disable Flush-To-Zero and Denormals-Are-Zero modes.
    #[inline]
    #[allow(dead_code)]
    pub unsafe fn disable_ftz_daz() {
        let mxcsr = get_mxcsr();
        set_mxcsr(mxcsr & !(FTZ_BIT | DAZ_BIT));
    }
}

#[cfg(target_arch = "x86_64")]
use x86::*;

/// Check if a value is denormal.
#[inline]
pub fn is_denormal(value: f32) -> bool {
    let bits = value.to_bits();
    let exponent = (bits >> 23) & 0xFF;
    let mantissa = bits & 0x7F_FFFF;
    exponent == 0 && mantissa != 0
}

/// Check if a value is denormal (f64 version).
#[inline]
pub fn is_denormal_f64(value: f64) -> bool {
    let bits = value.to_bits();
    let exponent = (bits >> 52) & 0x7FF;
    let mantissa = bits & 0xF_FFFF_FFFF_FFFF;
    exponent == 0 && mantissa != 0
}

/// Flush a single value to zero if denormal.
#[inline]
pub fn flush_denormal(value: f32) -> f32 {
    if is_denormal(value) {
        0.0
    } else {
        value
    }
}

/// Flush a single value to zero if denormal (f64 version).
#[inline]
pub fn flush_denormal_f64(value: f64) -> f64 {
    if is_denormal_f64(value) {
        0.0
    } else {
        value
    }
}

/// Count denormals in a slice.
pub fn count_denormals(values: &[f32]) -> usize {
    values.iter().filter(|&&v| is_denormal(v)).count()
}

/// Count denormals in a slice (f64 version).
pub fn count_denormals_f64(values: &[f64]) -> usize {
    values.iter().filter(|&&v| is_denormal_f64(v)).count()
}

/// Flush all denormals in a slice to zero.
pub fn flush_denormals_in_place(values: &mut [f32]) {
    for v in values.iter_mut() {
        if is_denormal(*v) {
            *v = 0.0;
        }
    }
}

/// Flush all denormals in a slice to zero (f64 version).
pub fn flush_denormals_in_place_f64(values: &mut [f64]) {
    for v in values.iter_mut() {
        if is_denormal_f64(*v) {
            *v = 0.0;
        }
    }
}

/// Denormal statistics for monitoring.
#[derive(Debug, Clone, Default)]
pub struct DenormalStats {
    /// Total values checked
    pub total_checked: usize,
    /// Number of denormals found
    pub denormals_found: usize,
    /// Number of times check was performed
    pub check_count: usize,
}

impl DenormalStats {
    /// Create new stats.
    pub fn new() -> Self {
        Self::default()
    }

    /// Update stats with a check result.
    pub fn update(&mut self, total: usize, denormals: usize) {
        self.total_checked += total;
        self.denormals_found += denormals;
        self.check_count += 1;
    }

    /// Get denormal ratio.
    pub fn denormal_ratio(&self) -> f64 {
        if self.total_checked > 0 {
            self.denormals_found as f64 / self.total_checked as f64
        } else {
            0.0
        }
    }

    /// Check if denormal ratio exceeds threshold.
    pub fn exceeds_threshold(&self, threshold: f64) -> bool {
        self.denormal_ratio() > threshold
    }

    /// Reset stats.
    pub fn reset(&mut self) {
        *self = Self::default();
    }
}

/// Denormal monitor for periodic checking during simulation.
pub struct DenormalMonitor {
    stats: DenormalStats,
    check_interval: u64,
    last_check: u64,
    warning_threshold: f64,
    warning_issued: bool,
}

impl DenormalMonitor {
    /// Create a new denormal monitor.
    pub fn new(check_interval: u64) -> Self {
        Self {
            stats: DenormalStats::new(),
            check_interval,
            last_check: 0,
            warning_threshold: 0.01, // 1% denormal ratio triggers warning
            warning_issued: false,
        }
    }

    /// Set warning threshold.
    pub fn set_warning_threshold(&mut self, threshold: f64) {
        self.warning_threshold = threshold;
    }

    /// Check if monitoring is due at this timestep.
    pub fn should_check(&self, timestep: u64) -> bool {
        // Check at timestep 0 always, then periodically based on interval
        self.check_interval > 0
            && (timestep == 0 || timestep - self.last_check >= self.check_interval)
    }

    /// Perform denormal check on a field.
    pub fn check_field(&mut self, values: &[f32], timestep: u64) -> bool {
        if !self.should_check(timestep) {
            return false;
        }

        self.last_check = timestep;
        let denormals = count_denormals(values);
        self.stats.update(values.len(), denormals);

        if !self.warning_issued && self.stats.exceeds_threshold(self.warning_threshold) {
            self.warning_issued = true;
            log::warn!(
                "High denormal ratio detected: {:.2}% ({} / {})",
                self.stats.denormal_ratio() * 100.0,
                self.stats.denormals_found,
                self.stats.total_checked
            );
            true
        } else {
            false
        }
    }

    /// Get current stats.
    pub fn stats(&self) -> &DenormalStats {
        &self.stats
    }

    /// Reset the monitor.
    pub fn reset(&mut self) {
        self.stats.reset();
        self.last_check = 0;
        self.warning_issued = false;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_denormal() {
        // Normal numbers
        assert!(!is_denormal(1.0));
        assert!(!is_denormal(-1.0));
        assert!(!is_denormal(1e-30));

        // Zero
        assert!(!is_denormal(0.0));
        assert!(!is_denormal(-0.0));

        // Denormal (very small but non-zero)
        let denormal = f32::from_bits(1); // Smallest positive denormal
        assert!(is_denormal(denormal));
    }

    #[test]
    fn test_is_denormal_f64() {
        assert!(!is_denormal_f64(1.0));
        assert!(!is_denormal_f64(0.0));

        let denormal = f64::from_bits(1); // Smallest positive denormal
        assert!(is_denormal_f64(denormal));
    }

    #[test]
    fn test_flush_denormal() {
        assert_eq!(flush_denormal(1.0), 1.0);
        assert_eq!(flush_denormal(0.0), 0.0);

        let denormal = f32::from_bits(1);
        assert_eq!(flush_denormal(denormal), 0.0);
    }

    #[test]
    fn test_count_denormals() {
        let normal_values = vec![1.0, 2.0, 3.0, 0.0];
        assert_eq!(count_denormals(&normal_values), 0);

        let denormal = f32::from_bits(1);
        let with_denormals = vec![1.0, denormal, 2.0, denormal];
        assert_eq!(count_denormals(&with_denormals), 2);
    }

    #[test]
    fn test_flush_in_place() {
        let denormal = f32::from_bits(1);
        let mut values = vec![1.0, denormal, 2.0, denormal, 3.0];

        flush_denormals_in_place(&mut values);

        assert_eq!(values[0], 1.0);
        assert_eq!(values[1], 0.0);
        assert_eq!(values[2], 2.0);
        assert_eq!(values[3], 0.0);
        assert_eq!(values[4], 3.0);
    }

    #[test]
    fn test_denormal_stats() {
        let mut stats = DenormalStats::new();

        stats.update(100, 5);
        assert_eq!(stats.total_checked, 100);
        assert_eq!(stats.denormals_found, 5);
        assert!((stats.denormal_ratio() - 0.05).abs() < 1e-10);

        stats.update(100, 10);
        assert_eq!(stats.total_checked, 200);
        assert_eq!(stats.denormals_found, 15);

        assert!(stats.exceeds_threshold(0.01));
        assert!(!stats.exceeds_threshold(0.1));
    }

    #[test]
    fn test_denormal_monitor() {
        let mut monitor = DenormalMonitor::new(10);

        assert!(monitor.should_check(0));
        assert!(!monitor.should_check(5));
        assert!(monitor.should_check(10));

        let values = vec![1.0; 100];
        monitor.check_field(&values, 0);

        assert_eq!(monitor.stats().total_checked, 100);
        assert_eq!(monitor.stats().denormals_found, 0);
    }

    #[cfg(target_arch = "x86_64")]
    #[test]
    fn test_denormal_guard() {
        // This test verifies the guard can be created and dropped without panic
        {
            let _guard = DenormalGuard::flush_to_zero();
            assert!(DenormalGuard::is_ftz_enabled());
        }
        // Guard dropped, FTZ should be disabled
        assert!(!DenormalGuard::is_ftz_enabled());
    }
}
