//! Post-processing and analysis module.
//!
//! Provides tools for analyzing simulation results including:
//! - Field probes (voltage, current, field values)
//! - S-parameter calculation
//! - Field dumps (time and frequency domain)
//! - SAR calculation
//! - Mode matching
//! - Processing array container

mod array;
mod fields_fd;
mod mode_match;
mod probes;
mod sar;

pub use array::{FieldRecorder, ProcessingArray, Processor};
pub use fields_fd::{FdDumpType, FrequencyDomainFieldDump, TimeDomainFieldDump};
pub use mode_match::{
    CircularTeMode, CustomModeFunction, FieldType, ModeFunction, ModeMatch, ModeMatchConfig,
    ModeMatchResult, RectangularTeMode, RectangularTmMode,
};
pub use probes::{CurrentProbe, FieldProbe, Processing, VoltageProbe};
pub use sar::{FrequencyDomainFields, SarAveragingMethod, SarCalculation, SarConfig, SarStats};

/// Probe types for extracting data from simulations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProbeType {
    /// Voltage probe (line integral of E-field)
    Voltage,
    /// Current probe (surface integral of H-field)
    Current,
    /// Point field probe
    Field,
    /// Box field dump
    FieldDump,
}

/// Frequency domain data.
#[derive(Debug, Clone)]
pub struct FrequencyData {
    /// Frequencies (Hz)
    pub frequencies: Vec<f64>,
    /// Complex values at each frequency
    pub values: Vec<num_complex::Complex64>,
}

impl FrequencyData {
    /// Create new frequency data.
    pub fn new(frequencies: Vec<f64>, values: Vec<num_complex::Complex64>) -> Self {
        Self {
            frequencies,
            values,
        }
    }

    /// Get magnitude in dB.
    pub fn magnitude_db(&self) -> Vec<f64> {
        self.values
            .iter()
            .map(|v| 20.0 * v.norm().log10())
            .collect()
    }

    /// Get phase in degrees.
    pub fn phase_deg(&self) -> Vec<f64> {
        self.values.iter().map(|v| v.arg().to_degrees()).collect()
    }
}
