//! # openEMS - High-Performance FDTD Electromagnetic Field Solver
//!
//! openEMS is a free and open electromagnetic field solver using the
//! Finite-Difference Time-Domain (FDTD) method.
//!
//! ## Features
//!
//! - Fast FDTD solver with SIMD acceleration (AVX2/AVX-512/SSE)
//! - Automatic parallelization using Rayon
//! - Support for Cartesian and cylindrical coordinate systems
//! - Various boundary conditions (PML, MUR ABC, periodic)
//! - Dispersive materials (Lorentz, Debye, Drude)
//! - Near-field to far-field (NF2FF) transformation
//! - HDF5 and VTK output support
//! - Python bindings via PyO3
//!
//! ## Example
//!
//! ```rust,no_run
//! use openems::{Simulation, Grid};
//! use openems::fdtd::EndCondition;
//!
//! let grid = Grid::uniform(100, 100, 100, 1e-3);
//! let mut sim = Simulation::new(grid);
//! sim.set_end_condition(EndCondition::Timesteps(10000));
//! sim.run().unwrap();
//! ```

#![warn(missing_docs)]
#![warn(clippy::all)]
#![allow(clippy::many_single_char_names)] // Common in physics code

pub mod arrays;
pub mod constants;
pub mod extensions;
pub mod fdtd;
pub mod geometry;
pub mod io;
pub mod nf2ff;
pub mod processing;
pub mod tools;

#[cfg(feature = "python")]
pub mod python;

// Re-exports for convenience
pub use arrays::{Field3D, FieldComponent};
pub use constants::*;
pub use fdtd::{Engine, Operator, Simulation};
pub use geometry::{CoordinateSystem, Grid};
pub use processing::Processing;

// Additional re-exports for the new features
pub use fdtd::{EngineInterface, CompressedCoefficients, CylindricalMultigrid};
pub use io::{Hdf5Reader, Hdf5Writer};
pub use processing::{ProcessingArray, FrequencyDomainFieldDump};
pub use tools::{SignalHandler, GlobalOptions, DenormalGuard};

/// Error types for openEMS
#[derive(Debug, thiserror::Error)]
pub enum Error {
    /// Invalid simulation configuration
    #[error("Configuration error: {0}")]
    Config(String),

    /// I/O error
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    /// Grid error
    #[error("Grid error: {0}")]
    Grid(String),

    /// Material error
    #[error("Material error: {0}")]
    Material(String),

    /// Numerical error
    #[error("Numerical error: {0}")]
    Numerical(String),

    /// XML parsing error
    #[error("XML error: {0}")]
    Xml(String),
}

/// Result type for openEMS operations
pub type Result<T> = std::result::Result<T, Error>;

/// Version information
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
