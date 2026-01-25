//! FDTD (Finite-Difference Time-Domain) solver implementation.
//!
//! This module provides the core FDTD solver with SIMD-accelerated
//! field updates and multi-threaded execution.
//!
//! ## Architecture
//!
//! The FDTD solver is split into two main components:
//!
//! - [`Operator`]: Pre-computes and stores material coefficients
//! - [`Engine`]: Performs the time-stepping field updates
//!
//! ## Usage
//!
//! ```rust,no_run
//! use openems::fdtd::{Simulation, EndCondition};
//! use openems::geometry::Grid;
//!
//! let grid = Grid::uniform(100, 100, 100, 1e-3);
//! let mut sim = Simulation::new(grid);
//! sim.set_end_condition(EndCondition::Timesteps(10000));
//! sim.run().unwrap();
//! ```

mod engine;
mod excitation;
mod operator;
mod simulation;
mod timestep;
mod cylindrical;

pub use engine::Engine;
pub use excitation::{Excitation, ExcitationType};
pub use operator::Operator;
pub use simulation::{EndCondition, Simulation, SimulationStats};
pub use timestep::TimestepInfo;
pub use cylindrical::{CylindricalEngine, CylindricalGrid, CylindricalOperator};

/// Engine type selection
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EngineType {
    /// Single-threaded reference implementation
    Basic,
    /// SIMD-accelerated with automatic dispatch
    Simd,
    /// Multi-threaded SIMD (default)
    Parallel,
}

impl Default for EngineType {
    fn default() -> Self {
        Self::Parallel
    }
}

/// Boundary condition types
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BoundaryCondition {
    /// Perfect Electric Conductor
    Pec,
    /// Perfect Magnetic Conductor
    Pmc,
    /// Mur's first-order absorbing BC
    MurAbc,
    /// Uniaxial Perfectly Matched Layer
    Pml {
        /// Number of PML layers
        layers: usize,
    },
    /// Periodic boundary
    Periodic,
}

impl Default for BoundaryCondition {
    fn default() -> Self {
        Self::Pec
    }
}

/// Boundary conditions for all six faces of the simulation domain.
#[derive(Debug, Clone)]
pub struct BoundaryConditions {
    /// Boundary at x_min
    pub x_min: BoundaryCondition,
    /// Boundary at x_max
    pub x_max: BoundaryCondition,
    /// Boundary at y_min
    pub y_min: BoundaryCondition,
    /// Boundary at y_max
    pub y_max: BoundaryCondition,
    /// Boundary at z_min
    pub z_min: BoundaryCondition,
    /// Boundary at z_max
    pub z_max: BoundaryCondition,
}

impl Default for BoundaryConditions {
    fn default() -> Self {
        Self {
            x_min: BoundaryCondition::Pml { layers: 8 },
            x_max: BoundaryCondition::Pml { layers: 8 },
            y_min: BoundaryCondition::Pml { layers: 8 },
            y_max: BoundaryCondition::Pml { layers: 8 },
            z_min: BoundaryCondition::Pml { layers: 8 },
            z_max: BoundaryCondition::Pml { layers: 8 },
        }
    }
}

impl BoundaryConditions {
    /// Create all-PEC boundaries
    pub fn all_pec() -> Self {
        Self {
            x_min: BoundaryCondition::Pec,
            x_max: BoundaryCondition::Pec,
            y_min: BoundaryCondition::Pec,
            y_max: BoundaryCondition::Pec,
            z_min: BoundaryCondition::Pec,
            z_max: BoundaryCondition::Pec,
        }
    }

    /// Create all-PML boundaries with specified number of layers
    pub fn all_pml(layers: usize) -> Self {
        let bc = BoundaryCondition::Pml { layers };
        Self {
            x_min: bc,
            x_max: bc,
            y_min: bc,
            y_max: bc,
            z_min: bc,
            z_max: bc,
        }
    }

    /// Set boundary from an array [x_min, x_max, y_min, y_max, z_min, z_max]
    pub fn from_array(bc: [BoundaryCondition; 6]) -> Self {
        Self {
            x_min: bc[0],
            x_max: bc[1],
            y_min: bc[2],
            y_max: bc[3],
            z_min: bc[4],
            z_max: bc[5],
        }
    }
}
