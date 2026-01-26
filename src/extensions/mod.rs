//! FDTD Extensions for boundary conditions and special materials.
//!
//! This module provides extensions to the base FDTD engine for:
//! - Perfectly Matched Layer (PML) absorbing boundaries
//! - Mur's absorbing boundary condition (domain boundaries)
//! - Local absorbing boundary condition (arbitrary sheet primitives)
//! - Dispersive materials (Lorentz, Debye, Drude)
//! - Lumped RLC elements
//! - Total-Field/Scattered-Field (TF/SF) boundaries
//! - Steady-state detection
//! - Conducting sheet model

mod pml;
mod mur_abc;
mod local_absorbing_bc;
mod dispersive;
mod lumped_rlc;
mod tfsf;
mod steady_state;
mod conducting_sheet;

// PML exports
pub use pml::{Pml, PmlConfig, PmlBoundaries, PmlBoundary, Upml};

// Mur ABC exports (domain boundaries)
pub use mur_abc::{MurAbc, MurAbcConfig};

// Local absorbing BC exports (arbitrary sheet primitives)
pub use local_absorbing_bc::{
    AbcType, LocalAbsorbingBc, LocalAbsorbingBcConfig, LocalAbsorbingSheet,
};

// Dispersive material exports
pub use dispersive::{
    DispersiveMaterial,
    LorentzMaterial, LorentzParams,
    DrudeMaterial, DrudeParams,
    DebyeMaterial, DebyeParams,
};

// Lumped RLC exports
pub use lumped_rlc::{LumpedRlc, RlcElement, RlcType};

// TF/SF exports
pub use tfsf::{TfsfBoundary, TfsfConfig, PropagationDirection, Polarization};

// Steady-state exports
pub use steady_state::{SteadyStateDetector, SteadyStateConfig, SteadyStateResult};

// Conducting sheet exports
pub use conducting_sheet::{
    ConductingSheet, ConductingSheetConfig, ConductingSheetManager, SheetLorentzParams,
};

/// Extension trait for FDTD operator modifications.
pub trait OperatorExtension {
    /// Apply extension to operator coefficients.
    fn apply_to_operator(&self, operator: &mut crate::fdtd::Operator);

    /// Get extension name for logging.
    fn name(&self) -> &str;
}

/// Extension trait for FDTD engine modifications.
pub trait EngineExtension {
    /// Pre-update hook (called before H-field update).
    fn pre_update_h(&mut self, engine: &mut crate::fdtd::Engine);

    /// Post-update H hook (called after H-field update).
    fn post_update_h(&mut self, engine: &mut crate::fdtd::Engine);

    /// Pre-update E hook (called before E-field update).
    fn pre_update_e(&mut self, engine: &mut crate::fdtd::Engine);

    /// Post-update E hook (called after E-field update).
    fn post_update_e(&mut self, engine: &mut crate::fdtd::Engine);

    /// Get extension name for logging.
    fn name(&self) -> &str;
}
