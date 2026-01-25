//! FDTD Extensions for boundary conditions and special materials.
//!
//! This module provides extensions to the base FDTD engine for:
//! - Perfectly Matched Layer (PML) absorbing boundaries
//! - Mur's absorbing boundary condition
//! - Dispersive materials (Lorentz, Debye, Drude)
//! - Lumped RLC elements
//! - Total-Field/Scattered-Field (TF/SF) boundaries

mod pml;

pub use pml::{Pml, PmlConfig};

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
