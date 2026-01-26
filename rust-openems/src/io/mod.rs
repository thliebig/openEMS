//! File I/O for simulation data.
//!
//! Supports reading and writing simulation data in various formats:
//! - HDF5 for scientific data storage
//! - VTK for visualization
//! - XML for configuration files

pub mod vtk;
pub mod xml;
pub mod hdf5;

use crate::Result;
use std::path::Path;

// Re-export HDF5 types for convenience
pub use hdf5::{Hdf5Reader, Hdf5Writer, MeshType};

/// Trait for types that can be written to files.
pub trait Writable {
    /// Write to HDF5 format.
    #[cfg(feature = "hdf5")]
    fn write_hdf5(&self, path: &Path) -> Result<()>;

    /// Write to VTK format.
    fn write_vtk(&self, path: &Path) -> Result<()>;
}

/// Trait for types that can be read from files.
pub trait Readable: Sized {
    /// Read from HDF5 format.
    #[cfg(feature = "hdf5")]
    fn read_hdf5(path: &Path) -> Result<Self>;
}
