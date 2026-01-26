//! Python bindings for openEMS using PyO3.
//!
//! This module provides Python bindings that expose the core openEMS
//! functionality to Python users.

#[cfg(feature = "python")]
use pyo3::exceptions::PyRuntimeError;
#[cfg(feature = "python")]
use pyo3::prelude::*;

#[cfg(feature = "python")]
use crate::fdtd::{BoundaryConditions, EndCondition, EngineType, Simulation as RustSimulation};
#[cfg(feature = "python")]
use crate::geometry::Grid as RustGrid;

/// Python wrapper for Grid
#[cfg(feature = "python")]
#[pyclass(name = "Grid")]
pub struct PyGrid {
    inner: RustGrid,
}

#[cfg(feature = "python")]
#[pymethods]
impl PyGrid {
    /// Create a uniform grid
    #[staticmethod]
    fn uniform(nx: usize, ny: usize, nz: usize, delta: f64) -> Self {
        Self {
            inner: RustGrid::uniform(nx, ny, nz, delta),
        }
    }

    /// Create a grid from mesh line arrays
    #[staticmethod]
    fn from_lines(x: Vec<f64>, y: Vec<f64>, z: Vec<f64>) -> Self {
        Self {
            inner: RustGrid::cartesian(x, y, z),
        }
    }

    /// Get the number of cells
    fn num_cells(&self) -> usize {
        self.inner.dimensions().total()
    }

    /// Get the minimum cell size
    fn cell_size(&self) -> (f64, f64, f64) {
        self.inner.cell_size()
    }
}

/// Python wrapper for openEMS simulation
#[cfg(feature = "python")]
#[pyclass(name = "OpenEMS")]
pub struct PyOpenEMS {
    grid: Option<RustGrid>,
    num_timesteps: u64,
    end_criteria_db: Option<f64>,
    engine_type: EngineType,
    excitations: Vec<ExcitationConfig>,
}

#[cfg(feature = "python")]
struct ExcitationConfig {
    freq: f64,
    bandwidth: f64,
    direction: usize,
    position: (usize, usize, usize),
}

#[cfg(feature = "python")]
#[pymethods]
impl PyOpenEMS {
    /// Create a new openEMS simulation
    #[new]
    #[pyo3(signature = (num_timesteps=10000))]
    fn new(num_timesteps: u64) -> Self {
        Self {
            grid: None,
            num_timesteps,
            end_criteria_db: None,
            engine_type: EngineType::Parallel,
            excitations: Vec::new(),
        }
    }

    /// Set the computational grid
    fn set_grid(&mut self, grid: &PyGrid) {
        self.grid = Some(grid.inner.clone());
    }

    /// Set Gaussian excitation
    #[pyo3(signature = (f0, fc))]
    fn set_gauss_excite(&mut self, f0: f64, fc: f64) {
        // f0 = center frequency, fc = cutoff frequency
        // Store for later use when creating simulation
        self.excitations.clear();
        // Will be properly configured when we have position info
    }

    /// Set the boundary conditions
    /// [x_min, x_max, y_min, y_max, z_min, z_max]
    /// 0 = PEC, 1 = PMC, 2 = MUR, 3 = PML
    fn set_boundary_cond(&mut self, _bc: Vec<i32>) {
        // Store boundary conditions for later
    }

    /// Run the simulation
    fn run(&self, sim_path: &str, cleanup: bool) -> PyResult<()> {
        let grid = self
            .grid
            .as_ref()
            .ok_or_else(|| PyRuntimeError::new_err("Grid not set. Call set_grid() first."))?;

        let mut sim = RustSimulation::new(grid.clone());
        sim.set_engine_type(self.engine_type);

        if let Some(db) = self.end_criteria_db {
            sim.set_end_condition(EndCondition::EnergyDecay(db));
        } else {
            sim.set_end_condition(EndCondition::Timesteps(self.num_timesteps));
        }

        sim.setup()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        sim.run()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;

        Ok(())
    }
}

/// Python module definition
#[cfg(feature = "python")]
#[pymodule]
fn openems(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyGrid>()?;
    m.add_class::<PyOpenEMS>()?;

    // Add constants
    m.add("C0", crate::constants::C0)?;
    m.add("EPS0", crate::constants::EPS0)?;
    m.add("MU0", crate::constants::MU0)?;
    m.add("Z0", crate::constants::Z0)?;
    m.add("VERSION", crate::VERSION)?;

    Ok(())
}
