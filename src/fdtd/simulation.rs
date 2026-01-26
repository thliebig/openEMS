//! High-level simulation control.
//!
//! The Simulation struct provides a convenient interface for setting up
//! and running FDTD simulations.

use super::engine::Engine;
use super::excitation::Excitation;
use super::operator::Operator;
use super::timestep::TimestepInfo;
use super::{BoundaryConditions, EngineType};
use crate::geometry::Grid;
use crate::{Error, Result};

use indicatif::{ProgressBar, ProgressStyle};
use instant::Instant;
use log::info;

/// Simulation state
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SimulationState {
    /// Initial state, not yet set up
    Created,
    /// Set up and ready to run
    Ready,
    /// Currently running
    Running,
    /// Completed
    Finished,
    /// Stopped early (by user or convergence)
    Stopped,
}

/// End condition for simulation
#[derive(Debug, Clone)]
pub enum EndCondition {
    /// Run for fixed number of timesteps
    Timesteps(u64),
    /// Run until energy decays below threshold (in dB relative to peak)
    EnergyDecay(f64),
    /// Run for fixed simulation time
    SimulationTime(f64),
}

impl Default for EndCondition {
    fn default() -> Self {
        Self::Timesteps(10000)
    }
}

/// Statistics from a simulation run
#[derive(Debug, Clone)]
pub struct SimulationStats {
    /// Total timesteps executed
    pub timesteps: u64,
    /// Total simulation time (seconds)
    pub sim_time: f64,
    /// Wall clock time (seconds)
    pub wall_time: f64,
    /// Peak energy during simulation
    pub peak_energy: f64,
    /// Final energy
    pub final_energy: f64,
    /// Average speed (cells/second)
    pub speed_mcells_per_sec: f64,
}

/// Main simulation controller.
pub struct Simulation {
    /// Grid definition
    grid: Grid,
    /// Boundary conditions
    boundaries: BoundaryConditions,
    /// Engine type selection
    engine_type: EngineType,
    /// Excitation sources
    excitations: Vec<Excitation>,
    /// End condition
    end_condition: EndCondition,
    /// Current state
    state: SimulationState,
    /// Operator (created during setup)
    operator: Option<Operator>,
    /// Engine (created during setup)
    engine: Option<Engine>,
    /// Verbosity level
    verbose: u8,
    /// Show progress bar
    show_progress: bool,
}

impl Simulation {
    /// Create a new simulation with the given grid.
    pub fn new(grid: Grid) -> Self {
        Self {
            grid,
            boundaries: BoundaryConditions::default(),
            engine_type: EngineType::default(),
            excitations: Vec::new(),
            end_condition: EndCondition::default(),
            state: SimulationState::Created,
            operator: None,
            engine: None,
            verbose: 1,
            show_progress: true,
        }
    }

    /// Set boundary conditions.
    pub fn set_boundaries(&mut self, boundaries: BoundaryConditions) -> &mut Self {
        self.boundaries = boundaries;
        self
    }

    /// Set engine type.
    pub fn set_engine_type(&mut self, engine_type: EngineType) -> &mut Self {
        self.engine_type = engine_type;
        self
    }

    /// Add an excitation source.
    pub fn add_excitation(&mut self, excitation: Excitation) -> &mut Self {
        self.excitations.push(excitation);
        self
    }

    /// Set end condition.
    pub fn set_end_condition(&mut self, condition: EndCondition) -> &mut Self {
        self.end_condition = condition;
        self
    }

    /// Set verbosity level (0=quiet, 1=normal, 2=verbose).
    pub fn set_verbose(&mut self, level: u8) -> &mut Self {
        self.verbose = level;
        self
    }

    /// Enable/disable progress bar.
    pub fn set_show_progress(&mut self, show: bool) -> &mut Self {
        self.show_progress = show;
        self
    }

    /// Get timestep information.
    pub fn timestep_info(&self) -> TimestepInfo {
        TimestepInfo::calculate(&self.grid)
    }

    /// Set up the simulation (create operator and engine).
    pub fn setup(&mut self) -> Result<()> {
        if self.state != SimulationState::Created {
            return Err(Error::Config("Simulation already set up".into()));
        }

        let info = self.timestep_info();

        if self.verbose >= 1 {
            info!(
                "FDTD simulation size: {}x{}x{} -> {} cells",
                self.grid.dimensions().nx,
                self.grid.dimensions().ny,
                self.grid.dimensions().nz,
                info.num_cells
            );
            info!(
                "FDTD timestep: {:.6e} s, Nyquist: {:.3e} Hz",
                info.dt, info.nyquist_freq
            );
            info!("Estimated memory: {}", info.memory_display());
        }

        // Create operator
        let operator = Operator::new(self.grid.clone(), self.boundaries.clone())?;

        // Create engine
        let engine = Engine::new(&operator, self.engine_type);

        self.operator = Some(operator);
        self.engine = Some(engine);
        self.state = SimulationState::Ready;

        Ok(())
    }

    /// Run the simulation.
    pub fn run(&mut self) -> Result<SimulationStats> {
        // Set up if not already done
        if self.state == SimulationState::Created {
            self.setup()?;
        }

        if self.state != SimulationState::Ready {
            return Err(Error::Config("Simulation not ready to run".into()));
        }

        self.state = SimulationState::Running;

        let operator = self.operator.as_ref().unwrap();
        let engine = self.engine.as_mut().unwrap();
        let dt = operator.timestep();

        // Determine number of timesteps
        let max_timesteps = match &self.end_condition {
            EndCondition::Timesteps(n) => *n,
            EndCondition::SimulationTime(t) => (t / dt).ceil() as u64,
            EndCondition::EnergyDecay(_) => 1_000_000, // Upper limit for energy decay
        };

        let energy_threshold = match &self.end_condition {
            EndCondition::EnergyDecay(db) => Some(*db),
            _ => None,
        };

        // Progress bar
        let progress = if self.show_progress {
            let pb = ProgressBar::new(max_timesteps);
            pb.set_style(
                ProgressStyle::default_bar()
                    .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} ({per_sec})")
                    .unwrap()
                    .progress_chars("##-"),
            );
            Some(pb)
        } else {
            None
        };

        let start_time = Instant::now();
        let mut peak_energy = 0.0f64;
        let mut timesteps_run = 0u64;

        // Main simulation loop
        while timesteps_run < max_timesteps {
            // Apply excitations
            let t = timesteps_run as f64 * dt;
            for exc in &self.excitations {
                if exc.is_active(t) {
                    let value = exc.evaluate(t) as f32;
                    let (i, j, k) = exc.position;
                    let field = engine.e_field_mut().component_mut(exc.direction);

                    if exc.soft_source {
                        field.add(i, j, k, value);
                    } else {
                        field.set(i, j, k, value);
                    }
                }
            }

            // Perform timestep
            engine.step(operator)?;
            timesteps_run += 1;

            // Check energy for convergence/decay
            if timesteps_run.is_multiple_of(100) {
                let energy = engine.total_energy(operator);
                if energy > peak_energy {
                    peak_energy = energy;
                }

                // Check energy decay condition
                if let Some(threshold_db) = energy_threshold {
                    if peak_energy > 0.0 {
                        let decay_db = 10.0 * (energy / peak_energy).log10();
                        if decay_db < threshold_db {
                            if self.verbose >= 1 {
                                info!(
                                    "Energy decay reached: {:.1} dB at timestep {}",
                                    decay_db, timesteps_run
                                );
                            }
                            break;
                        }
                    }
                }
            }

            // Update progress
            if let Some(ref pb) = progress {
                pb.set_position(timesteps_run);
            }
        }

        if let Some(pb) = progress {
            pb.finish_with_message("Simulation complete");
        }

        let wall_time = start_time.elapsed().as_secs_f64();
        let final_energy = engine.total_energy(operator);
        let num_cells = operator.dimensions().total();
        let speed = (timesteps_run as f64 * num_cells as f64) / wall_time / 1e6;

        self.state = SimulationState::Finished;

        let stats = SimulationStats {
            timesteps: timesteps_run,
            sim_time: timesteps_run as f64 * dt,
            wall_time,
            peak_energy,
            final_energy,
            speed_mcells_per_sec: speed,
        };

        if self.verbose >= 1 {
            info!(
                "Completed {} timesteps in {:.2}s ({:.2} MC/s)",
                stats.timesteps, stats.wall_time, stats.speed_mcells_per_sec
            );
        }

        Ok(stats)
    }

    /// Get reference to the engine (if set up).
    pub fn engine(&self) -> Option<&Engine> {
        self.engine.as_ref()
    }

    /// Get mutable reference to the engine (if set up).
    pub fn engine_mut(&mut self) -> Option<&mut Engine> {
        self.engine.as_mut()
    }

    /// Get reference to the operator (if set up).
    pub fn operator(&self) -> Option<&Operator> {
        self.operator.as_ref()
    }

    /// Get the current state.
    pub fn state(&self) -> SimulationState {
        self.state
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::CoordinateSystem;

    #[test]
    fn test_simulation_basic() {
        let grid = Grid::uniform(20, 20, 20, 1e-3);
        let mut sim = Simulation::new(grid);

        sim.set_end_condition(EndCondition::Timesteps(100))
            .set_verbose(0)
            .set_show_progress(false);

        // Add a simple excitation
        let exc = Excitation::gaussian(1e9, 0.5, 2, (10, 10, 10));
        sim.add_excitation(exc);

        let stats = sim.run().unwrap();
        assert_eq!(stats.timesteps, 100);
        assert!(stats.speed_mcells_per_sec > 0.0);
    }
}
