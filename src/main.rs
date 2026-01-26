//! openEMS command-line interface.
//!
//! High-performance FDTD electromagnetic field solver.

use clap::Parser;
use log::{error, info, warn, LevelFilter};
use std::path::PathBuf;
use std::process::ExitCode;

use openems::fdtd::{EndCondition, EngineType, Simulation};
use openems::geometry::Grid;
use openems::io::xml;

/// openEMS - High-performance FDTD electromagnetic field solver
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input XML configuration file
    #[arg(required = true)]
    input: PathBuf,

    /// Disable all field dumps for faster simulation
    #[arg(long)]
    disable_dumps: bool,

    /// Dump material distribution to VTK file for debugging
    #[arg(long)]
    debug_material: bool,

    /// Dump operator to VTK file for debugging
    #[arg(long)]
    debug_operator: bool,

    /// Engine type: basic, simd, parallel (default: parallel)
    #[arg(long, default_value = "parallel")]
    engine: String,

    /// Force number of threads for parallel engine
    #[arg(long, default_value = "0")]
    num_threads: usize,

    /// Only run preprocessing; do not simulate
    #[arg(long)]
    no_simulation: bool,

    /// Dump simulation statistics
    #[arg(long)]
    dump_statistics: bool,

    /// Verbose level (0=quiet, 1=normal, 2=verbose, 3=debug)
    #[arg(short, long, default_value = "1")]
    verbose: u8,
}

fn main() -> ExitCode {
    let args = Args::parse();

    // Initialize logging
    let log_level = match args.verbose {
        0 => LevelFilter::Warn,
        1 => LevelFilter::Info,
        2 => LevelFilter::Debug,
        _ => LevelFilter::Trace,
    };

    env_logger::Builder::new()
        .filter_level(log_level)
        .format_timestamp(None)
        .init();

    // Print banner
    println!(" ---------------------------------------------------------------------- ");
    println!(
        " | openEMS-Rust {} -- High-Performance FDTD Solver",
        openems::VERSION
    );
    println!(" | (c) 2024 openEMS contributors  GPL license");
    println!(" ---------------------------------------------------------------------- ");

    // Set number of threads if specified
    if args.num_threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.num_threads)
            .build_global()
            .unwrap();
        info!("Using {} threads", args.num_threads);
    }

    // Parse engine type
    let engine_type = match args.engine.as_str() {
        "basic" => EngineType::Basic,
        "simd" => EngineType::Simd,
        "parallel" => EngineType::Parallel,
        _ => {
            warn!("Unknown engine type '{}', using parallel", args.engine);
            EngineType::Parallel
        }
    };

    info!("Engine type: {:?}", engine_type);

    // Parse input file
    if !args.input.exists() {
        error!("Input file not found: {:?}", args.input);
        return ExitCode::FAILURE;
    }

    info!("Reading configuration from: {:?}", args.input);

    let config = match xml::parse_config(&args.input) {
        Ok(c) => c,
        Err(e) => {
            error!("Failed to parse configuration: {}", e);
            return ExitCode::FAILURE;
        }
    };

    // Create grid from configuration
    let grid = if let Some(ref grid_config) = config.grid {
        Grid::cartesian(
            grid_config.x_lines.clone(),
            grid_config.y_lines.clone(),
            grid_config.z_lines.clone(),
        )
    } else {
        // Default grid for testing
        warn!("No grid specified, using default 20x20x20 grid");
        Grid::uniform(20, 20, 20, 1e-3)
    };

    // Create simulation
    let mut sim = Simulation::new(grid);
    sim.set_engine_type(engine_type)
        .set_verbose(args.verbose)
        .set_show_progress(args.verbose > 0);

    // Set end condition
    if let Some(end_db) = config.end_criteria_db {
        sim.set_end_condition(EndCondition::EnergyDecay(end_db));
    } else if config.num_timesteps > 0 {
        sim.set_end_condition(EndCondition::Timesteps(config.num_timesteps));
    }

    // Setup simulation
    if let Err(e) = sim.setup() {
        error!("Failed to setup simulation: {}", e);
        return ExitCode::FAILURE;
    }

    if args.no_simulation {
        info!("Preprocessing complete (--no-simulation specified)");
        return ExitCode::SUCCESS;
    }

    // Run simulation
    match sim.run() {
        Ok(stats) => {
            info!(
                "Simulation completed: {} timesteps in {:.2}s ({:.2} MC/s)",
                stats.timesteps, stats.wall_time, stats.speed_mcells_per_sec
            );

            if args.dump_statistics {
                // Write statistics file
                let stats_path = args.input.with_extension("stats.json");
                let stats_json = serde_json::json!({
                    "timesteps": stats.timesteps,
                    "sim_time": stats.sim_time,
                    "wall_time": stats.wall_time,
                    "peak_energy": stats.peak_energy,
                    "final_energy": stats.final_energy,
                    "speed_mcells_per_sec": stats.speed_mcells_per_sec,
                });

                if let Ok(json_str) = serde_json::to_string_pretty(&stats_json) {
                    if std::fs::write(&stats_path, json_str).is_ok() {
                        info!("Statistics written to {:?}", stats_path);
                    }
                }
            }

            ExitCode::SUCCESS
        }
        Err(e) => {
            error!("Simulation failed: {}", e);
            ExitCode::FAILURE
        }
    }
}
