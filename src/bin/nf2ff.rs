//! NF2FF (Near-Field to Far-Field) command-line tool.
//!
//! Computes far-field radiation patterns from near-field simulation data.

use clap::Parser;
use log::{error, info, LevelFilter};
use std::path::PathBuf;
use std::process::ExitCode;

use openems::nf2ff::Nf2ff;

/// NF2FF - Near-Field to Far-Field transformation tool
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input HDF5 file with near-field data
    #[arg(required = true)]
    input: PathBuf,

    /// Output file for far-field results
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Frequency for calculation (Hz)
    #[arg(short, long)]
    frequency: f64,

    /// Number of theta angles
    #[arg(long, default_value = "181")]
    n_theta: usize,

    /// Number of phi angles
    #[arg(long, default_value = "361")]
    n_phi: usize,

    /// Verbose level
    #[arg(short, long, default_value = "1")]
    verbose: u8,
}

fn main() -> ExitCode {
    let args = Args::parse();

    // Initialize logging
    let log_level = match args.verbose {
        0 => LevelFilter::Warn,
        1 => LevelFilter::Info,
        _ => LevelFilter::Debug,
    };

    env_logger::Builder::new()
        .filter_level(log_level)
        .format_timestamp(None)
        .init();

    // Print banner
    println!(" NF2FF - Near-Field to Far-Field Transformation");
    println!(" openEMS-Rust {}", openems::VERSION);
    println!();

    // Check input file
    if !args.input.exists() {
        error!("Input file not found: {:?}", args.input);
        return ExitCode::FAILURE;
    }

    info!("Reading near-field data from: {:?}", args.input);
    info!("Frequency: {:.3e} Hz", args.frequency);

    // Create NF2FF calculator
    let nf2ff = Nf2ff::new(args.frequency);

    // TODO: Load near-field data from HDF5 file
    // For now, just compute empty result

    // Calculate far-field
    let theta: Vec<f64> = (0..args.n_theta)
        .map(|i| std::f64::consts::PI * i as f64 / (args.n_theta - 1) as f64)
        .collect();

    let phi: Vec<f64> = (0..args.n_phi)
        .map(|i| 2.0 * std::f64::consts::PI * i as f64 / (args.n_phi - 1) as f64)
        .collect();

    info!("Computing far-field at {} x {} angles...", args.n_theta, args.n_phi);

    let result = nf2ff.calculate(&theta, &phi);

    info!("Maximum directivity: {:.2} dB", 10.0 * result.max_directivity().log10());

    // Write output
    let output_path = args.output.unwrap_or_else(|| {
        args.input.with_extension("ff.json")
    });

    info!("Writing results to: {:?}", output_path);

    // Serialize result
    let output_data = serde_json::json!({
        "frequency": result.frequency,
        "theta": result.theta,
        "phi": result.phi,
        "directivity": result.directivity(),
    });

    match serde_json::to_string_pretty(&output_data) {
        Ok(json) => {
            if std::fs::write(&output_path, json).is_err() {
                error!("Failed to write output file");
                return ExitCode::FAILURE;
            }
        }
        Err(e) => {
            error!("Failed to serialize results: {}", e);
            return ExitCode::FAILURE;
        }
    }

    info!("NF2FF calculation complete");
    ExitCode::SUCCESS
}
