//! Performance benchmark for openEMS Rust implementation.
//!
//! Compares different engine implementations:
//! - Basic: Single-threaded reference
//! - SIMD: Single-threaded with SIMD acceleration
//! - Parallel: Multi-threaded SIMD
//! - Compressed: Memory-bandwidth optimized with coefficient compression

use openems::fdtd::{EndCondition, EngineType, Simulation};
use openems::geometry::Grid;
use std::time::Instant;

fn benchmark_engine(size: usize, timesteps: u64, engine_type: EngineType) -> f64 {
    let grid = Grid::uniform(size, size, size, 1e-3);
    let mut sim = Simulation::new(grid);

    sim.set_engine_type(engine_type)
        .set_end_condition(EndCondition::Timesteps(timesteps))
        .set_verbose(0)
        .set_show_progress(false);

    // Setup and warm up
    sim.setup().unwrap();

    let start = Instant::now();
    let _stats = sim.run().unwrap();
    let elapsed = start.elapsed().as_secs_f64();

    let cells_per_step = (size * size * size) as f64;
    let speed = (timesteps as f64 * cells_per_step) / elapsed / 1e6;

    // Print compression info for compressed engine
    let compression_info = if engine_type == EngineType::Compressed {
        if let Some(stats) = sim.engine_wrapper().and_then(|e| e.compression_stats()) {
            format!(" (E ratio: {:.2}, H ratio: {:.2})", stats.0, stats.1)
        } else {
            String::new()
        }
    } else {
        String::new()
    };

    println!(
        "  {:?}: {} cells, {} steps in {:.3}s = {:.2} MC/s{}",
        engine_type,
        size * size * size,
        timesteps,
        elapsed,
        speed,
        compression_info
    );

    speed
}

fn main() {
    println!("openEMS Rust Performance Benchmark");
    println!("===================================\n");

    // Warm up the CPU
    let _ = benchmark_engine(20, 10, EngineType::Parallel);

    // Run benchmarks for different grid sizes
    for &size in &[50, 100, 150] {
        println!(
            "\nGrid size: {}x{}x{} = {} cells",
            size,
            size,
            size,
            size * size * size
        );
        println!("-----------------------------------------");

        let timesteps = 100u64;

        let basic_speed = benchmark_engine(size, timesteps, EngineType::Basic);
        let simd_speed = benchmark_engine(size, timesteps, EngineType::Simd);
        let parallel_speed = benchmark_engine(size, timesteps, EngineType::Parallel);
        let compressed_speed = benchmark_engine(size, timesteps, EngineType::Compressed);

        println!();
        println!(
            "  SIMD speedup over Basic: {:.2}x",
            simd_speed / basic_speed
        );
        println!(
            "  Parallel speedup over Basic: {:.2}x",
            parallel_speed / basic_speed
        );
        println!(
            "  Compressed speedup over Basic: {:.2}x",
            compressed_speed / basic_speed
        );
        println!(
            "  Compressed speedup over Parallel: {:.2}x",
            compressed_speed / parallel_speed
        );
    }

    println!("\n===========================================");
    println!("Memory Bandwidth Analysis");
    println!("===========================================\n");

    // Test with a large grid to show memory bandwidth benefits
    let size = 200;
    let timesteps = 50u64;

    println!(
        "Large grid test: {}x{}x{} = {} cells ({:.1} MB field data)",
        size,
        size,
        size,
        size * size * size,
        (size * size * size * 6 * 4) as f64 / 1e6 // 6 field components * f32
    );
    println!("-----------------------------------------");

    let parallel_speed = benchmark_engine(size, timesteps, EngineType::Parallel);
    let compressed_speed = benchmark_engine(size, timesteps, EngineType::Compressed);

    println!();
    println!(
        "  Compressed vs Parallel: {:.2}x",
        compressed_speed / parallel_speed
    );

    println!("\nDone!");
}
