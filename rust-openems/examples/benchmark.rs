//! Performance benchmark for openEMS Rust implementation.

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

    // Warm up
    sim.setup().unwrap();

    let start = Instant::now();
    let stats = sim.run().unwrap();
    let elapsed = start.elapsed().as_secs_f64();

    let cells_per_step = (size * size * size) as f64;
    let speed = (timesteps as f64 * cells_per_step) / elapsed / 1e6;

    println!(
        "  {:?}: {} cells, {} steps in {:.3}s = {:.2} MC/s",
        engine_type,
        size * size * size,
        timesteps,
        elapsed,
        speed
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

        println!(
            "  SIMD speedup over Basic: {:.2}x",
            simd_speed / basic_speed
        );
        println!(
            "  Parallel speedup over Basic: {:.2}x",
            parallel_speed / basic_speed
        );
    }

    println!("\nDone!");
}
