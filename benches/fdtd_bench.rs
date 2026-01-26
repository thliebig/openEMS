//! Benchmarks for FDTD engine performance.

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};

use openems::arrays::Dimensions;
use openems::fdtd::{BoundaryConditions, Engine, EngineType, Operator};
use openems::geometry::Grid;

fn bench_fdtd_step(c: &mut Criterion) {
    // Sizes: Small, Medium, Large
    let sizes = [
        (50, 50, 50),
        (100, 100, 100),
        (200, 200, 200),
        (400, 400, 400),
    ];

    for (nx, ny, nz) in sizes {
        let grid = Grid::uniform(nx, ny, nz, 1e-3);
        let op = Operator::new(grid.clone(), BoundaryConditions::all_pec()).unwrap();
        let total_cells = nx * ny * nz;

        let mut group = c.benchmark_group(format!("fdtd_{}x{}x{}", nx, ny, nz));
        group.throughput(Throughput::Elements(total_cells as u64));
        group.sample_size(20); // Reduce sample size for slower/large benchmarks

        // Benchmark basic engine
        group.bench_function("basic", |b| {
            let mut engine = Engine::new(&op, EngineType::Basic);
            b.iter(|| {
                engine.step(&op).unwrap();
                // Use black_box to prevent optimization, but don't return the reference
                black_box(&engine);
            });
        });

        // Benchmark SIMD engine
        group.bench_function("simd", |b| {
            let mut engine = Engine::new(&op, EngineType::Simd);
            b.iter(|| {
                engine.step(&op).unwrap();
                black_box(&engine);
            });
        });

        // Benchmark parallel engine
        group.bench_function("parallel", |b| {
            let mut engine = Engine::new(&op, EngineType::Parallel);
            b.iter(|| {
                engine.step(&op).unwrap();
                black_box(&engine);
            });
        });

        // Benchmark GPU engine (if available)
        // We wrap this in a closure that might panic if no GPU is found,
        // but typically Engine::new will succeed if there's a fallback or actual GPU.
        // If it panics, the benchmark harness catches it.
        // Note: Engine::new for GPU might take time, so we do it outside iter.
        // If the machine has no GPU support, this might fail, but that's expected.
        let has_gpu = std::panic::catch_unwind(|| {
            let _ = Engine::new(&op, EngineType::Gpu);
        })
        .is_ok();

        if has_gpu {
            group.bench_function("gpu", |b| {
                let mut engine = Engine::new(&op, EngineType::Gpu);
                b.iter(|| {
                    engine.step(&op).unwrap();
                    // Crucial: Wait for GPU to finish!
                    engine.wait_idle();
                    black_box(&engine);
                });
            });
        }

        group.finish();
    }
}

fn bench_array_operations(c: &mut Criterion) {
    use openems::arrays::{Field3D, SimdOps};

    let dims = Dimensions::new(100, 100, 100);
    let total = dims.total();

    let mut group = c.benchmark_group("array_ops");
    group.throughput(Throughput::Elements(total as u64));

    group.bench_function("fill", |b| {
        let mut field = Field3D::new(dims);
        b.iter(|| {
            field.fill(1.0);
            black_box(&field);
        });
    });

    group.bench_function("simd_fmadd", |b| {
        let mut a = Field3D::new(dims);
        let b_field = Field3D::new(dims);
        a.fill(1.0);

        b.iter(|| {
            a.simd_fmadd(2.0, &b_field);
            black_box(&a);
        });
    });

    group.bench_function("energy", |b| {
        let field = Field3D::new(dims);
        b.iter(|| black_box(field.energy()));
    });

    group.finish();
}

criterion_group!(benches, bench_fdtd_step, bench_array_operations);
criterion_main!(benches);
