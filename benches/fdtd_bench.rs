//! Benchmarks for FDTD engine performance.

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};

use openems::arrays::Dimensions;
use openems::fdtd::{BoundaryConditions, Engine, EngineType, Operator};
use openems::geometry::Grid;

fn bench_fdtd_step(c: &mut Criterion) {
    let sizes = [(32, 32, 32), (64, 64, 64), (128, 128, 128)];

    for (nx, ny, nz) in sizes {
        let grid = Grid::uniform(nx, ny, nz, 1e-3);
        let op = Operator::new(grid.clone(), BoundaryConditions::all_pec()).unwrap();
        let total_cells = nx * ny * nz;

        let mut group = c.benchmark_group(format!("fdtd_{}x{}x{}", nx, ny, nz));
        group.throughput(Throughput::Elements(total_cells as u64));

        // Benchmark basic engine
        group.bench_function("basic", |b| {
            let mut engine = Engine::new(&op, EngineType::Basic);
            b.iter(|| {
                engine.step(&op).unwrap();
                black_box(())
            });
        });

        // Benchmark SIMD engine
        group.bench_function("simd", |b| {
            let mut engine = Engine::new(&op, EngineType::Simd);
            b.iter(|| {
                engine.step(&op).unwrap();
                black_box(())
            });
        });

        // Benchmark parallel engine
        group.bench_function("parallel", |b| {
            let mut engine = Engine::new(&op, EngineType::Parallel);
            b.iter(|| {
                engine.step(&op).unwrap();
                black_box(())
            });
        });

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
            black_box(())
        });
    });

    group.bench_function("simd_fmadd", |b| {
        let mut a = Field3D::new(dims);
        let b_field = Field3D::new(dims);
        a.fill(1.0);

        b.iter(|| {
            a.simd_fmadd(2.0, &b_field);
            black_box(())
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
