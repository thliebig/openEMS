# openEMS - Open Electromagnetic Field Solver

A high-performance FDTD (Finite-Difference Time-Domain) electromagnetic field solver written in Rust.

## Features

- **FDTD Engine**: Full 3D electromagnetic field simulation using the FDTD method
- **SIMD Acceleration**: AVX2/AVX-512 optimized field updates for maximum performance
- **Parallel Processing**: Multi-threaded simulation using Rayon
- **Boundary Conditions**:
  - Perfect Electric Conductor (PEC)
  - Perfectly Matched Layer (PML/UPML)
  - Mur Absorbing Boundary Conditions (1st and 2nd order)
  - Local Absorbing Boundary Conditions
- **Excitation Sources**:
  - Gaussian pulse
  - Sinusoidal (modulated)
  - Dirac impulse
  - Step function
  - Custom waveforms
- **Materials**:
  - Dielectrics (lossy and lossless)
  - Dispersive materials (Drude, Lorentz, Debye models)
  - Conducting sheets
  - Lumped RLC elements
- **Probes and Processing**:
  - Voltage and current probes
  - Field probes (E and H)
  - Near-field to far-field transformation (NF2FF)
  - SAR (Specific Absorption Rate) calculation
  - S-parameter extraction
- **Coordinate Systems**:
  - Cartesian
  - Cylindrical (with optional multi-grid support)

## Building

```bash
# Build the library
cargo build --release

# Run tests (289 tests)
cargo test

# Run benchmarks
cargo bench
```

## Usage

```rust
use openems::fdtd::{Simulation, BoundaryConditions};
use openems::geometry::Grid;

// Create a uniform grid (10x10x10 cells, 1mm cell size)
let grid = Grid::uniform(10, 10, 10, 0.001);

// Set up boundary conditions
let boundaries = BoundaryConditions::default();

// Create and run simulation
let mut sim = Simulation::new(grid, boundaries)?;
sim.run(1000)?; // Run for 1000 timesteps
```

## Test Coverage

The codebase maintains high test coverage on core modules:

| Module | Coverage |
|--------|----------|
| `operator.rs` | 100% |
| `local_absorbing_bc.rs` | 99% |
| `geometry/mod.rs` | 99% |
| `excitation.rs` | 98% |
| `engine_interface.rs` | 97% |
| `probes.rs` | 97% |
| **Overall** | **85%** |

## Python Bindings

Python bindings are available via PyO3/Maturin:

```bash
# Build Python wheel
maturin build --release

# Install in development mode
maturin develop
```

## License

This project is licensed under the GNU General Public License v3.0 - see the [COPYING](COPYING) file for details.

## Acknowledgments

This is a Rust port of the original [openEMS](https://openems.de) C++ project, maintaining API compatibility while providing improved performance and memory safety.
