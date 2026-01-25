# openEMS Rust Implementation - Developer Notes

## Project Overview

This is a complete Rust rewrite of the openEMS FDTD electromagnetic field solver.
The goal is to match or exceed C++ performance while providing better safety,
cross-platform support, and static linking.

## Key Documentation

- **PORTING_STATUS.md** - Tracks all C++ modules and their Rust porting status
- **LANGUAGE_DECISION.md** - Documents why Rust was chosen over Zig

## Architecture

### Core Modules
- `src/fdtd/` - FDTD engine, operator, simulation control
- `src/arrays/` - SIMD-optimized field arrays
- `src/geometry/` - Grid and coordinate systems
- `src/extensions/` - Boundary conditions, materials, excitations

### Extensions System
Extensions use pre/post hooks around the main engine updates:
1. `pre_update_h` - Before H-field update
2. `post_update_h` - After H-field update
3. `pre_update_e` - Before E-field update
4. `post_update_e` - After E-field update

### Performance Targets
- Engine should achieve 700+ MC/s (mega-cells per second)
- SIMD speedup: ~2x over scalar
- Parallel speedup: ~10x on 8+ core machines

## Testing

Run all tests:
```bash
cargo test --all-features
```

Run benchmarks:
```bash
cargo run --example benchmark --release
```

## Building

Debug build:
```bash
cargo build
```

Release build with all optimizations:
```bash
cargo build --release
```

Cross-platform builds (requires cross-rs):
```bash
./scripts/build-all.sh
```

## C++ Reference

When porting, reference the original C++ files:
- Extensions: `FDTD/extensions/`
- Processing: `Common/`
- Tools: `tools/`

See PORTING_STATUS.md for detailed file mappings.
