# Rust vs Zig Decision for openEMS Rewrite

## Decision: **Rust**

## Rationale

### Requirements Analysis

The openEMS rewrite has these key requirements:
1. **Python bindings** - Must maintain Python interface
2. **SIMD optimizations** - Runtime CPU feature detection for vectorization
3. **Static linking** - Final binary should be self-contained
4. **Cross-platform** - Windows, macOS, Linux
5. **Performance** - Must significantly beat existing C++ code
6. **No external dependencies** - Avoid cloning/building external repos

### Comparison Matrix

| Criterion | Rust | Zig | Winner |
|-----------|------|-----|--------|
| **Python Bindings** | PyO3 - mature, well-documented, automatic type conversions | Manual FFI through C ABI | **Rust** |
| **SIMD** | `std::simd` (nightly) or `simdeez`/`wide` crates with runtime detection | Built-in `@Vector` with runtime CPU detection | Tie |
| **HDF5 Support** | `hdf5` crate - mature bindings | Requires manual C FFI wrapping | **Rust** |
| **C++ Interop** | `cxx` crate - good but requires binding generation | Native C interop (no C++ directly) | Tie |
| **Static Linking** | Well-supported via `rustflags` | Native support | Tie |
| **Cross-Platform** | Excellent (cargo handles platform specifics) | Good (unified build system) | **Rust** |
| **Ecosystem/Libraries** | Large ecosystem (nalgebra, ndarray, rayon) | Smaller ecosystem | **Rust** |
| **Build System** | Cargo - mature, handles dependencies | Zig build - simpler, built-in | Tie |
| **Memory Safety** | Compile-time guaranteed | Runtime safety, manual control | **Rust** |
| **Performance Potential** | Competitive with C++ | Competitive with C++ | Tie |
| **Learning Resources** | Extensive documentation and community | Growing but smaller | **Rust** |

### Critical Factors

#### 1. Python Bindings (Critical - Rust wins decisively)
- **PyO3** provides seamless Python-Rust interop with automatic type conversions
- Can expose Rust structs directly to Python with `#[pyclass]`
- Supports async, numpy arrays via `numpy` crate
- Battle-tested in production (used by Polars, Ruff, cryptography)

For Zig:
- Would require manual C FFI layer
- No automatic type conversions
- Much more boilerplate code

#### 2. SIMD with Runtime Detection
Both languages can achieve this:

**Rust approach:**
```rust
use std::arch::x86_64::*;

fn compute() {
    if is_x86_feature_detected!("avx2") {
        unsafe { compute_avx2() }
    } else if is_x86_feature_detected!("sse4.1") {
        unsafe { compute_sse4() }
    } else {
        compute_scalar()
    }
}
```

Or use `wide` crate for portable SIMD with automatic dispatch.

**Zig approach:**
```zig
const builtin = @import("builtin");
const cpu_features = std.Target.x86.Feature;

pub fn compute() void {
    if (std.Target.current.cpu.features.contains(.avx2)) {
        compute_avx2();
    } else if (std.Target.current.cpu.features.contains(.sse4_1)) {
        compute_sse4();
    } else {
        compute_scalar();
    }
}
```

Both work well, but Rust has more mature libraries for this.

#### 3. Scientific Computing Ecosystem
- **Rust**: `nalgebra`, `ndarray`, `num-complex`, `fftw` bindings
- **Zig**: Would need to write more from scratch or wrap C libraries

#### 4. File Format Support
- **HDF5**: Rust has `hdf5` crate; Zig would need manual bindings
- **VTK**: Neither has great support, but Rust's `xml` ecosystem helps

### Decision Justification

**Rust is chosen because:**

1. **PyO3 is a game-changer** - Writing Python bindings in Zig would be 5-10x more work and error-prone

2. **Ecosystem maturity** - The `wide` crate provides excellent portable SIMD, `ndarray` provides numpy-like arrays, `rayon` provides effortless parallelism

3. **No external dependencies goal** - Cargo handles vendoring dependencies automatically, making static builds trivial

4. **Cross-platform reliability** - Rust's tier-1 support for all three major platforms is excellent

5. **Performance equivalence** - Both Rust and Zig can match C++ performance; Rust's abstractions don't add runtime cost

### Rust Libraries to Use

| Purpose | Library |
|---------|---------|
| SIMD | `wide` (portable SIMD with auto-dispatch) |
| Arrays | `ndarray` |
| Linear Algebra | `nalgebra` |
| Parallelism | `rayon` |
| HDF5 | `hdf5` |
| XML | `quick-xml` |
| Python bindings | `pyo3` + `numpy` |
| Complex numbers | `num-complex` |
| FFT | `rustfft` |

### Migration Strategy

1. Start with `tools/arraylib` - pure data structures, easy to port
2. Port utility modules (`tools/`)
3. Rewrite core FDTD engine with SIMD optimizations
4. Port extensions one by one
5. Port processing/analysis code
6. Create Python bindings
7. Remove C++ code, ensure static linking

### Performance Goals

The new Rust implementation should:
- Use AVX2/AVX-512 where available (via `wide` crate)
- Use rayon for parallel field updates
- Minimize memory allocations in hot loops
- Use cache-friendly data layouts (SoA where beneficial)
- Single optimized engine (no multiple implementations)

### Not Choosing Zig Because

While Zig has excellent properties:
- Better C interop (not C++)
- Simpler mental model
- Built-in SIMD with runtime detection

The Python bindings requirement makes Zig impractical. The effort to manually create and maintain Python bindings through C FFI would:
- Add thousands of lines of boilerplate
- Be error-prone
- Slow down development significantly
- Make the codebase harder to maintain

---

**Final Decision: Rust**

Date: 2026-01-25
