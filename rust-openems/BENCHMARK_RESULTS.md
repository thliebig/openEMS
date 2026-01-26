# FDTD Engine Benchmarks

Benchmark run on: aarch64-apple-darwin (Apple Silicon M2 Pro)

## Results Summary (Optimized with Thread Coarsening ILP=4)

| Grid Size | Engine | Time per Step | Throughput (Cells/s) |
|---|---|---|---|
| **50x50x50** | Basic (Ref) | 2.01 ms | ~62 M |
| (125k cells) | SIMD | 0.81 ms | ~154 M |
| | Parallel | 0.38 ms | ~329 M |
| | **GPU** | **1.38 ms** | **~91 M** |
| | | | |
| **100x100x100** | Basic (Ref) | 23.6 ms | ~42 M |
| (1M cells) | SIMD | 6.3 ms | ~159 M |
| | Parallel | 4.26 ms | ~235 M |
| | **GPU** | **1.39 ms** | **~718 M** |
| | | | |
| **200x200x200** | Basic (Ref) | 175.5 ms | ~45 M |
| (8M cells) | SIMD | 49.8 ms | ~160 M |
| | Parallel | 32.0 ms | ~250 M |
| | **GPU** | **5.27 ms** | **~1.52 G** |

## Analysis

1.  **Massive Speedup:** Implementing **Thread Coarsening (ILP=4)** resulted in a nearly **2x speedup** for large grids compared to the previous GPU optimization (9.6ms -> 5.27ms).
2.  **Comparison vs CPU:**
    - At 100^3, the GPU is **~3x faster** than the parallel CPU engine.
    - At 200^3, the GPU is **~6x faster** than the parallel CPU engine.
3.  **Throughput:** The GPU engine now achieves **~1.5 Gigacells/second**, which is a high-performance result for an integrated GPU like the M2 Pro, effectively saturating a significant portion of the available bandwidth with useful work (given the register reuse optimization).
4.  **Register Reuse:** By processing 4 elements per thread, we reduce the number of neighbor loads for `Ey` and `Hy` significantly (register caching), which directly translates to the observed performance gain.
