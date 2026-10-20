# HIP UPML allocation control (Linux/NVIDIA)

This optional synthetic fixture compares peak requested `cudaMalloc` bytes
before and after the UPML storage change. It also writes E/H probes and checks
that they are finite and nonzero. It is independent of antenna/material data.

Build and install each revision and its Python bindings into separate prefixes,
using the same dependencies and Release options. For example, in addition to
the usual dependency paths:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DENABLE_HIP=ON -DOPENEMS_HIP_AS_CUDA=ON \
  -DCMAKE_CUDA_RUNTIME_LIBRARY=Shared
cmake --build build --parallel 4
cmake --install build
```

Select the same GPU for both builds, using the backend's normal device selection.
The shared CUDA runtime is required for interception; this does not change the
solver sources. Compile the test-only tracer with the host compiler:

```sh
c++ -std=c++17 -shared -fPIC -O2 trace.cpp -o trace.so -ldl -pthread
```

With each build's Python environment and library path active, run each mode in a
fresh process and a new output directory (`before`/`after` identify the revision):

```sh
LD_PRELOAD="$PWD/trace.so" python pulse.py before-fused --mode fused > before-fused.log 2>&1
```

Repeat for `regions`, `separate`, and the no-PML `control`. If necessary, set
`OPENEMS_CUDART_LIBRARY` to the shared runtime's path or soname (e.g.
`libcudart.so.12`). Check each log for the HIP backend, all extensions on device,
the -60 dB energy criterion, and a final `CUDA_ALLOCATION_AUDIT` line. Require
positive `calls` and zero `remaining_bytes`; zero calls is an invalid trace.
The tracer is for this single-process fixture (one host allocation thread),
not a general CUDA profiler.

Compare the `et`/`ht` probe files with `openEMS.compare` (tolerance 1e-4 of peak)
and require matching stop steps. This 81³-node fixture has 281,394 PML locations;
omitting 18 FP32 coefficients per location should reduce the allocation peak by
20,260,368 bytes (19.32 MiB) in every PML mode, with no change in `control`.
This is requested buffer memory, excluding driver/context overhead. Other
solver revisions may change unrelated buffers; inspect those before changing
the expected difference. The usual GPU engine regression cases additionally
exercise full-array fallback and 16/32-bit UPML set indices.

Do not use interposition for timing comparisons. This fixture measures storage
and decay, not general antenna throughput or AMD/Metal compatibility.
