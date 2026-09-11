# Metal UPML boundary conditioning

`--engine=metal` now runs both voltage and current UPML pre/post conditioning
on the GPU. Operator construction and grading remain on the CPU. The kernels
use the existing UPML coefficients and persistent flux arrays through shared,
no-copy Metal buffers; no per-step field or coefficient copies are needed.

Each scalar UPML cell maps into the engine's SSE-compatible z-lane layout.
Only physical cells are touched, including boundary lines and corner regions;
SIMD padding is excluded. Extension execution keeps the original reverse-pre /
forward-post priority order. Separate compute encoders preserve dependencies,
including overlapping slabs. Pre-conditioning, field updates and
post-conditioning normally share one command submission per half-step. Pending
work completes before any CPU extension hook and before sources, probes or dumps
can read the fields.

Metal objects use ARC, including command buffers retained across autorelease
pools. The source flags are set in the CMake directory that creates the library
target (the previous subdirectory-only setting did not apply to that target).

## Controls

- Default: GPU UPML; startup prints `Metal: GPU UPML conditioning: N regions`.
- `OPENEMS_METAL_PML=0`: retain CPU UPML with Metal field updates for comparison.
- `OPENEMS_METAL_COMPRESS=0`: retain dense field-update coefficients, independently
  of UPML. UPML coefficients themselves are not dictionary-compressed.
- `OPENEMS_METAL_FP64_REFERENCE=1`: with UPML, checks each field stencil in FP64
  from the actual conditioned input. It does **not** simulate FP64 flux evolution.
  Evolving the old additive reference through flux swaps is invalid and can
  diverge. Without UPML the existing evolving diagnostic remains unchanged.

Other engines and non-UPML extensions retain their CPU implementations. This does
not add cylindrical or MPI support to the Cartesian Metal engine.

## Regression tests

Use a Python environment with the openEMS/CSXCAD bindings, NumPy and h5py:

```sh
python TESTSUITE/enginetests/metal_pml.py --openems /absolute/path/to/openEMS
python TESTSUITE/enginetests/metal_pml.py --openems /absolute/path/to/openEMS \
  --fp64-reference
MTL_DEBUG_LAYER=1 MTL_SHADER_VALIDATION=1 \
  python TESTSUITE/enginetests/metal_pml.py --openems /absolute/path/to/openEMS
```

The 12 cases cover each face separately, six-face edges/corners, overlapping
opposing slabs, unequal PML thicknesses, nonuniform spacing, mixed Mur/PEC/PML,
no PML, all four physical
z-line-count remainders, lossy dielectric, and a 1200-step post-pulse run.
Complete E/H dumps are compared between GPU UPML, CPU UPML on Metal, and SSE;
non-finite/missing/zero field data fails. Dense/compressed Metal outputs must be
bit-identical. CPU/GPU UPML uses pointwise and relative-L2 checks. SSE uses a
relative-L2 limit of 1e-4 and reports pointwise discrepancies separately because
existing Metal/SSE cancellation differences also occur without GPU UPML.

Validated on M4 Pro, Release:

- All 12 cases pass with fast math OFF (normal batching, Metal shader validation,
  and FP64 diagnostic modes) and fast math ON (normal batching).
- CPU/GPU UPML dumps had zero numeric differences in these cases;
  dense/compressed dumps were bit-identical.
- Fast-math-OFF SSE comparison: maximum E relative L2 1.15e-6, H 1.05e-5.
- Existing coefficient suite retains its known SSE long-run discrepancy
  (3895 E pointwise failures); dense/compressed output remains bit-identical.
- Existing PEC mapping suite passes. Builds with `WITH_METAL=OFF` and ON pass.

## Performance spot check

M4 Pro, Release, fast math OFF, compressed field coefficients, no FP64 diagnostic
or field dumps; uniform 128 x 127 x 126 cells, six `PML_8` faces, lossy dielectric,
600 steps. Three runs per variant, alternating order, medians:

| UPML execution | Full process | Stepping |
|---|---:|---:|
| CPU (Metal fields) | 18.14 s | 16.00 s |
| GPU (Metal fields) | 3.94 s | 1.81 s |

About 4.6x full-run and 8.9x stepping speedup on this case. Small grids can be
submission-bound; this is not a general speedup guarantee. Diagnostic FP64 mode
adds CPU work and extra synchronizations and is not suitable for benchmarking.
