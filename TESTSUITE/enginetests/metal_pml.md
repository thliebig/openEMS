# Metal UPML boundary conditioning

`--engine=metal` now runs both voltage and current UPML pre/post conditioning
on the GPU. Operator construction and grading remain on the CPU. The default
kernels use an indexed layout: existing UPML coefficients and initial
fluxes are reordered once into increasing SSE-compatible field addresses. A
persistent uint32 index buffer replaces per-thread coordinate division/modulo.
Adjacent threads traverse contiguous auxiliary data and monotonically increasing
field addresses, instead of sweeping the field buffer separately per component.
No per-step field or coefficient copies are needed.

Only physical cells are included in the mapping, including boundary lines and
corner regions; SIMD padding is excluded. CPU UPML arrays remain allocated for
the original implementation. Indexed storage adds 8 float arrays plus one uint32
index per scalar component (108 bytes per PML grid position; about 887 MB for the
14.42M-cell thin PCB). Runtime switching between layouts is not supported.
Extension execution keeps the original reverse-pre /
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
- `OPENEMS_METAL_PML_LAYOUT=scalar`: original no-copy scalar GPU UPML for
  lower-memory operation and A/B tests. Default is indexed; startup prints the
  selected layout. GPU allocation failure stops the run rather than silently
  selecting a different layout.
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

The 13 cases cover each face separately, six-face edges/corners, overlapping
opposing slabs, unequal PML thicknesses, nonuniform spacing, mixed Mur/PEC/PML,
no PML, all four physical z-line-count remainders, lossy dielectric, a 1200-step
post-pulse run, and a thin nonuniform board-shaped grid. Complete E/H dumps are
compared between GPU UPML, CPU UPML on Metal, and SSE; non-finite/missing/zero
field data fails. Dense/compressed Metal outputs and scalar/indexed GPU UPML
outputs must be bit-identical. CPU/GPU UPML uses pointwise and relative-L2 checks. SSE uses a
relative-L2 limit of 1e-4 and reports pointwise discrepancies separately because
existing Metal/SSE cancellation differences also occur without GPU UPML.

Indexed-layout validation on M4 Pro, Release, fast math OFF:

- All 13 cases pass normally, with Metal API/shader validation, and with the
  FP64 diagnostic. Scalar/indexed complete E/H datasets are bit-identical.
- PEC mapping suite passes; coefficient suite retains the existing 3895 SSE
  E pointwise failures, with dense/compressed Metal dumps still identical.
- Indexed layout has not been retested with fast math ON.

Historical scalar-layout validation on M4 Pro, Release:

- All 12 cases pass with fast math OFF (normal batching, Metal shader validation,
  and FP64 diagnostic modes) and fast math ON (normal batching).
- CPU/GPU UPML dumps had zero numeric differences in these cases;
  dense/compressed dumps were bit-identical.
- Fast-math-OFF SSE comparison: maximum E relative L2 1.15e-6, H 1.05e-5.
- Existing coefficient suite retains its known SSE long-run discrepancy
  (3895 E pointwise failures); dense/compressed output remains bit-identical.
- Existing PEC mapping suite passes. Builds with `WITH_METAL=OFF` and ON pass.

## Indexed-layout benchmark

M4 Pro, Release, fast math OFF; enlarged PCB, 658 x 664 x 33 = 14,418,096
solver cells, six PML_8 faces, dense field coefficients. Identical XML, 788
steps, order indexed/scalar/scalar/indexed (two runs per layout):

| Layout | Stepping | Throughput | Whole-process median |
|---|---:|---:|---:|
| Original scalar | 39.17–39.71 s | 286.1–290.1 MC/s | 139.06 s |
| Indexed | 21.98–22.06 s | 515.0–516.9 MC/s | 117.71 s |

Approximately **1.79x stepping**, **1.18x whole-process** speedup. All numeric
port probe outputs are bit-identical for these finite runs; 788 steps is before
the Gaussian pulse peak, not full waveform/convergence or SI validation. The
13-fixture suite separately exercises nonzero fields and post-pulse decay.

Specializing integer divisors without reordering arrays produced no meaningful
speedup (~283–291 MC/s); increasing threadgroup size from 32 to 256 likewise
made no meaningful difference after reordering. Neither change is retained.
The gain supports improving access order, not merely removing integer divides.
No hardware DRAM-bandwidth claim is made. Original CPU arrays plus reordered
copies cost ~887 MB extra for this board; further packing/fusion may improve
memory use and speed, but must preserve region ordering and boundary semantics.

## Historical scalar-layout performance spot check

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
