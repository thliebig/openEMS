# Metal engine

`openEMS model.xml --engine=metal` runs the FDTD field updates on Apple GPUs via
Metal. Build with `-DWITH_METAL=ON`; the option is off by default. It enables
every Metal feature; there are no per-feature switches. The same binary keeps
the SSE and multithreaded engines for comparison.

Operator construction, mesh grading, material EC sampling, UPML grading, the
coefficient build and the coefficient dictionaries stay on the CPU. The GPU
executes:

- the fused voltage/current field update,
- UPML pre/post conditioning,
- the PEC / `MATERIAL|METAL` geometry pass, whose winners are reused by the
  EC-consuming extensions (conducting sheets and dispersive materials),
- the conducting-sheet volt-ADE recurrence.

Cylindrical and MPI operators are not supported by the Metal engine.

## Field updates

The fused E/H pipeline submits one command buffer per timestep (extension hooks,
voltage update, excitation, current update), so the GPU is not drained between
half steps. `OPENEMS_METAL_FP64_REFERENCE` disables fusion and evolves a FP64
reference for diagnostics. Fast math is always off.

## PEC and geometry mapping

The Metal pass resolves the winning `MATERIAL|METAL` primitive per Yee component
and records it via `Operator::GetGeometryWinners()`; `Operator_Ext_ConductingSheet`
and `Operator_Ext_LorentzMaterial` consume those winners instead of re-querying
CSXCAD per cell. Unsupported primitives and genuine near-boundary ties fall back
to the original CPU query, so no approximate decision is used at an edge.

- Supports Cartesian boxes, polygons, linearly extruded polygons, and
  cylinders / cylindrical shells (via-like annuli).
- Yee coordinates and inclusive index ranges are computed in CPU FP64; polygon
  interiors use the CSXCAD winding rule with exact `orient2d` sign (double-float
  coordinates), with a conservative fallback band near the decision boundary.
- Transforms, curves, spheres and other unsupported primitives make the affected
  queries use the CSXCAD path; they are never silently dropped, and can remove
  most of the speedup for such geometries.
- The conducting-sheet extension keeps per-cell sheet state (`sigma`, thickness,
  tangent direction) only for resolved sheet cells instead of three full-grid
  lookup tables (~19 GB on a 696M-cell model).

## UPML

- **Indexed layout (default).** UPML coefficients and fluxes are permuted once
  into increasing packed-field addresses with a per-component `uint32` index.
  Stepping then does one indexed field access per lane and contiguous auxiliary
  I/O, with no per-step coordinate math. Reordering the access order is what
  produces the speedup; specializing integer divisors alone did not.
- **Scalar layout.** The original no-copy scalar kernel remains as a low-memory
  fallback; the two layouts must be bit-identical.
- **In-place reuse.** The operator's coefficient arrays are permuted in place and
  restored on teardown, so a later CPU or scalar engine sees valid data.
- **Lossless coefficient dictionaries.** UPML coefficient triples are
  deduplicated by exact 32-bit pattern; the dense CPU arrays are rebuilt on
  teardown for the FP64 diagnostic and later CPU engines.
- The mapping excludes SIMD padding lanes and skips zero-extent regions left by
  opposing slabs. Hooks keep CPU order (pre reverse priority, post forward) and
  pending GPU work is completed before any CPU hook, source, probe or dump.

## Coefficient dictionaries

The engine deduplicates each packed position's 48 FP32 values (VV, VI, II, IV ×
three components × four lanes) by exact 32-bit pattern, with a shared `uint16`
index buffer. This reduces the GPU coefficient working set, **not** total process
RAM; index and dictionary copies happen only at initialization. Construction
falls back to dense reads when unique records exceed `min(65536, positions/4)`
(the packed `uint16` index limit). Very nonuniform meshes can hit the limit and retain
dense reads automatically; initialization reports which path was taken.
Coefficients must remain immutable during stepping.

## Conducting-sheet ADE

The conducting-sheet model advances two ADE poles per active edge every step.
The Metal engine runs that recurrence in two kernels: `ade_advance` before the
voltage update and `ade_apply` from `Apply2Voltages` after it. One thread owns all
poles of one packed field edge, so the apply is race-free and matches the CPU
subtraction order. Previously the recurrence ran on the CPU and the engine
drained the GPU before each hook, serializing CPU and GPU.

Only the plain volt-ADE scheme is offloaded. Models that need Lorentz flux states
or ADE currents (Lorentz, Drude, Debye) and the FP64 reference mode keep the CPU
path.

## Diagnostic overrides

These environment variables exist for A/B testing and debugging. They all default
to the feature enabled and are not required to use the engine.

| Variable | Default | Effect |
|---|---|---|
| `OPENEMS_METAL_PEC` | GPU | `0` = CPU PEC mapping, `verify` = GPU + CPU compare |
| `OPENEMS_METAL_PML` | on | `0` = CPU UPML conditioning |
| `OPENEMS_METAL_PML_LAYOUT` | `indexed` | `scalar` selects the no-copy scalar kernel |
| `OPENEMS_METAL_PML_REUSE` | on | `0` stores separate packed copies |
| `OPENEMS_METAL_PML_COMPRESS` | on | `0` keeps dense UPML coefficients |
| `OPENEMS_METAL_COMPRESS` | on | `0` keeps dense operator coefficients |
| `OPENEMS_METAL_COEFF_RECORDS` | full range | lowers the dictionary limit |
| `OPENEMS_METAL_FUSED_PIPELINE` | on | `0` runs the unfused pipeline |
| `OPENEMS_METAL_SERIAL_COEFFICIENTS` | off | `1` builds coefficients single-threaded |
| `OPENEMS_METAL_EARLY_EC_FREE` | on | `0` keeps EC arrays before extensions |
| `OPENEMS_METAL_FP64_REFERENCE` | off | `1` enables the FP64 diagnostic |

## Validation

```sh
python TESTSUITE/enginetests/metal_fields.py --openems /absolute/path/to/openEMS --suite
python TESTSUITE/enginetests/metal_pec.py --openems /absolute/path/to/openEMS
python TESTSUITE/enginetests/metal_conductingsheet.py --openems /absolute/path/to/openEMS
python TESTSUITE/enginetests/metal_dispersive.py --openems /absolute/path/to/openEMS
python TESTSUITE/enginetests/metal_ade.py --openems /absolute/path/to/openEMS
python TESTSUITE/enginetests/metal_pml.py --openems /absolute/path/to/openEMS
```

`metal_fields.py` compares SSE against Metal with relative-L2 limits and requires
dense/compressed, fused/unfused and scalar/indexed variants to be bit-identical.
The PEC, conducting-sheet and dispersive suites compare CPU vs GPU winner
resolution and require bit-identical dumps. `metal_ade.py` compares the GPU
conducting-sheet ADE against the SSE CPU recurrence.

## Known limitations

- SSE and Metal long-run comparison: 3895 E values exceed the default pointwise
  tolerance on the 1000-step cavity case (max abs 7.06e-5). The dense Metal path
  shows the same, so it is not a compression or UPML artifact.
- Small or simple geometries are submission-bound and can be slower on GPU than
  on SSE.
- Performance figures are finite-run measurements, not convergence or SI
  validation. Regular grids compress exceptionally well; no speedup is claimed
  for real PCB models that hit the dense-coefficient fallback.

Representative M4 Pro measurements (Release, fast math off):

| Feature | Workload | Result |
|---|---|---|
| Coefficient compression | 16.8M cells, 1000 steps | ~1.54–1.58x stepping, ~1.16x process |
| Indexed UPML | 658×664×33 board, 788 steps | ~1.79x stepping, ~1.18x process |
| GPU PEC mapping | 60-pair CoSwitch pilot | ~8.7 s → 0.09 s PEC pass, ~10.9 s → 2.3 s setup |
| GPU ADE | 240×240×10, 100800 sheet edges | ~1.26x stepping |
