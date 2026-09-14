# Metal UPML boundary conditioning

`--engine=metal` runs the voltage and current UPML pre/post conditioning on the
GPU. Operator construction, grading and the EC/UPML coefficient build stay on the
CPU, as do all non-UPML extension hooks.

## Design

- **Indexed layout (default).** At engine creation the UPML coefficients and
  fluxes are permuted once into increasing SSE packed-field addresses and paired
  with a per-component uint32 index. Stepping then performs one indexed field
  access per lane and contiguous auxiliary I/O, with no per-step coordinate math
  and no per-step copies. Reordering the access order is what produces the
  speedup; specializing integer divisors alone did not.
- **Scalar layout (`OPENEMS_METAL_PML_LAYOUT=scalar`).** Keeps the original
  no-copy scalar kernel as a low-memory fallback and for A/B comparison. The two
  layouts must produce bit-identical fields.
- **In-place reuse (`OPENEMS_METAL_PML_REUSE`, default on).** The operator's
  coefficient arrays are permuted in place and restored on engine teardown, so a
  later CPU or scalar engine sees valid data. Disabling reuse stores separate
  packed copies instead (more memory, same results).
- **Lossless coefficient dictionaries (`OPENEMS_METAL_PML_COMPRESS`, default
  on).** UPML coefficient triples are deduplicated by exact 32-bit pattern. The
  dense CPU arrays are rebuilt on teardown so the FP64 diagnostic and a later CPU
  engine still work.
- **Physical cells only.** The mapping excludes SIMD padding lanes and skips
  zero-extent regions left by opposing slabs.
- **Hook ordering and synchronization.** Extension hooks keep the CPU order (pre
  in reverse priority, post forward). Pre-conditioning, the field update and
  post-conditioning normally share one command submission per half-step, and
  separate compute encoders preserve dependencies between overlapping slabs.
  Pending GPU work is completed before any CPU hook, source, probe or dump can
  read the fields.
- CPU UPML stays available for comparison; this adds no cylindrical or MPI
  support to the Metal engine.

## Controls

| Variable | Default | Effect |
|---|---|---|
| `OPENEMS_METAL_PML` | on | `0` runs CPU UPML with Metal field updates |
| `OPENEMS_METAL_PML_LAYOUT` | `indexed` | `scalar` selects the no-copy scalar kernel |
| `OPENEMS_METAL_PML_REUSE` | on | `0` stores separate packed copies instead of reordering in place |
| `OPENEMS_METAL_PML_COMPRESS` | on | `0` keeps dense UPML coefficients |
| `OPENEMS_METAL_FP64_REFERENCE` | off | `1` enables the FP64 diagnostic below |

`OPENEMS_METAL_FP64_REFERENCE=1` checks each update stencil in FP64 from the
actual conditioned input. It deliberately does not evolve a FP64 flux reference:
flux swaps are not additive field corrections, so an evolved reference has no
matching flux state and diverges. Without UPML the existing evolving reference is
unchanged.

## Validation

```sh
python TESTSUITE/enginetests/metal_pml.py --openems /absolute/path/to/openEMS
python TESTSUITE/enginetests/metal_pml.py --openems /absolute/path/to/openEMS \
  --fp64-reference
MTL_DEBUG_LAYER=1 MTL_SHADER_VALIDATION=1 \
  python TESTSUITE/enginetests/metal_pml.py --openems /absolute/path/to/openEMS
```

The suite covers each face, edges and corners, overlapping opposing slabs,
unequal thicknesses, nonuniform spacing, mixed Mur/PEC/PML, no PML, every
physical z-line-count remainder, a lossy dielectric and a post-pulse run. It
requires complete dense/compressed and scalar/indexed E/H dumps to be
bit-identical, verifies CPU UPML on Metal pointwise, and compares against SSE
with a relative-L2 limit because Metal/SSE cancellation differences also occur
without GPU UPML. The known SSE long-run discrepancy (3895 E values above the
default tolerance) is unrelated to UPML and present without it.

## Performance

Representative M4 Pro, Release, fast math OFF; 658 x 664 x 33 board, six
`PML_8` faces, 788 steps, medians of alternating runs:

| Layout | Stepping | Whole process |
|---|---:|---:|
| scalar | ~39.4 s | ~139 s |
| indexed | ~22.0 s | ~118 s |

Roughly **1.79x stepping** and **1.18x whole-process** speedup. These are finite
runs, not convergence or SI validation; small grids can be submission-bound.
