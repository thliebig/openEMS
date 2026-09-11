# Lossless Metal coefficient dictionaries

The Metal engine deduplicates each packed position's 48 FP32 coefficient values
(VV, VI, II, IV, three components, four lanes) by their exact 32-bit patterns.
A shared `uint32` index buffer maps packed positions to unique records. Only
coefficient addresses change in the shaders; field storage, stencil arithmetic,
fast-math settings, CPU extensions and half-step synchronization are unchanged.

The dictionary is built once after operator construction. Original dense CPU
coefficients remain available, including for the FP64 diagnostic. Consequently,
this reduces the GPU coefficient working set, **not total process RAM**. Index
and dictionary copies occur only during initialization, never per timestep.

Compression is enabled by default. `OPENEMS_METAL_COMPRESS=0` selects dense
coefficients for A/B comparisons. A Metal function constant removes the index
lookup entirely from the dense shader specialization.

Construction aborts if there are more than `min(65535, packed_positions / 4)`
unique records; `65535` is the packed `uint16` index limit and can be lowered
with `OPENEMS_METAL_COEFF_RECORDS` for A/B comparison. This bounds the dictionary
to 12 MiB and requires at least fourfold record reuse. The original 4096 cap was
too small for graded UPML coefficients, so any model with PML fell back to dense
reads for the whole operator; the full index range keeps those models
compressed. These are still heuristics, not a cache-size query or a guarantee of
speedup. Very nonuniform meshes can reach the limit and automatically retain
dense reads. Initialization reports compression counts or dictionary-limit
fallback. Coefficients must remain immutable during stepping, as in the current
solver; future time-varying coefficient support would need to refresh the
dictionary or disable it.

## Validation

```sh
python TESTSUITE/enginetests/metal_fields.py \
  --openems /absolute/path/to/openEMS --suite --compare-dense --fp64-reference
python TESTSUITE/enginetests/metal_fields.py \
  --openems /absolute/path/to/openEMS --compare-dense \
  --cells 128 128 122 --timesteps 300
python TESTSUITE/enginetests/metal_fields.py \
  --openems /absolute/path/to/openEMS --compare-dense --nonuniform \
  --cells 48 47 46 --timesteps 300
```

`--compare-dense` requires bit-identical numeric HDF5 datasets, including complete
E/H field dumps. With `--fp64-reference`, the reported diagnostic errors must
also match. The test checks that compression or its explicit fallback actually
ran. `--nonuniform` uses varying mesh spacings to exercise fallback.

On M4 Pro, fast math OFF:

- Tiny boundaries, odd dimensions, dielectric/conductive material, 1000-step
  long run, and the ~2-million-cell case: dense/compressed dumps bit-identical.
- FP64 diagnostic summaries: identical for all four suite cases.
- Nonuniform 48 x 47 x 46 case: bit-identical dumps; whether it compresses or
  keeps dense reads depends on the `OPENEMS_METAL_COEFF_RECORDS` cap.
- PML models: graded UPML coefficients stay compressed with the full index range
  instead of falling back to dense reads for the whole operator.
- Existing SSE versus Metal long-run discrepancies remain: 3895 E values exceed
  default tolerances, maximum absolute error 7.05719e-5. The dense Metal path has
  the same discrepancies; compression does not fix or worsen them. Other tested
  cases have no SSE tolerance failures. The original SSE comparison reports
  tolerance counts rather than failing the process; the bitwise check does fail.

## Performance measurements

M4 Pro, 20 GPU cores, Release, fast math OFF; uniform Cartesian grids with PEC
boundaries and the test script's lossy dielectric box, dumps removed,
EndCriteria=0. Three runs per variant, alternating dense/compressed order;
medians below include the full process, not just shaders. Dimensions are cell
counts supplied to the model generator; solver counts include boundary lines.

| Grid | Steps | Dense wall | Compressed wall | Dense stepping | Compressed stepping |
|---|---:|---:|---:|---:|---:|
| 256 x 256 x 254 (~16.84M solver cells) | 1000 | 23.01 s | 19.82 s | 8.77 s | 5.55 s |
| 1024 x 512 x 30 (~16.30M solver cells) | 1000 | 21.48 s | 18.47 s | 8.64 s | 5.60 s |

That is ~1.16x complete-run speedup and ~1.54–1.58x stepping speedup.
The first grid has 62 unique records; its GPU coefficient/index working set
shrinks from 811.61 MB to 16.92 MB. The thin grid has 42 unique records and
shrinks from 807.67 MB to 16.83 MB. Setup dominates shorter runs, so stepping
speedups must not be presented as whole-simulation speedups.

For CPU context, the first grid at **300** steps measured 33.14 s wall for SSE
and 19.03 s for multithreaded (three-run medians). Representative stepping
throughputs were 267 and 434 MCells/s respectively, versus roughly 1.9 GCells/s
dense Metal and 3.0 GCells/s compressed Metal. These regular grids compress
exceptionally well; no speedup is claimed for real PCB models that hit fallback.
