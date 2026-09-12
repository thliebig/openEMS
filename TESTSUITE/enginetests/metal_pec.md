# Experimental Metal PEC mapping

`--engine=metal` enables GPU-assisted PEC mapping as well as Metal field updates.
This accelerates only PEC mapping, not the earlier EC coefficient construction.
The field-update and coefficient-compression paths are unchanged.

```sh
# Enable Metal field updates and GPU-assisted PEC mapping together
openEMS model.xml --engine=metal
# Compare every winning primitive against the original CPU query before applying
OPENEMS_METAL_PEC=verify openEMS model.xml --engine=metal
# Original CPU PEC mapping, with PEC timing printed
OPENEMS_METAL_PEC=0 openEMS model.xml --engine=metal
```

No separate enable flag is needed. `OPENEMS_METAL_PEC` is retained only as a
diagnostic override (`0` for CPU, `verify` for full CPU comparison).
Full 17M-cell production-model validation remains pending. Small/simple
geometries are slower on GPU.

## Scope and semantics

- Supports Cartesian boxes, polygons, linearly extruded polygons, and
  cylinders / cylindrical shells (via-like annuli).
- Computes primal/dual Yee coordinates and inclusive bounding-box index ranges
  in CPU FP64. Nonuniform grids and zero-thickness sheets retain exact normal
  coordinate inclusion; no float tolerance thickens the copper.
- Flattens polygons and builds X-slab / XY-column candidates from exact ranges.
  Preserves CSXCAD's sorted primitive order (including equal-priority ordering)
  and `IsInsideBox` domain culling. A thread handles one electric component.
- Polygon interiors use the CSXCAD winding rule with fast math disabled. FP32
  coordinate/edge comparisons close to equality are sent to CSXCAD FP64, using
  conservative scale-dependent error margins. Large/nonfinite coordinates are
  rejected from the GPU path. No approximate decision is used at a near edge.
- Cylinders and cylindrical shells are tested on the GPU against the axis
  segment and radius (`dist <= radius`; for a shell `|dist - radius| <=
  ShellWidth/2`). Queries within a conservative FP32 margin of a wall, an end
  cap, or a degenerate axis fall back to the CPU, preserving CSXCAD's exact
  double-precision decision.
- Transforms, curves, spheres, and other unsupported primitives cause affected
  queries to use the original CPU candidate list. They are never silently
  ignored. This can remove most of the speedup for some geometries.
- Returns winning primitive IDs (a richer PEC mask), preserving primitive-used
  flags for both material and metal. The CPU sets VV/VI to zero for metal and
  maintains PEC counts. `CalcPEC_Curves()` still runs afterwards.
- One X slab at a time bounds output/candidate memory. Geometry copies occur
  during setup only. New Objective-C++ source uses ARC to release GPU resources.
- GPU runtime errors and verification mismatches stop the run. Missing device /
  kernel support falls back before coefficients are changed.

## Tests

```sh
python TESTSUITE/enginetests/metal_pec.py --openems /absolute/path/to/openEMS
python TESTSUITE/enginetests/metal_pec.py --openems /absolute/path/to/openEMS \
  --model /absolute/path/to/model.xml --keep
```

Each case runs CPU, verified GPU, and default GPU PEC mapping (environment
variable unset), all with the same
Metal stepping engine. Checks include every winning primitive in verify mode,
complete numeric HDF5 datasets, probe numeric output, and unused-primitive
warnings. Fixtures cover nonuniform sheets/boxes, diagonal traces/narrow gaps,
extrusion, exact and nearly coincident edges, equal-priority material overlap,
cylinders and cylindrical shells (exact and off-grid walls), and transforms.
Metal API validation also passed these fixtures.

For the local CoSwitch pilot XML (585 polygons, 33 boxes; 926,970 solver cells,
1000 steps), all 2,780,910 winning-primitive queries matched CPU. Only 671 queries
needed CPU refinement. All port waveforms and unused warnings matched. A fixed
50-ohm S[:,0] column over 1–10 GHz computed from MSL center voltage / averaged
current probes was bit-identical. This is a finite-run regression comparison,
not a converged full multi-excitation production S matrix.

Representative same-machine M4 Pro runs (not medians):

| Mapping | PEC pass | Complete process | Stepping |
|---|---:|---:|---:|
| CPU | 8.723 s | 27.440 s | 16.53 s |
| GPU + all-winner verification | 8.930 s | 27.608 s | 16.55 s |
| GPU | 0.085 s | 18.794 s | 16.52 s |

The speedup includes spatial candidate pruning as well as GPU execution; it
should not be attributed to shaders alone. Setup excluding stepping dropped
from approximately 10.91 s to 2.27 s.

The separate 17.3M-cell XML (30,845 polygons, 469 boxes, 280 extrusions) did not
complete within the bounded test attempts (300 s CPU / 180 s GPU mode), and did
not emit a PEC completion result. No setup speedup, complete mask equivalence,
or S-parameter validation is claimed for that model yet.
