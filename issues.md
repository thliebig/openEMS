# Engine issues found while extending the simulation tests

Found while writing the simulation tests in `python/Tests/` (see their
docstrings). Each test works around the issue it hit rather than asserting
the current behaviour.

AI disclosure: found and written up with Claude Opus 5 (Claude Code).

## 1. Debye material is unstable for large `EpsilonDelta`

**Where:** `FDTD/extensions/operator_ext_lorentzmaterial.cpp`, Debye branch
(`v_int`/`v_ext` for `R_D`, `C_L` without `L_D`, around line 330), and the
update in `engine_ext_lorentzmaterial.cpp`.

**What:** The polarisation current of the series R-C branch is updated as
`i[n+1] = (2 dT/tau - 1) i[n] + 2/R (v - v_C)`, with `v_C` already advanced
by `dT/C i[n]`. Substituting that gives `i[n+1] = -i[n] + 2/R (v - v_C[n])`.
That recursion has an eigenvalue of -1, and together with the explicit
coupling back into the E update it has a root with |lambda| > 1.

**Reproduce:** a 1D plane-wave channel (see `Dispersive_Materials.py`) with
a 10 mm Debye slab, dz = 0.25 mm, Gaussian pulse 1-10 GHz, 20000 timesteps:

| eps_inf | d_eps | f_relax | result |
|---|---|---|---|
| 2 | 3 | 4 GHz | diverges (1e33, then NaN) |
| 2 | 3 | 1 GHz | diverges |
| 1 | 1 | 4 GHz | diverges |
| 2 | 3 | 0.2 GHz | grows slowly at late time |
| 2 | 0.5 | 4 GHz | stable |
| 5 | 1 | 4 GHz | stable |

Typical media such as water (eps_inf ~5, d_eps ~75) are well inside the
unstable range. The stable cases never reach the -50 dB end criterion
either.

**Suggested fix:** use a standard Debye ADE (e.g. a semi-implicit
update of the polarisation current driven by `(E[n+1] + E[n])/2`) instead of
the trapezoidal update of an algebraic relation.

## 2. B-field dumps use the primary time and mesh — fixed

**Where:** `openems.cpp`, `SetupProcessing` (around line 636).

**What:** Only H-field dumps (types 1 and 11) were given
`SetDualTime(true)` and `SetDualMesh(true)`. B-field dumps (types 5 and 15)
read the same dual-grid values (`GetBField` -> `GetRawInterpolatedDualField`)
but did not get these settings, so:

- the TD time stamps of B dumps were half a timestep early (the E time
  instead of the H time), and FD B dumps had the matching phase error;
- with `DumpMode` 0 (no interpolation), B values were labelled with
  primary-mesh coordinates. A B dump of the plane z=0 then held the dual
  plane z=+dz/2 while the H dump held z=-dz/2, and x/y were shifted by half
  a cell, so B and mue0*H did not match anywhere.

With node interpolation (the default) only the time stamp was affected.

**Fix (9638fae):** treat types 5 and 15 like 1 and 11. `FieldProbes.py` now also
compares the non-interpolated B dump against the H dump.

## 3. NF2FF `Prad` and `Dmax` count the mirrored half-space

**Where:** `nf2ff/nf2ff_calc.cpp`, `AddSinglePlane` (line 414) called for
the mirrored planes from `AddPlane`/`AddMirrorPlane`.

**What:** `Prad` is the Poynting flux through the recording surfaces. With
a PEC/PMC mirror, the mirrored (image) surfaces are added too, so `Prad` is
the power of antenna plus image: twice the physical value for one mirror.
`Dmax = 4 pi r^2 P_max / Prad` is halved accordingly.

**Reproduce:** `NF2FF_Antennas.py`, quarter-wave monopole on a PEC ground
(z-min boundary PEC, `CreateNF2FFBox` sets the mirror): the nf2ff result
gives `Prad / P_acc = 1.92` and `Dmax = 1.65`. Integrating the far field
over the upper half-space gives `0.96` and `3.30` (analytic 3.28).

**Suggested fix:** exclude the mirrored planes from `m_radPower` (they
still have to contribute to the far field), or document that `Prad`/`Dmax`
refer to the equivalent free-space problem.

## 4. End criteria are only checked every ~4 s of wall-clock time — fixed

**Where:** `openems.cpp`, `RunFDTD` (loop at line 1425).

**What:** The loop condition tests `change > endCrit`, but `change` is only
updated inside the progress-report branch (`if (t_diff>4)`). This applies to
both the energy end criterion and the steady-state detection. A run
therefore continues for up to ~4 s after the criterion is met, and the last
timestep depends on the machine speed and load. Small simulations run many
times longer than needed. Output from the same input and build is also not
reproducible, which works against the bit-exact comparison of builds (#205).

**Reproduce:** `SteadyState.py`: the steady state is reached after a few
dozen periods, but the run ends after ~990 periods (the first report after
4 s), with a period-to-period change of -70 dB against a -60 dB criterion.

**Fix (e757648):** the energy is evaluated once per Nyquist period (at most
every 100 timesteps) and the steady-state result after every iteration
block, independent of the report. Across the 18 simulation tests the runs
got up to 99 % shorter (`SteadyState.py`: ~990 -> 5 periods).

## 5. Steady-state probes collapse to line 0 — fixed

**Where:** `openems.cpp`, steady-state detection setup (around line 1221).

**What:** `pos[n] *= 1/4;` and `pos[n] *= 3/4;` are integer divisions,
so both set `pos[n]` to 0. Of the nine intended probe positions, six end up
on line 0 (usually a boundary, where the field is zero or not updated), and
only the centre position is useful. The intent of 38ff1ce was the quarter
and three-quarter positions.

**Fix (334741d):** `pos[n] = N/4` and `pos[n] = 3*N/4` of the number of lines.

## 6. Progress report leaks `cout` formatting into later output

**Where:** `openems.cpp`, `RunFDTD` progress report (lines 1449-1462).

**What:** The report sets `std::fixed` / `std::scientific` and
`setprecision(1..3)` on `cout` and never restores them. Everything printed
afterwards in the same process uses that format. A second simulation run
from Python (e.g. a reference run followed by the device run) prints
`Create a steady state detection using a period of 0.00 s` and
`FDTD timestep is: 0.00 s`.

**Suggested fix:** save and restore the stream state around the report
(`std::ios_base::fmtflags` + precision), or format the report into an
`ostringstream`.

## 7. Extension hooks write fields that are documented as read-only

**Where:** `FDTD/extensions/engine_extension.h` (hook documentation) vs.
`engine_ext_upml.cpp`, `engine_ext_tfsf.cpp`, `engine_ext_cylinder.cpp`.

**What:** The documentation says `DoPre*Updates` and `DoPost*Updates` may
not change the engine voltages/currents; only `Apply2*` may. In practice:

- UPML writes the fields in all four pre/post hooks (it swaps its flux into
  the field before the main update and back after it);
- TF/SF and the cylinder extension write in `DoPostVoltageUpdates` and
  `DoPostCurrentUpdates`.

All other extensions only read in these hooks. The CPU engines work because
the hooks run in a fixed order, but the documented contract cannot be relied
on, e.g. by an engine that keeps the fields elsewhere (the GPU engine's host
fallback uploads both fields before every main update for this reason).

**Suggested fix:** document the actual contract (hooks may modify the field
in the region the extension owns), or move the writes to `Apply2*`.
