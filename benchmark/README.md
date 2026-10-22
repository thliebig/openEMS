# openEMS speed benchmark

A small, fixed set of cases that each isolate one cost of the solver, timed
with the `openEMS` binary. It answers two questions about a build: how fast
is it on this machine, and — against another build — is it faster, slower,
and does it still compute the same thing.

Needs the openEMS and CSXCAD Python modules (to write the cases and compare
the probes) and an `openEMS` binary.

## Run

```bash
python3 benchmark/bench.py run
```

times the `openEMS` in `PATH` on all cases at sizes S and M, prints a table
and writes `bench-<date>-<time>/result.json` with the XML, the logs and the
probes of every run next to it. Other binaries: `--openems path/to/openEMS`.

| Option | Default | |
|---|---|---|
| `--cases a,b` | all | `ideal`, `pml`, `inhomogeneous`, `dispersive`, `cylindrical` |
| `--size S,M` | `S,M` | `XS` (~50 k cells, smoke test), `S` (~1 M), `M` (~8 M), `L` (~30 M, 3–4 GB RAM) |
| `--engine a,b` | `multithreaded`, plus `gpu` with a device | also `basic`, `sse`, `sse-compressed` (single-threaded) |
| `--threads 1,4,8` | the physical cores | performance cores only on a hybrid CPU |
| `--gpu-device N` | the best device | |
| `--duration T` | 3 | seconds of time-stepping per run |
| `--repeat N` | 1, with two binaries 3 | timed runs per point |
| `--workdir DIR`, `-o FILE` | `./bench-<date>-<time>/`, `<workdir>/result.json` | |

The thread count is always passed explicitly: with the default
`--numThreads=0` openEMS ramps the threads up during the run and the reported
speed is an average over that ramp. Logical CPUs (SMT) are not used by
default; a memory-bound FDTD rarely gains from them.

## Compare two builds

Build the branch and master into two prefixes, then run them **interleaved**
on the same machine — each run of one binary followed by one of the other,
so that turbo, temperature and background load hit both alike:

```bash
python3 benchmark/bench.py run --openems ~/opt/master/bin/openEMS --label master \
                               --openems ~/opt/branch/bin/openEMS --label branch
```

The comparison table is printed at the end; it can be printed again with
`bench.py compare bench-*/result.json`. Two separately recorded results
compare with

```bash
python3 benchmark/bench.py compare base.json other.json   # or file.json:LABEL
```

— a stored baseline, another machine, a GPU against a CPU. Rows are matched by
case, size and engine; each side keeps its own thread count. When the machine
differs, the header says so: the ratio then is the hardware, not the code.

```
case           size engine         thr A/B    A MC/s    B MC/s    B/A  spread   change setup B/A  probes
ideal          M    multithreaded      8/8     442.1     501.7  1.135    1.2%     +13%      1.00  identical
pml            M    multithreaded      8/8     184.6     185.3  1.004    0.9%    noise      0.99  identical
```

- **MC/s** — median speed of the time-stepping in million cells per second,
  the `Speed:` openEMS prints at the end of a run. It includes the field
  processing of the run: probes and the energy estimate of the end criteria.
- **spread** — the larger of both sides' (max − min) / median, `-` for a
  single run; **change** is `noise` when the ratio lies within it or within
  3 %. A single run (the default with one binary) is a quick number; a
  comparison needs the spread, so with two binaries each runs 3 times.
- **setup B/A** — process wall time minus the time-stepping: reading the XML
  and building the operator, which is a speed target of its own.
- **probes** — the largest deviation between both sides' probes, relative to
  each probe's peak: `identical`, a value, or `MISMATCH` above 1e-4.
  A mismatch does not stop anything, a change in the last bits (fused
  multiply-add, another summation order) is legitimate; a broken update is
  not, and shows up here. `compare` exits with 1 on a mismatch.

## Cases

| Case | Mesh | Content | Isolates |
|---|---|---|---|
| `ideal` | uniform | vacuum, PEC and PMC walls | the bare update kernel, one coefficient set (both walls are only zeroed coefficients) |
| `pml` | uniform | vacuum, PML_8 on all sides | the UPML extension |
| `inhomogeneous` | graded 0.25–1 mm | microstrip with two lumped ports, dielectric, lossy and metal blocks, PML_8 | many distinct coefficients, the usual extensions |
| `dispersive` | uniform | a Lorentz and a Debye block, Mur | the dispersive material extensions |
| `cylindrical` | cylindrical incl. r=0 | a dielectric ring, PEC and Mur | the cylinder operator and extension |

Changing a case changes the results: bump `cases.CASES_VERSION` with it.

## How a run is timed

Every run is stopped after the same time, `--duration` seconds of
time-stepping (default 3), whatever the engine or the machine: `bench.py`
waits for openEMS to start its time-stepping, then creates the file `ABORT`
in the run's directory, which openEMS checks after every batch of timesteps
and stops gracefully. This works with every openEMS build, also an old
release to compare against. The speed is the one openEMS measured itself,
over the timesteps it did. The setup before comes on top and is reported
separately.

The cases are driven by a sinusoidal source, so the energy never decays and
no run ends on its own; one that does is reported as failed. A sinus also
enables the steady-state detection, which is part of every case. Since the
two sides of a comparison stop at different timesteps, the probes are
compared over the timesteps both did; a run is deterministic, so a timestep
has the same result however long the run went on.

A run of the defaults takes about 9 minutes per binary on an 8-core desktop
with a GPU, 18 minutes for an A/B comparison; the setup of the M cases is
almost half of it.

## In a pull request

A pull request that claims a speed-up includes the `compare` table of master
against the branch, run interleaved on one machine, together with the machine
line from its header. A speed-up that only shows on a case outside this suite
is shown with that case in the pull request, not by adding it here.
