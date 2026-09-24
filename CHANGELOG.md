# openEMS Changelog

Notable, user-visible changes in openEMS. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
`0.MINOR.PATCH`.

The rendered version of this file, together with the CSXCAD changelog, is part
of the [openEMS documentation](https://docs.openems.de/).

**When you change something user-visible, add a bullet under `Unreleased` in the
same commit**, creating that section if it is not there — it exists only while
there are unreleased changes. It becomes the next release entry.

## Unreleased

The version scheme changed with this release: the previous release was v0.0.36,
this one is v0.37.0. The leading `0.0.` was never going anywhere, so the minor
number moved up and patch releases now have somewhere to go.

### Added

- **SAR calculation reworked.** Averaging is done once for all frequencies
  instead of per frequency, and the calculation is multi-threaded, together
  giving a large speedup. Averaging follows IEEE/IEC 62704-1. The `--autorange`
  option restricts the calculation to the cells whose local SAR is within a
  given range of the peak, plus a padding of about one averaging cube; it is a
  speedup and not a guarantee to find the global peak, and a warning is printed
  when the averaged peak falls below the threshold. There is simple progress
  feedback. Available from Python (`sar_calculation`) and from Octave, where
  `CalcSAR.m` exposes `autoRange`, `numThreads` and `progress`. A tutorial
  demonstrating SAR averaging was added.
- **Waveguide mode excitation and probe from an HDF5 mode file**, complementing
  the analytic mode functions. `WaveguidePort`/`RectWGPort` accept a
  `mode_file` argument, the `Add*WaveGuidePort` methods accept a `local_origin`,
  and `matlab/h5writemode.m` (with an Octave oct-file) writes the mode files.
  Requires the matching CSXCAD support.
- **New Python port classes**: `CircWGPort`, `CoaxialPort`, `StripLinePort`,
  `CPWPort` and `CurvePort`, each with an integration test, plus a
  `check_mode_purity()` utility. Convenience methods `AddCircWaveGuidePort`,
  `AddCoaxialPort`, `AddStripLinePort`, `AddCPWPort` and `AddCurvePort` on the
  `openEMS` class mirror the existing `AddLumpedPort`/`AddMSLPort` pattern.
  Ports also store their time base, which makes plotting easier. On a
  cylindrical mesh `CircWGPort` uses the mode profile in its native
  (rho, a, z) form.
- **`CreateNF2FFBox()` (Python) accepts `directions` and `mirror`**, which
  overrule the settings derived from the boundary conditions, e.g. to leave
  the face an antenna feed passes through out of the Huygens surface.
- **`DelayFidelity()` (Python)**, the time delay from the source port to the
  antenna phase centre and the fidelity of the radiated pulse, for any
  polarisation via complex weights on E_theta and E_phi.
- **Localized Mur and SA-Mur absorbers** as stand-alone engine extensions, so an
  absorbing boundary can be placed inside the simulation domain rather than only
  at its edge.
- **Lumped RLC elements**: series and parallel R/L/C, implemented with an
  auxiliary differential equation.
- **`HDF5Dump` (Python)** for reading field dumps. A dump can be inspected
  before any field data is read (TD/FD, dump type, grid size, frequencies,
  timesteps); `SetPlane`/`SetLine`/`SetRange`/`SetSampling` push the selection
  down into HDF5 so only the requested part is read; and
  `GetFieldAtFrequency()` runs an on-the-fly DFT so time-domain and
  frequency-domain dumps are post-processed the same way. Stored values are
  plain attributes (`file`, `shape`, `frequencies`, `dump_type`) and everything
  that computes is a call (`GetNumTimesteps()`, `GetNumFrequencies()`,
  `GetTimes()`, `GetDumpTypeName()`, `IsTD()`, `IsFD()`, `IsVector()`).
- **`SetLibraryArguments()`**, which accepts almost every option of the
  `openEMS` executable as a string. The Python binding uses it, so engine
  selection and the debugging options are now reachable from Python.
- **Graceful abort on Ctrl-C.** SIGINT is handled explicitly, so a run can be
  stopped and its data written out instead of being killed. This also restores
  Ctrl-C when openEMS runs as a Python module, where it previously did nothing.
- FDTD setup and the CSX structure can be read from and written to XML.
- The FDTD object can be reset from Python, and an excitation can be enabled and
  disabled.
- Bundled head and body phantoms in `resources/phantoms/`, reachable from every
  language interface, installed to `share/openEMS/` and shipped inside the
  python wheel. `openEMS_resource_path()` (Octave) and
  `openEMS.utilities.get_resource_path()` (Python) return the path of such a
  file, so a tutorial keeps working from wherever it was copied to. The MRI
  tutorials fall back to the phantoms when the IT'IS Virtual Family dataset is
  not installed.
- New tutorials and examples: Python `Horn_Antenna` (coaxial pin feed with
  backshort), `StripLine2MSL`, `Dipole_SAR`, `MRI_Loop_Coil`, a rectangular
  resonant cavity example, and a SAR averaging tutorial. Ported from Octave:
  `Parallel_Plate_Waveguide`, `Circ_Waveguide`, `Conical_Horn_Antenna`,
  `CylindricalWave_CC` and `RadarUWBTutorial`.
- Python unit and integration tests, run in CI after each smoke test.
- **Optional oversampling for frequency-domain dumps and probes.** The running
  DFT accumulation for FD/SAR dumps and probes is sampled at exactly the
  Nyquist rate, which aliases the spectrum just above the highest excited
  frequency onto the upper band edge. A probe/dump can now oversample its own
  FD accumulation (e.g. by the same factor `OverSampling` gives the
  time-domain recording, default 4) via CSXCAD's new `OverSampling` property
  on that box — see the CSXCAD changelog. It defaults to the plain Nyquist
  rate, matching prior behavior, so existing simulations and their
  performance are unaffected unless a box opts in.
- **`--verbose`/`-vv` now reports the actual TD/FD sampling interval** for
  each probe and dump box during `SetupProcessing`.
- **`--exact-endcriteria`.** The energy end-criteria is normally re-evaluated
  every few seconds of wall-clock time, to keep its cost (a full-domain
  energy estimate) off the hot path; this makes the exact stopping timestep
  depend on machine speed/load. This option instead evaluates it every
  Nyquist period, for a stopping point that is reproducible across
  machines/builds, at the cost of performance — mainly useful for engine or
  code verification. The steady-state detection extension is unaffected: its
  diff estimate is cheap and is now always kept current every timestep
  rather than only at the wall-clock report interval.

### Changed

- **The MPI engine was removed.** It had not compiled for years, as it used
  the C++ MPI bindings that MPI-3 dropped, it had no tests, and several
  extensions never supported it (#260). The multithreaded engine is
  unaffected. `WITH_MPI`, `--engine=MPI`, `openEMS_MPI.sh`, `RunOpenEMS_MPI`
  and `SetupMPI` are gone. See *Upgrade notes*.
- **nf2ff result format.** The far field is written as one compound complex
  dataset per frequency, `/nf2ff/E_theta/FD/f{n}`, stored in (theta, phi)
  order — the format every other frequency-domain dump has used since HDF5
  version 0.3 — instead of a split `f{n}_real`/`f{n}_imag` pair in
  (phi, theta) order. `h5py` reads it as a native complex array, so no axis
  has to be swapped after reading. The Octave/Matlab interface keeps the old
  layout, which `CalcNF2FF` requests through the new `LegacyHDF5` attribute of
  the nf2ff XML file: Octave reads a compound complex dataset as zeros without
  any error. See *Upgrade notes*.
- **Simulation directory cleanup no longer deletes the directory.**
  `CleanupSimPath()` (Octave/Matlab) and `cleanup=True` in Python's `FDTD.Run()`
  now remove only known openEMS output files. Pointing `Sim_Path` at `$HOME`, or
  any other directory that matters, no longer destroys its contents. Generic
  `*.h5` files are only removed when they carry the `openEMS_HDF5_version` root
  attribute, so user-supplied HDF5 files such as mode files survive.
- Command-line argument parsing was rewritten. Options are declared by the
  module that uses them rather than in `openems.cpp`, which is what makes
  `SetLibraryArguments()` possible.
- The Gaussian excitation now ends at exactly zero, and the excitation types
  were renamed more descriptively.
- Frequency-domain dump files get more readable names, and HDF5 dumps carry more
  metadata attributes.
- `boost/program_options.hpp` was removed from the public `openems.h`.
- The `INVALID` lumped-element type was removed, following the same change in
  CSXCAD.
- Octave/Matlab docstrings were reformatted as Markdown so that the online
  function reference can be generated from them.
- **Python packaging modernised**: `pyproject.toml`, installable with `pip`,
  dynamic versioning via `setuptools_scm`, and more robust detection of an
  installed CSXCAD.
- Internally, the multi-dimensional field arrays were replaced by a new
  `ArrayLib`, and the engine, operator and their extensions converted to it.
- Default thread counts (multithreaded engine, nf2ff, SAR) now respect CPU
  affinity and cgroup CPU quotas (Linux only) instead of always using every
  CPU of the host, so a container or systemd unit with a CPU limit no longer
  oversubscribes it.

### Fixed

- Python: a `ContinuousStructure` handed to `SetCSX()` is no longer destroyed
  twice. `SetCSX()` takes ownership, which the binding now states, and
  `GetCSX()` no longer leaks a fresh structure per call. Requires a CSXCAD
  providing `CSObject`. See *Upgrade notes*.
- Python: `GetCSX()` returned an empty list after the structure had been read
  with `ReadFromXML()`, and a wrapper could dangle after `Reset()`.
- UPML: copy-paste errors in the update coefficients (#221).
- Tutorial `CRLH_LeakyWaveAnt.m`: the ground plane was added at the same
  priority as the substrate that spans z=0 as well, so it lost the tie and
  never made it into the operator.
- Modes higher than 0 in the parallel-plate direction were not excited
  correctly.
- The mode-match probe coordinates now match the excitation coordinates.
- nf2ff: `m_maxDir` was wrong for a radius other than 1.
- A steady-state engine extension could be freed twice on shutdown, and an
  operator extension that was never initialised could be freed invalidly. An
  extension that finds nothing to do — a conducting sheet without a primitive,
  for example — is now dropped instead of being kept and run empty.
- `Dmax` is a linear power quantity and was added to a dB value without
  conversion in the patch-antenna tutorials.
- Lumped RLC: the auxiliary-differential-equation update had several bugs, and
  `C = 0` was not handled.
- Python: `Run()` failed on a relative `sim_path`, and on one containing
  symlinks, with an assertion; both work now.
- Octave: an oct-file left over from an older Octave version is rebuilt
  instead of failing the run. It was still found by `exist()`, so `setup` was
  never re-run and the call died with "failed to load" or, on Windows, "the
  specified module could not be found" (#318).
- Octave/Matlab: paths containing spaces are quoted for the binary and the log
  file, and HDF5 detection in `setup.m` was improved and is now tested in CI.
- Python: `SetCustomExcite` encoding, and argument parsing with several
  `openEMS()` instances in one process.
- The excitation amplitude was ignored by `AddCoaxialPort` (`'ExciteAmp'`) in
  Octave/Matlab and by the waveguide ports (`excite`) in Python: any non-zero
  value excited with amplitude 1.
- Python: a port with the number of an existing port now raises `ValueError`.
  Both ports wrote their probes to the same files and corrupted them.
  Octave/Matlab already rejected this.
- Octave: `plotRefl` died with "vertical dimensions mismatch" instead of
  drawing the Smith chart, because it added the trace after the legend and
  Octave's legend autoupdate could not append it (#172).
- B-field dumps (`DumpType` 5/15) were not placed on the dual time/mesh like
  H-field dumps, despite reading the same dual-grid values: values were
  labelled half a cell and half a timestep off.
- Six of the nine steady-state detection probes sat on the first mesh line of
  their direction (usually a field-free boundary) instead of a quarter/three
  quarters across, due to an integer-division bug.

### Build

- C++11 is now required, and CMake 3.1 or newer.
- VTK 9 and newer are supported without deprecated names.
- Windows: builds via vcpkg manifest with MSVC and clang-cl under Visual
  Studio 2022; the `openEMS` and `nf2ff` import libraries are installed to
  `lib/`.
- Builds on ppc64le.
- Two new knobs for comparing the output of two builds bit by bit, both off the
  default path: debug builds compile with `-ffp-contract=off` on GCC and Clang,
  so multiply-add pairs are no longer contracted into FMA instructions, and the
  new `ENABLE_FLUSH_TO_ZERO` CMake option can be set to `OFF` to keep denormal
  values in the engines instead of flushing them to zero.
- The `WITH_MPI` CMake option and the `--with-MPI` option of
  `update_openEMS.sh` were removed.
- CI covers Linux, macOS, FreeBSD and Windows, and compiles with warnings
  enabled.

### Upgrade notes

- **Rebuild all components together.** The CSXCAD `CSObject` change alters the
  layout of every class deriving from it. The soname is unchanged and will not
  catch a partial rebuild. This release needs a CSXCAD that provides `CSObject`.
- A Python script that used a CSXCAD wrapper after its C++ object had been
  destroyed now stops with a `RuntimeError` instead of reading freed memory.
- Scripts that relied on `cleanup` wiping the whole simulation directory now
  keep any file that is not recognised openEMS output. This is deliberate.
- The bundled phantoms moved from `matlab/Tutorials/phantoms/` to
  `resources/phantoms/`, installed under `share/openEMS/`.
- A tool that reads nf2ff result files directly has to handle the compound
  complex datasets described above; `nf2ff_results` (Python) and `ReadNF2FF`
  (Matlab) read both formats, `ReadNF2FF` under Octave only the legacy one.
  Files written by the Octave/Matlab interface are unchanged.
- Octave/Matlab scripts calling `SetupMPI` fail, as the function is gone; drop
  the call. `RunOpenEMS` warns about a `Settings.MPI` field and runs the
  multithreaded engine instead.

## Older releases

Releases v0.0.32 (2013-11-27) through v0.0.36 (2023-10-22) were not recorded
here. See the [commit history](https://github.com/thliebig/openEMS/commits/master)
or compare two tags, for example
[`v0.0.35...v0.0.36`](https://github.com/thliebig/openEMS/compare/v0.0.35...v0.0.36).

The entries below are the original `NEWS` file, kept verbatim.

### v0.0.31

- nf2ff: calculate circular polarization
- improvements to calcPort
- allow 1D and 2D lumped ports
- improvements to SAR calculations
- curve primitives and port fixes & improvements
- FDTD operator now supports different material averaging methods
- support for full multi-polar Lorentz/Drude/Debye dispersive material types
- new ports for waveguides (rectangular and circular waveguides)
- improved PEC debugging
- improvements/simplifications in plotting far-fields (thanks to Stefan)
- new tutorials and examples
- many fixes and updates

### v0.0.30

- meshing improved with new detect edges and new smoothing capabilities
- new calcPort function for simplified port analysis
- cylindrical mesh improvement by considering 360° rotation symmetry
- support for harminv on all platforms
- update to auto-regressive model for voltage/current probes
- new SAR calculation options, incl. 1g/10g averaging
- support for a new primitive: polyhedron
- CAD import: STL/PLY surface solids supported (matlab: ImportSTL / ImportPLY)
- CAD export: STL/PLY export (using AppCSXCAD)
- lot of minor fixes and updates

### v0.0.29

- Cylindrical sub-grids now fully support alpha-graded meshes
- Property Electrode has been renamed to Excitation
    This doesn't have any effect on the Matlab/Octave interface,
    but old *.xml files cannot be run with a current openEMS/CSXCAD version.
- Overall memory usage reduced during pre-processing
- New excitation: Total-field/scattered field (TFSF)
- New tutorial on radar cross section on a metallic sphere using the TFSF excitation
- official support for 64-bit windows version
- check for engine extensions MPI compatibility
- CSXCAD: support for new CSXGeomPlot export options (see help for more infos)
