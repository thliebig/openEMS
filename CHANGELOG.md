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

### Added

- Opt-in `--dft-snapshots` requests export accumulated E/H HDF5 data at a
  timestep boundary without stopping or resetting DFT sums. Snapshot directories
  are published atomically and labelled unconverged with sample timing. Python
  and Matlab request helpers are provided; final field dumps are unchanged.

- **A GPU engine**, `--engine=gpu` (`engine='gpu'` from Python), with a Metal
  backend on macOS and a HIP backend for AMD (ROCm) and NVIDIA (CUDA) GPUs. It
  is selected like the other engines, runs Cartesian and cylindrical meshes
  (including multi-grid) with the usual extensions on the device, and writes
  the same probes and dumps. Without a usable GPU the setup fails with an
  error instead of running something slower. `--engine=gpu-reference` runs it
  on a CPU reference backend, to check the device backends against.
- `python -m openEMS.compare <sim_path_a> <sim_path_b>` compares the probes
  and field dumps of two runs, e.g. of two engines, and
  `openEMS.compare.compare_sim_results()` does the same from Python.

### Fixed

- A run switched flush-to-zero on in the calling process and left it on: in
  Python, numpy computations after `Run()` lost their denormal (subnormal)
  values, and the results of a following run with the basic engine depended
  on it. Every engine now flushes denormals to zero during the run only.
- Python: an option passed as `False` to `Run()`, e.g. `debug_PEC=False`, made
  openEMS ignore all options after it without a message, so that for example
  a following `engine=...` or `numThreads=...` had no effect.
- A node-interpolated B-field dump of a closed cylindrical mesh was zero on
  its first alpha line: the B field was not evaluated across the seam of the
  mesh as the H field is.
- A simulation that was aborted (`SetAbort()`, an `ABORT` file, Ctrl-C) no
  longer ends with the warning that the max. number of timesteps was reached;
  it reports the abort.
- An error during the setup or the run of a simulation, e.g. running out of
  memory, no longer terminates the process: openEMS ends with an error message,
  Python raises a `RuntimeError`, and Ctrl-C works again afterwards.
- A simulation with more than about 1.4 billion mesh nodes no longer
  overflows the size of its field and coefficient arrays, which corrupted
  memory without a message. The operator statistics show the right memory
  size above about 89 million nodes, and the cell count, the speed, frequency
  domain and SAR dumps, the SAR calculation, VTK field dumps and the node
  counts of materials, excitations and lumped elements no longer wrap around
  at 2 to 4 billion nodes.
- A SAR dump kept a second, unused copy of its frequency domain fields:
  24 bytes per dumped cell and frequency less memory now.
- `Write2XML()` crashed when no excitation or no structure was set. It now
  reports the missing part and returns false, as the file could not be read
  back anyway.

### Changed

- The Python bindings need Python 3.6 or newer (was 3.4), the oldest version
  CI still tests; CentOS 7 and Ubuntu 14.04 are no longer tested.
- Time-domain HDF5 field dumps are written in the background instead of
  stalling the timestepping, with the file kept open for the whole run; with
  the GPU engine from snapshots of the fields.
- Field dumps are evaluated on several threads, with every engine; large
  dumps, like the surfaces of an NF2FF box, slow the simulation down less.
- The end criteria is checked every Nyquist period during the first minute of
  a simulation, and only after that with the progress report every few
  seconds. A short simulation used to run on until the first report, on a fast
  machine several times longer than necessary, with correspondingly larger
  dumps; it now stops where the end criteria is met, at the same timestep on
  every machine. `--exact-endcriteria` keeps this for the whole run.
- The memory of the operator and the engine is freed at the end of the
  simulation, not when the `openEMS` object is destroyed. A post-processing in
  the same script (NF2FF, SAR) has it available, and running one `openEMS`
  object several times no longer keeps the memory of every run.
- Reading a time-domain field dump with many timesteps is much faster: its
  datasets are listed once instead of looked up one by one, which took
  quadratic time. A far-field (NF2FF) calculation from a recording of 3200
  timesteps went from 79 s to 2.4 s per frequency.
- An unknown `--engine` value (`engine=` from Python) is an error; it used to
  run the default engine without a word, e.g. for `--engine=GPU` or `cuda`.
- A mesh line within 1/1000 of the smallest mesh step of a primitive's
  surface counts as lying on it, and so as metal or material there (CSXCAD's
  drawing tolerance): for STL files, coordinates rounded in the XML or
  computed ones. A few edges more can become metal, e.g. on sloped walls; a
  slit at x=0 in the split ground plane of the StripLine2MSL tutorial is
  closed.

### Build

- `ENABLE_HIP` and `ENABLE_METAL` CMake options for the GPU backends. HIP
  (CMake 3.28) is on when a HIP compiler is found. `OPENEMS_HIP_AS_CUDA`
  compiles the same sources as CUDA, for NVIDIA GPUs without a HIP compiler
  (CMake 3.24, a CUDA compiler and the HIP headers, looked for in `ROCM_PATH`
  and `/opt/rocm`); it is the default on Windows. Metal (macOS, CMake 3.16)
  is on by default. The GPU code is built for the GPU of the build machine
  unless `CMAKE_HIP_ARCHITECTURES` / `CMAKE_CUDA_ARCHITECTURES` name others.
  Whatever a backend is missing switches it off with a warning; openEMS
  builds without it.
- C++17 is now required (was C++11), e.g. GCC 8 or newer, also for the
  Python bindings and the GPU sources (HIP/CUDA, were C++14).
- VTK is no longer needed: openEMS writes its VTK files (`.vtr`, `.vts`,
  `.vtp`) itself, in the same format as before. zlib is needed instead,
  which HDF5 usually brings along.

### Upgrade notes

- **Rebuild the Python bindings** together with openEMS: the layout of its
  array classes changed, the soname did not.

## 0.37.0 — 2026-10-09

The version scheme changed with this release: the previous release was v0.0.36,
this one is v0.37.0. The leading `0.0.` was never going anywhere, so the minor
number moved up and patch releases now have somewhere to go.

### Added

- **SAR calculation reworked**, following IEEE/IEC 62704-1: averaging is done
  once for all frequencies instead of per frequency, and is multi-threaded,
  together giving a large speedup. The new `--autorange` option restricts the
  calculation to the cells near the peak, a further speedup that does not
  guarantee finding the global peak. Available from Python
  (`sar_calculation`) and from Octave (`CalcSAR.m`, with `autoRange`,
  `numThreads` and `progress`), with a new SAR averaging tutorial.
- **Waveguide mode excitation and probe from an HDF5 mode file**, complementing
  the analytic mode functions. `WaveguidePort`/`RectWGPort` accept a
  `mode_file` argument, the `Add*WaveGuidePort` methods accept a `local_origin`,
  and `matlab/h5writemode.m` (with an Octave oct-file) writes the mode files.
  Requires the matching CSXCAD support.
- **New Python port classes**: `CircWGPort`, `CoaxialPort`, `StripLinePort`,
  `CPWPort` and `CurvePort`, with matching `Add*Port` convenience methods on
  the `openEMS` class following the existing `AddLumpedPort`/`AddMSLPort`
  pattern, plus a `check_mode_purity()` utility. Ports also store their time
  base, which makes plotting easier.
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
- **`HDF5Dump` (Python)** for reading field dumps: a dump can be inspected
  before any field data is read, `SetPlane`/`SetLine`/`SetRange`/`SetSampling`
  push the selection down into HDF5 so only the requested part is read, and
  `GetFieldAtFrequency()` runs an on-the-fly DFT so time-domain and
  frequency-domain dumps are post-processed the same way.
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
- **Optional oversampling for frequency-domain dumps and probes**, to avoid
  the aliasing onto the upper band edge that sampling the DFT accumulation at
  exactly the Nyquist rate causes. A box opts in through CSXCAD's new
  `OverSampling` property; the default is unchanged, so existing simulations
  and their performance are unaffected.
- **`--verbose`/`-vv` now reports the actual TD/FD sampling interval** for
  each probe and dump box during `SetupProcessing`.
- **`--exact-endcriteria`**, which evaluates the energy end-criteria every
  Nyquist period instead of every few seconds of wall-clock time. The run
  then stops at the same timestep regardless of machine speed and load, at
  the cost of performance — mainly useful for engine or code verification.

### Changed

- **The MPI engine was removed.** It had not compiled for years, as it used
  the C++ MPI bindings that MPI-3 dropped, it had no tests, and several
  extensions never supported it (#260). The multithreaded engine is
  unaffected. `WITH_MPI`, `--engine=MPI`, `openEMS_MPI.sh`, `RunOpenEMS_MPI`
  and `SetupMPI` are gone. See *Upgrade notes*.
- **nf2ff result format.** The far field is written as one compound complex
  dataset per frequency, `/nf2ff/E_theta/FD/f{n}` in (theta, phi) order — the
  format every other frequency-domain dump has used since HDF5 version 0.3 —
  instead of a split `f{n}_real`/`f{n}_imag` pair in (phi, theta) order.
  Octave cannot read compound complex datasets, so the Octave/Matlab interface
  keeps the old layout, requested by `CalcNF2FF` through the new `LegacyHDF5`
  attribute of the nf2ff XML file. See *Upgrade notes*.
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
  oversubscribes it. They also follow the `numThreads` option of a run: the
  field dumps, and the NF2FF and SAR calculations of the same Python process.
- Python `CalcNF2FF()` reads an existing result file by default
  (`read_cached=True`), as its docstring already stated and as the Octave
  `CalcNF2FF` (`'Mode', 0`) always did. Pass `read_cached=False` to force
  recalculation.

### Fixed

- **`DebyeMaterial` diverged** once `sum(eps_delta)` approached `epsilon`, and did
  so silently -- the NaN energy satisfied the end criteria, so the run reported a
  normal finish. Its branch carries no inductance, so integrating it like the
  Drude/Lorentz branch put the update's root on -1. It now integrates the pole's
  capacitor voltage, in its own extension pair; `eps_delta` is no longer limited,
  and a pole shorter than three timesteps is refused with a warning (#229).
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
  instead of failing the run with "failed to load" or, on Windows, "the
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
- Two knobs for comparing the output of two builds bit by bit, both off the
  default path: debug builds compile with `-ffp-contract=off` on GCC and
  Clang, and the new `ENABLE_FLUSH_TO_ZERO` CMake option can be set to `OFF`
  to keep denormals in the engines.
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
