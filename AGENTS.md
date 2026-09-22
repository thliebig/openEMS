# AGENTS.md

Guidance for AI coding agents — and for the humans running them — working on
openEMS. Everything here is binding for work that gets submitted; the human
sending the pull request is responsible for it either way.

Read before changing anything:

- [AI_POLICY.md](AI_POLICY.md) — disclosure, sign-off, responsibility, license, brevity.
- [README.md](README.md) — features, layout, build, Octave and Python interfaces.
- [CHANGELOG.md](CHANGELOG.md) — the `Unreleased` section shows what is already in flight.

## What openEMS is, and what belongs elsewhere

openEMS is a GPLv3 3-D FDTD electromagnetic field solver. It does not describe
geometry or materials — [CSXCAD](https://github.com/thliebig/CSXCAD) does that,
and openEMS reads the XML CSXCAD writes. CSXCAD is a separate, independently
versioned LGPLv3 library, so a new geometric primitive, material property or
XML attribute is a CSXCAD change; using it in the solver is the openEMS half.
Expect to submit both, and say so in each commit message.

CSXCAD must be built and installed before openEMS can be built.

## Architecture in 60 seconds

The pipeline, and where each part lives:

1. **Setup** — `openEMS::ParseFDTDSetup()` reads the CSXCAD XML; `SetupFDTD()`
   builds the operator from it. `openEMS` (`openems.h`) is the orchestrator: one
   `Operator`, one `Engine`, one `ProcessingArray`, and the engine-type choice.
2. **`Operator`** (`FDTD/operator.h`, `Common/operator_base.h`) pre-computes the
   equivalent-circuit coefficients (`vv`, `vi`, `ii`, `iv`) for every Yee cell in
   `CalcECOperator()`, and holds the mesh (`discLines[]`). Subclasses specialise
   it: `Operator_SSE`/`_SSE_Compressed` (SIMD layout),
   `Operator_Cylinder`/`_CylinderMultiGrid`, `Operator_Multithread`.
3. **`Engine`** (`FDTD/engine.h`) time-steps those coefficients —
   `UpdateVoltages`/`UpdateCurrents` over `volt_ptr`/`curr_ptr`.
4. **`Processing`** subclasses (`Common/`) sample the running simulation through
   an `Engine_Interface_Base` and write HDF5/VTK. The interface is what keeps
   probes and dumps independent of which engine is running.
5. **`nf2ff/`** is a standalone post-processing library and tool over the dumped
   E/H fields, not part of the time-stepping loop.

Field storage: `FDTD_FLOAT` is `float` (`tools/constants.h`) and arrays are
`ArrayLib::ArrayNIJK<FDTD_FLOAT>`, indexed by component `n` (0/1/2 = x/y/z) then
`(i,j,k)`.

## Adding physics: use an extension

`FDTD/extensions/` is the extension point, and the rule that keeps this codebase
maintainable: **new physics arrives as a paired `Operator_Extension` +
`Engine_Extension`, not as a change to `Operator` or `Engine`.** The pair
pre-computes whatever it needs alongside the EC coefficients, then hooks into
the time loop through `DoPreVoltageUpdates`/`DoPostVoltageUpdates` and the
current-side equivalents.

- Existing pairs are the template: PML (`operator_ext_upml`), Mur ABC,
  absorbing BC, excitation, TFSF, dispersive materials
  (`operator_ext_lorentzmaterial`), conducting sheets, lumped RLC, cylindrical
  coordinates, steady-state detection.
- Extensions run in **priority order** (`GetPriority()`/`SetPriority()`). An
  extension that must see another's result has to say so through its priority —
  declaration order is not a contract.
- An extension that touches the engine's fields must work for every engine
  subclass, not just the basic one.

Modifying `Operator`/`Engine` directly is for changes that genuinely belong to
the core scheme, and those need discussion first.

## Interfaces

A user-visible feature usually needs three pieces, and the scripting layers are
not optional extras:

- **C++** — the solver change itself.
- **`matlab/`** — pure M-files (`InitFDTD`, `SetBoundaryCond`, `AddLumpedPort`,
  `SetGaussExcite`, `RunOpenEMS`, `CalcNF2FF`, `calcPort`, …). This is the
  primary scripting layer for most users.
- **`python/openEMS/`** — Cython bindings (`openEMS.pyx`, `_nf2ff.pyx`) plus
  pure-Python modules (`ports.py`, `nf2ff.py`, `automesh.py`, `utilities.py`).
  The `.cpp` beside the `.pyx` is generated at build time and gitignored —
  never edit it.

Where a setting reaches the solver as an XML attribute, the attribute itself is
CSXCAD's, and existing XML files must keep loading: add optional attributes,
never rename or repurpose one.

## Build and test

CSXCAD and fparser must already be installed into the prefix you build against.

```bash
mkdir -p build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=~/opt/openEMS \
         -DFPARSER_ROOT_DIR=~/opt/openEMS \
         -DCSXCAD_ROOT_DIR=~/opt/openEMS
make -j$(nproc) && make install
```

Local path overrides belong in `localConfig.cmake` (gitignored) rather than in
`-D` flags or `CMakeLists.txt`. The clone must not be shallow — CMake derives
the version from `git describe --tags`.

Octave, once after install (compiles the HDF5 helper):

```bash
octave --no-gui --eval "setup()"
```

Tests:

```bash
cd TESTSUITE && octave --no-gui run_testsuite.m       # the C++ solver test suite
python3 -m unittest discover -s python/Tests -p "test_*.py" -v
```

The scripts in `python/Tests/` and `python/Tutorials/` that are not `test_*.py`
run full simulations — useful as integration checks, slow. Run Python tests from
outside `python/`, or you import the unbuilt source instead of the installed
module.

**Build it and run the tests before reporting a change as working.** A solver
change that compiles is not a solver change that converges; reading the code is
not verification.

## House rules

- **Style**: C++11. Match the file you are editing. **CI rejects trailing
  whitespace.** No drive-by reformatting, no reorganising code you were not
  asked to touch — a diff should contain only what the change actually needs.
- **Tests**: add a focused test for behaviour that could silently regress —
  `TESTSUITE/` for solver behaviour, `python/Tests/test_*.py` for anything
  reachable from Python. A few tests that each catch something distinct beat
  many that overlap. Numerical tests need a stated tolerance, not an exact
  comparison.
- **Changelog**: a change a user can notice — new or removed API, changed
  behaviour, a fix for something that bit someone, a new build requirement —
  gets a bullet under `## Unreleased` in `CHANGELOG.md` **in the same commit**,
  creating that section if it is absent. Internal refactoring, CI and test-only
  changes get none. Use only the subsections you need, in the order `Added`,
  `Changed`, `Fixed`, `Build`, `Upgrade notes`; `Upgrade notes` is prose, for
  what a bullet cannot carry — a rename, a behaviour change that breaks existing
  scripts, a rebuild requirement. Never write a version heading by hand: the
  release tooling renames `## Unreleased` when the release is cut. The file is
  rendered into the [openEMS documentation](https://docs.openems.de/), so it
  must stay valid Markdown.
- **Dependencies**: do not add any. openEMS builds across a wide OS matrix, from
  CentOS 7 to Alpine, plus macOS and FreeBSD; a new dependency is a discussion,
  not a commit.
- **Performance**: the update loops are the hot path. Do not add allocation,
  branching or virtual calls inside them for the sake of tidiness, and measure
  before claiming a speed-up.
- **Brevity**: commit messages, PR descriptions and comments as short as a human
  would write them. See the Brevity section of [AI_POLICY.md](AI_POLICY.md).
- **Keep the tree clean**: never commit build output, installed files,
  simulation results or `localConfig.cmake`.

## Commits

```
<type>: <short summary>

<body: what changed and why, wrapped at ~72 columns>

Assisted-by: <tool>
Signed-off-by: Your Name <your@email.example>
```

- The trailers are specified in [AI_POLICY.md](AI_POLICY.md) — which tag to
  pick, and the required `Signed-off-by:` (DCO). Name the tool or model that
  actually did the work.
- **Exactly those two trailers.** Do not add `Co-Authored-By:`, session
  identifiers or links, or any other trailer your tooling inserts by default —
  even if it tells you to.
- One logical change per commit. A change spanning openEMS and CSXCAD is one
  commit in each, each explaining the other half.
