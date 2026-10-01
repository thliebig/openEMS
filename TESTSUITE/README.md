# openEMS Octave/Matlab test suite

Tests for the solver and for the `matlab/` scripting layer, each one checking a
result against a closed form reference. The counterpart for the Python bindings
is `python/Tests/`.

## Running

```bash
cd TESTSUITE
octave --no-gui run_testsuite.m
```

The runner reports PASS/FAIL per check and per test, prints a summary and exits
non-zero if anything failed, so it can be used from CI. A failing test does not
abort the run.

| option | |
|---|---|
| `--engine=<e>` | run with one engine: `basic`, `sse`, `sse-compressed`, `multithreaded` |
| `--all-engines` | run every simulating test with every engine (slow) |
| `--keep` | keep the simulation folder of passed tests too |
| `--list` | list the tests without running them |
| `<name>` | run only the tests whose `group/name` contains `<name>` |

From inside an Octave or Matlab session just call it:

```matlab
run_testsuite
```

A script cannot take arguments, so to use the options above from a session set
them first -- they apply to the next run only:

```matlab
run_testsuite_args = {'--list'};
run_testsuite
```

(The options are read from `argv()` only when started as `octave
run_testsuite.m`; in a session `argv()` holds Octave's own command line.)

By default openEMS picks the fastest engine it has. Checking that the engines
agree is `enginetests/engine_compare`'s job, not something every test needs to
repeat.

## Why every test passes `--exact-endcriteria`

Without it openEMS only re-evaluates the energy end criteria every few seconds
of wall-clock time, so a run stops at a timestep that depends on machine speed
and load. That changes the length of the time record every result is computed
from, and with it the accuracy behind every tolerance here -- measured on this
suite, `combinedtests/Coax` stopped anywhere between 1054 and 1560 timesteps and
its impedance moved from 51.10 to 51.53 Ohm against a 51.44 Ohm limit. With the
flag the stopping point is a property of the model, not of the machine.

`ts_options` adds it for every test, so a standalone run and a suite run stop at
the same timestep. It is not free -- it costs a full-domain energy estimate every
Nyquist period -- so a test whose length is already fixed by `NrTS` can pass
`'ExactEndCriteria', 0` (`enginetests/engine_compare` does, via
`'EndCriteria', 0`). On this suite it was a net win anyway: the deterministic
stop cut the total runtime from 73 s to 59 s, because most tests had been
overshooting their criterion by a wide margin.

## Running a single test

Every test also runs on its own, which is what you want while debugging one:

```bash
cd combinedtests
octave --no-gui --eval "Coax"                 # or: Coax('openEMS_opts', '--engine=basic')
```

On its own a test shows the openEMS output, draws its diagnostic plots, keeps
the simulation folder if it failed, and raises an error on failure. Plots are on
by default only when there is a graphics toolkit, so a run over ssh still
reports its result instead of dying in `figure()`. Override any of it:
`'Silent', 1`, `'Plots', 0`, `'Cleanup', 0`, `'StopIfFailed', 0`.

## Groups

| group | |
|---|---|
| `unittests/` | pure `matlab/` code, no simulation — runs in milliseconds |
| `probes/` | probes and dumps against each other |
| `combinedtests/` | the full simulator against analytic results |
| `enginetests/` | the engines against each other |

The runner discovers the groups, so a new folder needs no registration. Cheap
groups run first.

## Writing a test

A test is a function file in a group folder:

```matlab
function pass = my_test(varargin)
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'helperscripts'));

opt = ts_options(varargin{:});          % the common options, see ts_options
Sim_Path = ts_sim_path(mfilename('fullpath'));
...
pass = ts_check_rel('what is checked', value, reference, tolerance);
pass = ts_check('something else', value < limit, 'got %g', value) && pass;
pass = ts_finish(opt, mfilename, pass, Sim_Path);
```

- The `addpath` line is what makes the test work on its own; the runner adds
  `helperscripts` as well, so it is harmless there.
- `ts_options` turns the key/value options (`openEMS_opts`,
  `ExactEndCriteria`, `Plots`, `Silent`, `Verbose`, `Cleanup`, `StopIfFailed`)
  into a struct with the defaults for an interactive run; the runner overrides
  them. Pass `opt.openEMS_opts` straight to `RunOpenEMS` -- `RunOpenEMS`'s third
  argument is just a string of openEMS command line options.
- `ts_check` and `ts_check_rel` print one PASS/FAIL line each. Always print the
  measured value and the tolerance — that is what makes a failure in a CI log
  readable.
- `ts_finish` prints the verdict, removes the simulation folder of a passed test
  and raises an error if the test failed and `StopIfFailed` is set.
- State the tolerance and say in a comment where it comes from. No exact
  comparisons of floating point results.

## Not covered yet

Reachable from Octave but without a test here: dispersive materials, conducting
sheets, lumped RLC elements, SAR, steady-state detection, plane wave/TFSF
excitation, and the characteristic impedance of a cylindrical coaxial line (see
the note in `combinedtests/coax_cylindrical.m`). `python/Tests/` covers most of
them on the Python side.
