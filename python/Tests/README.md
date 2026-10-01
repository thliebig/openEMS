# openEMS Python test suite

Tests for the Python bindings and, through them, for the solver. Two kinds live
side by side here:

| group | |
|---|---|
| `test_*.py` | `unittest` modules for the `python/openEMS/` layer — no simulation, milliseconds |
| `test_sim_scripts.py` | the two properties the runner needs from the scripts below, since a passing run cannot show them |
| everything else | full simulations checked against a closed-form reference — minutes |

The counterpart for the Octave/Matlab interface is `TESTSUITE/`.

## Running

```bash
cd python/Tests
python3 run_testsuite.py                    # the full-simulation tests
python3 run_testsuite.py --unittests        # those plus the unittest modules
python3 run_testsuite.py --only-unittests   # just the fast ones
```

The runner reports PASS/FAIL/ERROR per test, prints a summary and exits non-zero
if anything failed, so it can be used from CI. A failing test does not abort the
run.

| option | |
|---|---|
| `--unittests` | also run the `test_*.py` unit tests, before the simulations |
| `--only-unittests` | run only those |
| `--keep` | keep the simulation files of passed tests too |
| `--timeout=<s>` | kill a test running longer than this (default 1800, `0` disables) |
| `-v`, `--verbose` | stream the test output instead of showing it on failure only |
| `--list` | list the tests without running them |
| `<name>` | run only the tests whose `group/name` contains `<name>` |

Anything but `python/` works as the current directory — from inside it `import
openEMS` picks up the unbuilt source tree instead of the installed module.

What the runner does to a test it starts:

- It points `TMPDIR` at an empty per-test folder. Every test builds its
  `Sim_Path` below `tempfile.gettempdir()`, so this is what hands each run a
  clean folder and keeps it from reading the files of the previous one. The
  folder is removed after a pass and kept, with its path printed, after a
  failure.
- It sets `MPLBACKEND=Agg`, so a test left with its debugging plots switched on
  cannot block the whole run in `plt.show()`. It also removes
  `DISPLAY`/`WAYLAND_DISPLAY`, which stops a switched-on AppCSXCAD block from
  waiting for a window nobody is there to close -- on X11 and Wayland. Windows
  and macOS have no such variable, so there the source-level guard in
  `test_sim_scripts.py` is what keeps a viewer out of a batch run. A single test
  started by hand keeps its display and its viewer either way.
- The verdict comes from the exit status: a failed `assert` is a FAIL, any other
  exception or a crash is an ERROR.

## Running a single test

Every test also runs on its own, which is what you want while debugging one:

```bash
python3 Coax.py
```

On its own it shows the openEMS output and leaves its simulation folder in the
system temp directory. Set the `if 0:` block at the end to `1` for the
diagnostic plots.

## Writing a test

A full-simulation test is a plain script, not a `unittest` case — the setup
reads like any other openEMS script, which is half the point of having it:

```python
Sim_Path = os.path.join(tempfile.gettempdir(), 'My_Test')
...
FDTD.Run(Sim_Path, cleanup=True, exact_endcriteria=True)
...
print('max(dB(S11)) = {:.1f} dB'.format(value))
assert value < -20, 'FAIL: max(dB(S11)) = {:.1f} dB, expected < -20 dB'.format(value)
print('PASS')
```

- Put `Sim_Path` below `tempfile.gettempdir()`, so the runner can place and
  clear it.
- Print the measured value and the tolerance in the assertion message — that is
  what makes a failure in a CI log readable.
- State the tolerance and say in a comment where it comes from. No exact
  comparisons of floating point results.
- Pass `exact_endcriteria=True` to `FDTD.Run`, so the run stops at a timestep
  that is a property of the model and not of the machine. `TESTSUITE/README.md`
  has the measurements behind that.
- Keep plotting code behind `if 0:` and import `matplotlib` inside that block,
  so the test does not need it installed. The same goes for a geometry check:
  `Write2XML` plus `AppCSXCAD_BIN` behind `if 0:` is one character away when you
  need it, and never leaves a suite run waiting for a window
  (`MSL_With_Local_Absorbers.py` has both).
- A unit test for the Python layer alone goes into a `test_*.py` module instead;
  it runs in milliseconds and in CI on every platform.

`test_sim_scripts.py` enforces the last two points: it fails on a script that
has an `if 1:` block committed, or that contains no assertion at all. Neither
shows up in a passing run — a viewer only stalls the suite on a machine that has
a display, and a script without a check passes by not crashing.

## Other platforms

The runner itself is portable: `TMPDIR`/`TEMP`/`TMP` is how `tempfile` finds its
directory on every platform, the subprocesses are started from an argument list
rather than a shell, and captured output is decoded with replacement so a
console encoding cannot end a run. Of the behaviour above, only the
`DISPLAY`/`WAYLAND_DISPLAY` removal is specific to X11 and Wayland.

The full-simulation tests themselves have only ever been run on Linux, though --
CI runs `test_*.py` everywhere but no simulation, so expect to adjust a
tolerance when you first run them on Windows or macOS.

## No engine selection

The Octave suite can sweep a test over all four FDTD engines; this one cannot,
because every script here passes its own options to `FDTD.Run` and the
`openEMS` extension type cannot be patched from the outside. Checking that the
engines agree is `TESTSUITE/enginetests/engine_compare`'s job anyway, and it is
not specific to the language bindings.
