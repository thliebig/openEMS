#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""run_testsuite -- run the openEMS Python test suite

The counterpart of ``TESTSUITE/run_testsuite.m`` for the Python bindings. It
runs the full-simulation tests in this folder, reports PASS/FAIL/ERROR per test
plus a summary and exits non-zero, so it can be used from CI. A failing test
neither aborts the run nor hides the tests behind it.

usage:

    python3 run_testsuite.py [<options>] [<name>...]

options:

  --unittests     also run the test_*.py unit tests (fast, no simulation)
  --only-unittests  run only those
  --keep          keep the simulation files of passed tests too
  --timeout=<s>   kill a test that runs longer than this (default 1800,
                  0 disables)
  -v, --verbose   stream the test output instead of showing it on failure only
  --list          list the tests without running them
  <name>          run only the tests whose group/name contains <name>

Every test also runs on its own, which is what you want while debugging one:

    python3 Coax.py

Start either from anywhere but ``python/`` itself -- from there ``import
openEMS`` picks up the unbuilt source tree instead of the installed module.
``cd python/Tests`` and run it from here is fine.

openEMS testsuite
-----------------

See also the test contract in README.md
"""

import argparse
import locale
import os
import shutil
import subprocess
import sys
import tempfile
import time
from datetime import datetime

TS_DIR = os.path.dirname(os.path.abspath(__file__))

# Scripts in this folder that are not tests this runner can judge, because they
# check nothing: "it exited 0" would be the only verdict they could get. Empty,
# and meant to stay that way -- a demonstration script belongs in Examples/ or
# Tutorials/. Listing one here with its reason keeps it out of the run without
# pretending it passed.
NOT_A_TEST = {}

# groups in the order they run -- the cheap one first, so a broken installation
# shows up in seconds instead of after the first simulation
GROUPS = ('unittests', 'simulations')


class Test(object):
    """One entry of the suite: a group, a name and the command that runs it."""

    def __init__(self, group, name, cmd):
        self.group = group
        self.name  = name
        self.cmd   = cmd

    @property
    def id(self):
        return self.group + '/' + self.name


def collect(groups, patterns):
    """Find the tests of the given groups whose id matches one of patterns."""
    tests   = []
    skipped = []

    for group in GROUPS:
        if group not in groups:
            continue
        for entry in sorted(os.listdir(TS_DIR)):
            if not entry.endswith('.py') or entry.startswith('_'):
                continue
            if entry == os.path.basename(__file__):
                continue
            is_unittest = entry.startswith('test_')
            if is_unittest != (group == 'unittests'):
                continue

            name = entry[:-3]
            test_id = group + '/' + name
            if patterns and not any(p in test_id for p in patterns):
                continue

            if entry in NOT_A_TEST:
                skipped.append((test_id, NOT_A_TEST[entry]))
            elif is_unittest:
                # discover rather than the module name: it puts this folder on
                # sys.path itself, so the test runs from a cwd outside python/
                # just like the simulation tests. -v names the failing method.
                tests.append(Test(group, name,
                                  [sys.executable, '-m', 'unittest', 'discover',
                                   '-v', '-s', TS_DIR, '-p', entry]))
            else:
                tests.append(Test(group, name,
                                  [sys.executable, os.path.join(TS_DIR, entry)]))

    return tests, skipped


def child_env(work_dir):
    """The environment a test is run in."""
    env = os.environ.copy()

    # Every test builds its simulation folder below tempfile.gettempdir(), so
    # pointing that at a per-test directory is what lets the runner hand each
    # test an empty folder and drop it again afterwards -- no test run ever
    # sees the files of the previous one. (TMP/TEMP for Windows.)
    for var in ('TMPDIR', 'TMP', 'TEMP'):
        env[var] = work_dir

    # A test left with its debugging plots switched on would otherwise block the
    # whole run in plt.show().
    env['MPLBACKEND'] = 'Agg'

    # And one left with its AppCSXCAD block switched on would block it in
    # os.system(AppCSXCAD_BIN), which waits for the window to be closed. Without
    # a display Qt gives up in a fraction of a second and the test carries on, so
    # a batch run cannot be stopped by a window nobody is there to close. A
    # single test started by hand keeps its display and its viewer.
    for var in ('DISPLAY', 'WAYLAND_DISPLAY', 'QT_QPA_PLATFORM'):
        env.pop(var, None)

    env['PYTHONUNBUFFERED'] = '1'

    return env


def decode(raw):
    """A test's captured output as text, never raising over its encoding."""
    if not raw:
        return ''
    return raw.decode(locale.getpreferredencoding(False), 'replace')


def run_one(test, work_dir, timeout, verbose):
    """Run one test and return (status, seconds, output)."""
    os.makedirs(work_dir)
    env = child_env(work_dir)

    out = ''
    t0  = time.time()

    # cwd is the test's own folder, so that anything a test writes relative to
    # the current directory stays out of the source tree. A verbose run lets the
    # test write to the terminal instead of capturing its output.
    if verbose:
        proc = subprocess.Popen(test.cmd, cwd=work_dir, env=env)
    else:
        # bytes, decoded below: text mode would decode with the locale
        # encoding and raise on the first byte that does not fit, which on a
        # Windows console is a real possibility and no reason to lose the run
        proc = subprocess.Popen(test.cmd, cwd=work_dir, env=env,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT)
    try:
        out = decode(proc.communicate(timeout=timeout or None)[0])
        rc  = proc.returncode
    except subprocess.TimeoutExpired:
        proc.kill()
        out = decode(proc.communicate()[0]) + \
              '\n*** killed after {:.0f} s\n'.format(timeout)
        rc  = None

    secs = time.time() - t0

    if rc == 0 and 'Ran 0 tests' in out:
        # unittest discover is happy with a pattern that matches nothing; a
        # test file that is never collected must not pass by default
        status = 'ERROR'
    elif rc == 0:
        status = 'PASS'
    elif rc is None:
        status = 'TIMEOUT'
    elif rc == 1 and (verbose or 'AssertionError' in out):
        # a check said no -- as opposed to the test falling over on the way
        status = 'FAIL'
    else:
        status = 'ERROR'

    return status, secs, out


def tail(text, lines=20):
    """The last lines of text, indented, for the log of a failed test."""
    kept = [l for l in text.splitlines() if l.strip()][-lines:]
    return '\n'.join('    | ' + l for l in kept)


def main(argv):
    parser = argparse.ArgumentParser(
        description='run the openEMS Python test suite',
        epilog='With no <name> every full-simulation test in this folder runs.')
    parser.add_argument('patterns', metavar='name', nargs='*',
                        help='run only the tests whose group/name contains this')
    parser.add_argument('--unittests', action='store_true',
                        help='also run the test_*.py unit tests')
    parser.add_argument('--only-unittests', action='store_true',
                        help='run only the test_*.py unit tests')
    parser.add_argument('--keep', action='store_true',
                        help='keep the simulation files of passed tests too')
    parser.add_argument('--timeout', type=float, default=1800, metavar='S',
                        help='kill a test running longer than S seconds '
                             '(default 1800, 0 disables)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='stream the test output while it runs')
    parser.add_argument('--list', action='store_true', dest='list_only',
                        help='list the tests without running them')
    args = parser.parse_args(argv)

    if args.only_unittests:
        groups = ('unittests',)
    elif args.unittests:
        groups = ('unittests', 'simulations')
    else:
        groups = ('simulations',)

    tests, skipped = collect(groups, args.patterns)

    if not tests and not skipped:
        parser.error('no test matches the given name(s)')

    if args.list_only:
        print('*** openEMS Python testsuite -- available tests:')
        for test in tests:
            print('  ' + test.id)
        for test_id, why in skipped:
            print('  {}   (not run: {})'.format(test_id, why))
        return 0

    run_dir = tempfile.mkdtemp(prefix='openEMS_testsuite_')

    print('\n*** {} openEMS Python testsuite started -- {} test(s)\n'.format(
        datetime.now().strftime('%Y-%m-%d %H:%M:%S'), len(tests)))

    results = []
    total   = time.time()
    for n, test in enumerate(tests, 1):
        label = '[{:2d}/{:2d}] {:<44s} '.format(n, len(tests), test.id)
        if args.verbose:
            print(label)
        else:
            sys.stdout.write(label)
        sys.stdout.flush()

        work_dir = os.path.join(run_dir, test.name)
        status, secs, out = run_one(test, work_dir, args.timeout, args.verbose)

        if args.verbose:
            print('  ==> {}: {} ({:.1f} s)'.format(test.id, status, secs))
        else:
            print('{:<7s} {:7.1f} s'.format(status, secs))
        if status != 'PASS' and out:
            print(tail(out))

        if status == 'PASS' and not args.keep:
            shutil.rmtree(work_dir, ignore_errors=True)
        elif os.path.isdir(work_dir) and os.listdir(work_dir):
            # what the test left behind is what one needs to look at next --
            # unless it got no further than writing nothing at all
            print('    simulation files kept in {}'.format(work_dir))
        else:
            shutil.rmtree(work_dir, ignore_errors=True)

        results.append((test.id, status, secs))
        sys.stdout.flush()

    # ------------------------------------------------------------- summary
    counts = {}
    for test_id, status, secs in results:
        counts[status] = counts.get(status, 0) + 1

    # Only what the lines above do not already say: which tests went wrong, so
    # they can be read off without scrolling back through the output of the ones
    # that failed, plus the tests that never ran at all.
    problems = [r for r in results if r[1] != 'PASS']
    if problems or skipped:
        print('')
        for test_id, status, secs in problems:
            print('  {:<7s} {:<44s}  {:7.1f} s'.format(status, test_id, secs))
        for test_id, why in skipped:
            print('  {:<7s} {:<44s}  {}'.format('SKIP', test_id, '(' + why + ')'))

    failed = len(results) - counts.get('PASS', 0)
    print('\n  {} passed, {} failed, {} errors, {} skipped in {:.1f} s'.format(
        counts.get('PASS', 0),
        counts.get('FAIL', 0),
        counts.get('ERROR', 0) + counts.get('TIMEOUT', 0),
        len(skipped),
        time.time() - total))

    if failed:
        print('\n*** TESTSUITE FAILED')
    else:
        print('\n*** ALL TESTS PASSED')

    # empty unless a failed test (or --keep) left something behind
    try:
        os.rmdir(run_dir)
    except OSError:
        pass

    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
