"""Check the estimator against the measurements it comes from.

    python3 test_estimate.py

The rate model has to reproduce the two points that fixed it, and the estimate
of a whole run has to land near the runs of benchmarks.md. The horn of the
benchmarks ran 14900 timesteps on every machine, so its rows are 32 test cases.
"""

import os
import re
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from estimate import MACHINES, Scope, estimate, rate_mcells_s  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
BENCHMARKS = os.path.join(os.path.dirname(HERE), 'benchmarks.md')
HORN_CELLS = 2416581
HORN_TIMESTEPS = 14900
FREE_CELLS = 300 ** 3

# the horn of the benchmarks, as openEMS --dry-run reports it
HORN_FD = {'cells': HORN_CELLS, 'max_timesteps': 1000000000, 'end_criteria': 1e-4,
           'time_domain_dump_bytes_per_timestep': 0,
           'frequency_domain_dump_bytes': 3960000}
HORN_TD = {'cells': HORN_CELLS, 'max_timesteps': 1000000000, 'end_criteria': 1e-4,
           'time_domain_dump_bytes_per_timestep': 329976,
           'frequency_domain_dump_bytes': 0}


def table(header):
    text = open(BENCHMARKS).read()
    body = text[text.index(header):]
    body = re.search(r'\|---[^\n]*\n((?:\|[^\n]*\n)+)', body).group(1)
    return {r.split('|')[1].strip(): [c.strip().replace('**', '')
                                      for c in r.split('|')[1:-1]]
            for r in body.strip().split('\n')}


def test_rate_reproduces_the_measurements():
    """the fitted curve has to go through both measured points"""
    for name, m in MACHINES.items():
        for cells, want in ((HORN_CELLS, m['measured']['horn_fd_mcells_s']),
                            (FREE_CELLS, m['measured']['free_space_pml_mcells_s'])):
            got = rate_mcells_s(name, cells)
            assert abs(got / want - 1) < 0.01, \
                '%s at %d cells: %.0f MCells/s, measured %.0f' % (name, cells, got, want)
    print('PASS rate model reproduces %d machines x 2 points' % len(MACHINES))


def test_rate_grows_with_the_mesh():
    """a bigger mesh cannot be slower per cell, and the curve has to saturate"""
    for name in MACHINES:
        rates = [rate_mcells_s(name, n) for n in (1e5, 1e6, 1e7, 1e8)]
        assert rates == sorted(rates), '%s: %s' % (name, rates)
        assert rates[-1] <= MACHINES[name]['rate_inf_mcells_s'] * 1.001, name
    print('PASS rate grows with the mesh and stays under its limit')


def test_against_the_benchmarks():
    """the estimate of a whole run against the runs of benchmarks.md"""
    # A run that dumps the time domain fields is as slow as the machine writes
    # them, and the rented machines differed by a factor of three there (239 to
    # 1171 MB/s implied by their runs). Thus the bound of that case is wider: the
    # estimator cannot know the write speed of a machine it has not run on.
    for header, data, what, tol, want in (
            ('### Horn antenna, frequency-domain NF2FF', HORN_FD, 'frequency domain',
             0.30, 0.85),
            ('### Horn antenna, time-domain NF2FF', HORN_TD, 'time domain',
             0.40, 0.80)):
        scope = Scope(data)
        errors = []
        for name, row in table(header).items():
            if name not in MACHINES:
                continue
            got = estimate(scope, name, timesteps=HORN_TIMESTEPS)['total_s']
            measured = float(row[2].split()[0])
            errors.append((got / measured - 1, name))
        errors.sort()
        median = statistics.median(e for e, _ in errors)
        near = sum(1 for e, _ in errors if abs(e) < tol)
        print('%-17s: median %+3.0f%%, %d of %d within %.0f%% (worst %+.0f%% %s, '
              '%+.0f%% %s)' % (what, 100 * median, near, len(errors), 100 * tol,
                               100 * errors[0][0], errors[0][1],
                               100 * errors[-1][0], errors[-1][1]))
        assert abs(median) < 0.20, '%s: median error %+.0f%%' % (what, 100 * median)
        assert near >= want * len(errors), \
            '%s: only %d of %d within %.0f%%' % (what, near, len(errors), 100 * tol)
    print('PASS estimates agree with the benchmarks')


def test_price():
    """the price follows the time, plus what the machine bills before it runs"""
    scope = Scope(HORN_FD)
    r = estimate(scope, 'RTX 5080', timesteps=HORN_TIMESTEPS, price_per_hour=0.36,
                 rent_overhead_s=300)
    assert abs(r['cost_usd'] - 0.36 * (r['total_s'] + 300) / 3600) < 1e-3, r
    two = estimate(scope, 'RTX 5080', timesteps=HORN_TIMESTEPS, excitations=2)
    one = estimate(scope, 'RTX 5080', timesteps=HORN_TIMESTEPS)
    assert abs(two["total_s"] - 2 * one["total_s"]) < 0.2, (one, two)
    print('PASS price and excitations')


if __name__ == '__main__':
    test_rate_reproduces_the_measurements()
    test_rate_grows_with_the_mesh()
    test_against_the_benchmarks()
    test_price()
    print('PASS')
