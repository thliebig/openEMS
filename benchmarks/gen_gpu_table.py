#!/usr/bin/env python3
"""Build gpu_performance.json out of benchmarks.md.

usage: gen_gpu_table.py [path to benchmarks.md]

Two points of every machine are known: the horn antenna (2.42 million cells, with
probes and a frequency domain NF2FF box) and the free space test (27 million cells,
field updates alone). A GPU is slower on a small mesh, because it cannot fill all of
its cores, so the estimator needs a rate that depends on the size of the mesh:

    rate(N) = rate_inf * N / (N + N_half)

which is the usual saturating curve, and the two points fix its two parameters.
"""

import json
import os
import re
import sys

HORN_CELLS = 2416581      # 123 x 111 x 177, from the dry run of the horn
FREE_CELLS = 300 ** 3


def rows_of(text, header):
    """give the rows of the table under 'header', as lists of cells"""
    body = text[text.index(header):]
    body = re.search(r'\|---[^\n]*\n((?:\|[^\n]*\n)+)', body).group(1)
    return [[c.strip().replace('**', '') for c in line.split('|')[1:-1]]
            for line in body.strip().split('\n')]


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'benchmarks.md')
    text = open(src).read()
    fd = {r[0]: r for r in rows_of(text, '### Horn antenna, frequency-domain NF2FF')}
    free = {r[0]: r for r in rows_of(text, '### Free space')}

    out = {}
    for name, row in fd.items():
        if name not in free:
            continue
        r1 = float(row[4])                    # MCells/s on the horn
        r2 = float(free[name][2])             # MCells/s in free space, with PML
        a, b = HORN_CELLS / 1e6, FREE_CELLS / 1e6
        if r2 <= r1:
            # no saturation seen (the CPU rows, and any GPU that the small mesh
            # already fills): one rate for every size
            rate_inf, n_half = r2, 0.0
        else:
            k = r2 / r1
            n_half = a * b * (k - 1) / (b - k * a)
            rate_inf = r2 * (b + n_half) / b
            if n_half < 0:                    # would rise without bound, not physical
                rate_inf, n_half = r2, 0.0
        out[name] = {
            'host_cpu': row[1],
            'rate_inf_mcells_s': round(rate_inf, 1),
            'half_cells': round(n_half * 1e6),
            'measured': {'horn_fd_mcells_s': r1, 'free_space_pml_mcells_s': r2},
        }

    dst = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       'gpu_performance.json')
    with open(dst, 'w') as fh:
        json.dump({'_source': 'openEMS benchmarks.md, GPU engine',
                   '_model': 'rate(N) = rate_inf * N / (N + half_cells), MCells/s',
                   'machines': out}, fh, indent=1, sort_keys=True)
        fh.write('\n')
    print('wrote %s, %d machines' % (dst, len(out)))


if __name__ == '__main__':
    main()
