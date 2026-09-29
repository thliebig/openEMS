#!/usr/bin/env python3
"""Validate the conducting-sheet model (and its GPU volt-ADE update) against SSE.

The Metal engine runs the conducting-sheet ADE recurrence on the GPU; the SSE
engine runs the original CPU recurrence. Both must agree to a small relative L2
(the pointwise counts are reported, as in metal_fields.py), and the Metal run
must report the offloaded edge count. No Metal feature flag is involved:
`--engine=metal` enables it.
"""
import argparse
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import time

from metal_conductingsheet import fixture
from metal_fields import compare

# SSE and Metal use different stencil evaluation order, so pointwise tolerance
# counts are reported (as in metal_fields.py) rather than fatal. A gross ADE
# error diverges far past this gate.
MAX_RELATIVE_L2 = 1e-4
TOLERANCE_RTOL = 2e-4
TOLERANCE_ATOL = 1e-6


def run(binary, model, engine, output, expect_offload):
    output.mkdir()
    start = time.perf_counter()
    env = os.environ.copy()
    if engine == 'metal':
        env['OPENEMS_METAL_FUSED_PIPELINE'] = '0'  # ADE migration is pending; legacy is explicit.
    p = subprocess.run([binary, str(model), '--engine=' + engine], cwd=output,
                       env=env, text=True,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    wall = time.perf_counter() - start
    (output / 'solver.log').write_text(p.stdout)
    if p.returncode:
        raise RuntimeError(p.stdout)
    offload = [s for s in p.stdout.splitlines() if s.startswith('Metal: ADE offload:')]
    if engine == 'metal' and expect_offload and not offload:
        raise AssertionError('conducting sheet did not offload its ADE update')
    if offload and not (engine == 'metal' and expect_offload):
        raise AssertionError('unexpected ADE offload')
    return p.stdout, wall


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--openems', required=True)
    parser.add_argument('--keep', action='store_true')
    args = parser.parse_args()
    binary = str(Path(args.openems).resolve())
    root = Path(tempfile.mkdtemp(prefix='openems-metal-ade-'))
    failures = 0
    try:
        for kind in ('polygon', 'fallback', 'pml', 'thick'):
            model = root / (kind + '.xml')
            fixture(model, kind)
            # The fallback fixture's only sheet primitive is a 3D box, so every
            # component falls back to PEC and no ADE edges remain.
            expect_offload = kind != 'fallback'
            print(kind, flush=True)
            _, sse_wall = run(binary, model, 'sse', root / (kind + '-sse'), False)
            metal_log, metal_wall = run(binary, model, 'metal',
                                        root / (kind + '-metal'), expect_offload)
            if expect_offload:
                edges = re.search(r'Metal: ADE offload: \d+ region\(s\), (\d+) active edges',
                                  metal_log).group(1)
            else:
                edges = '0'
            print('  SSE {:.3f}s, Metal {:.3f}s, {} active edges'.format(
                sse_wall, metal_wall, edges), flush=True)
            for field in ('Et.h5', 'Ht.h5'):
                stats = compare(root / (kind + '-sse') / field,
                                root / (kind + '-metal') / field,
                                TOLERANCE_RTOL, TOLERANCE_ATOL)
                print('    {}: max abs {:.3g}, relative L2 {:.3g}, pointwise over '
                      'tolerance {}'.format(field, stats['max_abs'],
                                            stats['relative_l2'], stats['failures']),
                      flush=True)
                if stats['relative_l2'] > MAX_RELATIVE_L2:
                    failures += 1
        if failures:
            raise AssertionError('SSE/Metal relative L2 exceeded {:.0e}'.format(
                MAX_RELATIVE_L2))
    finally:
        if args.keep:
            print('Kept', root)
        else:
            shutil.rmtree(root)


if __name__ == '__main__':
    main()
