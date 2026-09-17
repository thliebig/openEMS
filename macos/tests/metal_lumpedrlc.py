#!/usr/bin/env python3
"""Validate the in-kernel lumped RLC update against the CPU engine.

The Metal diamond wavefront folds the lumped RLC recurrence into the fused
E/H kernel; the SSE engine runs the original CPU extension. Both must agree to
a small relative L2, and the Metal run must stay on the diamond path and report
the offloaded edge count. No Metal feature flag is involved: ``--engine=metal``
enables the offload.
"""
import argparse
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import time
import xml.etree.ElementTree as ET

from metal_fields import compare, make_model

# SSE and Metal evaluate the stencil in a different order, so this is a
# correctness gate rather than a bit-identity check (as in metal_fields.py).
MAX_RELATIVE_L2 = 1e-4
TOLERANCE_RTOL = 2e-4
TOLERANCE_ATOL = 1e-6


def fixture(path):
    """A PEC cavity with series RLC, parallel RLC and series RC elements.

    PEC boundaries keep the run on the primary diamond path (no UPML, no ADE).
    """
    make_model(path, (24, 20, 18), 400, nonuniform=False, boundaries=['PEC'] * 6,
               frequency=1e9)
    tree = ET.parse(path)
    props = tree.find('.//Properties')

    elements = [
        # name, x0,y0,z0, x1,y1,z1, R, L, C, LEtype
        ('ser_rlc', 6.0, 9.5, 6.0, 6.5, 10.5, 12.0, 25.0, 8e-9, 3.2e-12, 1),
        ('par_rlc', 17.0, 9.5, 6.0, 17.5, 10.5, 12.0, 400.0, 12e-9, 8.0e-12, 0),
        ('ser_rc', 11.0, 5.0, 6.0, 11.5, 6.0, 12.0, 30.0, 'nan', 2.0e-12, 1),
    ]
    for name, x0, y0, z0, x1, y1, z1, r, l, c, letype in elements:
        le = ET.SubElement(props, 'LumpedElement', Name=name, Direction='2',
                           R=str(r), L=str(l), C=str(c), LEtype=str(letype))
        prims = ET.SubElement(le, 'Primitives')
        box = ET.SubElement(prims, 'Box', Priority='10')
        ET.SubElement(box, 'P1', X=str(x0), Y=str(y0), Z=str(z0))
        ET.SubElement(box, 'P2', X=str(x1), Y=str(y1), Z=str(z1))
    tree.write(path)


def run(binary, model, engine, output):
    output.mkdir()
    env = os.environ.copy()
    env.pop('OPENEMS_METAL_FUSED_PIPELINE', None)
    start = time.perf_counter()
    p = subprocess.run([binary, str(model), '--engine=' + engine], cwd=output,
                       env=env, text=True, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT)
    wall = time.perf_counter() - start
    (output / 'solver.log').write_text(p.stdout)
    if p.returncode:
        raise RuntimeError(p.stdout)
    return p.stdout, wall


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--openems', required=True)
    parser.add_argument('--keep', action='store_true')
    args = parser.parse_args()
    binary = str(Path(args.openems).resolve())
    root = Path(tempfile.mkdtemp(prefix='openems-metal-rlc-'))
    try:
        model = root / 'model.xml'
        fixture(model)
        _, sse_wall = run(binary, model, 'sse', root / 'sse')
        metal_log, metal_wall = run(binary, model, 'metal', root / 'metal')
        if 'Metal: in-place diamond E/H pipeline: enabled' not in metal_log:
            raise AssertionError('lumped RLC did not use the diamond wavefront')
        match = re.search(r'Metal: lumped RLC offload: (\d+) active edges', metal_log)
        if not match or int(match.group(1)) == 0:
            raise AssertionError('lumped RLC was not offloaded to the GPU')
        print('SSE {:.3f}s, Metal diamond {:.3f}s, {} active edges'.format(
            sse_wall, metal_wall, match.group(1)), flush=True)
        for field in ('Et.h5', 'Ht.h5'):
            stats = compare(root / 'sse' / field, root / 'metal' / field,
                            TOLERANCE_RTOL, TOLERANCE_ATOL)
            print('  {}: max abs {:.3g}, relative L2 {:.3g}, pointwise over '
                  'tolerance {}'.format(field, stats['max_abs'],
                                        stats['relative_l2'], stats['failures']),
                  flush=True)
            if stats['relative_l2'] > MAX_RELATIVE_L2:
                raise AssertionError('SSE/Metal relative L2 exceeded {:.0e}'.format(
                    MAX_RELATIVE_L2))
    finally:
        if args.keep:
            print('Kept', root)
        else:
            shutil.rmtree(root)


if __name__ == '__main__':
    main()
