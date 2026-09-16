#!/usr/bin/env python3
"""Verify the conducting-sheet (lossy conductor) setup reuses the Metal winners.

Mode '0' sets OPENEMS_METAL_PEC=0, forcing both the CPU PEC mapping and the
CSXCAD per-row fallback in Operator_Ext_ConductingSheet. Mode 'default' uses the
Metal geometry pass, where the same MATERIAL|METAL winners are recorded and
consumed by the extension. The two must produce bit-identical output.

Use --model /absolute/path/model.xml to run an existing CoSwitch model.
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

from metal_fields import h5_arrays, make_model


def fixture(path, kind):
    boundaries = ['PML_4', 'PML_4', 'PEC', 'PEC', 'PEC', 'PEC'] if kind == 'pml' else None
    frequency = 2e10 if kind == 'thick' else 1e9
    make_model(path, (20, 18, 16), 400, nonuniform=False, boundaries=boundaries,
               frequency=frequency)
    tree = ET.parse(path)
    props = tree.find('.//Properties')

    if kind == 'fallback':
        # A 3D primitive assigned to a conducting sheet must fall back to PEC,
        # exercising the fast-path fallback branch.
        cs = ET.SubElement(props, 'ConductingSheet', Name='cs',
                           Conductivity='5.8e7', Thickness='3.5e-5')
        prims = ET.SubElement(cs, 'Primitives')
        box = ET.SubElement(prims, 'Box', Priority='3')
        ET.SubElement(box, 'P1', X='4', Y='4', Z='2')
        ET.SubElement(box, 'P2', X='8', Y='8', Z='6')
    elif kind == 'thick':
        # 30 mm "copper" at 20 GHz saturates the ADE optimization table. This
        # must be reported once with a count, not per Yee component.
        cs = ET.SubElement(props, 'ConductingSheet', Name='cs',
                           Conductivity='5.8e7', Thickness='3.0e-2')
        prims = ET.SubElement(cs, 'Primitives')
        poly = ET.SubElement(prims, 'Polygon', Priority='3', NormDir='2',
                             Elevation='5')
        for x, y in ((2, 2), (17, 2), (17, 15), (2, 15)):
            ET.SubElement(poly, 'Vertex', X1=str(x), X2=str(y))
    else:
        cs = ET.SubElement(props, 'ConductingSheet', Name='cs',
                           Conductivity='5.8e7', Thickness='3.5e-5')
        prims = ET.SubElement(cs, 'Primitives')
        # Two coincident-priority polygons with diagonal edges (edge cases for
        # the winner boundaries), plus an extruded line-polygon sheet.
        for _ in range(2):
            poly = ET.SubElement(prims, 'Polygon', Priority='3', NormDir='2',
                                 Elevation='5')
            for x, y in ((2, 2), (7, 3), (15, 14), (10, 16)):
                ET.SubElement(poly, 'Vertex', X1=str(x), X2=str(y))
        lp = ET.SubElement(prims, 'LinPoly', Priority='4', NormDir='2',
                           Elevation='9', Length='6')
        for x, y in ((2, 2), (17, 2)):
            ET.SubElement(lp, 'Vertex', X1=str(x), X2=str(y))
    tree.write(path)


def run(binary, model, mode, output):
    output.mkdir()
    env = os.environ.copy()
    if mode == 'default':
        env.pop('OPENEMS_METAL_PEC', None)
    else:
        env['OPENEMS_METAL_PEC'] = mode
    start = time.perf_counter()
    p = subprocess.run([binary, str(model), '--engine=metal'], cwd=output,
                       env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    wall = time.perf_counter() - start
    (output / 'solver.log').write_text(p.stdout)
    if p.returncode:
        raise RuntimeError(p.stdout)
    fast = 'resolved geometry winners' in p.stdout
    if mode == 'default':
        if 'Metal PEC:' not in p.stdout:
            raise AssertionError('Metal geometry pass was not exercised')
        if not fast:
            raise AssertionError('conducting sheet did not use resolved winners')
        count = int(re.search(r'Metal conducting sheet: (\d+) resolved', p.stdout).group(1))
        if count == 0:
            raise AssertionError('no conducting-sheet winners were resolved')
    elif fast:
        raise AssertionError('fallback mode unexpectedly used resolved winners')
    timings = [s for s in p.stdout.splitlines()
               if s.startswith(('Metal PEC:', 'CPU PEC:', 'Time for ', 'Metal conducting sheet:'))]
    print(mode, 'wall {:.3f}s'.format(wall), '; '.join(timings), flush=True)
    return p.stdout


def compare(a, b):
    files = {p.relative_to(a) for p in a.rglob('*') if p.is_file() and p.name != 'solver.log'}
    other = {p.relative_to(b) for p in b.rglob('*') if p.is_file() and p.name != 'solver.log'}
    if files != other:
        raise AssertionError('Output file sets differ')
    for name in files:
        if name.suffix == '.h5':
            x, y = h5_arrays(a / name), h5_arrays(b / name)
            if x.keys() != y.keys():
                raise AssertionError('HDF5 datasets differ')
            for key in x:
                if (x[key].dtype != y[key].dtype or x[key].shape != y[key].shape
                        or x[key].tobytes() != y[key].tobytes()):
                    raise AssertionError(str(name) + '/' + key)
        else:
            rows = lambda p: [s for s in p.read_text().splitlines()
                              if s.strip() and not s.startswith(('#', '%'))]
            if rows(a / name) != rows(b / name):
                raise AssertionError(str(name))
    print('  Bit-identical field datasets / probe numeric outputs')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--openems', required=True)
    p.add_argument('--model', type=Path)
    p.add_argument('--keep', action='store_true')
    args = p.parse_args()
    root = Path(tempfile.mkdtemp(prefix='openems-metal-cs-'))
    try:
        cases = ['external'] if args.model else ['polygon', 'fallback', 'pml', 'thick']
        for kind in cases:
            model = args.model.resolve() if args.model else root / (kind + '.xml')
            if not args.model:
                fixture(model, kind)
            print(kind, flush=True)
            logs = {}
            for mode in ('0', 'default'):
                logs[mode] = run(str(Path(args.openems).resolve()), model, mode,
                                 root / (kind + '-' + mode))
            if kind == 'fallback':
                for mode in ('0', 'default'):
                    if 'fell back to PEC' not in logs[mode]:
                        raise AssertionError('fallback warning summary missing in ' + mode)
            if kind == 'thick':
                for mode in ('0', 'default'):
                    if 'exceed the ADE optimization table' not in logs[mode]:
                        raise AssertionError('ADE overflow summary missing in ' + mode)
            compare(root / (kind + '-0'), root / (kind + '-default'))
    finally:
        if args.keep:
            print('Kept', root)
        else:
            shutil.rmtree(root)


if __name__ == '__main__':
    main()
