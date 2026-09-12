#!/usr/bin/env python3
"""Verify Metal PEC winners, used flags, and complete output against CPU mapping.

Use --model /absolute/path/model.xml for an existing CoSwitch simulation.
All modes use Metal field updates, isolating geometry mapping differences.
"""
import argparse
import os
import re
from pathlib import Path
import shutil
import subprocess
import tempfile
import time
import xml.etree.ElementTree as ET

from metal_fields import h5_arrays, make_model
import numpy as np
from openEMS.ports import UI_data


def fixture(path, kind):
    make_model(path, (24, 23, 22), 300, nonuniform=True)
    tree = ET.parse(path)
    props = tree.find('.//Properties')
    # Include sheets on actual mesh lines, diagonal edges, overlapping equal
    # priorities, and tiny gaps. XML avoids dependence on binding signatures.
    sheet = ET.SubElement(props, 'Metal', Name='sheet')
    prims = ET.SubElement(sheet, 'Primitives')
    if kind == 'edge-cases':
        # Exact horizontal/vertical/diagonal hits, plus nearly coincident edges.
        for offset in (0, 1e-8):
            poly = ET.SubElement(prims, 'LinPoly', Priority='3', NormDir='2',
                                 Elevation='0', Length='22')
            for x, y in ((0, 0), (24, 0), (24, 22), (0, 22)):
                ET.SubElement(poly, 'Vertex', X1=str(x+offset), X2=str(y))
        mesh = tree.find('.//RectilinearGrid')
        for axis, count in zip('XYZ', (24, 23, 22)):
            mesh.find(axis+'Lines').text = ','.join(str(i) for i in range(count+1))
    elif kind == 'boxes':
        for z in ('0', '10'):
            box = ET.SubElement(prims, 'Box', Priority='3')
            ET.SubElement(box, 'P1', X='2', Y='2', Z=z)
            ET.SubElement(box, 'P2', X='20', Y='20', Z=z)
    elif kind == 'cylinders':
        # Via-like z cylinders with exact and off-grid walls, plus a shell.
        for x, y, r in ((6, 6, 1), (12, 12, 2), (18, 6, 0.5), (6.5, 18.5, 1.5)):
            cyl = ET.SubElement(prims, 'Cylinder', Priority='5', Radius=str(r))
            ET.SubElement(cyl, 'P1', X=str(x), Y=str(y), Z='0')
            ET.SubElement(cyl, 'P2', X=str(x), Y=str(y), Z='22')
        # Shell keeps only points with |dist - radius| <= ShellWidth/2.
        shell = ET.SubElement(prims, 'CylindricalShell', Priority='5',
                              Radius='3', ShellWidth='1')
        ET.SubElement(shell, 'P1', X='18', Y='18', Z='0')
        ET.SubElement(shell, 'P2', X='18', Y='18', Z='22')
        mesh = tree.find('.//RectilinearGrid')
        for axis, count in zip('XYZ', (24, 23, 22)):
            mesh.find(axis+'Lines').text = ','.join(str(i) for i in range(count+1))
    else:
        for offset in (0, 5.000001):
            poly = ET.SubElement(prims, 'Polygon' if kind == 'sheets' else 'LinPoly',
                                 Priority='3', NormDir='2', Elevation='0', Length='16')
            for x, y in ((2, 2), (7, 2), (17, 18), (12, 18)):
                ET.SubElement(poly, 'Vertex', X1=str(x+offset), X2=str(y))
        box = ET.SubElement(prims, 'Box', Priority='3')
        ET.SubElement(box, 'P1', X='8', Y='8', Z='0')
        ET.SubElement(box, 'P2', X='12', Y='12', Z='12')
    material = ET.SubElement(props, 'Material', Name='overlap')
    ET.SubElement(material, 'Property', Epsilon='2')
    mp = ET.SubElement(material, 'Primitives')
    box = ET.SubElement(mp, 'Box', Priority='3')
    ET.SubElement(box, 'P1', X='9', Y='9', Z='0')
    ET.SubElement(box, 'P2', X='15', Y='15', Z='15')
    if kind == 'fallback':
        cyl = ET.SubElement(prims, 'Cylinder', Priority='5', Radius='1')
        ET.SubElement(cyl, 'P1', X='10', Y='10', Z='0')
        ET.SubElement(cyl, 'P2', X='10', Y='10', Z='20')
        transform = ET.SubElement(box, 'Transformation')
        ET.SubElement(transform, 'Translate', Argument='1,0,0')
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
    wall = time.perf_counter()-start
    (output/'solver.log').write_text(p.stdout)
    if p.returncode:
        raise RuntimeError(p.stdout)
    if mode != '0' and 'Metal PEC:' not in p.stdout:
        raise AssertionError('Metal PEC was not exercised')
    if mode == 'verify' and 'all winners CPU-verified' not in p.stdout:
        raise AssertionError('CPU winner verification did not complete')
    warnings = [s for s in p.stdout.splitlines() if 'Unused primitive' in s]
    timings = [s for s in p.stdout.splitlines() if s.startswith(('Metal PEC:', 'CPU PEC:', 'Time for '))]
    print(mode, 'wall {:.3f}s'.format(wall), '; '.join(timings), flush=True)
    return warnings


def compare(a, b):
    files = {p.relative_to(a) for p in a.rglob('*') if p.is_file() and p.name != 'solver.log'}
    other = {p.relative_to(b) for p in b.rglob('*') if p.is_file() and p.name != 'solver.log'}
    if files != other:
        raise AssertionError('Output file sets differ')
    for name in files:
        if name.suffix == '.h5':
            x, y = h5_arrays(a/name), h5_arrays(b/name)
            if x.keys() != y.keys():
                raise AssertionError('HDF5 datasets differ')
            for key in x:
                if x[key].dtype != y[key].dtype or x[key].shape != y[key].shape or x[key].tobytes() != y[key].tobytes():
                    raise AssertionError(str(name) + '/' + key)
        else:
            # Probe headers include creation timestamps; compare all numeric rows.
            rows = lambda p: [s for s in p.read_text().splitlines() if s.strip() and not s.startswith(('#', '%'))]
            if rows(a/name) != rows(b/name):
                raise AssertionError(str(name))
    print('  Bit-identical field datasets / probe numeric outputs; identical unused warnings')
    # MSL voltage center / spatially averaged current, fixed 50-ohm reference.
    # This compares the available excitation column, not a full multi-run S matrix.
    ports = sorted(int(m.group(1)) for p in a.glob('port_ut_*B')
                   if (m := re.fullmatch(r'port_ut_(\d+)B', p.name)))
    if 0 in ports:
        def s_column(path):
            freq = np.linspace(1e9, 10e9, 101)
            waves = []
            for port in ports:
                u = UI_data(['port_ut_{}B'.format(port)], str(path), freq).ui_f_val[0]
                i = UI_data(['port_it_{}A'.format(port), 'port_it_{}B'.format(port)],
                            str(path), freq).ui_f_val
                current = 0.5*(i[0]+i[1])
                waves.append((0.5*(u+50*current), 0.5*(u-50*current)))
            incident = waves[ports.index(0)][0]
            if not np.all(np.isfinite(incident)) or np.any(incident == 0):
                raise AssertionError('Invalid incident spectrum')
            return np.asarray([w[1]/incident for w in waves])
        x, y = s_column(a), s_column(b)
        if not np.all(np.isfinite(x)) or x.tobytes() != y.tobytes():
            raise AssertionError('50-ohm S-parameter column differs')
        print('  50-ohm S[:,0], 1–10 GHz: bit-identical (finite-run comparison)')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--openems', required=True)
    p.add_argument('--model', type=Path)
    p.add_argument('--keep', action='store_true')
    args = p.parse_args()
    root = Path(tempfile.mkdtemp(prefix='openems-metal-pec-'))
    try:
        cases = ['external'] if args.model else ['boxes', 'sheets', 'extruded',
                                                 'edge-cases', 'fallback', 'cylinders']
        for kind in cases:
            model = args.model.resolve() if args.model else root/(kind+'.xml')
            if not args.model:
                fixture(model, kind)
            print(kind, flush=True)
            warnings = []
            for mode in ('0', 'verify', 'default'):
                warnings.append(run(str(Path(args.openems).resolve()), model, mode, root/(kind+'-'+mode)))
            if warnings[0] != warnings[1] or warnings[0] != warnings[2]:
                raise AssertionError('Primitive-used warnings differ')
            if kind == 'cylinders':
                # Cylinders and shells must be flattened for the GPU, not left
                # as unsupported primitives that force CPU refinement.
                for mode in ('verify', 'default'):
                    log = (root/(kind+'-'+mode)/'solver.log').read_text()
                    if 'unsupported primitives' in log:
                        raise AssertionError('Cylinders/shells not mapped on the GPU')
            compare(root/(kind+'-0'), root/(kind+'-verify'))
            compare(root/(kind+'-0'), root/(kind+'-default'))
    finally:
        if args.keep:
            print('Kept', root)
        else:
            shutil.rmtree(root)


if __name__ == '__main__':
    main()
