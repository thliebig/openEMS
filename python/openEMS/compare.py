# -*- coding: utf-8 -*-
#
# Copyright (C) 2026 Thorsten Liebig (Thorsten.Liebig@gmx.de)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
"""Compare the results of two openEMS runs.

Typical use is checking one engine against another: run the same simulation
with, e.g., ``engine='multithreaded'`` and ``engine='gpu'`` into two
simulation folders, then::

    python -m openEMS.compare <sim_path_a> <sim_path_b>

Run both with ``exact_endcriteria=True``, otherwise they may stop after a
different number of timesteps and every probe differs in length.
"""

import os
import re
import sys
import argparse

import numpy as np
import h5py


def deviation(a, b, peak=None):
    """Maximum deviation of `b` from `a`, relative to `peak`.

    `peak` defaults to the peak absolute value of `a`.  Returns 0 for identical
    arrays and ``inf`` if the shapes differ.
    """
    a = np.asarray(a)
    b = np.asarray(b)
    if a.shape != b.shape:
        return np.inf
    if np.array_equal(a, b):
        return 0.0
    if peak is None:
        peak = np.max(np.abs(a))
    return np.max(np.abs(a - b)) / (peak if peak > 0 else 1.0)


def _peak(a):
    return np.max(np.abs(a)) if np.size(a) else 0.0


def _peak_group(name):
    """The datasets of a HDF5 file that are compared relative to a common peak.

    These are the timesteps or frequencies of a field dump, and both
    polarisations of a far field: the cross-polarisation is often only
    numerical noise, which has no meaningful peak of its own.
    """
    if name.startswith('FieldData'):
        return 'FieldData'
    if name.startswith(('nf2ff/E_theta', 'nf2ff/E_phi')):
        return 'nf2ff/E'
    if name.startswith('nf2ff/P_rad'):
        return 'nf2ff/P_rad'
    return None


def _read_probe(fn):
    """The axis (time or frequency) and the value columns of a probe file."""
    data = np.loadtxt(fn, comments='%', ndmin=2)
    return data[:, 0], data[:, 1:]


def _read_vtr(fn):
    """The point data arrays of a VTK rectilinear grid dump, by name."""
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy
    reader = vtk.vtkXMLRectilinearGridReader()
    reader.SetFileName(fn)
    reader.Update()
    pd = reader.GetOutput().GetPointData()
    return {pd.GetArrayName(n): vtk_to_numpy(pd.GetArray(n)) for n in range(pd.GetNumberOfArrays())}


def _compare(path_a, path_b):
    """List of ``(group, name, deviation)``, see `compare_sim_results`."""
    files = sorted(f for f in os.listdir(path_a) if os.path.isfile(os.path.join(path_a, f)))
    result = []

    # the VTK dumps are one file per timestep, a series shares its peak
    series = {}
    for f in files:
        if f.endswith('.vtr'):
            series.setdefault(re.sub(r'_\d+\.vtr$', '', f), []).append(f)
    if series:
        try:
            import vtk  # noqa: F401
        except ImportError:
            print('compare: vtk not available, skipping {} *.vtr files'.format(sum(map(len, series.values()))),
                  file=sys.stderr)
            series = {}
    for group, names in series.items():
        diffs = []
        peak = 0.0
        for f in names:
            fb = os.path.join(path_b, f)
            if not os.path.isfile(fb):
                diffs.append((f, np.inf))
                continue
            a, b = _read_vtr(os.path.join(path_a, f)), _read_vtr(fb)
            for n, arr in a.items():
                peak = max(peak, _peak(arr))
                diffs.append((f, deviation(arr, b[n], 1.0) if n in b else np.inf))
        result += [(group + '_*.vtr', f, d / (peak if peak > 0 else 1.0)) for f, d in diffs]

    for f in files:
        fa, fb = os.path.join(path_a, f), os.path.join(path_b, f)
        if '.' not in f:
            if not os.path.isfile(fb):
                result.append((f, f, np.inf))
                continue
            try:
                (xa, a), (xb, b) = _read_probe(fa), _read_probe(fb)
            except ValueError:
                continue    # not a probe file
            # the time or frequency axis has to agree exactly
            result.append((f, f, deviation(a, b) if np.array_equal(xa, xb) else np.inf))
        elif f.endswith('.h5'):
            if not os.path.isfile(fb):
                result.append((f, f, np.inf))
                continue
            with h5py.File(fa, 'r') as a, h5py.File(fb, 'r') as b:
                names = []
                a.visititems(lambda name, obj: names.append(name) if isinstance(obj, h5py.Dataset) else None)
                # datasets of one kind share their peak
                peak = {}
                for n in names:
                    key = _peak_group(n)
                    if key is not None:
                        peak[key] = max(peak.get(key, 0.0), _peak(a[n][()]))
                for n in names:
                    if n not in b:
                        result.append((f, f'{f}:{n}', np.inf))
                    else:
                        result.append((f, f'{f}:{n}', deviation(a[n][()], b[n][()], peak.get(_peak_group(n)))))

    # an output that only the second run wrote is no agreement either
    for f in sorted(os.listdir(path_b)):
        fb = os.path.join(path_b, f)
        if (f in files) or not os.path.isfile(fb):
            continue
        if f.endswith('.h5') or f.endswith('.vtr'):
            result.append((f, f, np.inf))
        elif '.' not in f:
            try:
                _read_probe(fb)
            except ValueError:
                continue    # not a probe file
            result.append((f, f, np.inf))
    return result


def compare_sim_results(path_a, path_b):
    """Compare all probe files and field dumps of two simulation folders.

    Probe and port files (the files without an extension, time and frequency
    domain) have to agree exactly in their time or frequency axis; the values
    are compared relative to their peak.  In a field dump, all field datasets
    are compared relative to the peak of the whole dump, so that late, decayed
    timesteps are not judged against their own tiny peak.  That holds for the
    datasets of a HDF5 dump and for the files of a VTK (``*.vtr``) dump
    series; the latter are only compared if the `vtk` module is available.
    In a far field result (NF2FF), both polarisations share their peak.

    Parameters
    ----------
    path_a, path_b : str
        The simulation folders; `path_a` is the reference.

    Returns
    -------
    list of (str, float)
        One ``(name, deviation)`` per probe file, per HDF5 dump dataset
        (named ``file.h5:dataset``) and per array of a VTK dump file.  A file
        or dataset missing in `path_b`, and a probe file or field dump that only
        `path_b` has, has deviation ``inf``.
    """
    return [(name, dev) for _, name, dev in _compare(path_a, path_b)]


def main(argv=None):
    parser = argparse.ArgumentParser(prog='python -m openEMS.compare',
                                     description='Compare the probe files and field dumps of two openEMS runs.')
    parser.add_argument('path_a', help='reference simulation folder')
    parser.add_argument('path_b', help='simulation folder to check')
    parser.add_argument('--rtol', type=float, default=1e-4,
                        help='max. deviation relative to the peak value (default: %(default)g)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='list every dump dataset and file, not only the worst per dump')
    args = parser.parse_args(argv)

    result = _compare(args.path_a, args.path_b)
    if not result:
        print(f'no probe files or field dumps in {args.path_a}')
        return 1

    if args.verbose:
        rows = [(name, dev) for _, name, dev in result]
    else:
        worst = {}
        for group, _, dev in result:
            worst[group] = max(worst.get(group, 0.0), dev)
        rows = list(worst.items())
    width = max(len(n) for n, _ in rows)
    for name, dev in rows:
        print(f'{name:{width}s}  {dev:.2e}{"  FAIL" if dev > args.rtol else ""}')

    failed = [n for _, n, d in result if d > args.rtol]
    print(f'{len(result) - len(failed)}/{len(result)} within rtol={args.rtol:g}')
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
