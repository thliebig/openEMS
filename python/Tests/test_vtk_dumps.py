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
"""
The VTK files openEMS writes: field dumps against the HDF5 dump of the same
box, and the material, operator and PEC debug dumps, on a Cartesian and a
cylindrical mesh. Two short runs of 60 timesteps. The files are read with the
standard library, and also with VTK where it is installed.
"""

import base64
import glob
import os
import shutil
import tempfile
import unittest
import xml.etree.ElementTree as ET
import zlib

import h5py
import numpy as np

from CSXCAD import ContinuousStructure
from openEMS.openEMS import openEMS

try:
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy
    HAS_VTK = True
except ImportError:
    HAS_VTK = False

DEBUG = dict(debug_material=True, debug_PEC=True, debug_operator=True, verbose=0)


def _run_cartesian(path):
    # a lossy substrate, a metal plate on it, fed by a lumped port
    fdtd = openEMS(NrTS=60, EndCriteria=0)
    fdtd.SetGaussExcite(2e9, 1e9)
    fdtd.SetBoundaryCond(['PEC', 'PEC', 'PEC', 'PEC', 'MUR', 'MUR'])
    csx = ContinuousStructure()
    fdtd.SetCSX(csx)
    mesh = csx.GetGrid()
    mesh.SetDeltaUnit(1e-3)
    mesh.AddLine('x', np.linspace(-20, 20, 29))
    mesh.AddLine('y', np.linspace(-15, 15, 11))
    mesh.AddLine('z', np.concatenate([np.linspace(-20, 0, 7), [1.5, 4, 8, 13, 20]]))
    csx.AddMaterial('sub', epsilon=3.6, kappa=0.05, density=1000).AddBox([-20, -15, -5], [20, 15, 0])
    csx.AddMetal('plate').AddBox([-10, -8, 0], [10, 8, 0], priority=10)
    fdtd.AddLumpedPort(1, 50, [0, 0, -5], [0, 0, 0], 'z', excite=1)
    csx.AddDump('Et').AddBox([-20, -15, -20], [20, 15, 20])
    csx.AddDump('Et_h5', file_type=1).AddBox([-20, -15, -20], [20, 15, 20])
    csx.AddDump('Hf', dump_type=11, frequency=[2e9]).AddBox([-20, -15, -2], [20, 15, -2])
    csx.AddDump('SAR', dump_type=20, frequency=[2e9]).AddBox([-20, -15, -5], [20, 15, 0])
    csx.AddDump('SAR_h5', dump_type=20, frequency=[2e9], file_type=1).AddBox([-20, -15, -5], [20, 15, 0])
    fdtd.Run(path, **DEBUG)


def _run_cylindrical(path):
    # a dielectric ring and a metal wall in a radial waveguide
    fdtd = openEMS(NrTS=60, EndCriteria=0, CoordSystem=1)
    fdtd.SetGaussExcite(2e9, 1e9)
    fdtd.SetBoundaryCond(['PEC', 'PEC', 'PEC', 'PEC', 'MUR', 'MUR'])
    csx = ContinuousStructure(CoordSystem=1)
    fdtd.SetCSX(csx)
    mesh = csx.GetGrid()
    mesh.SetDeltaUnit(1e-3)
    mesh.AddLine('r', np.linspace(2, 20, 10))
    mesh.AddLine('a', np.linspace(0, 2 * np.pi, 17))
    mesh.AddLine('z', np.linspace(0, 40, 11))
    csx.AddMaterial('ring', epsilon=2.2).AddCylindricalShell([0, 0, 10], [0, 0, 20], 10, 4)
    csx.AddMetal('wall').AddBox([8, 0, 32], [14, 1, 32], priority=10)
    csx.AddExcitation('exc', exc_type=0, exc_val=[1, 0, 0]).AddBox([2, 0, 0], [20, 2 * np.pi, 0])
    csx.AddDump('Et').AddBox([2, 0, 0], [20, 2 * np.pi, 40])
    csx.AddDump('Et_h5', file_type=1).AddBox([2, 0, 0], [20, 2 * np.pi, 40])
    fdtd.Run(path, **DEBUG)


_TYPES = {'Float32': '<f4', 'Float64': '<f8', 'UInt8': 'u1', 'Int64': '<i8', 'UInt32': '<u4'}


def _decode(element):
    """ The values of a DataArray of VTK's zlib compressed binary format: a base64 header of
    UInt32 (blocks, block size, last block size, compressed size of each block), then the
    base64 of the zlib compressed blocks. """
    text = element.text.strip()
    nblocks = int(np.frombuffer(base64.b64decode(text[:16]), '<u4')[0])
    head_len = 4 * -(-(3 + nblocks) * 4 // 3)   # base64 length of the header
    header = np.frombuffer(base64.b64decode(text[:head_len]), '<u4')
    data = base64.b64decode(text[head_len:])
    raw, pos = b'', 0
    for size in header[3:]:
        raw += zlib.decompress(data[pos:pos + size])
        pos += size
    assert pos == len(data) and len(raw) == max(nblocks - 1, 0) * header[1] + header[2], 'wrong block sizes'
    values = np.frombuffer(raw, _TYPES[element.get('type')])
    components = int(element.get('NumberOfComponents', 1))
    return values.reshape(-1, components) if components > 1 else values


def _read(fn):
    """ A VTK XML file as dict, read with the standard library: the arrays by section and name,
    the extent, the coordinates or points, and the lines. """
    root = ET.parse(fn).getroot()
    assert root.tag == 'VTKFile', fn
    assert root.get('byte_order') == 'LittleEndian' and root.get('header_type') == 'UInt32', fn
    assert root.get('compressor') == 'vtkZLibDataCompressor', fn
    data = root.find(root.get('type'))
    piece = data.find('Piece')
    out = {'type': root.get('type')}
    if piece.get('Extent'):
        out['extent'] = [int(v) for v in piece.get('Extent').split()]
    for section in ('PointData', 'CellData'):
        out[section] = {a.get('Name'): _decode(a) for a in piece.find(section)}
    out['FieldData'] = {a.get('Name'): _decode(a) for a in data.findall('FieldData/DataArray')}
    if piece.find('Coordinates') is not None:
        out['Coordinates'] = [_decode(a) for a in piece.find('Coordinates')]
    if piece.find('Points') is not None:
        out['Points'] = _decode(piece.find('Points/DataArray'))
    if piece.find('Lines') is not None:
        lines = {a.get('Name'): _decode(a) for a in piece.find('Lines')}
        assert np.array_equal(lines['offsets'], 2 * np.arange(1, len(lines['offsets']) + 1)), fn
        out['Lines'] = lines['connectivity'].reshape(-1, 2)
        assert int(piece.get('NumberOfLines')) == len(out['Lines']), fn
    return out


def _as_nxyz(values, extent):
    """ A VTK point array (x fastest) in the order of the HDF5 dump: component, x, y, z. """
    nx, ny, nz = (extent[1] - extent[0] + 1, extent[3] - extent[2] + 1, extent[5] - extent[4] + 1)
    return values.reshape(nz, ny, nx, -1).transpose(3, 2, 1, 0)


class Test_VTK_Dumps(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cwd = os.getcwd()   # Run() changes into the simulation path
        cls.dir = tempfile.mkdtemp()
        try:
            _run_cartesian(os.path.join(cls.dir, 'cart'))
            _run_cylindrical(os.path.join(cls.dir, 'cyl'))
        finally:
            os.chdir(cls.cwd)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dir)

    def path(self, *names):
        return os.path.join(self.dir, *names)

    def files(self):
        return sorted(glob.glob(self.path('*', '*.vt[rsp]')))

    def last_step(self, sim, name, ext):
        """ The VTK and HDF5 field dump of the last timestep both have. """
        with h5py.File(self.path(sim, name + '_h5.h5'), 'r') as f:
            td = f['FieldData/TD']
            steps = {int(k): td[k][()] for k in td}
            mesh = {k: f['Mesh'][k][()] for k in f['Mesh']}   # in metres
        vtk_steps = {int(os.path.basename(fn)[len(name) + 1:-len(ext) - 1]): fn
                     for fn in glob.glob(self.path(sim, name + '_*.' + ext))}
        ts = max(steps.keys() & vtk_steps.keys())
        self.assertGreater(np.abs(steps[ts]).max(), 0, 'the compared timestep has no field')
        return _read(vtk_steps[ts]), steps[ts], mesh

    def test_files(self):
        # the files stay what they have been: VTK XML with zlib compressed arrays
        self.assertEqual({os.path.splitext(f)[1] for f in self.files()}, {'.vtr', '.vts', '.vtp'})
        # the vector fields of the Cartesian dumps are compressed in more than one block of 32 KiB
        self.assertGreater(_read(self.path('cart', 'material_dump.vtr'))['PointData']['epsilon'].nbytes, 32768)
        for fn in self.files():
            with self.subTest(os.path.relpath(fn, self.dir)):
                _read(fn)

    @unittest.skipUnless(HAS_VTK, 'vtk not available')
    def test_vtk_reads_the_same(self):
        errors = vtk.vtkFileOutputWindow()   # VTK's complaints into a file, not the test log
        errors.SetFileName(self.path('vtk_errors.txt'))
        vtk.vtkOutputWindow.SetInstance(errors)
        for fn in self.files():
            with self.subTest(os.path.relpath(fn, self.dir)):
                reader = vtk.vtkXMLGenericDataObjectReader()
                reader.SetFileName(fn)
                reader.Update()
                self.assertEqual(reader.GetErrorCode(), 0)
                data, ours = reader.GetOutput(), _read(fn)
                for section in ('PointData', 'CellData', 'FieldData'):
                    arrays = getattr(data, 'Get' + section)()
                    self.assertEqual(arrays.GetNumberOfArrays(), len(ours[section]))
                    for name, values in ours[section].items():
                        np.testing.assert_array_equal(vtk_to_numpy(arrays.GetAbstractArray(name)), values)
                if 'Coordinates' in ours:
                    for c, values in zip((data.GetXCoordinates(), data.GetYCoordinates(), data.GetZCoordinates()),
                                         ours['Coordinates']):
                        np.testing.assert_array_equal(vtk_to_numpy(c), values)
                    self.assertEqual(list(data.GetExtent()), ours['extent'])
                if 'Points' in ours:
                    np.testing.assert_array_equal(vtk_to_numpy(data.GetPoints().GetData()).ravel(), ours['Points'].ravel())
                if 'Lines' in ours:
                    np.testing.assert_array_equal(vtk_to_numpy(data.GetLines().GetConnectivityArray()), ours['Lines'].ravel())

    def test_field_dump_cartesian(self):
        data, h5, mesh = self.last_step('cart', 'Et', 'vtr')
        self.assertEqual(data['type'], 'RectilinearGrid')
        for c, axis in zip(data['Coordinates'], 'xyz'):
            np.testing.assert_allclose(c, mesh[axis], rtol=1e-12)
        E = data['PointData']['E-Field']
        self.assertEqual(E.dtype, np.float32)
        np.testing.assert_array_equal(_as_nxyz(E, data['extent']), h5)

    def test_field_dump_cylindrical(self):
        # in Cartesian coordinates: points and field components
        data, h5, mesh = self.last_step('cyl', 'Et', 'vts')
        self.assertEqual(data['type'], 'StructuredGrid')
        r, a, z = np.meshgrid(mesh['rho'], mesh['alpha'], mesh['z'], indexing='ij')
        pts = _as_nxyz(data['Points'], data['extent'])
        expected = np.array([r * np.cos(a), r * np.sin(a), z])
        np.testing.assert_allclose(pts, expected, rtol=0, atol=1e-7 * np.abs(expected).max())
        E = _as_nxyz(data['PointData']['E-Field'], data['extent'])
        cos, sin = np.cos(a), np.sin(a)
        expected = np.array([h5[0] * cos - h5[1] * sin, h5[0] * sin + h5[1] * cos, h5[2]])
        np.testing.assert_allclose(E, expected, rtol=0, atol=1e-6 * np.abs(expected).max())

    def test_frequency_and_sar_dumps(self):
        for fn, name in (('Hf_f=2000000000_abs.vtr', 'H-Field'), ('SAR_f=2000000000.000000.vtr', 'SAR')):
            with self.subTest(fn):
                arrays = _read(self.path('cart', fn))['PointData']
                self.assertEqual(len(arrays), 1)
                (key, values), = arrays.items()
                self.assertTrue(key.startswith(name), key)
                self.assertGreater(np.abs(values).max(), 0)
        # the SAR, a scalar field, against its HDF5 dump
        data = _read(self.path('cart', 'SAR_f=2000000000.000000.vtr'))
        with h5py.File(self.path('cart', 'SAR_h5.h5'), 'r') as f:
            (sar,) = data['PointData'].values()
            np.testing.assert_array_equal(_as_nxyz(sar, data['extent'])[0], f['FieldData/FD/f0'][()])

    def test_material_dump(self):
        data = _read(self.path('cart', 'material_dump.vtr'))
        z = data['Coordinates'][2]
        eps = _as_nxyz(data['PointData']['epsilon'], data['extent'])
        kappa = _as_nxyz(data['PointData']['kappa'], data['extent'])
        inside, air = np.argmin(np.abs(z + 3.3e-3)), np.argmin(np.abs(z - 13e-3))
        np.testing.assert_allclose(eps[:, 3:-3, 3:-3, inside], 3.6, rtol=1e-6)
        np.testing.assert_allclose(kappa[:, 3:-3, 3:-3, inside], 0.05, rtol=1e-6)
        np.testing.assert_allclose(eps[:, 3:-3, 3:-3, air], 1, rtol=1e-6)
        for name in ('mue', 'sigma'):
            self.assertEqual(data['PointData'][name].shape[1], 3)
        data = _read(self.path('cyl', 'material_dump.vts'))
        self.assertAlmostEqual(float(data['PointData']['epsilon'].max()), 2.2, places=5)

    def test_operator_dump(self):
        for sim, ext in (('cart', 'vtr'), ('cyl', 'vts')):
            with self.subTest(sim):
                arrays = _read(self.path(sim, 'operator_dump.' + ext))['PointData']
                for name in ('vv', 'vi', 'iv', 'ii', 'exc_volt'):
                    self.assertEqual(arrays[name].shape[1], 3, name)
                self.assertGreater(np.abs(arrays['exc_volt']).max(), 0)

    def test_pec_dump(self):
        # the metal edges: in the plane of the plate (Cartesian) or the wall (cylindrical)
        data = _read(self.path('cart', 'PEC_dump.vtp'))
        self.assertEqual(data['type'], 'PolyData')
        pts, lines = data['Points'], data['Lines']
        self.assertGreater(len(lines), 50)
        np.testing.assert_allclose(pts[:, 2], 0, atol=1e-9)
        mid = (pts[lines[:, 0]] + pts[lines[:, 1]]) / 2
        self.assertTrue(np.all(np.abs(mid[:, 0]) <= 10e-3) and np.all(np.abs(mid[:, 1]) <= 8e-3))
        # each edge runs along one axis, from one mesh line to the next
        step = np.abs(pts[lines[:, 1]] - pts[lines[:, 0]])
        np.testing.assert_array_equal(np.count_nonzero(step > 1e-12, axis=1), 1)
        data = _read(self.path('cyl', 'PEC_dump.vtp'))
        self.assertGreater(len(data['Lines']), 0)
        np.testing.assert_allclose(data['Points'][:, 2], 32e-3, rtol=1e-6)


if __name__ == '__main__':
    unittest.main()
