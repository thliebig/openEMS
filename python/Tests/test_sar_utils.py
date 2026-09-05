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

import os
import tempfile
import unittest
import numpy as np
import h5py

from openEMS.sar_utils import readSAR


def _make_sar_result_hdf5(path, nx=4, ny=3, nz=2, version=0.3):
    """Create a minimal valid SAR result HDF5 file at *path*.

    Parameters
    ----------
    path : str
        Destination file path (will be overwritten).
    nx, ny, nz : int
        Grid dimensions along x, y, z.
    version : float
        Value written to the ``openEMS_HDF5_version`` root attribute.
        Use 0.1 or 0.2 to exercise the legacy swapaxes path.
    """
    with h5py.File(path, 'w') as h5:
        h5.attrs['openEMS_HDF5_version'] = np.float64(version)
        h5.attrs['mass'] = np.float64(10.0)

        h5['Mesh/x'] = np.linspace(0.0, 1.0, nx)
        h5['Mesh/y'] = np.linspace(0.0, 1.0, ny)
        h5['Mesh/z'] = np.linspace(0.0, 1.0, nz)

        data = np.arange(1, nx * ny * nz + 1, dtype=np.float32).reshape(nx, ny, nz)
        ds = h5.create_dataset('FieldData/FD/f0', data=data)
        ds.attrs['frequency'] = np.float64(1e9)
        ds.attrs['power'] = np.float64(1.0)


class Test_readSAR_MissingFile(unittest.TestCase):
    def test_nonexistent_file_raises(self):
        """readSAR on a non-existent file raises an OSError / FileNotFoundError."""
        with self.assertRaises(OSError):
            readSAR('/nonexistent/path/sar_result.h5')


class Test_readSAR_NotADump(unittest.TestCase):
    def setUp(self):
        fd, self.path = tempfile.mkstemp(suffix='.h5')
        os.close(fd)
        with h5py.File(self.path, 'w') as h5:
            h5.create_dataset('dummy', data=[1, 2, 3])

    def tearDown(self):
        os.unlink(self.path)

    def test_raises(self):
        """A file that is not an openEMS dump raises, rather than returning None."""
        with self.assertRaises(KeyError):
            readSAR(self.path)


class Test_readSAR_ValidFile(unittest.TestCase):
    NX, NY, NZ = 4, 3, 2

    def setUp(self):
        fd, self.path = tempfile.mkstemp(suffix='.h5')
        os.close(fd)
        _make_sar_result_hdf5(self.path, nx=self.NX, ny=self.NY, nz=self.NZ)

    def tearDown(self):
        os.unlink(self.path)

    def test_shapes(self):
        sar, mesh, data = readSAR(self.path)
        self.assertEqual(sar.shape, (self.NX, self.NY, self.NZ))
        self.assertEqual([len(mesh[n]) for n in range(3)], [self.NX, self.NY, self.NZ])

    def test_sar_values_float_nonzero(self):
        sar, _, _ = readSAR(self.path)
        self.assertTrue(np.issubdtype(sar.dtype, np.floating))
        self.assertTrue(np.any(sar != 0))

    def test_metadata(self):
        _, _, data = readSAR(self.path)
        self.assertAlmostEqual(float(data['mass']), 10.0)
        self.assertAlmostEqual(float(data['frequency']), 1e9)
        self.assertAlmostEqual(float(data['power']), 1.0)

    def test_default_f_idx_zero(self):
        sar0_default, _, _ = readSAR(self.path)
        sar0_explicit, _, _ = readSAR(self.path, f_idx=0)
        np.testing.assert_array_equal(sar0_default, sar0_explicit)


class Test_readSAR_LegacyVersion(unittest.TestCase):
    """Version <= 0.2 triggers swapaxes(0, 2) on the returned array."""
    NX, NY, NZ = 4, 3, 2  # nx != nz so shape change is detectable

    def setUp(self):
        fd, self.path = tempfile.mkstemp(suffix='.h5')
        os.close(fd)
        _make_sar_result_hdf5(self.path, nx=self.NX, ny=self.NY, nz=self.NZ,
                               version=0.1)

    def tearDown(self):
        os.unlink(self.path)

    def test_shape_is_swapped(self):
        """After swapaxes(0, 2) a (NX, NY, NZ) array becomes (NZ, NY, NX)."""
        sar, mesh, data = readSAR(self.path)
        self.assertEqual(sar.shape[0], self.NZ)
        self.assertEqual(sar.shape[2], self.NX)

    def test_middle_axis_unchanged(self):
        sar, mesh, data = readSAR(self.path)
        self.assertEqual(sar.shape[1], self.NY)

    def test_version_02_also_swapped(self):
        fd, path = tempfile.mkstemp(suffix='.h5')
        os.close(fd)
        try:
            _make_sar_result_hdf5(path, nx=self.NX, ny=self.NY, nz=self.NZ,
                                   version=0.2)
            sar, mesh, data = readSAR(path)
            self.assertEqual(sar.shape[0], self.NZ)
            self.assertEqual(sar.shape[2], self.NX)
        finally:
            os.unlink(path)


class Test_readSAR_OutOfRangeFIdx(unittest.TestCase):
    """Requesting an f_idx that does not exist in the file raises IndexError."""

    def setUp(self):
        fd, self.path = tempfile.mkstemp(suffix='.h5')
        os.close(fd)
        _make_sar_result_hdf5(self.path)

    def tearDown(self):
        os.unlink(self.path)

    def test_f_idx_1_raises_index_error(self):
        """f_idx=1 on a file that only has f0 is out of range."""
        with self.assertRaises(IndexError):
            readSAR(self.path, f_idx=1)

    def test_f_idx_999_raises_index_error(self):
        with self.assertRaises(IndexError):
            readSAR(self.path, f_idx=999)


if __name__ == '__main__':
    unittest.main()
