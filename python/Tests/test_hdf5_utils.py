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

"""Tests for the openEMS HDF5 dump reader in openEMS.utilities."""

import os
import tempfile
import unittest
import numpy as np
import h5py

from openEMS.utilities import HDF5Dump

COMPLEX_T = np.dtype([('r', np.float32), ('i', np.float32)])

NX, NY, NZ = 5, 4, 3   # all different, so a wrong axis order is detectable


def _write_mesh(h5, mesh_type=0, scaling=1e-3, names=('x', 'y', 'z')):
    grp = h5.require_group('Mesh')
    grp.attrs['mesh_type'] = np.int32(mesh_type)
    grp.attrs['mesh_scaling'] = np.float64(scaling)
    for name, n in zip(names, (NX, NY, NZ)):
        grp[name] = np.arange(n, dtype=np.float64) * scaling


def _ref_vector_field():
    """Reference vector field in the natural (3, NX, NY, NZ) order."""
    n = 3 * NX * NY * NZ
    re = np.arange(n, dtype=np.float32).reshape(3, NX, NY, NZ)
    return re + 1j * (re + 0.5)


def _ref_scalar_field():
    return np.arange(NX * NY * NZ, dtype=np.float32).reshape(NX, NY, NZ)


def _write_fd_vector(path, legacy=False, version=0.3, d_order=True):
    """Write a frequency domain vector field dump in the given format."""
    data = _ref_vector_field()
    with h5py.File(path, 'w') as h5:
        h5.attrs['openEMS_HDF5_version'] = np.float64(version)
        h5.attrs['dump_type'] = np.int32(10)
        if legacy:
            h5.attrs['legacy_fmt'] = True
        _write_mesh(h5)
        grp = h5.require_group('FieldData/FD')
        grp.attrs['frequency'] = np.array([1e9], dtype=np.float64)
        if legacy:
            swapped = np.swapaxes(data, 1, 3)   # (3, NZ, NY, NX)
            for part, val in (('real', swapped.real), ('imag', swapped.imag)):
                ds = grp.create_dataset('f0_' + part, data=val.astype(np.float32))
                ds.attrs['frequency'] = np.float64(1e9)
                if d_order:
                    ds.attrs['d_order'] = np.bytes_(b'NZYX')
        else:
            buf = np.empty(data.shape, dtype=COMPLEX_T)
            buf['r'] = data.real
            buf['i'] = data.imag
            ds = grp.create_dataset('f0', data=buf)
            ds.attrs['frequency'] = np.float64(1e9)
            if d_order:
                ds.attrs['d_order'] = np.bytes_(b'NXYZ')
    return data


class _TempFile(unittest.TestCase):
    def setUp(self):
        fd, self.path = tempfile.mkstemp(suffix='.h5')
        os.close(fd)

    def tearDown(self):
        os.unlink(self.path)


class Test_Mesh(_TempFile):
    def test_cartesian(self):
        _write_fd_vector(self.path)
        with HDF5Dump(self.path) as dump:
            mesh = dump.GetMesh()
        self.assertEqual(mesh['type'], 0)
        self.assertEqual(mesh['names'], ['x', 'y', 'z'])
        self.assertEqual([len(l) for l in mesh['lines']], [NX, NY, NZ])
        self.assertAlmostEqual(mesh['scaling'], 1e-3)

    def test_lines_are_si_units(self):
        """Lines are returned as stored, i.e. already scaled to metres."""
        _write_fd_vector(self.path)
        with HDF5Dump(self.path) as dump:
            mesh = dump.GetMesh()
        np.testing.assert_allclose(mesh['lines'][0] / mesh['scaling'],
                                   np.arange(NX))

    def test_cylindrical(self):
        with h5py.File(self.path, 'w') as h5:
            h5.attrs['openEMS_HDF5_version'] = np.float64(0.3)
            _write_mesh(h5, mesh_type=1, names=('rho', 'alpha', 'z'))
            ds = h5.create_dataset('FieldData/FD/f0', data=_ref_scalar_field())
            ds.attrs['d_order'] = np.bytes_(b'XYZ')
        with HDF5Dump(self.path) as dump:
            mesh = dump.GetMesh()
        self.assertEqual(mesh['type'], 1)
        self.assertEqual(mesh['names'], ['rho', 'alpha', 'z'])

    def test_missing_mesh_type_attribute(self):
        """Without a mesh_type attribute the type is derived from the datasets."""
        with h5py.File(self.path, 'w') as h5:
            h5.attrs['openEMS_HDF5_version'] = np.float64(0.3)
            grp = h5.require_group('Mesh')
            for name, n in zip(('rho', 'alpha', 'z'), (NX, NY, NZ)):
                grp[name] = np.arange(n, dtype=np.float64)
            ds = h5.create_dataset('FieldData/FD/f0', data=_ref_scalar_field())
            ds.attrs['d_order'] = np.bytes_(b'XYZ')
        with HDF5Dump(self.path) as dump:
            mesh = dump.GetMesh()
        self.assertEqual(mesh['type'], 1)
        self.assertAlmostEqual(mesh['scaling'], 1.0)

    def test_no_mesh_group_raises(self):
        with h5py.File(self.path, 'w') as h5:
            h5['dummy'] = [1, 2, 3]
        with self.assertRaises(KeyError):
            HDF5Dump(self.path)

    def test_no_field_data_raises(self):
        with h5py.File(self.path, 'w') as h5:
            _write_mesh(h5)
        with self.assertRaises(KeyError):
            HDF5Dump(self.path)


class Test_Metadata(_TempFile):
    def test_fd_metadata(self):
        _write_fd_vector(self.path)
        with HDF5Dump(self.path) as dump:
            self.assertTrue(dump.IsFD)
            self.assertFalse(dump.IsTD)
            self.assertTrue(dump.IsVector)
            self.assertEqual(dump.Shape, (NX, NY, NZ))
            self.assertEqual(dump.NumFrequencies, 1)
            self.assertEqual(dump.NumTimesteps, 0)
            self.assertEqual(dump.DumpType, 10)
            self.assertEqual(dump.DumpTypeName, 'E-field (FD)')
            np.testing.assert_allclose(dump.Frequencies, [1e9])

    def test_legacy_shape_is_reported_in_logical_order(self):
        _write_fd_vector(self.path, legacy=True)
        with HDF5Dump(self.path) as dump:
            self.assertEqual(dump.Shape, (NX, NY, NZ))

    def test_repr_does_not_raise(self):
        _write_fd_vector(self.path)
        with HDF5Dump(self.path) as dump:
            self.assertIn('E-field (FD)', repr(dump))

    def test_file_property_exposes_handle(self):
        _write_fd_vector(self.path)
        with HDF5Dump(self.path) as dump:
            self.assertIn('Mesh', dump.File)

    def test_open_file_handle_is_not_closed(self):
        _write_fd_vector(self.path)
        with h5py.File(self.path, 'r') as h5:
            with HDF5Dump(h5) as dump:
                dump.GetFieldAtIndex(f_idx=0)
            self.assertTrue(bool(h5))       # still usable after the class closed


class Test_FD(_TempFile):
    def test_complex_vector_field(self):
        ref = _write_fd_vector(self.path)
        with HDF5Dump(self.path) as dump:
            data = dump.GetFieldAtIndex(f_idx=0)
        self.assertEqual(data.shape, (3, NX, NY, NZ))
        self.assertTrue(np.iscomplexobj(data))
        np.testing.assert_allclose(data, ref)

    def test_complex_dtype_is_not_upcast(self):
        """float32 on disk stays complex64, in both storage formats."""
        for legacy in (False, True):
            _write_fd_vector(self.path, legacy=legacy)
            with HDF5Dump(self.path) as dump:
                self.assertEqual(dump.GetFieldAtIndex(f_idx=0).dtype,
                                 np.complex64, 'legacy={}'.format(legacy))

    def test_attributes_merged(self):
        """Attributes come from the file root, the FD group and the dataset."""
        _write_fd_vector(self.path)
        with HDF5Dump(self.path) as dump:
            attrs = dump.GetAttributes(f_idx=0)
        self.assertEqual(int(attrs['dump_type']), 10)             # root
        self.assertAlmostEqual(float(attrs['openEMS_HDF5_version']), 0.3)
        self.assertAlmostEqual(float(attrs['frequency']), 1e9)    # dataset wins

    def test_missing_index_raises(self):
        _write_fd_vector(self.path)
        with HDF5Dump(self.path) as dump:
            with self.assertRaises(IndexError):
                dump.GetFieldAtIndex(f_idx=1)

    def test_get_at_frequency(self):
        ref = _write_fd_vector(self.path)
        with HDF5Dump(self.path) as dump:
            np.testing.assert_allclose(dump.GetFieldAtFrequency(1e9), ref)
            # a tiny relative offset still matches
            np.testing.assert_allclose(dump.GetFieldAtFrequency(1e9*(1+1e-9)), ref)

    def test_unknown_frequency_raises(self):
        _write_fd_vector(self.path)
        with HDF5Dump(self.path) as dump:
            with self.assertRaises(ValueError):
                dump.GetFieldAtFrequency(2e9)

    def test_iter_fd(self):
        ref = _write_fd_vector(self.path)
        with HDF5Dump(self.path) as dump:
            got = list(dump.IterFD())
        self.assertEqual(len(got), 1)
        self.assertAlmostEqual(got[0][0], 1e9)
        np.testing.assert_allclose(got[0][1], ref)

    def test_scalar_field(self):
        """A real valued FD dump (e.g. local SAR) is read as (NX, NY, NZ)."""
        ref = _ref_scalar_field()
        with h5py.File(self.path, 'w') as h5:
            h5.attrs['openEMS_HDF5_version'] = np.float64(0.3)
            h5.attrs['dump_type'] = np.int32(20)
            _write_mesh(h5)
            ds = h5.create_dataset('FieldData/FD/f0', data=ref)
            ds.attrs['d_order'] = np.bytes_(b'XYZ')
        with HDF5Dump(self.path) as dump:
            self.assertFalse(dump.IsVector)
            self.assertEqual(dump.DumpTypeName, 'local SAR')
            data = dump.GetFieldAtIndex(f_idx=0)
        self.assertEqual(data.shape, (NX, NY, NZ))
        np.testing.assert_allclose(data, ref)


class Test_Legacy(_TempFile):
    """Legacy dumps are transposed back into the natural x-inner order."""

    def test_legacy_by_d_order(self):
        ref = _write_fd_vector(self.path, legacy=True, version=0.3)
        with HDF5Dump(self.path) as dump:
            data = dump.GetFieldAtIndex(f_idx=0)
        self.assertEqual(data.shape, (3, NX, NY, NZ))
        np.testing.assert_allclose(data, ref)

    def test_legacy_by_version(self):
        """Old files carry no d_order attribute; the version decides."""
        ref = _write_fd_vector(self.path, legacy=True, version=0.2, d_order=False)
        with HDF5Dump(self.path) as dump:
            data = dump.GetFieldAtIndex(f_idx=0)
        self.assertEqual(data.shape, (3, NX, NY, NZ))
        np.testing.assert_allclose(data, ref)

    def test_legacy_scalar_field(self):
        ref = _ref_scalar_field()
        with h5py.File(self.path, 'w') as h5:
            h5.attrs['openEMS_HDF5_version'] = np.float64(0.3)
            h5.attrs['legacy_fmt'] = True
            _write_mesh(h5)
            ds = h5.create_dataset('FieldData/FD/f0', data=np.swapaxes(ref, 0, 2))
            ds.attrs['d_order'] = np.bytes_(b'ZYX')
        with HDF5Dump(self.path) as dump:
            data = dump.GetFieldAtIndex(f_idx=0)
        self.assertEqual(data.shape, (NX, NY, NZ))
        np.testing.assert_allclose(data, ref)

    def test_legacy_region_matches_new_format(self):
        """A region selection gives identical results in both storage formats."""
        out = []
        for legacy in (False, True):
            ref = _write_fd_vector(self.path, legacy=legacy)
            with HDF5Dump(self.path) as dump:
                dump.SetPlane('z', idx=1)
                dump.SetRange('x', idx_start=1, idx_stop=4)
                dump.SetSampling(2, 1, 1)
                out.append(dump.GetFieldAtIndex(f_idx=0, component='y'))
            np.testing.assert_allclose(out[-1], ref[1, 1:4:2, :, 1])
        np.testing.assert_allclose(out[0], out[1])


class _TDFile(_TempFile):
    TIMES = (1e-12, 2e-12, 3e-12)

    def setUp(self):
        super().setUp()
        self.ref = []
        with h5py.File(self.path, 'w') as h5:
            h5.attrs['openEMS_HDF5_version'] = np.float64(0.3)
            h5.attrs['dump_type'] = np.int32(0)
            _write_mesh(h5)
            grp = h5.require_group('FieldData/TD')
            for n, time in enumerate(self.TIMES):
                val = _ref_vector_field().real + n
                self.ref.append(val)
                ds = grp.create_dataset('{:06d}'.format(100 * (n + 1)), data=val)
                ds.attrs['time'] = np.float64(time)
                ds.attrs['d_order'] = np.bytes_(b'NXYZ')


class Test_TD(_TDFile):
    def test_metadata(self):
        with HDF5Dump(self.path) as dump:
            self.assertTrue(dump.IsTD)
            self.assertFalse(dump.IsFD)
            self.assertEqual(dump.NumTimesteps, 3)
            self.assertEqual(dump.DumpTypeName, 'E-field (TD)')
            np.testing.assert_allclose(dump.Times, self.TIMES)

    def test_index_is_the_position_in_the_file(self):
        with HDF5Dump(self.path) as dump:
            data = dump.GetFieldAtIndex(t_idx=2)
            attrs = dump.GetAttributes(t_idx=2)
        self.assertAlmostEqual(float(attrs['time']), self.TIMES[2])
        np.testing.assert_allclose(data, self.ref[2])

    def test_index_by_dataset_name(self):
        with HDF5Dump(self.path) as dump:
            data = dump.GetFieldAtIndex(t_idx='000200')
            attrs = dump.GetAttributes(t_idx='000200')
        self.assertAlmostEqual(float(attrs['time']), self.TIMES[1])
        np.testing.assert_allclose(data, self.ref[1])

    def test_requesting_fd_raises(self):
        with HDF5Dump(self.path) as dump:
            with self.assertRaises(KeyError):
                dump.GetFieldAtIndex(f_idx=0)

    def test_iter_td(self):
        with HDF5Dump(self.path) as dump:
            got = list(dump.IterTD())
        self.assertEqual(len(got), 3)
        for n, (time, data) in enumerate(got):
            self.assertAlmostEqual(time, self.TIMES[n])
            np.testing.assert_allclose(data, self.ref[n])

    def test_on_the_fly_dft(self):
        """Requesting a frequency from a TD dump transforms it on the fly."""
        freq = 0.25e12
        dt = self.TIMES[1] - self.TIMES[0]
        want = 2 * dt * sum(self.ref[n] * np.exp(-2j*np.pi*freq*self.TIMES[n])
                            for n in range(len(self.TIMES)))
        with HDF5Dump(self.path) as dump:
            got = dump.GetFieldAtFrequency(freq)
        np.testing.assert_allclose(got, want, rtol=1e-5)

    def test_dft_honours_the_region(self):
        freq = 0.25e12
        with HDF5Dump(self.path) as dump:
            full = dump.GetFieldAtFrequency(freq)
            dump.SetPlane('z', idx=2)
            plane = dump.GetFieldAtFrequency(freq)
        self.assertEqual(plane.shape, (3, NX, NY))
        np.testing.assert_allclose(plane, full[:, :, :, 2], rtol=1e-5)


class Test_Region(_TempFile):
    def setUp(self):
        super().setUp()
        self.ref = _write_fd_vector(self.path)
        self.dump = HDF5Dump(self.path)

    def tearDown(self):
        self.dump.Close()
        super().tearDown()

    def test_plane_by_index(self):
        self.dump.SetPlane('z', idx=1)
        data = self.dump.GetFieldAtIndex(f_idx=0)
        self.assertEqual(data.shape, (3, NX, NY))
        np.testing.assert_allclose(data, self.ref[:, :, :, 1])

    def test_plane_by_position(self):
        lines = self.dump.GetMesh()['lines'][2]
        self.dump.SetPlane('z', pos=lines[2] + 1e-9)
        np.testing.assert_allclose(self.dump.GetFieldAtIndex(f_idx=0),
                                   self.ref[:, :, :, 2])

    def test_plane_by_axis_number(self):
        self.dump.SetPlane(1, idx=2)
        np.testing.assert_allclose(self.dump.GetFieldAtIndex(f_idx=0),
                                   self.ref[:, :, 2, :])

    def test_range_by_index_is_half_open(self):
        self.dump.SetRange('x', idx_start=1, idx_stop=4)
        np.testing.assert_allclose(self.dump.GetFieldAtIndex(f_idx=0),
                                   self.ref[:, 1:4])

    def test_range_by_position_is_inclusive(self):
        lines = self.dump.GetMesh()['lines'][0]
        self.dump.SetRange('x', start=lines[1], stop=lines[3])
        np.testing.assert_allclose(self.dump.GetFieldAtIndex(f_idx=0),
                                   self.ref[:, 1:4])

    def test_sampling(self):
        self.dump.SetSampling(2, 2, 1)
        np.testing.assert_allclose(self.dump.GetFieldAtIndex(f_idx=0),
                                   self.ref[:, ::2, ::2, :])

    def test_single_sampling_factor_applies_to_all(self):
        self.dump.SetSampling(2)
        np.testing.assert_allclose(self.dump.GetFieldAtIndex(f_idx=0),
                                   self.ref[:, ::2, ::2, ::2])

    def test_component_selection(self):
        data = self.dump.GetFieldAtIndex(f_idx=0, component='z')
        self.assertEqual(data.shape, (NX, NY, NZ))
        np.testing.assert_allclose(data, self.ref[2])

    def test_combined(self):
        self.dump.SetPlane('z', idx=0)
        self.dump.SetRange('x', idx_start=1, idx_stop=5)
        self.dump.SetSampling(2, 1, 1)
        np.testing.assert_allclose(
            self.dump.GetFieldAtIndex(f_idx=0, component=0),
            self.ref[0, 1:5:2, :, 0])

    def test_region_mesh_matches_the_data(self):
        self.dump.SetPlane('z', idx=1)
        self.dump.SetSampling(2, 1, 1)
        data = self.dump.GetFieldAtIndex(f_idx=0)
        mesh = self.dump.GetMesh(region=True)
        self.assertEqual(len(mesh['lines'][0]), data.shape[1])
        self.assertEqual(len(mesh['lines'][1]), data.shape[2])
        self.assertEqual(len(mesh['lines'][2]), 1)

    def test_plane_then_range_same_direction_raises(self):
        """A plane and a range on one direction contradict each other."""
        self.dump.SetPlane('z', idx=1)
        with self.assertRaises(ValueError):
            self.dump.SetRange('z', idx_start=0, idx_stop=2)

    def test_range_then_plane_same_direction_raises(self):
        self.dump.SetRange('x', idx_start=0, idx_stop=2)
        with self.assertRaises(ValueError):
            self.dump.SetPlane('x', idx=1)

    def test_conflict_is_per_direction(self):
        """A plane on one direction and a range on another is fine."""
        self.dump.SetPlane('z', idx=1)
        self.dump.SetRange('x', idx_start=1, idx_stop=4)
        np.testing.assert_allclose(self.dump.GetFieldAtIndex(f_idx=0),
                                   self.ref[:, 1:4, :, 1])

    def test_replacing_a_setting_of_the_same_kind_is_allowed(self):
        self.dump.SetPlane('z', idx=1)
        self.dump.SetPlane('z', idx=2)          # move the plane
        np.testing.assert_allclose(self.dump.GetFieldAtIndex(f_idx=0),
                                   self.ref[:, :, :, 2])
        self.dump.SetRange('x', idx_start=0, idx_stop=2)
        self.dump.SetRange('x', idx_start=1, idx_stop=4)   # redefine the range
        np.testing.assert_allclose(self.dump.GetFieldAtIndex(f_idx=0),
                                   self.ref[:, 1:4, :, 2])

    def test_plane_on_another_direction_replaces_the_first(self):
        """There is only one plane; switching orientation needs no reset."""
        self.dump.SetPlane('z', idx=1)
        self.dump.SetPlane('y', idx=2)
        data = self.dump.GetFieldAtIndex(f_idx=0)
        self.assertEqual(data.shape, (3, NX, NZ))       # an xz-plane, not a line
        np.testing.assert_allclose(data, self.ref[:, :, 2, :])

    def test_plane_replacement_keeps_ranges_on_other_directions(self):
        self.dump.SetRange('x', idx_start=1, idx_stop=4)
        self.dump.SetPlane('z', idx=1)
        self.dump.SetPlane('y', idx=2)                  # replaces the z-plane
        np.testing.assert_allclose(self.dump.GetFieldAtIndex(f_idx=0),
                                   self.ref[:, 1:4, 2, :])

    def test_plane_replacement_keeps_sampling(self):
        self.dump.SetSampling(2, 1, 1)
        self.dump.SetPlane('z', idx=1)
        self.dump.SetPlane('y', idx=2)
        np.testing.assert_allclose(self.dump.GetFieldAtIndex(f_idx=0),
                                   self.ref[:, ::2, 2, :])

    def test_reset_single_direction_clears_the_conflict(self):
        self.dump.SetPlane('z', idx=1)
        self.dump.ResetRegion('z')
        self.dump.SetRange('z', idx_start=0, idx_stop=2)    # no longer a conflict
        np.testing.assert_allclose(self.dump.GetFieldAtIndex(f_idx=0),
                                   self.ref[:, :, :, 0:2])

    def test_reset_single_direction_keeps_the_others(self):
        self.dump.SetPlane('z', idx=1)
        self.dump.SetRange('x', idx_start=1, idx_stop=4)
        self.dump.ResetRegion('x')
        np.testing.assert_allclose(self.dump.GetFieldAtIndex(f_idx=0),
                                   self.ref[:, :, :, 1])

    def test_line_by_index(self):
        """A line along x collapses y and z."""
        self.dump.SetLine('x', idx=(None, 2, 1))
        data = self.dump.GetFieldAtIndex(f_idx=0)
        self.assertEqual(data.shape, (3, NX))
        np.testing.assert_allclose(data, self.ref[:, :, 2, 1])

    def test_line_by_position(self):
        lines = self.dump.GetMesh()['lines']
        self.dump.SetLine('z', pos=(lines[0][1], lines[1][2], None))
        data = self.dump.GetFieldAtIndex(f_idx=0)
        self.assertEqual(data.shape, (3, NZ))
        np.testing.assert_allclose(data, self.ref[:, 1, 2, :])

    def test_line_entry_order_is_by_axis_position(self):
        """The entry for direction n is always at position n."""
        self.dump.SetLine('y', idx=(3, None, 1))    # x=3, z=1
        np.testing.assert_allclose(self.dump.GetFieldAtIndex(f_idx=0),
                                   self.ref[:, 3, :, 1])

    def test_line_overrides_perpendicular_ranges(self):
        self.dump.SetRange('y', idx_start=0, idx_stop=3)
        self.dump.SetLine('x', idx=(None, 2, 1))
        np.testing.assert_allclose(self.dump.GetFieldAtIndex(f_idx=0),
                                   self.ref[:, :, 2, 1])

    def test_line_resets_a_plane_on_the_line_direction(self):
        """A former plane normal cannot stay collapsed for the line axis."""
        self.dump.SetPlane('x', idx=4)
        self.dump.SetLine('x', idx=(None, 2, 1))
        np.testing.assert_allclose(self.dump.GetFieldAtIndex(f_idx=0),
                                   self.ref[:, :, 2, 1])

    def test_line_replaces_a_previous_line(self):
        self.dump.SetLine('x', idx=(None, 2, 1))
        self.dump.SetLine('z', idx=(1, 2, None))
        data = self.dump.GetFieldAtIndex(f_idx=0)
        self.assertEqual(data.shape, (3, NZ))
        np.testing.assert_allclose(data, self.ref[:, 1, 2, :])

    def test_line_keeps_a_range_on_the_line_direction(self):
        """A range of two or more lines along the line direction survives."""
        self.dump.SetRange('x', idx_start=1, idx_stop=4)
        self.dump.SetLine('x', idx=(None, 2, 1))
        np.testing.assert_allclose(self.dump.GetFieldAtIndex(f_idx=0),
                                   self.ref[:, 1:4, 2, 1])

    def test_line_resets_a_single_valued_range_on_the_line_direction(self):
        self.dump.SetRange('x', idx_start=2, idx_stop=3)   # only one line
        self.dump.SetLine('x', idx=(None, 2, 1))
        np.testing.assert_allclose(self.dump.GetFieldAtIndex(f_idx=0),
                                   self.ref[:, :, 2, 1])

    def test_line_keeps_sampling(self):
        self.dump.SetSampling(2, 1, 1)
        self.dump.SetLine('x', idx=(None, 2, 1))
        np.testing.assert_allclose(self.dump.GetFieldAtIndex(f_idx=0),
                                   self.ref[:, ::2, 2, 1])

    def test_plane_after_line_gives_a_plane(self):
        self.dump.SetLine('x', idx=(None, 2, 1))
        self.dump.SetPlane('z', idx=1)
        data = self.dump.GetFieldAtIndex(f_idx=0)
        self.assertEqual(data.shape, (3, NX, NY))
        np.testing.assert_allclose(data, self.ref[:, :, :, 1])

    def test_line_region_mesh(self):
        self.dump.SetLine('x', idx=(None, 2, 1))
        mesh = self.dump.GetMesh(region=True)
        self.assertEqual([len(l) for l in mesh['lines']], [NX, 1, 1])

    def test_line_needs_exactly_one_of_pos_or_idx(self):
        with self.assertRaises(ValueError):
            self.dump.SetLine('x')
        with self.assertRaises(ValueError):
            self.dump.SetLine('x', pos=(None, 0.0, 0.0), idx=(None, 0, 0))

    def test_line_direction_entry_must_be_none(self):
        with self.assertRaises(ValueError):
            self.dump.SetLine('x', idx=(0, 2, 1))

    def test_line_perpendicular_entry_must_not_be_none(self):
        with self.assertRaises(ValueError):
            self.dump.SetLine('x', idx=(None, None, 1))

    def test_line_needs_three_entries(self):
        with self.assertRaises(ValueError):
            self.dump.SetLine('x', idx=(2, 1))

    def test_line_index_out_of_range_raises(self):
        with self.assertRaises(IndexError):
            self.dump.SetLine('x', idx=(None, NY, 1))

    def test_reset_region(self):
        self.dump.SetPlane('z', idx=1)
        self.dump.ResetRegion()
        np.testing.assert_allclose(self.dump.GetFieldAtIndex(f_idx=0), self.ref)

    def test_nearest_index(self):
        lines = self.dump.GetMesh()['lines'][1]
        self.assertEqual(self.dump.NearestIndex('y', lines[2] + 1e-9), 2)
        self.assertEqual(self.dump.NearestIndex(1, lines[2]), 2)

    def test_invalid_direction_raises(self):
        with self.assertRaises(ValueError):
            self.dump.SetPlane('q', idx=0)

    def test_plane_index_out_of_range_raises(self):
        with self.assertRaises(IndexError):
            self.dump.SetPlane('z', idx=NZ)

    def test_plane_needs_exactly_one_of_pos_or_idx(self):
        with self.assertRaises(ValueError):
            self.dump.SetPlane('z')
        with self.assertRaises(ValueError):
            self.dump.SetPlane('z', pos=0.0, idx=0)

    def test_range_cannot_mix_coord_and_index(self):
        with self.assertRaises(ValueError):
            self.dump.SetRange('x', start=0.0, idx_start=0)

    def test_sampling_must_be_positive(self):
        with self.assertRaises(ValueError):
            self.dump.SetSampling(0)

    def test_index_arguments_are_exclusive(self):
        with self.assertRaises(ValueError):
            self.dump.GetFieldAtIndex()
        with self.assertRaises(ValueError):
            self.dump.GetFieldAtIndex(f_idx=0, t_idx=0)


if __name__ == '__main__':
    unittest.main()
