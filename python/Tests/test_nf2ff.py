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

"""Tests for the nf2ff result reader in openEMS.nf2ff."""

import os
import tempfile
import unittest
import numpy as np
import h5py

from openEMS.nf2ff import nf2ff_results

N_THETA, N_PHI, N_FREQ = 5, 3, 2   # all different, so a wrong axis order shows

THETA = np.linspace(0, np.pi, N_THETA)
PHI   = np.linspace(0, 2*np.pi, N_PHI)
FREQ  = np.array([1e9, 2e9])


def _ref_field(n):
    """Reference far field of frequency *n* in the natural (theta, phi) order."""
    re = np.arange(N_THETA*N_PHI, dtype=np.float64).reshape(N_THETA, N_PHI)
    return (re + n) + 1j*(re + n + 0.5)


def _ref_p_rad(n):
    return np.arange(N_THETA*N_PHI, dtype=np.float64).reshape(N_THETA, N_PHI) + 10*n


def _write_result(path, legacy=False):
    """Write an nf2ff result file in the current or in the legacy format."""
    with h5py.File(path, 'w') as h5:
        h5.attrs['openEMS_HDF5_version'] = np.float64(0.3)
        if legacy:
            h5.attrs['legacy_fmt'] = True

        mesh = h5.create_group('Mesh')
        mesh['theta'] = THETA.astype(np.float32)
        mesh['phi']   = PHI.astype(np.float32)
        mesh['r']     = np.array([1.0], dtype=np.float32)
        mesh.attrs['MeshType'] = np.float32(2)

        nf2ff = h5.create_group('nf2ff')
        nf2ff.attrs['Frequency'] = FREQ
        nf2ff.attrs['Prad'] = np.ones(N_FREQ)
        nf2ff.attrs['Dmax'] = np.ones(N_FREQ)

        for name in ('E_theta', 'E_phi'):
            grp = h5.create_group('/nf2ff/{}/FD'.format(name))
            for n in range(N_FREQ):
                data = _ref_field(n)
                if name == 'E_phi':
                    data = data.conj()
                if legacy:
                    grp['f{}_real'.format(n)] = np.real(data).T
                    grp['f{}_imag'.format(n)] = np.imag(data).T
                else:
                    grp['f{}'.format(n)] = data

        grp = h5.create_group('/nf2ff/P_rad/FD')
        for n in range(N_FREQ):
            data = _ref_p_rad(n)
            grp['f{}'.format(n)] = data.T if legacy else data


class TestNF2FFResults(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.path = os.path.join(self.tmp, 'nf2ff.h5')

    def tearDown(self):
        for fn in os.listdir(self.tmp):
            os.remove(os.path.join(self.tmp, fn))
        os.rmdir(self.tmp)

    def test_mesh_and_attributes(self):
        _write_result(self.path)
        res = nf2ff_results(self.path)
        np.testing.assert_allclose(res.theta, THETA, rtol=1e-6)
        np.testing.assert_allclose(res.phi, PHI, rtol=1e-6)
        np.testing.assert_allclose(res.freq, FREQ)

    def test_both_formats_read_identically(self):
        """The legacy and the current format must decode to the same arrays."""
        for legacy in (False, True):
            # a file of its own per format, so that neither run depends on
            # the other one being closed again
            path = os.path.join(self.tmp, 'nf2ff_legacy_{}.h5'.format(legacy))
            _write_result(path, legacy=legacy)
            res = nf2ff_results(path)
            for n in range(N_FREQ):
                msg = 'legacy={}, f{}'.format(legacy, n)
                self.assertEqual(res.E_theta[n].shape, (N_THETA, N_PHI), msg)
                np.testing.assert_allclose(res.E_theta[n], _ref_field(n),
                                           err_msg=msg)
                np.testing.assert_allclose(res.E_phi[n], _ref_field(n).conj(),
                                           err_msg=msg)
                np.testing.assert_allclose(res.P_rad[n], _ref_p_rad(n),
                                           err_msg=msg)

    def test_reading_does_not_keep_the_file_open(self):
        """A result file must be writable again right after it was read.

        Older HDF5 refuses to create a file that the same process still holds
        open, which is what a leaked reader handle would cause.
        """
        _write_result(self.path)
        res = nf2ff_results(self.path)
        _write_result(self.path, legacy=True)
        self.assertEqual(res.E_theta[0].shape, (N_THETA, N_PHI))

    def test_derived_quantities(self):
        _write_result(self.path)
        res = nf2ff_results(self.path)
        for n in range(N_FREQ):
            E_norm = np.sqrt(np.abs(res.E_theta[n])**2 + np.abs(res.E_phi[n])**2)
            np.testing.assert_allclose(res.E_norm[n], E_norm)
            self.assertEqual(res.E_cprh[n].shape, (N_THETA, N_PHI))


if __name__ == '__main__':
    unittest.main()
