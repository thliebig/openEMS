# -*- coding: utf-8 -*-
#
# Copyright (C) 2026 openEMS contributors
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

"""Tests for the nf2ff backends (CPU, GPU, auto).

The recording box is filled with the analytic field of a z-directed Hertzian
dipole, so the transform has a known answer: a sin(theta) pattern and a
directivity of 1.5. The GPU tests are skipped if openEMS was built without
NF2FF_HIP or no GPU is usable.
"""

import os
import shutil
import tempfile
import unittest
import numpy as np
import h5py

from openEMS import _nf2ff
from openEMS.nf2ff import nf2ff_results

C0 = 299792458.0
Z0 = 376.730313668
FREQ = 3e9                      # exact in float32, the recorded frequency is matched exactly
K = 2*np.pi*FREQ/C0
HALF = 0.06                     # half edge of the recording box [m], ~0.6 wavelengths
N_LINES = 25

THETA = np.deg2rad(np.arange(0, 181, 15.0))
PHI = np.deg2rad(np.arange(0, 360, 30.0))


def _dipole_field(x, y, z):
    """E and H of a z-directed Hertzian dipole (I*l = 1) at the origin, exp(+j*w*t)."""
    r = np.sqrt(x**2 + y**2 + z**2)
    cos_t = z/r
    sin_t = np.sqrt(1 - cos_t**2)
    phase = np.exp(-1j*K*r)
    jkr = 1j*K*r

    E_r = Z0*cos_t/(2*np.pi*r**2)*(1 + 1/jkr)*phase
    E_t = 1j*Z0*K*sin_t/(4*np.pi*r)*(1 + 1/jkr - 1/(K*r)**2)*phase
    H_p = 1j*K*sin_t/(4*np.pi*r)*(1 + 1/jkr)*phase

    rho = np.sqrt(x**2 + y**2)
    with np.errstate(invalid='ignore', divide='ignore'):
        cos_p, sin_p = np.where(rho > 0, x/rho, 1.0), np.where(rho > 0, y/rho, 0.0)
    # spherical -> Cartesian
    E = np.array([E_r*sin_t*cos_p + E_t*cos_t*cos_p,
                  E_r*sin_t*sin_p + E_t*cos_t*sin_p,
                  E_r*cos_t - E_t*sin_t])
    H = np.array([-H_p*sin_p, H_p*cos_p, np.zeros_like(H_p)])
    return E, H


def _write_plane(path, lines, field_index):
    """One recorded plane, laid out as openEMS dumps it: (3, nz, ny, nx), float32."""
    x, y, z = np.meshgrid(*lines, indexing='ij')
    data = _dipole_field(x, y, z)[field_index]            # (3, nx, ny, nz)
    data = np.transpose(data, (0, 3, 2, 1))
    with h5py.File(path, 'w') as h5:
        h5.attrs['openEMS_HDF5_version'] = np.float64(0.2)
        mesh = h5.create_group('Mesh')
        for ax, name in enumerate('xyz'):
            mesh[name] = lines[ax].astype(np.float32)
        fd = h5.create_group('FieldData/FD')
        fd.attrs['frequency'] = np.array([FREQ], dtype=np.float32)
        for part, val in (('real', data.real), ('imag', data.imag)):
            ds = fd.create_dataset('f0_' + part, data=val.astype(np.float32))
            ds.attrs['frequency'] = np.array([FREQ], dtype=np.float32)


def _write_box(path):
    """The six faces x-,x+,y-,y+,z-,z+ as nf2ff_E_<n>.h5 / nf2ff_H_<n>.h5."""
    line = np.linspace(-HALF, HALF, N_LINES)
    for n in range(6):
        lines = [line, line, line]
        lines[n//2] = np.array([-HALF if n % 2 == 0 else HALF])
        for kind, idx in (('E', 0), ('H', 1)):
            _write_plane(os.path.join(path, 'nf2ff_%s_%d.h5' % (kind, n)), lines, idx)


def _run(path, backend, out):
    """nf2ff on the recorded box with the given backend, returns the result."""
    nfc = _nf2ff._nf2ff([FREQ], THETA, PHI, [0, 0, 0])
    nfc.SetBackend(backend)
    for n in range(6):
        ok = nfc.AnalyseFile(os.path.join(path, 'nf2ff_E_%d.h5' % n),
                             os.path.join(path, 'nf2ff_H_%d.h5' % n))
        if not ok:
            raise RuntimeError('AnalyseFile failed')
    fn = os.path.join(path, out)
    nfc.Write2HDF5(fn)
    return nf2ff_results(fn)


class TestNF2FFBackend(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        _write_box(cls.tmp)
        cls.gpu = _nf2ff._nf2ff.GetGpuDevice() != ''

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def test_cpu_hertzian_dipole(self):
        """The CPU path reproduces the analytic dipole pattern and directivity."""
        res = _run(self.tmp, 'cpu', 'cpu.h5')
        E = res.E_norm[0]
        E = E/E.max()
        sin_t = np.sin(THETA)[:, None]*np.ones((1, len(PHI)))
        self.assertLess(np.abs(E - sin_t).max(), 0.02)
        self.assertAlmostEqual(res.Dmax[0], 1.5, delta=0.03)

    def test_auto_equals_cpu_or_gpu(self):
        """'auto' agrees with the CPU path whichever one it ends up using."""
        cpu = _run(self.tmp, 'cpu', 'cpu.h5')
        auto = _run(self.tmp, 'auto', 'auto.h5')
        peak = np.abs(cpu.E_theta[0]).max()
        np.testing.assert_allclose(auto.E_theta[0], cpu.E_theta[0], atol=1e-4*peak)
        np.testing.assert_allclose(auto.E_phi[0], cpu.E_phi[0], atol=1e-4*peak)
        np.testing.assert_allclose(auto.Prad, cpu.Prad, rtol=1e-6)

    def test_gpu_matches_cpu(self):
        if not self.gpu:
            self.skipTest('no GPU backend')
        cpu = _run(self.tmp, 'cpu', 'cpu.h5')
        gpu = _run(self.tmp, 'gpu', 'gpu.h5')
        peak = np.abs(cpu.E_theta[0]).max()
        np.testing.assert_allclose(gpu.E_theta[0], cpu.E_theta[0], atol=1e-4*peak)
        np.testing.assert_allclose(gpu.E_phi[0], cpu.E_phi[0], atol=1e-4*peak)
        np.testing.assert_allclose(gpu.Dmax, cpu.Dmax, rtol=1e-4)

    def test_gpu_unavailable_is_an_error(self):
        """Asking for the GPU must fail, not silently run on the CPU."""
        if self.gpu:
            self.skipTest('a GPU is available')
        with self.assertRaises(RuntimeError):
            _run(self.tmp, 'gpu', 'gpu.h5')

    def test_unknown_backend(self):
        nfc = _nf2ff._nf2ff([FREQ], THETA, PHI, [0, 0, 0])
        with self.assertRaises(ValueError):
            nfc.SetBackend('tpu')


if __name__ == '__main__':
    unittest.main()
