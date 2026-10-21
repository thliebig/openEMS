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

"""Tests for the field items of openEMS.Viewer, on synthetic data."""

import os
import shutil
import tempfile
import types
import unittest

import numpy as np
import h5py

try:
    import vtk  # noqa: F401
    from CSXCAD.Viewer import Viewer
    HAS_VIEWER = True
except ImportError:
    HAS_VIEWER = False

from CSXCAD import ContinuousStructure

COMPLEX_T = np.dtype([('r', np.float32), ('i', np.float32)])


def _HasGL():
    if not HAS_VIEWER:
        return False
    from CSXCAD.Viewer.Core import OffscreenBackend
    try:
        OffscreenBackend()
        return True
    except RuntimeError:
        return False


def WriteDump(path, lines, fields, freqs=None, times=None, mesh_type=0, scaling=1e-3, dump_type=None):
    """ Write a dump in the openEMS HDF5 layout: `lines` in drawing units
    (alpha in radians), `fields` a list of (3, Nx, Ny, Nz) / (Nx, Ny, Nz)
    arrays, one per frequency or time step. """
    names = ('rho', 'alpha', 'z') if mesh_type == 1 else ('x', 'y', 'z')
    with h5py.File(path, 'w') as h5:
        h5.attrs['openEMS_HDF5_version'] = np.float64(0.3)
        h5.attrs['dump_type'] = np.int32(dump_type if dump_type is not None else (10 if freqs is not None else 0))
        grp = h5.require_group('Mesh')
        grp.attrs['mesh_type'] = np.int32(mesh_type)
        grp.attrs['mesh_scaling'] = np.float64(scaling)
        for n, (name, l) in enumerate(zip(names, lines)):
            grp[name] = np.asarray(l, float) * (1.0 if (mesh_type == 1 and n == 1) else scaling)
        if freqs is not None:
            fd = h5.require_group('FieldData/FD')
            fd.attrs['frequency'] = np.asarray(freqs, np.float64)
            for n, data in enumerate(fields):
                if np.iscomplexobj(data):
                    buf = np.empty(data.shape, dtype=COMPLEX_T)
                    buf['r'], buf['i'] = data.real, data.imag
                else:
                    buf = data.astype(np.float32)
                ds = fd.create_dataset('f{}'.format(n), data=buf)
                ds.attrs['frequency'] = np.float64(freqs[n])
                ds.attrs['d_order'] = np.bytes_(b'NXYZ')
        else:
            td = h5.require_group('FieldData/TD')
            for n, (t, data) in enumerate(zip(times, fields)):
                ds = td.create_dataset('{:06d}'.format(10 * n), data=data.astype(np.float32))
                ds.attrs['time'] = np.float64(t)
                ds.attrs['d_order'] = np.bytes_(b'NXYZ')


LINES = (np.linspace(-10, 10, 11), np.linspace(-5, 5, 6), np.linspace(0, 8, 9))


def _Field():
    X, Y, Z = np.meshgrid(*LINES, indexing='ij')
    E = np.zeros((3,) + X.shape, complex)
    E[0] = 1.0 + 0j
    E[2] = (1 + X / 10) * np.exp(-1j * Z)
    return E


def _Structure():
    csx = ContinuousStructure()
    g = csx.GetGrid()
    g.SetDeltaUnit(1e-3)
    for d, l in zip('xyz', LINES):
        g.SetLines(d, l)
    csx.AddMetal('plate').AddBox([-8, -4, 0], [8, 4, 1])
    return csx


def _Pattern(n_theta=37, n_phi=73):
    th = np.linspace(0, np.pi, n_theta)
    ph = np.linspace(0, 2 * np.pi, n_phi)
    TH, PH = np.meshgrid(th, ph, indexing='ij')
    E = np.abs(np.sin(TH)) + 1e-6          # a short dipole along z
    P = E ** 2 / (2 * 377)
    prad = np.sum(P * np.sin(TH)) * (th[1] - th[0]) * (ph[1] - ph[0])
    return types.SimpleNamespace(theta=th, phi=ph, r=1.0, freq=np.array([1e9]), Prad=np.array([prad]),
                                 Dmax=np.array([4 * np.pi * P.max() / prad]), P_rad=[P], E_norm=[E],
                                 E_theta=[E], E_phi=[0 * E], E_cprh=[E], E_cplh=[E], fn=None)


@unittest.skipUnless(HAS_VIEWER, 'vtk or the CSXCAD viewer not available')
class Test_Viewer(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_source(self):
        # plane reads, drawing units and snapping, Cartesian and cylindrical
        from openEMS.Viewer import FieldSource, FieldSlices
        E = _Field()
        fn = os.path.join(self.tmp, 'Ef.h5')
        WriteDump(fn, LINES, [E], freqs=[1e9])
        src = FieldSource(fn)
        self.assertEqual((src.name, src.kind, src.unit, src.mesh_type), ('Ef', 'E', 'V/m', 'cartesian'))
        np.testing.assert_allclose(src.Lines(1e-3)[0], LINES[0])
        np.testing.assert_allclose(src.Read(2, 3, f_idx=0), E[:, :, :, 3], rtol=1e-6)
        sl = FieldSlices(src, planes=[('z', 2.6), ('x', None)])
        self.assertEqual(sl.GetPlanes(), [('z', 3.0), ('x', 0.0)])

        rho, alpha, z = np.linspace(1, 5, 5), np.linspace(0, np.pi, 7), np.linspace(0, 2, 3)
        H = np.random.RandomState(1).normal(size=(2, 3, 5, 7, 3)).astype(np.float32)
        fn = os.path.join(self.tmp, 'Ht.h5')
        WriteDump(fn, (rho, alpha, z), list(H), times=[1e-9, 2e-9], mesh_type=1, dump_type=1)
        cyl = FieldSlices(fn, planes=[('rho', 3)], quantity='x', time_step=0)
        self.assertEqual((cyl.source.kind, cyl.domain), ('H', 'time'))
        pts, _ = cyl._Points(0)
        np.testing.assert_allclose(np.linalg.norm(pts[..., :2], axis=-1), 3.0)   # a cylinder surface
        a = alpha[:, None]
        Hx = H[0, 0, 2] * np.cos(a) - H[0, 1, 2] * np.sin(a)
        frame = cyl._Frame(0)
        np.testing.assert_allclose(np.asarray(frame[0][0]).reshape(7, 3, order='F'), Hx, rtol=1e-5, atol=1e-6)

    def test_quantities(self):
        from openEMS.Viewer import Quantity as Q
        F = np.array([[1j], [0], [2]])
        self.assertAlmostEqual(float(Q.Evaluate(F, 'x', phase=90)[0]), -1.0, places=6)
        self.assertAlmostEqual(float(Q.Evaluate(F, 'x', amplitude=True)[0]), 1.0)
        self.assertAlmostEqual(float(Q.Evaluate(F, 'magnitude', amplitude=True)[0]), np.sqrt(5), places=6)
        self.assertAlmostEqual(float(Q.Evaluate(F, 'magnitude', phase=0)[0]), 2.0)
        np.testing.assert_allclose(Q.ToDB(np.array([10.0]), False, -100), [20.0])
        np.testing.assert_allclose(Q.ToDB(np.array([10.0]), True, -100), [10.0])
        self.assertEqual(Q.AutoRange(3.0, True), (-3.0, 3.0))
        self.assertEqual(Q.AutoRange(100.0, False, 'log', dynamic_range=30), (10.0, 40.0))
        self.assertEqual(Q.ParseRange('-1, 2'), [-1.0, 2.0])
        self.assertEqual(Q.ParseNumbers('[0, 0, 1e2]'), [0.0, 0.0, 100.0])
        self.assertEqual((Q.FormatTime(0), Q.FormatTime(2.5e-9)), ('0 s', '2.5 ns'))
        # arrows: the whole field for the magnitude, else only the component
        V = np.array([[1.0], [2.0], [3.0]])
        np.testing.assert_allclose(Q.ArrowVector(V, 'z'), [[0], [0], [3]])
        np.testing.assert_allclose(Q.ArrowVector(V, 'magnitude'), V)
        np.testing.assert_allclose(Q.ArrowVector(V, 'rho', 'cylindrical', np.array([np.pi / 2])),
                                   [[0], [1], [0]], atol=1e-12)

    def test_farfield(self):
        from openEMS.Viewer import FarField
        csx = _Structure()
        v = Viewer(csx)
        ff = v.AddItem(FarField(_Pattern(), scale='log', dynamic_range=20, size='50%'))
        self.assertEqual(ff.GetPath(), 'results/farfield/far field')
        b = v.GetStructureBounds()
        R = 0.5 * max(b[1] - b[0], b[3] - b[2], b[5] - b[4])
        pd, xyz, shown = ff._BuildPattern(0)
        rel = xyz - ff._Origin()
        self.assertAlmostEqual(float(np.max(np.linalg.norm(rel, axis=-1))), R, places=6)
        # the dipole null along z is at the -20 dB floor: radius 0
        self.assertAlmostEqual(float(np.linalg.norm(rel[0, 0])), 0.0, places=6)
        self.assertAlmostEqual(float(np.max(shown)), 10 * np.log10(1.5), delta=0.1)   # 1.76 dBi
        # the main lobe is in the xy plane (theta=90): the auto offset clears the plate along it
        th, ph, _ = ff.MainLobe()
        u = ff._Dir(th, ph)
        corners = np.array([[x, y, z] for x in b[:2] for y in b[2:4] for z in b[4:]])
        self.assertGreater(float(np.min(xyz.reshape(-1, 3) @ u)), float(np.max(corners @ u)))
        self.assertEqual(len(v.GetState()['items'][ff.GetPath()]['origin']), 3)
        # the -3 dB beam width of a short dipole: 90 degree (sin^2 = 1/2 at 45 and 135 degree)
        self.assertEqual(ff._BeamWidths(0), {'φ=0°': 90.0, 'φ=90°': 90.0})
        with self.assertRaisesRegex(ValueError, 'CalcNF2FF'):
            FarField(_Pattern(n_phi=2))

    def test_animation(self):
        from openEMS.Viewer import FieldSlices, FieldSource
        E = _Field()
        v = Viewer(_Structure())
        fd = v.AddItem(FieldSlices(FieldSource.FromArrays(LINES, E[None], freqs=[1e9], delta_unit=1e-3, name='Ef'),
                                   planes=[('y', 0)], quantity='z'))
        self.assertEqual(v.GetAnimationSteps(), 36)
        v.SetAnimationStep(9)    # 90 degree
        vals = np.asarray(fd._Frame(9)[0][0]).reshape(11, 9, order='F')
        np.testing.assert_allclose(vals, np.real(E[2, :, 2] * 1j), atol=1e-6)
        # the values divided by the normalisation; a complex one shifts the phase
        v.SetItemParameter(fd, 'norm', '2j')
        vals = np.asarray(fd._Frame(9)[0][0]).reshape(11, 9, order='F')
        np.testing.assert_allclose(vals, np.real(E[2, :, 2] * 1j / 2j), atol=1e-6)
        self.assertEqual(v.GetState()['items'][fd.GetPath()]['norm'], '2j')
        v.SetItemParameter(fd, 'norm', 1)
        # iso-surfaces of a component's instantaneous value: both signs
        from openEMS.Viewer import FieldIsoSurface
        iso = v.AddItem(FieldIsoSurface(FieldSource.FromArrays(LINES, E[None], freqs=[1e9], delta_unit=1e-3,
                                                               name='Ef'), quantity='z', amplitude=False))
        self.assertEqual((iso.Levels(), iso._Label()), ([-0.5, 0.5], 'E_z'))
        lo, hi = iso.actor.GetMapper().GetInput().GetPointData().GetScalars().GetRange()
        self.assertTrue(lo < 0 < hi, (lo, hi))
        v.RemoveItem(iso)

        # TD: the auto range scans every animation step once (when created);
        # the clock starts at the largest step, playback reads the others
        times = np.arange(5) * 1e-10
        fields = [np.sin(n) * np.ones((3, 11, 6, 9), np.float32) for n in range(5)]
        fn = os.path.join(self.tmp, 'Et.h5')
        WriteDump(fn, LINES, fields, times=times)
        td = FieldSlices(fn, planes=[('z', 4)], quantity='x', time_sample=2)
        reads = []
        read = td.source.Read
        td.source.Read = lambda *a, **kw: (reads.append(kw.get('t_idx')), read(*a, **kw))[1]
        v = Viewer(_Structure())
        v.AddItem(td)
        self.assertEqual(reads, [2])
        self.assertAlmostEqual(td._GlobalPeaks()[0], np.sin(2), places=6)
        self.assertEqual((td.time_step, v.GetAnimationStep()), (2, 1))
        reads.clear()
        v.Play()   # without an event loop: prepares, then stays at the step
        self.assertEqual(sorted(reads), [0, 4])
        v.SetAnimationStep(2)
        self.assertEqual(td.time_step, 4)
        v.SetItemParameter(td, 'time_step', 1)   # snapped to an animated step
        self.assertEqual((td.time_step, v.GetAnimationStep()), (0, 0))
        reads.clear()
        v.SetItemParameter(td, 'planes', [['z', 2]])
        self.assertEqual(sorted(reads), [0, 2, 4])
        # iso-surfaces start at the largest time step too
        self.assertEqual(FieldIsoSurface(fn, quantity='x', time_sample=2).time_step, 2)

    def test_state(self):
        from openEMS.Viewer import FieldSlices
        fn = os.path.join(self.tmp, 'Ef.h5')
        WriteDump(fn, LINES, [_Field()], freqs=[1e9])
        v = Viewer(_Structure())
        sl = v.AddItem(FieldSlices(fn))
        v.SetItemParameter(sl, 'plane0_pos', 5.2)
        v.SetItemParameter(sl, 'quantity', 'E_x')
        v.SetItemParameter(sl, 'contours', True)
        v.SetItemParameter(sl, 'value_range', '-2, 2')
        v.SetItemParameter(sl, 'legend', False)
        v.SetVisible(sl.GetPath(), False)
        v.SetVisible(sl.GetPath(), True)
        self.assertFalse(sl._bar.GetVisibility())   # stays off when shown again
        self.assertTrue(all(p.get('help') for p in v.GetItemParameters(sl)))   # every setting explained
        st = v.GetState()
        self.assertEqual(st['items']['results/fields/Ef']['planes'], [['z', 5.0]])
        v2 = Viewer(_Structure(), items=st['items'])
        sl2 = v2.AddItem(FieldSlices(fn))
        self.assertEqual(v2.GetState()['items'], st['items'])
        self.assertEqual(sl2.GetDescription(), 'E_x, 1 GHz')   # animated: no phase
        v2.SetItemParameter(sl2, 'animate', False)
        self.assertEqual(sl2.GetDescription(), 'E_x, 1 GHz, φ=0°')

    def test_plane_parts(self):
        # each plane of a FieldSlices is a tree row to show or hide
        from openEMS.Viewer import FieldSlices
        fn = os.path.join(self.tmp, 'Ef.h5')
        WriteDump(fn, LINES, [_Field()], freqs=[1e9])
        v = Viewer(_Structure())
        sl = v.AddItem(FieldSlices(fn, planes=[('z', 4), ('x', 0)]))
        rows = v.GetTree(sl.GetPath())['children']
        self.assertEqual([(r['kind'], r['label']) for r in rows], [('part', 'plane 1: z = 4'), ('part', 'plane 2: x = 0')])
        part = rows[1]['path']
        v.SetVisible(part, False)
        self.assertFalse(sl._planes[1].color.GetVisibility())
        self.assertTrue(sl._planes[0].color.GetVisibility())
        v.SetVisible(sl.GetPath(), False)
        v.SetVisible(sl.GetPath(), True)            # the item shown again keeps the plane hidden
        self.assertFalse(sl._planes[1].color.GetVisibility())
        st = v.GetState()
        self.assertFalse(st['visible'][part])
        v2 = Viewer(_Structure(), visible=st['visible'])   # before the item exists
        sl2 = v2.AddItem(FieldSlices(fn, planes=[('z', 4), ('x', 0)]))
        self.assertFalse(v2.GetVisible(part))
        self.assertFalse(sl2._planes[1].color.GetVisibility())
        sl.AddPlane()                                # from the context menu
        sl.Update()
        self.assertEqual(len(v.GetTree(sl.GetPath())['children']), 3)
        v.SetVisible(sl.GetPath() + '/plane 2', True)
        v.SetVisible(sl.GetPath() + '/plane 3', False)   # the copy of x = 0
        sl.RemovePlane(0)                                # the others move up, hidden stays hidden
        sl.Update()
        rows = v.GetTree(sl.GetPath())['children']
        self.assertEqual([(r['label'], r['visible']) for r in rows], [('plane 1: x = 0', True),
                                                                     ('plane 2: x = 0', False)])
        sl.SetPlanes([('z', 4)])
        sl.Update()
        self.assertEqual(v.GetTree(sl.GetPath())['children'], [])   # one plane: no rows

    def test_file_closed_when_removed(self):
        from openEMS.Viewer import FieldSlices
        fn = os.path.join(self.tmp, 'Ef.h5')
        WriteDump(fn, LINES, [_Field()], freqs=[1e9])
        v = Viewer(_Structure())
        sl = v.AddItem(FieldSlices(fn))
        h5 = sl.source._dump.file
        v.RemoveItem(sl)
        self.assertFalse(h5.id.valid)            # e.g. the simulation folder can be deleted
        v.AddItem(sl)                            # shown again, from the read data
        self.assertTrue(sl._planes[0].color.GetVisibility())
        v.SetItemParameter(sl, 'sample', 2)      # reads again: the file opens again
        self.assertTrue(sl.source._dump.file.id.valid)
        self.assertEqual(sl._planes[0].grid.GetNumberOfPoints(), 6 * 3)   # every 2nd line of 11 x 6

    def test_drawing_unit(self):
        # dumps are in m: the items place them in the structure's drawing unit
        from openEMS.Viewer import FieldIsoSurface, FieldSlices
        fn = os.path.join(self.tmp, 'Ef.h5')
        WriteDump(fn, LINES, [_Field()], freqs=[1e9])   # mm
        csx = _Structure()
        csx.GetGrid().SetDeltaUnit(1e-2)                 # cm
        v = Viewer(csx)
        iso = v.AddItem(FieldIsoSurface(fn))
        np.testing.assert_allclose(iso._grid.GetBounds()[:2], (-1, 1))
        sl = v.AddItem(FieldSlices(fn, planes=[('z', 0.4)]))
        self.assertEqual(sl.GetPlanes(), [('z', 0.4)])

    @unittest.skipUnless(_HasGL(), 'no offscreen OpenGL')
    def test_screenshot(self):
        from openEMS.Viewer import FieldSlices, FieldSource, FarField
        v = Viewer(_Structure(), offscreen=True, size=(400, 300))
        src = FieldSource.FromArrays(LINES, _Field()[None], freqs=[1e9], delta_unit=1e-3, name='Ef')
        sl = v.AddItem(FieldSlices(src, planes=[('z', 4)], quantity='z', glyphs=True, contours=True, warp=True))
        v.AddItem(FarField(_Pattern(), cuts=True))
        v.SetView('+z')
        img = v.Capture()
        self.assertGreater(len(np.unique(img.reshape(-1, img.shape[-1]), axis=0)), 50)
        # the value under the mouse: E_z = Re((1 + x/10) exp(-j z)) at z = 4, linear in x
        v.SetVisible('results/farfield', False)
        for name in ('glyphs', 'contours'):   # on the coloured plane itself
            v.SetItemParameter(sl, name, False)
        v.Capture()
        ren = v.GetRenderer()
        ren.SetWorldPoint(3, 1, 4, 1)   # a point of the plane
        ren.WorldToDisplay()
        px, py, _ = ren.GetDisplayPoint()
        hit = v.Pick(px, py)
        self.assertEqual(hit['item'], sl.GetPath())
        x = hit['world_pos'][0]
        text = v.GetItemValueText(px, py)
        self.assertRegex(text, r'^E_z = ')
        self.assertAlmostEqual(float(text.split()[2]), (1 + x / 10) * np.cos(4), places=3)
        for name in ('glyphs', 'contours'):
            v.SetItemParameter(sl, name, True)
        n = sl._planes[0].points.GetNumberOfPoints()
        v.SetItemParameter(sl, 'glyph_density', 2)
        self.assertLess(sl._planes[0].points.GetNumberOfPoints(), n)
        f = sl._planes[0].glyph.GetScaleFactor()
        v.SetItemParameter(sl, 'glyph_size', 3)
        self.assertAlmostEqual(sl._planes[0].glyph.GetScaleFactor(), 3 * f)


if __name__ == '__main__':
    unittest.main()
