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
openEMS.results on synthetic probe files and dumps, no simulation. The
far-field quantities need real near-field dumps; they are checked against
the tutorials' own formulas outside the unit tests.
"""

import json
import os
import shutil
import tempfile
import types
import unittest

import numpy as np
import matplotlib
matplotlib.use('Agg')

from openEMS.openEMS import openEMS
from openEMS.physical_constants import Z0
from openEMS.results import Results, Format

from test_setup_file import _make_fdtd, PORTS
from test_viewer import WriteDump, LINES, _Field


def _probe_files(path, port):
    t = np.linspace(0, 10e-9, 400)
    for k, fn in enumerate(port.U_filenames + port.I_filenames):
        sig = np.exp(-((t - 2e-9 - k * 0.02e-9) / 0.5e-9) ** 2) * (1 + 0.1 * k)
        np.savetxt(os.path.join(path, fn), np.c_[t, sig], header='t/s voltage', comments='% ')


def _ShortDipole(name, freq, theta, phi, center=None):
    """ Stands in for Results._NF2FF: the far field (r = 1 m) of a short
    dipole along z, its amplitude growing with the frequency (f / 1 GHz). """
    freq = np.atleast_1d(np.asarray(freq, float))
    th, ph = np.deg2rad(np.asarray(theta, float)), np.deg2rad(np.asarray(phi, float))
    TH, _ = np.meshgrid(th, ph, indexing='ij')
    amp = freq / 1e9
    E = [a * np.abs(np.sin(TH)) for a in amp]
    P = [e ** 2 / (2 * Z0) for e in E]
    prad = amp ** 2 / (2 * Z0) * 8 * np.pi / 3   # the integral of sin^2 over the sphere
    return types.SimpleNamespace(theta=th, phi=ph, r=1.0, freq=freq, Dmax=np.full(len(freq), 1.5), Prad=prad,
                                 E_norm=E, E_theta=E, E_phi=[0 * e for e in E], E_cprh=E, E_cplh=E, P_rad=P)


class Test_Results(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp()
        self.fdtd, csx = _make_fdtd()   # Gauss excitation 1 GHz +- 0.5 GHz
        self.port = PORTS['MSL'](self.fdtd, csx)   # port 7
        _probe_files(self.dir, self.port)

    def tearDown(self):
        shutil.rmtree(self.dir)

    def test_default_frequency(self):
        res = Results(self.fdtd, self.dir)
        self.assertAlmostEqual(res.f[0], 0.5e9)
        self.assertAlmostEqual(res.f[-1], 1.5e9)
        self.assertEqual(len(res.f), 401)
        with self.assertRaises(ValueError):
            Results(openEMS(), self.dir)   # no excitation, no f

    def test_port_quantities(self):
        f = np.linspace(0.5e9, 1.5e9, 11)
        res = Results(self.fdtd, self.dir, f=f)
        self.port.CalcPort(self.dir, f)
        x, s = res.S(7, 7).GetXY()
        np.testing.assert_allclose(x, f)
        np.testing.assert_allclose(s, self.port.uf_ref / self.port.uf_inc)
        np.testing.assert_allclose(res.Zin(7).GetXY()[1], self.port.uf_tot / self.port.if_tot)
        np.testing.assert_allclose(res.P_acc(7).GetXY()[1], self.port.P_acc)
        np.testing.assert_allclose(res.Yin(7).GetXY()[1], self.port.if_tot / self.port.uf_tot)
        np.testing.assert_allclose(res.Beta(7).GetXY()[1], self.port.beta)
        t, u = res.U(7, domain='t').GetXY()
        np.testing.assert_allclose(u, self.port.ut_tot)
        self.assertAlmostEqual(res.P_acc(7).At(0.95e9), np.interp(0.95e9, f, self.port.P_acc))

        res2 = Results(self.fdtd, self.dir, f=f)
        res2.P_acc(7).GetXY()                 # Results has done its CalcPort ...
        self.port.CalcPort(self.dir, f[:3])   # ... the script's own one, later, does not change its data
        np.testing.assert_allclose(res2.S(7, 7).GetXY()[1], s)
        self.port.CalcPort(self.dir, f)
        res.SetPort(7, ref_impedance=50)   # a scalar reference: also incident/reflected in time
        self.port.CalcPort(self.dir, f, ref_impedance=50)
        np.testing.assert_allclose(res.S(7, 7).GetXY()[1], self.port.uf_ref / self.port.uf_inc)
        np.testing.assert_allclose(res.U(7, 'inc', 't').GetXY()[1], self.port.ut_inc)
        with self.assertRaises(ValueError):
            res.SetPort(7, Z_ref=50)

    def test_format(self):
        v = np.array([10 + 0j, -1j])
        np.testing.assert_allclose(Format(v, 'dB'), [20, 0])
        np.testing.assert_allclose(Format(v, 'dB', power=True), [10, 0])
        np.testing.assert_allclose(Format(v, 'phase'), [0, -90])
        np.testing.assert_allclose(Format(v, 'imag'), [0, -1])
        np.testing.assert_allclose(Format([4, 9], 'amplitude', power=True), [2, 3])
        np.testing.assert_allclose(Format([4, -9], 'amplitude'), [4, 9])
        with self.assertRaises(ValueError):
            Format(v, 'dBm')

    def test_rcs_kinds(self):
        res = Results(self.fdtd, self.dir)
        q = res.RCS('nf2ff', 1e9, 'et', theta=90)
        self.assertEqual((q.KIND, q.X, q.XName()), ('curve', 'angle', 'phi'))
        q = res.RCS('nf2ff', [0.5e9, 1e9], 'et', theta=90, phi=180)   # backscatter over f
        self.assertEqual((q.KIND, q.X), ('curve', 'f'))
        self.assertEqual(res.RCS('nf2ff', 1e9, 'et', theta=90, phi=180).KIND, 'value')
        with self.assertRaises(ValueError):
            res.RCS('nf2ff', 1e9, 'et')
        with self.assertRaises(ValueError):
            res.RCS('nf2ff', [0.5e9, 1e9], 'et', phi=0)   # a cut needs one frequency

    def test_farfield_quantities(self):
        # the formulas, on the analytic far field of a short dipole
        res = Results(self.fdtd, self.dir)
        res._NF2FF = _ShortDipole
        ang, D = res.FarFieldCut('nf2ff', 1e9, phi=0, angles=np.arange(0, 181, 15.0)).GetXY()
        np.testing.assert_allclose(D, 1.5 * np.sin(np.deg2rad(ang)) ** 2, atol=1e-12)   # 1.5 sin^2
        _, D = res.FarFieldCut('nf2ff', 1e9, theta=90, component='phi').GetXY()
        np.testing.assert_allclose(D, 0)                     # no phi component
        self.assertAlmostEqual(res.Dmax('nf2ff', 1e9).GetValue(), 1.5)
        f = [1e9, 2e9]
        _, prad = res.Prad('nf2ff', f).GetXY()
        np.testing.assert_allclose(prad, [1, 4] / (2 * Z0) * 8 * np.pi / 3)
        self.assertAlmostEqual(res.Prad('nf2ff', 2e9).GetValue(), prad[1])
        res.SetFrequency(f)
        p_acc = res.P_acc(7).GetXY()[1]
        _, eff = res.Efficiency('nf2ff', f, 7).GetXY()
        np.testing.assert_allclose(eff, prad / p_acc)

        # RCS: 4 pi r^2 S_scattered / S_incident, the incident amplitude from the probe
        from openEMS.ports import UI_data
        t = np.linspace(0, 10e-9, 400)
        np.savetxt(os.path.join(self.dir, 'et'), np.c_[t, np.exp(-((t - 2e-9) / 0.5e-9) ** 2)])
        et = np.abs(np.asarray(UI_data(['et'], self.dir, f).ui_f_val[0]))
        _, rcs = res.RCS('nf2ff', f, 'et', theta=60, phi=0, E_amp=2.0).GetXY()
        np.testing.assert_allclose(rcs, 4 * np.pi * np.array([1, 4]) * np.sin(np.pi / 3) ** 2 / (4 * et ** 2))
        ang, cut = res.RCS('nf2ff', 1e9, 'et', phi=0, angles=[30, 90], E_amp=2.0).GetXY()
        np.testing.assert_allclose(cut, 4 * np.pi * np.sin(np.deg2rad(ang)) ** 2 / (4 * et[0] ** 2))

    def test_probe(self):
        # a probe file as recorded, and its spectrum at the result frequencies
        from openEMS.ports import UI_data
        t = np.linspace(0, 10e-9, 400)
        sig = np.exp(-((t - 2e-9) / 0.5e-9) ** 2)
        np.savetxt(os.path.join(self.dir, 'ut_probe'), np.c_[t, sig])
        res = Results(self.fdtd, self.dir, f=[0.5e9, 1e9])
        x, y = res.Probe('ut_probe', domain='t').GetXY()
        np.testing.assert_allclose(x, t)
        np.testing.assert_allclose(y, sig)
        q = res.Probe('ut_probe')
        self.assertEqual((q.Label(), q.XName()), ('ut_probe', 'frequency'))
        np.testing.assert_allclose(q.GetXY()[1], UI_data(['ut_probe'], self.dir, res.f).ui_f_val[0])

    def test_dump_slice(self):
        E = _Field()
        WriteDump(os.path.join(self.dir, 'E.h5'), LINES, [E], freqs=[1e9])
        res = Results(self.fdtd, self.dir, f=[0.9e9, 1.1e9])
        q = res.DumpSlice('E.h5', ('z', 4), component='z', freq=1e9)
        x, y, c = q.GetMap()
        np.testing.assert_allclose(x, LINES[0])
        np.testing.assert_allclose(y, LINES[1])
        np.testing.assert_allclose(c, np.abs(E[2][:, :, 4]), rtol=1e-6)
        fig = res.AddFigure('E', kind='map')
        e = fig.AddMap(q, norm=res.Array(res.f, [2.0, 2.0]))
        np.testing.assert_allclose(fig.MapValues(e)[2], np.abs(E[2][:, :, 4]) / 2, rtol=1e-6)
        self.assertEqual(q.AxisNames(), ['x', 'y'])
        e = fig.AddMap(res.DumpSlice('E.h5', ('z', 4), component='z'), norm=res.P_acc(7))
        with self.assertRaisesRegex(ValueError, 'frequency'):
            fig.MapValues(e)   # normalising needs the map's frequency

    def test_kind_checks(self):
        res = Results(self.fdtd, self.dir)
        fig = res.AddFigure('S')
        with self.assertRaises(TypeError):
            fig.AddMap(res.S(7, 7))
        with self.assertRaises(TypeError):
            res.S(7, 7).GetMap()
        with self.assertRaises(ValueError):
            fig.AddCurve(res.S(7, 7), axes=1)   # a 1x1 grid
        with self.assertRaises(ValueError):
            res.FarFieldCut('nf2ff', 1e9, phi=0, theta=90)
        with self.assertRaises(ValueError):
            res.AddFigure('x', kind='smith')

    def _figures(self, res):
        fig = res.AddFigure('S-Parameter', ylabel='S (dB)')
        fig.AddCurve(res.S(7, 7), label='S77')
        fig.AddCurve(res.S(7, 7), fmt='phase', right=True, style='r--')
        fig.AddHLine(-10)
        z = res.AddFigure('Zin', grid=(1, 2), ylim=(-100, 100))
        z.AddCurve(res.Zin(7), fmt='real')
        z.AddCurve(res.Zin(7), fmt='imag', style='--', axes=1, color='r')
        z.SetAxes(1, title='imaginary', ylim=[-50, 50])
        a = res.AddFigure('Pattern', kind='polar', theta_zero='N', theta_direction=-1)
        a.AddCurve(res.Array(np.arange(0, 360, 10.0), np.ones(36), label='iso', domain='angle', power=True),
                   fmt='dBi')
        WriteDump(os.path.join(self.dir, 'E.h5'), LINES, [_Field()], freqs=[1e9])
        m = res.AddFigure('E', grid=(1, 2), shared_clim=True)
        m.AddMap(res.DumpSlice('E.h5', ('z', 4), freq=1e9), norm=res.P_acc(7))
        m.AddMap(res.ArrayMap([0, 1], [0, 1, 2], np.ones((2, 3)), label='ones'), axes=1)
        R, A = np.meshgrid([1, 2], np.linspace(0, np.pi, 4), indexing='ij')   # curvilinear
        c = res.AddFigure('ring', kind='map')
        c.AddMap(res.ArrayMap(R*np.cos(A), R*np.sin(A), R, label='r'), shading='gouraud')

    def test_write_read(self):
        xml = os.path.join(self.dir, 'model.xml')
        self.fdtd.Write2XML(xml)
        res = Results(self.fdtd, self.dir)
        res.SetPort(7, ref_impedance=50)
        self._figures(res)
        sidecar = xml[:-4] + '.json'
        res.Write(sidecar)
        with open(sidecar) as fh:
            data = json.load(fh)
        self.assertEqual(len(data['ports']), 1)   # the setup part is kept
        self.assertEqual(data['results']['f'], {'start': 0.5e9, 'stop': 1.5e9, 'num': 401})
        self.assertTrue(os.path.exists(xml[:-4] + '.npz'))

        fdtd = openEMS()
        fdtd.ReadFromXML(xml)
        fdtd.ReadSetup(sidecar)
        res2 = Results.Load(fdtd, self.dir, sidecar)
        self.assertEqual(res2.ToDict(), res.ToDict())
        np.testing.assert_allclose(res2.S(7, 7).GetXY()[1], res.S(7, 7).GetXY()[1])
        arr = res2.GetFigure('Pattern').entries[0]['q']
        self.assertTrue(arr.POWER)
        np.testing.assert_allclose(arr.GetXY()[1], np.ones(36))
        np.testing.assert_allclose(res2.GetFigure('E').entries[1]['q'].GetMap()[2], np.ones((2, 3)))

    def test_draw(self):
        res = Results(self.fdtd, self.dir)
        self._figures(res)
        figs = res.Show(block=False)
        self.assertEqual(len(figs), 5)
        s = figs[0]
        left, right = s.axes[0], s.axes[1]
        self.assertEqual(left.get_ylabel(), 'S (dB)')
        self.assertEqual(right.get_ylabel(), 'S-parameter phase (deg)')
        self.assertEqual(right.yaxis.label.get_color(), 'r')            # in the colour of its curve
        self.assertEqual(len(left.get_legend().get_texts()), 2)         # one legend for both axes
        self.assertEqual(left.get_lines()[-1].get_ydata()[0], -10)      # the reference line
        p = figs[2].axes[0]
        self.assertEqual(p.get_theta_direction(), -1)
        self.assertAlmostEqual(p.get_theta_offset(), np.pi / 2)        # zero at the top
        z = figs[1]
        self.assertEqual([l.get_label() for l in z.axes[0].get_lines()], ['Re Zin7'])
        self.assertEqual(z.axes[0].get_xlabel(), 'frequency (GHz)')
        self.assertEqual(z.axes[0].get_ylabel(), 'impedance (Ohm)')
        self.assertEqual(z.axes[0].get_ylim(), (-100, 100))   # the figure's
        self.assertEqual(z.axes[1].get_ylim(), (-50, 50))     # the axes' own
        self.assertEqual(z.axes[1].get_title(), 'imaginary')
        with self.assertRaises(ValueError):
            res.GetFigure('Zin').SetAxes(0, colour='r')
        self.assertEqual(len(figs[3].axes), 4)   # two maps with colour bars
        import matplotlib.pyplot as plt
        plt.close('all')

    def test_sharex(self):
        res = Results(self.fdtd, self.dir)
        fig = res.AddFigure('Feed', grid=(2, 1), sharex=True)
        fig.AddCurve(res.S(7, 7), fmt='dB')
        fig.AddCurve(res.Zin(7), fmt='real', axes=1)
        self.assertTrue(res.ToDict()['figures'][0]['sharex'])
        top, bottom = fig.Draw().axes
        self.assertTrue(top.get_shared_x_axes().joined(top, bottom))
        self.assertEqual(top.get_xlabel(), '')                    # labelled at the bottom only
        # the stored tick setting: get_tick_params() is missing before matplotlib 3.7 and
        # names the x axis labels as y axis ones in 3.9
        labelbottom = lambda ax: ax.xaxis._major_tick_kw.get('label1On', True)
        self.assertFalse(labelbottom(top))
        self.assertTrue(labelbottom(bottom))
        self.assertEqual(bottom.get_xlabel(), 'frequency (GHz)')
        sidecar = os.path.join(self.dir, 'model.json')
        res.Write(sidecar)
        self.assertTrue(Results.Load(self.fdtd, self.dir, sidecar).GetFigure('Feed').sharex)
        import matplotlib.pyplot as plt
        plt.close('all')

    def test_legend(self):
        # a legend for several curves, or for a single one with a label of its own
        res = Results(self.fdtd, self.dir)
        none = res.AddFigure('none')
        none.AddCurve(res.S(7, 7))
        self.assertIsNone(none.Draw().axes[0].get_legend())
        one = res.AddFigure('one')
        one.AddCurve(res.P_acc(7), label='P')
        self.assertEqual([t.get_text() for t in one.Draw().axes[0].get_legend().get_texts()], ['P'])
        import matplotlib.pyplot as plt
        plt.close('all')

    def test_curve_names(self):
        # the tree name of a curve: its name=, else the label
        res = Results(self.fdtd, self.dir)
        fig = res.AddFigure('P')
        fig.AddCurve(res.P_acc(7), label='P')
        fig.AddCurve(res.P_inc(7), label='$P_{inc}$', name='incident')
        fig.AddCurve(res.P_ref(7))   # neither: the quantity's label
        fig.Draw()
        self.assertEqual(list(fig.curves), ['P', 'incident', 'P_ref7'])
        sidecar = os.path.join(self.dir, 'model.json')
        res.Write(sidecar)
        self.assertEqual([e.get('name') for e in Results.Load(self.fdtd, self.dir, sidecar).GetFigure('P').entries],
                         [None, 'incident', None])
        import matplotlib.pyplot as plt
        plt.close('all')

    def test_viewer_settings(self):
        # figure settings changed in the viewer change the description, which
        # ToDict exports, and draw the figure anew
        try:
            import vtk  # noqa: F401
        except ImportError:
            self.skipTest('vtk not available')
        from CSXCAD.Viewer import Viewer
        res = Results(self.fdtd, self.dir)
        self._figures(res)
        v = Viewer(self.fdtd.GetCSX(), offscreen=True)
        vfs = {vf.GetName(): vf for vf in res.AddToViewer(v, store=False)}
        s, z = vfs['S-Parameter'], vfs['Zin']
        names = {p['name'] for p in v.GetItemParameters(s)}
        self.assertTrue({'title', 'ylabel', 'ylim', 'ylabel_right', 'grid', 'legend', 'S77/label'} <= names, names)
        v.SetItemParameter(s, 'ylim', '-30, 0')
        v.SetItemParameter(s, 'grid', False)
        v.SetItemParameter(s, 'S77/label', 'S$_{77}$')
        v.SetItemParameter(s, 'S77/color', '#00ff00')
        v.SetItemParameter(z, 'title', 'Input impedance')
        v.SetItemParameter(z, 'ylabel.1', 'X (Ohm)')
        d = {f['name']: f for f in res.ToDict()['figures']}
        self.assertEqual(d['S-Parameter']['axes']['0'], dict(ylim=[-30.0, 0.0], grid=False))
        e = d['S-Parameter']['entries'][0]
        self.assertEqual((e['label'], e['name'], e['kw']), ('S$_{77}$', 'S77', dict(color='#00ff00')))
        self.assertEqual(d['Zin']['opts']['title'], 'Input impedance')
        self.assertEqual(d['Zin']['axes']['1']['ylabel'], 'X (Ohm)')
        ax = s.GetFigure().axes[0]   # drawn anew from the description
        self.assertEqual(ax.get_ylim(), (-30, 0))
        self.assertEqual(ax.get_legend().get_texts()[0].get_text(), 'S$_{77}$')
        self.assertEqual(s.GetCurves()[0], 'S77')   # its tree row keeps its name
        self.assertEqual(z.GetFigure()._suptitle.get_text(), 'Input impedance')
        v.SetItemParameter(s, 'ylim', '')   # automatic again
        self.assertNotIn('ylim', res.ToDict()['figures'][0]['axes']['0'])
        v.SetVisible(s.GetPath() + '/S77', False)   # hidden curves stay hidden
        v.SetItemParameter(s, 'ylabel', 'dB')
        self.assertFalse(s.GetCurveVisible('S77'))
        with self.assertRaises(ValueError):
            v.SetItemParameter(s, 'xscale', 'cubic')
        self.assertNotIn('xscale', res.ToDict()['figures'][0]['axes']['0'])
        self.assertEqual(s.GetFigure().axes[0].get_ylabel(), 'dB')   # drawn as before
        # the view state carries them, e.g. to a viewer opened later
        v2 = Viewer(self.fdtd.GetCSX(), offscreen=True, items={s.GetPath(): dict(ylim=[-20, 5])})
        res.AddToViewer(v2, store=False)
        self.assertEqual(v2.GetFigures()['S-Parameter'].GetFigure().axes[0].get_ylim(), (-20, 5))

    def test_show_in_viewer(self):
        # the output mode of the CSXCAD viewer: structure and figures as images
        try:
            import vtk  # noqa: F401
        except ImportError:
            self.skipTest('vtk not available')
        res = Results(self.fdtd, self.dir)
        self._figures(res)
        out = os.path.join(self.dir, 'out')
        os.mkdir(out)
        old = os.environ.get('CSXCAD_VIEWER_OUTPUT')
        os.environ['CSXCAD_VIEWER_OUTPUT'] = out
        try:
            from CSXCAD.Viewer import Viewer
            v = Viewer(self.fdtd.GetCSX(), offscreen=True, name='res')
            res.AddToViewer(v)
            v.Show()
        finally:
            if old is None:
                os.environ.pop('CSXCAD_VIEWER_OUTPUT')
            else:
                os.environ['CSXCAD_VIEWER_OUTPUT'] = old
        files = sorted(os.listdir(out))
        for name in ('res_view.png', 'res_S-Parameter.png', 'res_Zin.png', 'res_Pattern.png', 'res_E.png'):
            self.assertIn(name, files)
        self.assertEqual(v.GetFigures()['Zin'].GetCurves(), ['Re Zin7', 'Im Zin7'])


if __name__ == '__main__':
    unittest.main()
