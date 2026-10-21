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
openems-view on a simulation folder as Run() leaves it: the setup Run()
writes, the figures Results stores there, and the field dumps. Synthetic
probe and dump files, no time stepping.
"""

import contextlib
import io
import json
import os
import shutil
import tempfile
import unittest

import numpy as np
import matplotlib
matplotlib.use('Agg')

from openEMS.openEMS import openEMS
from openEMS.results import Results

from test_setup_file import _make_fdtd, PORTS
from test_results import _probe_files
from test_viewer import WriteDump, LINES, _Field

try:
    import vtk  # noqa: F401
    HAS_VTK = True
except ImportError:
    HAS_VTK = False


class Test_ViewCLI(unittest.TestCase):
    def setUp(self):
        self.cwd = os.getcwd()   # Run() changes into the simulation path
        self.dir = tempfile.mkdtemp()
        self.fdtd, csx = _make_fdtd()
        self.port = PORTS['MSL'](self.fdtd, csx)   # port 7
        csx.AddDump('Et', file_type=1).AddBox([0, 0, 50], [100, 100, 50])
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(self.fdtd.Run(self.dir, setup_only=True), 0)
        _probe_files(self.dir, self.port)
        WriteDump(os.path.join(self.dir, 'Et.h5'), LINES, [_Field()], freqs=[1e9])

    def tearDown(self):
        os.chdir(self.cwd)
        shutil.rmtree(self.dir)

    def _results(self):
        res = Results(self.fdtd, self.dir)
        fig = res.AddFigure('Feed point', grid=(2, 1), sharex=True)
        fig.AddCurve(res.S(7, 7), label='$S_{77}$')
        fig.AddCurve(res.Zin(7), fmt='real', axes=1)
        return res

    def test_run_writes_setup(self):
        # the structure, ports and boxes, for post-processing without the script
        for name in ('openEMS_setup.xml', 'openEMS_setup.json'):
            self.assertTrue(os.path.isfile(os.path.join(self.dir, name)), name)
        fdtd = openEMS()
        fdtd.ReadFromXML(os.path.join(self.dir, 'openEMS_setup.xml'))
        fdtd.ReadSetup(os.path.join(self.dir, 'openEMS_setup.json'))
        self.assertEqual(fdtd.GetPort(7).U_filenames, self.port.U_filenames)

    def test_results_stored(self):
        # shown figures are stored next to that setup, the setup kept
        res = self._results()
        res.Show(block=False)
        with open(res.SetupFile()) as f:
            data = json.load(f)
        self.assertEqual([f['name'] for f in data['results']['figures']], ['Feed point'])
        self.assertEqual(len(data['ports']), 1)
        # without a setup of Run() nothing is written
        other = tempfile.mkdtemp()
        try:
            Results(self.fdtd, other, f=res.f).Show(block=False)
            self.assertEqual(os.listdir(other), [])
        finally:
            shutil.rmtree(other)
        import matplotlib.pyplot as plt
        plt.close('all')

    @unittest.skipUnless(HAS_VTK, 'vtk not available')
    def test_openems_view(self):
        from openEMS.Viewer import CLI
        self._results().Write()

        def run(*args):
            out, err = io.StringIO(), io.StringIO()
            with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
                code = CLI.Main([self.dir] + list(args))
            return code, out.getvalue(), err.getvalue()
        code, out, err = run('--describe', 'json')
        self.assertEqual(code, 0, err)
        d = json.loads(out)
        self.assertEqual(list(d['figures']), ['figures/Feed point'])
        self.assertEqual(d['figures']['figures/Feed point']['curves'], {'S₇₇': True, 'Re Zin7': True})
        et = d['items']['results/fields/Et']
        self.assertFalse(et['visible'])                       # the dumps start hidden
        self.assertEqual(et['settings']['quantity']['choices'], ['x', 'y', 'z', 'magnitude'])
        self.assertTrue(et['settings']['quantity']['help'])

        # the view's settings and visibility reach the items and figures added later
        code, out, err = run('--print-view', '--set', 'visible={"results/fields/Et": true, '
                             '"figures/Feed point/Re Zin7": false}',
                             '--set', 'items={"results/fields/Et": {"quantity": "z", "scale": "log"}}')
        self.assertEqual(code, 0, err)
        st = json.loads(out)
        self.assertTrue(st['visible']['results/fields/Et'])
        self.assertFalse(st['visible']['figures/Feed point/Re Zin7'])
        self.assertEqual((st['items']['results/fields/Et']['quantity'], st['items']['results/fields/Et']['scale']),
                         ('z', 'log'))

        figs = os.path.join(self.dir, 'figs')
        self.assertEqual(run('--figures', figs)[0], 0)
        self.assertEqual(os.listdir(figs), ['Feed_point.png'])

        # mistakes stop it, naming what is valid
        code, _, err = run('--print-view', '--set', 'items={"results/fields/Et": {"quantity": "q"}}')
        self.assertEqual(code, 1)
        self.assertIn('invalid quantity "q", use one of', err)
        code, _, err = run('--print-view', '--set', 'visible={"results/fields/E": true}')
        self.assertEqual(code, 1)
        self.assertIn('no tree path results/fields/E; the paths:', err)
        self.assertIn('results/fields/Et', err)
        code, _, err = run('--print-view', '--set', 'items={"results/E": {"scale": "log"}}')
        self.assertIn('no item or figure results/E; the items and figures: results/fields/Et', err)
        # a folder without the setup of Run() says how to get one
        empty = tempfile.mkdtemp()
        try:
            err = io.StringIO()
            with contextlib.redirect_stderr(err):
                self.assertEqual(CLI.Main([empty, '--describe']), 1)
            self.assertIn('no openEMS_setup.xml; Run() writes it', err.getvalue())
        finally:
            shutil.rmtree(empty)


    @unittest.skipUnless(HAS_VTK, 'vtk not available')
    def test_debug_pec(self):
        # the metal edges of debug_PEC: on the mesh lines of the strip, with
        # their unit and direction, shown hidden by openems-view
        from CSXCAD.Viewer import Viewer
        from CSXCAD.Viewer.Files import ReadFile, FileUnit
        from openEMS.Viewer import CLI
        sim = os.path.join(self.dir, 'pec')
        fdtd, csx = _make_fdtd()
        PORTS['MSL'](fdtd, csx)
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(fdtd.Run(sim, setup_only=True, debug_PEC=True), 0)
        self.assertEqual(CLI.DebugFiles(sim), [os.path.join(sim, 'PEC_dump.vtp')])
        data = ReadFile(CLI.DebugFiles(sim)[0])
        self.assertEqual(FileUnit(data), 1.0)                 # metres
        from vtkmodules.util.numpy_support import vtk_to_numpy
        direction = vtk_to_numpy(data.GetCellData().GetArray('direction'))
        self.assertEqual(len(direction), data.GetNumberOfLines())
        self.assertEqual(set(direction), {0, 1})              # a sheet in z: x and y edges
        np.testing.assert_allclose(data.GetBounds(), (0, 0.030, 0.040, 0.044, 0.010, 0.010), atol=1e-7)

        v = Viewer(fdtd.GetCSX(), offscreen=True)
        item = v.AddFile(CLI.DebugFiles(sim)[0])
        np.testing.assert_allclose(item.GetBounds(), v._scene.Entry('msl').GetBounds(), atol=1e-5)

        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            code = CLI.Main([sim, '--describe', 'json', '--set', 'isolate=type:Metal'])
        self.assertEqual(code, 0, err.getvalue())
        d = json.loads(out.getvalue())
        pec = d['items']['results/debug/PEC_dump.vtp']
        self.assertFalse(pec['visible'])
        self.assertEqual(pec['settings']['color_by']['choices'], ['solid', 'direction'])
        self.assertEqual([p['name'] for p in d['properties'] if p['shown']], ['msl'])

if __name__ == '__main__':
    unittest.main()
