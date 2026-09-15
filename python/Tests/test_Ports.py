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

import unittest
import numpy as np

from CSXCAD import ContinuousStructure
from CSXCAD.CSProperties import CSPropExcitation, CSPropProbeBox
from CSXCAD.CSRectGrid import CoordinateSystem
from openEMS.ports import (LumpedPort, RectWGPort, WaveguidePort,
                            CircWGPort, CoaxialPort, StripLinePort, CPWPort, CurvePort)

C0 = 299792458.0  # speed of light


def _make_csx():
    csx = ContinuousStructure()
    grid = csx.GetGrid()
    grid.SetDeltaUnit(1e-3)
    grid.SetLines('x', np.linspace(-50, 50, 11))
    grid.SetLines('y', np.linspace(-50, 50, 11))
    grid.SetLines('z', np.linspace(-5, 5, 5))
    return csx


def _make_csx_circ():
    """Grid large enough for a 320 mm radius circular waveguide (z-propagation)."""
    csx = ContinuousStructure()
    grid = csx.GetGrid()
    grid.SetDeltaUnit(1e-3)
    grid.SetLines('x', np.linspace(-350, 350, 15))
    grid.SetLines('y', np.linspace(-350, 350, 15))
    grid.SetLines('z', np.linspace(0, 200, 15))
    return csx


def _make_csx_circ_3d():
    """Symmetric grid large enough for a 320 mm radius waveguide in any direction."""
    csx = ContinuousStructure()
    grid = csx.GetGrid()
    grid.SetDeltaUnit(1e-3)
    for ax in ('x', 'y', 'z'):
        grid.SetLines(ax, np.linspace(-350, 350, 15))
    return csx


def _make_csx_cylindrical():
    """Cylindrical (rho, a, z) grid for a circular waveguide of radius 350 mm."""
    csx = ContinuousStructure(CoordSystem=CoordinateSystem.CYLINDRICAL)
    grid = csx.GetGrid()
    grid.SetDeltaUnit(1e-3)
    grid.SetLines('r', np.linspace(0, 350, 15))
    grid.SetLines('a', np.linspace(0, 2 * np.pi, 25))
    grid.SetLines('z', np.linspace(0, 2000, 15))
    return csx


def _make_csx_coax():
    """Grid for a z-directed coaxial port (r_i=2, r_o=6, r_os=7 mm)."""
    csx = ContinuousStructure()
    grid = csx.GetGrid()
    grid.SetDeltaUnit(1e-3)
    grid.SetLines('x', np.linspace(-10, 10, 25))
    grid.SetLines('y', np.linspace(-10, 10, 25))
    grid.SetLines('z', np.linspace(0, 100, 15))
    return csx


def _make_csx_tl():
    """Grid for an x-propagating stripline / CPW port."""
    csx = ContinuousStructure()
    grid = csx.GetGrid()
    grid.SetDeltaUnit(1e-3)
    grid.SetLines('x', np.linspace(0, 100, 15))
    grid.SetLines('y', np.linspace(-8, 8, 20))
    grid.SetLines('z', np.linspace(-10, 10, 20))
    return csx


class Test_LumpedPort(unittest.TestCase):
    def setUp(self):
        self.csx = _make_csx()

    def test_passive_port(self):
        port = LumpedPort(self.csx, port_nr=1, R=50, start=[0, 0, -1], stop=[0, 0, 1], exc_dir='z', excite=0)
        self.assertEqual(port.R, 50)
        self.assertEqual(port.number, 1)

    def test_active_port(self):
        port = LumpedPort(self.csx, port_nr=1, R=50, start=[0, 0, -1], stop=[0, 0, 1], exc_dir='z', excite=1)
        self.assertEqual(port.excite, 1)

    def test_short_r_zero(self):
        port = LumpedPort(self.csx, port_nr=2, R=0, start=[-5, 0, -1], stop=[-5, 0, 1], exc_dir='z', excite=0)
        self.assertEqual(port.R, 0)

    def test_start_stop_same_raises(self):
        with self.assertRaises(Exception):
            LumpedPort(self.csx, port_nr=1, R=50, start=[0, 0, 0], stop=[0, 0, 0], exc_dir='z', excite=0)

    def test_x_direction(self):
        port = LumpedPort(self.csx, port_nr=1, R=50, start=[-1, 0, 0], stop=[1, 0, 0], exc_dir='x', excite=0)
        self.assertIsNotNone(port)

    def test_y_direction(self):
        port = LumpedPort(self.csx, port_nr=1, R=50, start=[0, -1, 0], stop=[0, 1, 0], exc_dir='y', excite=0)
        self.assertIsNotNone(port)

    def test_probe_filenames_set(self):
        port = LumpedPort(self.csx, port_nr=3, R=50, start=[0, 0, -1], stop=[0, 0, 1], exc_dir='z', excite=0)
        self.assertEqual(len(port.U_filenames), 1)
        self.assertEqual(len(port.I_filenames), 1)

    def test_port_number_in_probe_names(self):
        port = LumpedPort(self.csx, port_nr=7, R=50, start=[0, 0, -1], stop=[0, 0, 1], exc_dir='z', excite=0)
        self.assertIn('7', port.U_filenames[0])

    def test_name_prefix(self):
        port = LumpedPort(self.csx, port_nr=1, R=50, start=[0, 0, -1], stop=[0, 0, 1],
                          exc_dir='z', excite=0, PortNamePrefix='test_')
        self.assertTrue(port.U_filenames[0].startswith('test_'))

    def test_port_props_created(self):
        port = LumpedPort(self.csx, port_nr=1, R=50, start=[0, 0, -1], stop=[0, 0, 1], exc_dir='z', excite=0)
        # passive port: lumped element + voltage probe + current probe
        self.assertGreaterEqual(len(port.port_props), 3)

    def test_excitation_adds_prop(self):
        port_passive = LumpedPort(self.csx, port_nr=1, R=50, start=[0, 0, -1], stop=[0, 0, 1], exc_dir='z', excite=0)
        csx2 = _make_csx()
        port_active  = LumpedPort(csx2,       port_nr=1, R=50, start=[0, 0, -1], stop=[0, 0, 1], exc_dir='z', excite=1)
        self.assertGreater(len(port_active.port_props), len(port_passive.port_props))


class Test_RectWGPort(unittest.TestCase):
    def setUp(self):
        self.csx = _make_csx()

    def test_te10_mode(self):
        port = RectWGPort(self.csx, port_nr=1,
                          start=[0, 0, 0], stop=[0, 20, 10], exc_dir='x',
                          a=20e-3, b=10e-3, mode_name='TE10', excite=0)
        self.assertTrue(port.TE)
        self.assertFalse(port.TM)

    def test_te11_mode(self):
        port = RectWGPort(self.csx, port_nr=1,
                          start=[0, 0, 0], stop=[0, 20, 10], exc_dir='x',
                          a=20e-3, b=10e-3, mode_name='TE11', excite=0)
        self.assertTrue(port.TE)

    def test_mode_name_too_short_raises(self):
        with self.assertRaises(Exception):
            RectWGPort(self.csx, port_nr=1,
                       start=[0, 0, 0], stop=[0, 20, 10], exc_dir='x',
                       a=20e-3, b=10e-3, mode_name='TE1', excite=0)

    def test_mode_name_too_long_raises(self):
        with self.assertRaises(Exception):
            RectWGPort(self.csx, port_nr=1,
                       start=[0, 0, 0], stop=[0, 20, 10], exc_dir='x',
                       a=20e-3, b=10e-3, mode_name='TE101', excite=0)

    def test_tm_mode_raises(self):
        with self.assertRaises(Exception):
            RectWGPort(self.csx, port_nr=1,
                       start=[0, 0, 0], stop=[0, 20, 10], exc_dir='x',
                       a=20e-3, b=10e-3, mode_name='TM11', excite=0)

    def test_kc_te10(self):
        port = RectWGPort(self.csx, port_nr=1,
                          start=[0, 0, 0], stop=[0, 20, 10], exc_dir='x',
                          a=20e-3, b=10e-3, mode_name='TE10', excite=0)
        expected_kc = np.pi / 20e-3
        self.assertAlmostEqual(port.kc, expected_kc, places=5)

    def test_probe_filenames_set(self):
        port = RectWGPort(self.csx, port_nr=1,
                          start=[0, 0, 0], stop=[0, 20, 10], exc_dir='x',
                          a=20e-3, b=10e-3, mode_name='TE10', excite=0)
        self.assertEqual(len(port.U_filenames), 1)
        self.assertEqual(len(port.I_filenames), 1)


class Test_CircWGPort(unittest.TestCase):
    def setUp(self):
        self.csx    = _make_csx_circ()
        self.radius = 320e-3  # 320 mm

    def _make_port(self, mode='TE11', pol_ang=0, excite=0):
        return CircWGPort(self.csx, port_nr=1,
                          start=[0, 0, 0], stop=[0, 0, 200],
                          exc_dir='z', radius=self.radius,
                          mode_name=mode, pol_ang=pol_ang, excite=excite)

    def test_te11_kc(self):
        port = self._make_port('TE11')
        expected_kc = CircWGPort._pnm[(1, 1)] / self.radius
        self.assertAlmostEqual(port.kc, expected_kc, places=6)

    def test_te01_kc(self):
        port = self._make_port('TE01')
        expected_kc = CircWGPort._pnm[(0, 1)] / self.radius
        self.assertAlmostEqual(port.kc, expected_kc, places=6)

    def test_te21_kc(self):
        port = self._make_port('TE21')
        expected_kc = CircWGPort._pnm[(2, 1)] / self.radius
        self.assertAlmostEqual(port.kc, expected_kc, places=6)

    def test_te11_cutoff_frequency(self):
        port = self._make_port('TE11')
        fc = C0 * port.kc / (2 * np.pi)
        expected_fc = C0 * CircWGPort._pnm[(1, 1)] / (2 * np.pi * self.radius)
        self.assertAlmostEqual(fc / expected_fc, 1.0, places=10)

    def test_probe_filenames(self):
        port = self._make_port()
        self.assertEqual(len(port.U_filenames), 1)
        self.assertEqual(len(port.I_filenames), 1)

    def test_passive_port_no_excitation(self):
        port_passive = self._make_port(excite=0)
        port_active  = CircWGPort(self.csx, port_nr=2,
                                  start=[0, 0, 0], stop=[0, 0, 200],
                                  exc_dir='z', radius=self.radius,
                                  mode_name='TE11', excite=1)
        self.assertGreater(len(port_active.port_props), len(port_passive.port_props))

    def test_tm_mode_raises(self):
        with self.assertRaises(Exception):
            self._make_port('TM11')

    def test_unknown_mode_raises(self):
        with self.assertRaises(Exception):
            self._make_port('TE99')

    def test_pol_ang_accepted(self):
        port = self._make_port('TE11', pol_ang=np.pi / 2)
        self.assertIsNotNone(port)

    def test_port_number_in_probe_names(self):
        port = CircWGPort(self.csx, port_nr=7,
                          start=[0, 0, 0], stop=[0, 0, 200],
                          exc_dir='z', radius=self.radius, mode_name='TE11')
        self.assertIn('7', port.U_filenames[0])

    def test_name_prefix(self):
        port = CircWGPort(self.csx, port_nr=1,
                          start=[0, 0, 0], stop=[0, 0, 200],
                          exc_dir='z', radius=self.radius, mode_name='TE11',
                          PortNamePrefix='wg_')
        self.assertTrue(port.U_filenames[0].startswith('wg_'))

    def test_all_pnm_modes_accepted(self):
        for (n, m) in CircWGPort._pnm:
            mode = 'TE{}{}'.format(n, m)
            port = CircWGPort(_make_csx_circ(), port_nr=1,
                              start=[0, 0, 0], stop=[0, 0, 200],
                              exc_dir='z', radius=self.radius, mode_name=mode)
            self.assertAlmostEqual(port.kc, CircWGPort._pnm[(n, m)] / self.radius, places=6)

    def test_x_direction(self):
        csx = _make_csx_circ_3d()
        port = CircWGPort(csx, port_nr=1,
                          start=[0, 0, 0], stop=[200, 0, 0],
                          exc_dir='x', radius=self.radius, mode_name='TE11')
        self.assertAlmostEqual(port.kc, CircWGPort._pnm[(1, 1)] / self.radius, places=6)

    def test_y_direction(self):
        csx = _make_csx_circ_3d()
        port = CircWGPort(csx, port_nr=1,
                          start=[0, 0, 0], stop=[0, 200, 0],
                          exc_dir='y', radius=self.radius, mode_name='TE11')
        self.assertAlmostEqual(port.kc, CircWGPort._pnm[(1, 1)] / self.radius, places=6)

    def test_mode_funcs_use_transverse_coords(self):
        """E-field functions must reference the two transverse coordinates, not exc axis."""
        stops = {'x': [200, 0, 0], 'y': [0, 200, 0], 'z': [0, 0, 200]}
        transverse = {'x': ('y', 'z'), 'y': ('z', 'x'), 'z': ('x', 'y')}
        for exc_dir in ('x', 'y', 'z'):
            csx = _make_csx_circ_3d()
            port = CircWGPort(csx, port_nr=1,
                              start=[0, 0, 0], stop=stops[exc_dir],
                              exc_dir=exc_dir, radius=self.radius, mode_name='TE11')
            for func in port.E_func:
                if func == '0':
                    continue
                for t in transverse[exc_dir]:
                    self.assertIn(t, func, msg='coord {} missing in E_func for dir {}'.format(t, exc_dir))


class Test_CoaxialPort(unittest.TestCase):
    def setUp(self):
        self.csx = _make_csx_coax()
        self.pec  = self.csx.AddMetal('pec')
        self.mat  = self.csx.AddMaterial('fill', epsilon=2.1)
        self.kw   = dict(r_i=2, r_o=6, r_os=7)

    def _make_port(self, excite=0, feed_R=np.inf, extra_kw=None):
        kw = dict(**self.kw)
        kw['Feed_R'] = feed_R
        if extra_kw:
            kw.update(extra_kw)
        return CoaxialPort(self.csx, port_nr=1, pec_prop=self.pec, mat_prop=self.mat,
                           start=[0, 0, 0], stop=[0, 0, 100],
                           prop_dir='z', excite=excite, **kw)

    def test_passive_port(self):
        port = self._make_port()
        self.assertEqual(port.excite, 0)

    def test_active_port_has_more_props(self):
        passive = self._make_port(excite=0)
        active  = CoaxialPort(_make_csx_coax(), port_nr=1,
                               pec_prop=_make_csx_coax().AddMetal('p'),
                               mat_prop=None,
                               start=[0, 0, 0], stop=[0, 0, 100],
                               prop_dir='z', excite=1, **self.kw)
        self.assertGreater(len(active.port_props), len(passive.port_props))

    def test_three_voltage_probes(self):
        port = self._make_port()
        self.assertEqual(len(port.U_filenames), 3)

    def test_two_current_probes(self):
        port = self._make_port()
        self.assertEqual(len(port.I_filenames), 2)

    def test_radii_stored(self):
        port = self._make_port()
        self.assertEqual(port.r_i, 2)
        self.assertEqual(port.r_o, 6)

    def test_measplane_shift_set(self):
        port = self._make_port()
        self.assertGreater(port.measplane_shift, 0)

    def test_u_delta_positive(self):
        port = self._make_port()
        self.assertTrue(np.all(port.U_delta > 0))

    def test_feed_r_inf_open(self):
        port = self._make_port(feed_R=np.inf)
        self.assertIsNotNone(port)

    def test_feed_r_zero_short(self):
        port = self._make_port(feed_R=0)
        self.assertIsNotNone(port)

    def test_feed_r_positive_raises(self):
        with self.assertRaises(NotImplementedError):
            self._make_port(feed_R=50)

    def test_feed_r_negative_raises(self):
        with self.assertRaises(Exception):
            self._make_port(feed_R=-1)

    def test_no_material_fill(self):
        port = CoaxialPort(self.csx, port_nr=2, pec_prop=self.pec, mat_prop=None,
                           start=[0, 0, 0], stop=[0, 0, 100],
                           prop_dir='z', **self.kw)
        self.assertIsNotNone(port)

    def test_meas_plane_shift_kwarg(self):
        port = self._make_port(extra_kw={'MeasPlaneShift': 20})
        self.assertGreater(port.measplane_shift, 0)

    def test_port_number_in_probe_names(self):
        port = CoaxialPort(self.csx, port_nr=5, pec_prop=self.pec, mat_prop=None,
                           start=[0, 0, 0], stop=[0, 0, 100],
                           prop_dir='z', **self.kw)
        self.assertIn('5', port.U_filenames[0])


class Test_StripLinePort(unittest.TestCase):
    def setUp(self):
        self.csx   = _make_csx_tl()
        self.metal = self.csx.AddMetal('strip')

    def _make_port(self, excite=False, feed_R=np.inf, extra_kw=None):
        kw = {'Feed_R': feed_R}
        if extra_kw:
            kw.update(extra_kw)
        return StripLinePort(self.csx, port_nr=1, metal_prop=self.metal,
                             start=[0, -3, 0], stop=[100, 3, 0],
                             prop_dir='x', exc_dir='z', height=8,
                             excite=excite, **kw)

    def test_passive_port(self):
        port = self._make_port()
        self.assertIsNotNone(port)

    def test_six_voltage_probes(self):
        port = self._make_port()
        self.assertEqual(len(port.U_filenames), 6)

    def test_two_current_probes(self):
        port = self._make_port()
        self.assertEqual(len(port.I_filenames), 2)

    def test_probe_suffixes(self):
        port = self._make_port()
        self.assertTrue(any('A1' in fn for fn in port.U_filenames))
        self.assertTrue(any('A2' in fn for fn in port.U_filenames))
        self.assertTrue(any('B1' in fn for fn in port.U_filenames))
        self.assertTrue(any('C2' in fn for fn in port.U_filenames))

    def test_active_port_has_more_props(self):
        passive = self._make_port(excite=False)
        active  = StripLinePort(_make_csx_tl(), port_nr=1,
                                metal_prop=_make_csx_tl().AddMetal('s2'),
                                start=[0, -3, 0], stop=[100, 3, 0],
                                prop_dir='x', exc_dir='z', height=8,
                                excite=True)
        self.assertGreater(len(active.port_props), len(passive.port_props))

    def test_height_direction_mismatch_raises(self):
        with self.assertRaises(Exception):
            StripLinePort(self.csx, port_nr=1, metal_prop=self.metal,
                          start=[0, -3, 0], stop=[100, 3, 5],
                          prop_dir='x', exc_dir='z', height=8, excite=False)

    def test_measplane_shift_set(self):
        port = self._make_port()
        self.assertGreater(port.measplane_shift, 0)

    def test_u_delta_positive(self):
        port = self._make_port()
        self.assertTrue(np.all(port.U_delta > 0))

    def test_feed_r_finite(self):
        port = self._make_port(feed_R=50)
        self.assertIsNotNone(port)

    def test_feed_r_zero_short(self):
        port = self._make_port(feed_R=0)
        self.assertIsNotNone(port)

    def test_meas_plane_shift_kwarg(self):
        port = self._make_port(extra_kw={'MeasPlaneShift': 20})
        self.assertGreater(port.measplane_shift, 0)

    def test_port_number_in_probe_names(self):
        port = StripLinePort(self.csx, port_nr=9, metal_prop=self.metal,
                             start=[0, -3, 0], stop=[100, 3, 0],
                             prop_dir='x', exc_dir='z', height=8)
        self.assertIn('9', port.U_filenames[0])


class Test_CPWPort(unittest.TestCase):
    def setUp(self):
        self.csx   = _make_csx_tl()
        self.metal = self.csx.AddMetal('cpw')

    def _make_port(self, excite=False, feed_R=np.inf, extra_kw=None):
        kw = {'Feed_R': feed_R}
        if extra_kw:
            kw.update(extra_kw)
        return CPWPort(self.csx, port_nr=1, metal_prop=self.metal,
                       start=[0, -3, 0], stop=[100, 3, 0],
                       prop_dir='x', exc_dir='y', gap_width=1,
                       excite=excite, **kw)

    def test_passive_port(self):
        port = self._make_port()
        self.assertIsNotNone(port)

    def test_six_voltage_probes(self):
        port = self._make_port()
        self.assertEqual(len(port.U_filenames), 6)

    def test_two_current_probes(self):
        port = self._make_port()
        self.assertEqual(len(port.I_filenames), 2)

    def test_probe_suffixes(self):
        port = self._make_port()
        self.assertTrue(any('A1' in fn for fn in port.U_filenames))
        self.assertTrue(any('B2' in fn for fn in port.U_filenames))

    def test_active_port_has_more_props(self):
        passive = self._make_port(excite=False)
        active  = CPWPort(_make_csx_tl(), port_nr=1,
                          metal_prop=_make_csx_tl().AddMetal('c2'),
                          start=[0, -3, 0], stop=[100, 3, 0],
                          prop_dir='x', exc_dir='y', gap_width=1,
                          excite=True)
        self.assertGreater(len(active.port_props), len(passive.port_props))

    def test_exc_dir_is_the_field_across_the_gaps(self):
        port = self._make_port(excite=1)
        excitations = [prop.GetExcitation() for prop in port.port_props
                       if isinstance(prop, CSPropExcitation)]
        self.assertEqual(len(excitations), 2)
        for val in excitations:
            self.assertNotEqual(val[1], 0)
            self.assertEqual(val[0], 0)
            self.assertEqual(val[2], 0)

    def test_exc_dir_along_prop_dir_raises(self):
        with self.assertRaises(Exception):
            CPWPort(self.csx, port_nr=1, metal_prop=self.metal,
                    start=[0, -3, 0], stop=[100, 3, 0],
                    prop_dir='x', exc_dir='x', gap_width=1)

    def test_height_direction_mismatch_raises(self):
        with self.assertRaises(Exception):
            CPWPort(self.csx, port_nr=1, metal_prop=self.metal,
                    start=[0, -3, 0], stop=[100, 3, 5],
                    prop_dir='x', exc_dir='y', gap_width=1, excite=False)

    def test_measplane_shift_set(self):
        port = self._make_port()
        self.assertGreater(port.measplane_shift, 0)

    def test_u_delta_positive(self):
        port = self._make_port()
        self.assertTrue(np.all(port.U_delta > 0))

    def test_feed_r_finite(self):
        port = self._make_port(feed_R=50)
        self.assertIsNotNone(port)

    def test_feed_r_zero_short(self):
        port = self._make_port(feed_R=0)
        self.assertIsNotNone(port)

    def test_port_number_in_probe_names(self):
        port = CPWPort(self.csx, port_nr=4, metal_prop=self.metal,
                       start=[0, -3, 0], stop=[100, 3, 0],
                       prop_dir='x', exc_dir='y', gap_width=1)
        self.assertIn('4', port.U_filenames[0])


class Test_CurvePort(unittest.TestCase):
    def setUp(self):
        self.csx = _make_csx()

    def test_z_ref_stored(self):
        port = CurvePort(self.csx, port_nr=1, R=75,
                         start=[0, 0, -5], stop=[0, 0, 5])
        self.assertEqual(port.Z_ref, 75)
        self.assertEqual(port.R, 75)

    def test_single_voltage_probe(self):
        port = CurvePort(self.csx, port_nr=1, R=50,
                         start=[0, 0, -5], stop=[0, 0, 5])
        self.assertEqual(len(port.U_filenames), 1)

    def test_single_current_probe(self):
        port = CurvePort(self.csx, port_nr=1, R=50,
                         start=[0, 0, -5], stop=[0, 0, 5])
        self.assertEqual(len(port.I_filenames), 1)

    def test_port_direction_z(self):
        port = CurvePort(self.csx, port_nr=1, R=50,
                         start=[0, 0, -5], stop=[0, 0, 5])
        self.assertEqual(port.port_dir, 2)

    def test_port_direction_x(self):
        port = CurvePort(self.csx, port_nr=1, R=50,
                         start=[-50, 0, 0], stop=[50, 0, 0])
        self.assertEqual(port.port_dir, 0)

    def test_active_port_has_more_props(self):
        passive = CurvePort(self.csx, port_nr=1, R=50,
                            start=[0, 0, -5], stop=[0, 0, 5], excite=False)
        active  = CurvePort(self.csx, port_nr=2, R=50,
                            start=[0, 0, -5], stop=[0, 0, 5], excite=True)
        self.assertGreater(len(active.port_props), len(passive.port_props))

    def test_r_zero_short(self):
        port = CurvePort(self.csx, port_nr=1, R=0,
                         start=[0, 0, -5], stop=[0, 0, 5])
        self.assertIsNotNone(port)

    def test_r_inf_open(self):
        port = CurvePort(self.csx, port_nr=1, R=np.inf,
                         start=[0, 0, -5], stop=[0, 0, 5])
        self.assertIsNotNone(port)

    def test_multi_cell_port_constructs(self):
        # Port spanning many cells - triggers PEC curve creation
        port = CurvePort(self.csx, port_nr=1, R=50,
                         start=[0, 0, -5], stop=[0, 0, 4], excite=False)
        self.assertIsNotNone(port)

    def test_port_number_in_probe_names(self):
        port = CurvePort(self.csx, port_nr=6, R=50,
                         start=[0, 0, -5], stop=[0, 0, 5])
        self.assertIn('6', port.U_filenames[0])

    def test_name_prefix(self):
        port = CurvePort(self.csx, port_nr=1, R=50,
                         start=[0, 0, -5], stop=[0, 0, 5],
                         PortNamePrefix='cv_')
        self.assertTrue(port.U_filenames[0].startswith('cv_'))


class Test_WaveguidePort_LocalOrigin(unittest.TestCase):
    """The mode origin must reach the excitation and both mode-match probes
    identically. If the two ever diverge, the excited and the probed mode are
    no longer the same field and the mode match silently degrades -- which is
    exactly what a 'center' origin resolved in cylindrical mesh coordinates
    used to do to Circ_Waveguide.
    """

    def setUp(self):
        self.csx = _make_csx_circ()
        self.start = [0, 0, 0]
        self.stop  = [0, 0, 200]

    @staticmethod
    def _origins(port):
        """Collect (weight origins, mode origins) from a port's properties."""
        weight, mode = [], []
        for prop in port.port_props:
            if isinstance(prop, CSPropExcitation):
                weight.append(prop.GetWeightOrigin())
            elif isinstance(prop, CSPropProbeBox):
                mode.append(prop.GetModeOrigin())
        return weight, mode

    def test_excitation_and_probes_share_the_same_origin(self):
        port = CircWGPort(self.csx, port_nr=1, start=self.start, stop=self.stop,
                          exc_dir='z', radius=320e-3, mode_name='TE11', excite=1)
        weight, mode = self._origins(port)
        self.assertEqual(len(weight), 1)   # one excitation
        self.assertEqual(len(mode), 2)     # U and I mode-match probes
        for origin in weight + mode:
            np.testing.assert_allclose(origin, [0, 0, 100])

    def test_center_shorthand_is_the_box_midpoint(self):
        port = CircWGPort(self.csx, port_nr=1, start=[-10, 20, 0], stop=[10, 40, 200],
                          exc_dir='z', radius=320e-3, mode_name='TE11', excite=1,
                          local_origin='center')
        weight, mode = self._origins(port)
        for origin in weight + mode:
            np.testing.assert_allclose(origin, [0, 30, 100])

    def test_corner_shorthand_is_the_per_axis_minimum(self):
        port = RectWGPort(self.csx, port_nr=1, start=[10, -20, 200], stop=[-10, 20, 0],
                          exc_dir='z', a=100e-3, b=50e-3, mode_name='TE10', excite=1,
                          local_origin='corner')
        weight, mode = self._origins(port)
        for origin in weight + mode:
            np.testing.assert_allclose(origin, [-10, -20, 0])

    def test_explicit_vector_is_passed_through(self):
        port = CircWGPort(self.csx, port_nr=1, start=self.start, stop=self.stop,
                          exc_dir='z', radius=320e-3, mode_name='TE11', excite=1,
                          local_origin=[1.5, -2.5, 3.0])
        weight, mode = self._origins(port)
        for origin in weight + mode:
            np.testing.assert_allclose(origin, [1.5, -2.5, 3.0])

    def test_no_origin_leaves_everything_at_zero(self):
        port = WaveguidePort(self.csx, port_nr=1, start=self.start, stop=self.stop,
                             exc_dir='z', E_WG_func=['0', 'sin(x)', '0'],
                             H_WG_func=['sin(x)', '0', '0'], kc=1.0, excite=1)
        weight, mode = self._origins(port)
        for origin in weight + mode:
            np.testing.assert_allclose(origin, [0, 0, 0])

    def test_passive_port_still_gets_probe_origins(self):
        port = CircWGPort(self.csx, port_nr=2, start=self.start, stop=self.stop,
                          exc_dir='z', radius=320e-3, mode_name='TE11', excite=0)
        weight, mode = self._origins(port)
        self.assertEqual(weight, [])       # passive: no excitation
        self.assertEqual(len(mode), 2)
        for origin in mode:
            np.testing.assert_allclose(origin, [0, 0, 100])

    def test_unknown_shorthand_raises(self):
        with self.assertRaises(Exception):
            CircWGPort(self.csx, port_nr=1, start=self.start, stop=self.stop,
                       exc_dir='z', radius=320e-3, mode_name='TE11',
                       local_origin='middle')

    def test_shorthand_rejected_on_cylindrical_mesh(self):
        # start/stop are mesh coordinates there, so their midpoint is a
        # (rho,a,z) triple and not the Cartesian point the origin must be
        csx = _make_csx_cylindrical()
        for shorthand in ('center', 'corner'):
            with self.assertRaises(Exception):
                WaveguidePort(csx, port_nr=1, start=[0, 0, 100], stop=[350, 2*np.pi, 300],
                              exc_dir='z', E_WG_func=['0', 'sin(a)', '0'],
                              H_WG_func=['sin(a)', '0', '0'], kc=1.0, excite=1,
                              local_origin=shorthand)

    def test_explicit_vector_allowed_on_cylindrical_mesh(self):
        csx = _make_csx_cylindrical()
        port = WaveguidePort(csx, port_nr=1, start=[0, 0, 100], stop=[350, 2*np.pi, 300],
                             exc_dir='z', E_WG_func=['0', 'sin(a)', '0'],
                             H_WG_func=['sin(a)', '0', '0'], kc=1.0, excite=1,
                             local_origin=[0, 0, 200])
        weight, mode = self._origins(port)
        for origin in weight + mode:
            np.testing.assert_allclose(origin, [0, 0, 200])


class Test_UniquePortNumber(unittest.TestCase):
    """openEMS writes each port probe to a file named after the port number, so
    two ports with the same number would corrupt each other's files."""

    def _lumped(self, csx, port_nr, x, **kw):
        return LumpedPort(csx, port_nr=port_nr, R=50, start=[x, 0, -1], stop=[x, 0, 1],
                          exc_dir='z', **kw)

    def test_duplicate_lumped_port_raises(self):
        csx = _make_csx()
        self._lumped(csx, 1, -20, excite=1)
        with self.assertRaises(ValueError):
            self._lumped(csx, 1, 20)

    def test_duplicate_transmission_line_port_raises(self):
        csx = _make_csx_tl()
        metal = csx.AddMetal('strip')
        StripLinePort(csx, 1, metal, [0, -3, 0], [100, 3, 0], 'x', 'z', height=8, excite=1)
        with self.assertRaises(ValueError):
            StripLinePort(csx, 1, metal, [100, -3, 0], [0, 3, 0], 'x', 'z', height=8)

    def test_different_number_or_prefix_allowed(self):
        csx = _make_csx()
        self._lumped(csx, 1, -20, excite=1)
        self._lumped(csx, 2, 20)
        self._lumped(csx, 1, 0, PortNamePrefix='other_')


class Test_ExciteAmplitude(unittest.TestCase):
    """excite is an amplitude for every port type, not just an on/off switch."""

    @staticmethod
    def _excitations(port):
        return [prop.GetExcitation() for prop in port.port_props
                if isinstance(prop, CSPropExcitation)]

    def _check(self, make_port):
        unit   = self._excitations(make_port(1))
        scaled = self._excitations(make_port(-2))
        self.assertGreater(len(unit), 0)
        self.assertEqual(len(unit), len(scaled))
        for u, s in zip(unit, scaled):
            self.assertTrue(np.any(u != 0))
            np.testing.assert_allclose(s, -2*u)

    def test_lumped_port(self):
        self._check(lambda excite: LumpedPort(_make_csx(), port_nr=1, R=50,
                    start=[0, 0, -1], stop=[0, 0, 1], exc_dir='z', excite=excite))

    def test_curve_port(self):
        self._check(lambda excite: CurvePort(_make_csx(), port_nr=1, R=50,
                    start=[0, 0, -5], stop=[0, 0, 5], excite=excite))

    def test_waveguide_port(self):
        self._check(lambda excite: CircWGPort(_make_csx_circ(), port_nr=1,
                    start=[0, 0, 0], stop=[0, 0, 200], exc_dir='z',
                    radius=320e-3, mode_name='TE11', excite=excite))

    def test_coaxial_port(self):
        def make_port(excite):
            csx = _make_csx_coax()
            return CoaxialPort(csx, port_nr=1, pec_prop=csx.AddMetal('pec'), mat_prop=None,
                               start=[0, 0, 0], stop=[0, 0, 100], prop_dir='z',
                               r_i=2, r_o=6, r_os=7, excite=excite)
        self._check(make_port)

    def test_stripline_port(self):
        def make_port(excite):
            csx = _make_csx_tl()
            return StripLinePort(csx, port_nr=1, metal_prop=csx.AddMetal('strip'),
                                 start=[0, -3, 0], stop=[100, 3, 0],
                                 prop_dir='x', exc_dir='z', height=8, excite=excite)
        self._check(make_port)

    def test_cpw_port(self):
        def make_port(excite):
            csx = _make_csx_tl()
            return CPWPort(csx, port_nr=1, metal_prop=csx.AddMetal('cpw'),
                           start=[0, -3, 0], stop=[100, 3, 0],
                           prop_dir='x', exc_dir='y', gap_width=1, excite=excite)
        self._check(make_port)


if __name__ == '__main__':
    unittest.main()
