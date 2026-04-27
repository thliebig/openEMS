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
from openEMS.ports import LumpedPort, RectWGPort, WaveguidePort


def _make_csx():
    csx = ContinuousStructure()
    grid = csx.GetGrid()
    grid.SetDeltaUnit(1e-3)
    grid.SetLines('x', np.linspace(-50, 50, 11))
    grid.SetLines('y', np.linspace(-50, 50, 11))
    grid.SetLines('z', np.linspace(-5, 5, 5))
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


if __name__ == '__main__':
    unittest.main()
