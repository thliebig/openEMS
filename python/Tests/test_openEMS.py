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
from openEMS.openEMS import openEMS


def _make_csx_with_grid():
    csx = ContinuousStructure()
    grid = csx.GetGrid()
    grid.SetDeltaUnit(1e-3)
    grid.SetLines('x', np.linspace(-50, 50, 11))
    grid.SetLines('y', np.linspace(-50, 50, 11))
    grid.SetLines('z', np.linspace(-5, 5, 5))
    return csx


class Test_Constructor(unittest.TestCase):
    def test_default(self):
        fdtd = openEMS()
        self.assertIsNotNone(fdtd)

    def test_nrts(self):
        fdtd = openEMS(NrTS=1e4)
        self.assertIsNotNone(fdtd)

    def test_end_criteria(self):
        fdtd = openEMS(EndCriteria=1e-6)
        self.assertIsNotNone(fdtd)

    def test_max_time(self):
        fdtd = openEMS(MaxTime=1e-9)
        self.assertIsNotNone(fdtd)

    def test_oversampling(self):
        fdtd = openEMS(OverSampling=10)
        self.assertIsNotNone(fdtd)

    def test_coord_system_cartesian(self):
        fdtd = openEMS(CoordSystem=0)
        self.assertIsNotNone(fdtd)

    def test_coord_system_cylindrical(self):
        fdtd = openEMS(CoordSystem=1)
        self.assertIsNotNone(fdtd)

    def test_timestep_factor(self):
        fdtd = openEMS(TimeStepFactor=0.9)
        self.assertIsNotNone(fdtd)

    def test_timestep_method(self):
        fdtd = openEMS(TimeStepMethod=1)
        self.assertIsNotNone(fdtd)

    def test_cell_constant_material(self):
        fdtd = openEMS(CellConstantMaterial=1)
        self.assertIsNotNone(fdtd)

    def test_multigrid(self):
        fdtd = openEMS(MultiGrid=[10.0])
        self.assertIsNotNone(fdtd)

    def test_unknown_kwarg_raises(self):
        with self.assertRaises(AssertionError):
            openEMS(UnknownOption=42)

    def test_multiple_kwargs(self):
        fdtd = openEMS(NrTS=5e4, EndCriteria=1e-5, TimeStepFactor=0.95)
        self.assertIsNotNone(fdtd)


class Test_CoordSystem(unittest.TestCase):
    def setUp(self):
        self.fdtd = openEMS()

    def test_cartesian(self):
        self.fdtd.SetCoordSystem(0)

    def test_cylindrical(self):
        self.fdtd.SetCoordSystem(1)

    def test_invalid_raises(self):
        with self.assertRaises(Exception):
            self.fdtd.SetCoordSystem(2)

    def test_negative_raises(self):
        with self.assertRaises(Exception):
            self.fdtd.SetCoordSystem(-1)


class Test_MultiGrid(unittest.TestCase):
    def setUp(self):
        self.fdtd = openEMS()

    def test_valid(self):
        self.fdtd.SetMultiGrid([10.0, 20.0])

    def test_empty_raises(self):
        with self.assertRaises(Exception):
            self.fdtd.SetMultiGrid([])


class Test_BoundaryConditions(unittest.TestCase):
    def setUp(self):
        self.fdtd = openEMS()

    def test_int_all_pec(self):
        self.fdtd.SetBoundaryCond([0, 0, 0, 0, 0, 0])

    def test_int_mixed(self):
        self.fdtd.SetBoundaryCond([0, 1, 2, 3, 0, 1])

    def test_string_pec(self):
        self.fdtd.SetBoundaryCond(['PEC', 'PEC', 'PEC', 'PEC', 'PEC', 'PEC'])

    def test_string_pmc(self):
        self.fdtd.SetBoundaryCond(['PMC', 'PMC', 'PMC', 'PMC', 'PMC', 'PMC'])

    def test_string_mur(self):
        self.fdtd.SetBoundaryCond(['MUR', 'MUR', 'MUR', 'MUR', 'MUR', 'MUR'])

    def test_string_pml(self):
        self.fdtd.SetBoundaryCond(['PML_8', 'PML_8', 'PML_8', 'PML_8', 'PML_8', 'PML_8'])

    def test_pml_custom_size(self):
        self.fdtd.SetBoundaryCond(['PML_16', 'PML_16', 'PML_8', 'PML_8', 'PEC', 'PEC'])

    def test_mixed_int_string(self):
        self.fdtd.SetBoundaryCond([0, 'PMC', 'MUR', 'PML_8', 'PEC', 1])

    def test_wrong_length_short_raises(self):
        with self.assertRaises(Exception):
            self.fdtd.SetBoundaryCond([0, 0, 0, 0, 0])

    def test_wrong_length_long_raises(self):
        with self.assertRaises(Exception):
            self.fdtd.SetBoundaryCond([0, 0, 0, 0, 0, 0, 0])

    def test_unknown_string_raises(self):
        with self.assertRaises(Exception):
            self.fdtd.SetBoundaryCond(['PEC', 'PEC', 'PEC', 'PEC', 'PEC', 'UNKNOWN'])


class Test_Excitation(unittest.TestCase):
    def setUp(self):
        self.fdtd = openEMS()

    def test_gauss(self):
        self.fdtd.SetGaussExcite(f0=1e9, fc=500e6)

    def test_gauss_zero_center(self):
        self.fdtd.SetGaussExcite(f0=0, fc=10e9)

    def test_sinus(self):
        self.fdtd.SetSinusExcite(2.4e9)

    def test_dirac(self):
        self.fdtd.SetDiracExcite(10e9)

    def test_step(self):
        self.fdtd.SetStepExcite(10e9)

    def test_custom(self):
        self.fdtd.SetCustomExcite('sin(2*pi*1e9*t)', f0=1e9, fmax=2e9)


class Test_SetCSX(unittest.TestCase):
    def setUp(self):
        self.fdtd = openEMS()
        self.csx  = ContinuousStructure()

    def test_set_and_get(self):
        self.fdtd.SetCSX(self.csx)
        self.assertIs(self.fdtd.GetCSX(), self.csx)

    def test_get_without_set_returns_none(self):
        self.assertIsNone(self.fdtd.GetCSX())

    def test_lumped_port_without_csx_raises(self):
        with self.assertRaises(Exception):
            self.fdtd.AddLumpedPort(1, 50, [0, 0, -1], [0, 0, 1], 'z', excite=1)

    def test_rect_wg_port_without_csx_raises(self):
        with self.assertRaises(Exception):
            self.fdtd.AddRectWaveGuidePort(1, [0, 0, 0], [0, 0, 10], 'z', 20e-3, 10e-3, 'TE10', excite=0)

    def test_msl_port_without_csx_raises(self):
        with self.assertRaises(Exception):
            metal = ContinuousStructure().AddMetal('metal')
            self.fdtd.AddMSLPort(1, metal, [0, 0, 0], [100, 0, 10], 'x', 'z', excite=0)

    def test_nf2ff_without_csx_raises(self):
        with self.assertRaises(Exception):
            self.fdtd.CreateNF2FFBox()


class Test_AddLumpedPort(unittest.TestCase):
    def setUp(self):
        self.fdtd = openEMS()
        self.csx  = _make_csx_with_grid()
        self.fdtd.SetCSX(self.csx)

    def test_passive_port(self):
        port = self.fdtd.AddLumpedPort(1, 50, [0, 0, -1], [0, 0, 1], 'z', excite=0)
        self.assertIsNotNone(port)

    def test_active_port(self):
        port = self.fdtd.AddLumpedPort(1, 50, [0, 0, -1], [0, 0, 1], 'z', excite=1)
        self.assertIsNotNone(port)

    def test_short_port_r_zero(self):
        port = self.fdtd.AddLumpedPort(2, 0, [-5, 0, -1], [-5, 0, 1], 'z', excite=0)
        self.assertIsNotNone(port)

    def test_edges2grid(self):
        n_lines_before = self.csx.GetGrid().GetQtyLines('z')
        self.fdtd.AddLumpedPort(1, 50, [0, 0, -2], [0, 0, 2], 'z', excite=0, edges2grid='z')
        n_lines_after = self.csx.GetGrid().GetQtyLines('z')
        self.assertGreaterEqual(n_lines_after, n_lines_before)


if __name__ == '__main__':
    unittest.main()
