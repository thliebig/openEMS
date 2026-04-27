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
import xml.etree.ElementTree as ET

import numpy as np

from CSXCAD import ContinuousStructure
from openEMS.openEMS import openEMS


def _make_fdtd(csx=None):
    """Return a fully configured openEMS object ready for Write2XML."""
    fdtd = openEMS(NrTS=1234, EndCriteria=1e-5)
    fdtd.SetBoundaryCond(['PML_8', 'PML_8', 'PEC', 'PEC', 'PMC', 'MUR'])
    fdtd.SetGaussExcite(f0=1e9, fc=500e6)
    if csx is None:
        csx = _make_csx()
    fdtd.SetCSX(csx)
    return fdtd


def _make_csx():
    csx = ContinuousStructure()
    grid = csx.GetGrid()
    grid.SetDeltaUnit(1e-3)
    grid.SetLines('x', np.linspace(-50, 50, 11))
    grid.SetLines('y', np.linspace(-50, 50, 11))
    grid.SetLines('z', np.linspace(-5, 5, 5))
    return csx


class Test_Write2XML(unittest.TestCase):
    def setUp(self):
        self.fn = os.path.join(tempfile.gettempdir(), 'test_openEMS_write.xml')
        self.fdtd = _make_fdtd()

    def tearDown(self):
        if os.path.exists(self.fn):
            os.remove(self.fn)

    def test_returns_true(self):
        self.assertTrue(self.fdtd.Write2XML(self.fn))

    def test_file_created(self):
        self.fdtd.Write2XML(self.fn)
        self.assertTrue(os.path.exists(self.fn))

    def test_file_nonempty(self):
        self.fdtd.Write2XML(self.fn)
        self.assertGreater(os.path.getsize(self.fn), 0)

    def test_valid_xml(self):
        self.fdtd.Write2XML(self.fn)
        tree = ET.parse(self.fn)
        self.assertIsNotNone(tree.getroot())

    def test_root_element_is_openEMS(self):
        self.fdtd.Write2XML(self.fn)
        root = ET.parse(self.fn).getroot()
        self.assertEqual(root.tag, 'openEMS')

    def test_fdtd_element_present(self):
        self.fdtd.Write2XML(self.fn)
        root = ET.parse(self.fn).getroot()
        self.assertIsNotNone(root.find('FDTD'))

    def test_nrts_in_xml(self):
        self.fdtd.Write2XML(self.fn)
        fdtd_el = ET.parse(self.fn).getroot().find('FDTD')
        self.assertEqual(int(fdtd_el.get('NumberOfTimesteps')), 1234)

    def test_end_criteria_in_xml(self):
        self.fdtd.Write2XML(self.fn)
        fdtd_el = ET.parse(self.fn).getroot().find('FDTD')
        self.assertAlmostEqual(float(fdtd_el.get('endCriteria')), 1e-5)

    def test_excitation_element_present(self):
        self.fdtd.Write2XML(self.fn)
        fdtd_el = ET.parse(self.fn).getroot().find('FDTD')
        self.assertIsNotNone(fdtd_el.find('Excitation'))

    def test_gauss_excitation_type(self):
        self.fdtd.Write2XML(self.fn)
        exc_el = ET.parse(self.fn).getroot().find('FDTD/Excitation')
        # Gaussian = type 0
        self.assertEqual(int(exc_el.get('Type')), 0)

    def test_gauss_f0_in_xml(self):
        self.fdtd.Write2XML(self.fn)
        exc_el = ET.parse(self.fn).getroot().find('FDTD/Excitation')
        self.assertAlmostEqual(float(exc_el.get('f0')), 1e9)

    def test_gauss_fc_in_xml(self):
        self.fdtd.Write2XML(self.fn)
        exc_el = ET.parse(self.fn).getroot().find('FDTD/Excitation')
        self.assertAlmostEqual(float(exc_el.get('fc')), 500e6)

    def test_boundary_cond_element_present(self):
        self.fdtd.Write2XML(self.fn)
        fdtd_el = ET.parse(self.fn).getroot().find('FDTD')
        self.assertIsNotNone(fdtd_el.find('BoundaryCond'))

    def test_boundary_conditions_in_xml(self):
        self.fdtd.Write2XML(self.fn)
        bc_el = ET.parse(self.fn).getroot().find('FDTD/BoundaryCond')
        self.assertIn('PML_', bc_el.get('xmin'))
        self.assertEqual(bc_el.get('ymin'), 'PEC')
        self.assertEqual(bc_el.get('zmin'), 'PMC')
        self.assertEqual(bc_el.get('zmax'), 'MUR')

    def test_sinus_excitation_type(self):
        fdtd = openEMS(NrTS=100)
        fdtd.SetCSX(_make_csx())
        fdtd.SetBoundaryCond([0]*6)
        fdtd.SetSinusExcite(2.4e9)
        fdtd.Write2XML(self.fn)
        exc_el = ET.parse(self.fn).getroot().find('FDTD/Excitation')
        self.assertEqual(int(exc_el.get('Type')), 1)

    def test_csx_structure_in_xml(self):
        self.fdtd.Write2XML(self.fn)
        root = ET.parse(self.fn).getroot()
        # CSXCAD structure is embedded as a child of the root
        self.assertIsNotNone(root.find('ContinuousStructure'))


class Test_ReadFromXML(unittest.TestCase):
    def setUp(self):
        self.fn = os.path.join(tempfile.gettempdir(), 'test_openEMS_read.xml')
        _make_fdtd().Write2XML(self.fn)

    def tearDown(self):
        if os.path.exists(self.fn):
            os.remove(self.fn)

    def test_returns_true_on_success(self):
        fdtd = openEMS()
        self.assertTrue(fdtd.ReadFromXML(self.fn))

    def test_read_resets_state(self):
        # ReadFromXML calls Reset() first; the read must still succeed
        fdtd = openEMS(NrTS=9999)
        self.assertTrue(fdtd.ReadFromXML(self.fn))

    def test_getcsx_not_none_after_read(self):
        fdtd = openEMS()
        fdtd.ReadFromXML(self.fn)
        self.assertIsNotNone(fdtd.GetCSX())

    def test_getcsx_replaced_after_read(self):
        csx_before = _make_csx()
        csx_before.AddMetal('pre_existing_metal')
        mat = csx_before.AddMaterial('pre_existing_dielectric')
        mat.SetMaterialProperty(epsilon=99)
        csx_before.GetGrid().SetDeltaUnit(1.0)  # intentionally wrong unit
        fdtd = openEMS()
        fdtd.SetCSX(csx_before)

        fdtd.ReadFromXML(self.fn)  # self.fn has no metals/dielectrics, DeltaUnit=1e-3

        csx_after = fdtd.GetCSX()
        self.assertIsNotNone(csx_after)

        # Grid delta unit must reflect the XML, not csx_before
        self.assertAlmostEqual(csx_after.GetGrid().GetDeltaUnit(), 1e-3)

        fn_out = os.path.join(tempfile.gettempdir(), 'test_getcsx_replaced.xml')
        try:
            fdtd.Write2XML(fn_out)
            with open(fn_out) as f:
                xml_str = f.read()
            self.assertNotIn('pre_existing_metal', xml_str)
            self.assertNotIn('pre_existing_dielectric', xml_str)
        finally:
            if os.path.exists(fn_out):
                os.remove(fn_out)


class Test_XML_RoundTrip(unittest.TestCase):
    """Write → Read → Write cycle: second file should be identical to first."""

    def setUp(self):
        self.fn1 = os.path.join(tempfile.gettempdir(), 'test_openEMS_rt1.xml')
        self.fn2 = os.path.join(tempfile.gettempdir(), 'test_openEMS_rt2.xml')

    def tearDown(self):
        for fn in (self.fn1, self.fn2):
            if os.path.exists(fn):
                os.remove(fn)

    def _xml_text(self, fn):
        with open(fn) as f:
            return f.read()

    def test_write_read_write_identical(self):
        _make_fdtd().Write2XML(self.fn1)

        fdtd2 = openEMS()
        fdtd2.ReadFromXML(self.fn1)
        fdtd2.Write2XML(self.fn2)

        self.assertEqual(self._xml_text(self.fn1), self._xml_text(self.fn2))

    def test_roundtrip_nrts(self):
        _make_fdtd().Write2XML(self.fn1)
        fdtd2 = openEMS()
        fdtd2.ReadFromXML(self.fn1)
        fdtd2.Write2XML(self.fn2)

        fdtd_el = ET.parse(self.fn2).getroot().find('FDTD')
        self.assertEqual(int(fdtd_el.get('NumberOfTimesteps')), 1234)

    def test_roundtrip_excitation(self):
        _make_fdtd().Write2XML(self.fn1)
        fdtd2 = openEMS()
        fdtd2.ReadFromXML(self.fn1)
        fdtd2.Write2XML(self.fn2)

        exc_el = ET.parse(self.fn2).getroot().find('FDTD/Excitation')
        self.assertEqual(int(exc_el.get('Type')), 0)
        self.assertAlmostEqual(float(exc_el.get('f0')), 1e9)
        self.assertAlmostEqual(float(exc_el.get('fc')), 500e6)

    def test_roundtrip_boundary_conditions(self):
        _make_fdtd().Write2XML(self.fn1)
        fdtd2 = openEMS()
        fdtd2.ReadFromXML(self.fn1)
        fdtd2.Write2XML(self.fn2)

        bc_el = ET.parse(self.fn2).getroot().find('FDTD/BoundaryCond')
        self.assertIn('PML_', bc_el.get('xmin'))
        self.assertIn('PML_', bc_el.get('xmax'))
        self.assertEqual(bc_el.get('ymin'), 'PEC')
        self.assertEqual(bc_el.get('ymax'), 'PEC')
        self.assertEqual(bc_el.get('zmin'), 'PMC')
        self.assertEqual(bc_el.get('zmax'), 'MUR')

    def test_roundtrip_csx_geometry(self):
        csx = _make_csx()
        metal = csx.AddMetal('test_metal')
        metal.AddBox([0, 0, 0], [10, 5, 2])
        _make_fdtd(csx).Write2XML(self.fn1)

        fdtd2 = openEMS()
        fdtd2.ReadFromXML(self.fn1)
        fdtd2.Write2XML(self.fn2)

        root2 = ET.parse(self.fn2).getroot()
        # Verify the metal property name survived
        xml_str = self._xml_text(self.fn2)
        self.assertIn('test_metal', xml_str)

    def test_roundtrip_cylindrical_coords(self):
        fdtd = openEMS(NrTS=100, CoordSystem=1)
        csx = ContinuousStructure()
        csx.SetMeshType(1)
        grid = csx.GetGrid()
        grid.SetDeltaUnit(1e-3)
        grid.SetLines('x', np.linspace(0, 50, 6))
        grid.SetLines('y', np.linspace(0, np.pi, 6))
        grid.SetLines('z', np.linspace(-5, 5, 5))
        fdtd.SetCSX(csx)
        fdtd.SetBoundaryCond([0]*6)
        fdtd.SetGaussExcite(0, 1e9)
        fdtd.Write2XML(self.fn1)

        fdtd2 = openEMS()
        fdtd2.ReadFromXML(self.fn1)
        fdtd2.Write2XML(self.fn2)

        fdtd_el = ET.parse(self.fn2).getroot().find('FDTD')
        self.assertEqual(int(fdtd_el.get('CylinderCoords')), 1)


if __name__ == '__main__':
    unittest.main()
