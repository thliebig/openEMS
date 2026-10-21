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
Ports and NF2FF boxes written by WriteSetup and recreated by ReadSetup
against the structure read from the XML (openEMS#234).
"""

import json
import os
import shutil
import tempfile
import unittest

import numpy as np

from CSXCAD import ContinuousStructure
from openEMS import ports, nf2ff
from openEMS.openEMS import openEMS
from openEMS.setup_file import Setup


def _make_fdtd():
    fdtd = openEMS(NrTS=100)
    fdtd.SetGaussExcite(1e9, 0.5e9)
    fdtd.SetBoundaryCond(['MUR'] * 6)
    csx = ContinuousStructure()
    grid = csx.GetGrid()
    grid.SetDeltaUnit(1e-3)
    for n in range(3):
        grid.SetLines(n, np.arange(0, 101, 2.0))
    fdtd.SetCSX(csx)
    return fdtd, csx


# one port per type: name -> function(fdtd, csx) adding it
PORTS = {
    'Lumped':    lambda F, C: F.AddLumpedPort(1, 50, [10, 10, 0], [10, 10, 4], 'z', excite=1, priority=5),
    'Curve':     lambda F, C: F.AddCurvePort(2, 50, [20, 21, 0], [20, 21, 10], excite=1),
    'WaveGuide': lambda F, C: F.AddWaveGuidePort(3, [0, 0, 10], [40, 20, 14], 'z',
                                                 E_func=['0', 'sin(x/40*pi)', '0'], H_func=['sin(x/40*pi)', '0', '0'],
                                                 kc=np.pi / 40e-3, excite=1),
    'RectWG':    lambda F, C: F.AddRectWaveGuidePort(4, [0, 0, 20], [40, 20, 24], 'z', 40e-3, 20e-3, 'TE10', excite=1),
    'CircWG':    lambda F, C: F.AddCircWaveGuidePort(5, [40, 40, 20], [60, 60, 24], 'z', 10e-3, 'TE11', excite=1),
    'Coaxial':   lambda F, C: F.AddCoaxialPort(6, C.AddMetal('pec'), C.AddMaterial('teflon', epsilon=2.1),
                                               [50, 50, 0], [50, 50, 40], 'z', 2, 6, 8, excite=1),
    'MSL':       lambda F, C: F.AddMSLPort(7, C.AddMetal('msl'), [0, 40, 10], [30, 44, 14], 'x', 'z', excite=1,
                                           Feed_R=50, MeasPlaneShift=10),
    'StripLine': lambda F, C: F.AddStripLinePort(8, C.AddMetal('sl'), [0, 60, 20], [30, 64, 20], 'x', 'z', 6,
                                                 excite=1, Feed_R=np.inf),
    'CPW':       lambda F, C: F.AddCPWPort(9, C.AddMetal('cpw'), [0, 80, 30], [30, 84, 30], 'x', 'y', 2, excite=1,
                                           PortNamePrefix='a_'),
}

# what CalcPort and SetEnabled use
ATTRS = ['number', 'prefix', 'excite', 'delay', 'U_filenames', 'I_filenames', 'Z_ref', 'R', 'direction',
         'measplane_shift', 'U_delta', 'I_delta', 'kc', 'ref_index', 'feed_R']


class _Direct:
    """ Stands in for the openEMS object in PORTS: creates the port from its
    class, with the structure only, and adds it with AddPort. """
    CLASSES = {'AddLumpedPort': ports.LumpedPort, 'AddCurvePort': ports.CurvePort,
               'AddWaveGuidePort': ports.WaveguidePort, 'AddRectWaveGuidePort': ports.RectWGPort,
               'AddCircWaveGuidePort': ports.CircWGPort, 'AddCoaxialPort': ports.CoaxialPort,
               'AddMSLPort': ports.MSLPort, 'AddStripLinePort': ports.StripLinePort, 'AddCPWPort': ports.CPWPort}
    RENAME = {'E_func': 'E_WG_func', 'H_func': 'H_WG_func'}

    def __init__(self, fdtd, csx):
        self.fdtd, self.csx = fdtd, csx

    def __getattr__(self, name):
        cls = self.CLASSES[name]
        def add(*args, **kw):
            kw = {self.RENAME.get(k, k): v for k, v in kw.items()}
            return self.fdtd.AddPort(cls(self.csx, *args, **kw))
        return add


class Test_SetupFile(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp()
        self.xml = os.path.join(self.dir, 'model.xml')
        self.json = os.path.join(self.dir, 'model.json')

    def tearDown(self):
        shutil.rmtree(self.dir)

    def _write(self, *names, direct=False):
        fdtd, csx = _make_fdtd()
        adder = _Direct(fdtd, csx) if direct else fdtd
        orig = [PORTS[n](adder, csx) for n in names]
        fdtd.Write2XML(self.xml)
        return fdtd, orig

    def _read(self):
        fdtd = openEMS()
        fdtd.ReadFromXML(self.xml)
        fdtd.ReadSetup(self.json)
        return fdtd

    def test_round_trip(self):
        self._round_trip(direct=False)

    def test_round_trip_direct(self):
        # created from the port classes, added with AddPort
        self._round_trip(direct=True)

    def _round_trip(self, direct):
        for name in PORTS:
            with self.subTest(port=name):
                _, (orig,) = self._write(name, direct=direct)
                fdtd = self._read()
                csx = fdtd.GetCSX()
                n_props = len(csx.GetAllProperties())
                port = fdtd.GetPort(orig.number, orig.prefix)
                self.assertIs(type(port), type(orig))
                for attr in ATTRS:
                    if hasattr(orig, attr):
                        np.testing.assert_equal(getattr(port, attr), getattr(orig, attr), err_msg=attr)
                self.assertIs(port.CSX, csx)
                names = [p.GetName() for p in csx.GetAllProperties()]
                for prop in port.port_props:
                    self.assertIn(prop.GetName(), names)
                self.assertFalse(any(p.GetName().startswith('__openEMS_setup') for p in port.port_props))
                self.assertEqual(len(csx.GetAllProperties()), n_props)
                port.SetEnabled(False)   # acts on the loaded excitation
                exc = [p for p in csx.GetAllProperties() if 'excite' in p.GetName()]
                self.assertTrue(exc and not any(p.GetEnabled() for p in exc))

    def test_several_ports_and_order(self):
        _, orig = self._write('MSL', 'Lumped', 'CPW')
        fdtd = self._read()
        self.assertEqual([(p.number, p.prefix) for p in fdtd.GetPorts()],
                         [(p.number, p.prefix) for p in orig])

    def test_add_twice(self):
        fdtd, csx = _make_fdtd()
        port = PORTS['Lumped'](fdtd, csx)
        self.assertIs(fdtd.AddPort(port), port)
        self.assertEqual(len(fdtd.GetPorts()), 1)

    def test_own_subclass(self):
        class MyPort(ports.MSLPort):
            def __init__(self, CSX, metal, length, **kw):
                super().__init__(CSX, 7, metal, [0, 40, 10], [length, 44, 14], 'x', 'z', **kw)
        fdtd, csx = _make_fdtd()
        orig = fdtd.AddPort(MyPort(csx, csx.AddMetal('msl'), 30, excite=1))
        fdtd.Write2XML(self.xml)
        port = self._read().GetPort(7)
        self.assertIs(type(port), ports.MSLPort)   # stored as the openEMS class
        self.assertEqual(port.U_filenames, orig.U_filenames)

    def test_not_a_port(self):
        with self.assertRaises(TypeError):
            Setup().AddPort(object())

    def test_without_openEMS(self):
        _, csx = _make_fdtd()
        setup = Setup()
        port = setup.AddPort(ports.MSLPort(csx, 1, csx.AddMetal('msl'), [0, 40, 10], [30, 44, 14], 'x', 'z', excite=1))
        setup.AddNF2FF(nf2ff.nf2ff(csx, 'ff', [10, 10, 10], [90, 90, 90]))
        csx.Write2XML(self.xml)
        setup.Write(self.json)

        csx2 = ContinuousStructure()
        csx2.ReadFromXML(self.xml)
        setup2 = Setup()
        setup2.Read(self.json, csx2)
        self.assertEqual(setup2.GetPort(1).U_filenames, port.U_filenames)
        self.assertIs(setup2.GetPort(1).CSX, csx2)
        self.assertEqual(setup2.GetNF2FFBoxes()[0].name, 'ff')

    def test_file_contents(self):
        self._write('MSL', 'StripLine')
        with open(self.json) as f:
            text = f.read()
        data = json.loads(text)
        self.assertEqual(data['format'], 'openEMS-setup')
        self.assertEqual(data['version'], 1)
        msl, sl = data['ports']
        self.assertEqual(msl['type'], 'MSL')
        self.assertNotIn('metal_prop', msl['args'])   # user properties are not stored
        self.assertEqual(sl['args']['Feed_R'], 'inf')  # strict JSON
        self.assertNotIn('Infinity', text)

    def test_calc_port_same_result(self):
        _, (orig,) = self._write('MSL')
        t = np.linspace(0, 10e-9, 400)
        for k, fn in enumerate(orig.U_filenames + orig.I_filenames):
            sig = np.exp(-((t - 2e-9 - k * 0.02e-9) / 0.5e-9) ** 2) * (1 + 0.1 * k)
            np.savetxt(os.path.join(self.dir, fn), np.c_[t, sig], header='t/s voltage', comments='% ')
        f = np.linspace(0.5e9, 1.5e9, 11)
        orig.CalcPort(self.dir, f, ref_plane_shift=5)
        fdtd = self._read()   # owns the structure the port uses
        port = fdtd.GetPort(7)
        port.CalcPort(self.dir, f, ref_plane_shift=5)
        for attr in ['uf_inc', 'uf_ref', 'if_inc', 'Z_ref', 'beta', 'P_acc']:
            np.testing.assert_allclose(getattr(port, attr), getattr(orig, attr), rtol=1e-12, err_msg=attr)

    def test_nf2ff(self):
        fdtd, csx = _make_fdtd()
        box = fdtd.CreateNF2FFBox(frequency=[1e9])
        fdtd.Write2XML(self.xml)
        fdtd = self._read()
        (new,) = fdtd.GetNF2FFBoxes()
        self.assertEqual(new.name, box.name)
        np.testing.assert_equal(new.start, box.start)
        np.testing.assert_equal(new.stop, box.stop)
        self.assertEqual(new.directions, box.directions)
        self.assertEqual(list(new.freq), list(box.freq))
        self.assertIs(new.CSX, fdtd.GetCSX())
        self.assertEqual(len(fdtd.GetCSX().GetPropertiesByName('nf2ff_E')), 1)

    def _edit_xml(self, old, new):
        with open(self.xml) as f:
            text = f.read()
        self.assertIn(old, text)
        with open(self.xml, 'w') as f:
            f.write(text.replace(old, new, 1))

    def test_missing_property(self):
        self._write('Lumped')
        self._edit_xml('Name="port_ut_1"', 'Name="renamed"')
        with self.assertRaisesRegex(ValueError, 'Lumped port 1: property "port_ut_1".*not found'):
            self._read()

    def test_changed_property(self):
        self._write('Lumped')
        self._edit_xml('Name="port_ut_1" Visible="0" Number="0" Type="0" Weight="-1"',
                       'Name="port_ut_1" Visible="0" Number="0" Type="0" Weight="1"')
        with self.assertRaisesRegex(ValueError, 'port_ut_1" differs .*Weight'):
            self._read()

    def test_ambiguous_property(self):
        fdtd, csx = _make_fdtd()
        PORTS['Lumped'](fdtd, csx)
        csx.AddProbe('port_ut_1', p_type=0)   # same name and type, outside the port
        fdtd.Write2XML(self.xml)
        with self.assertRaisesRegex(ValueError, 'more than once'):
            self._read()

    def test_failed_read_keeps_ports(self):
        self._write('Lumped')
        fdtd = self._read()
        with open(self.json) as f:
            data = json.load(f)
        data['ports'][0]['args']['stop'] = [10, 10, 6]
        with open(self.json, 'w') as f:
            json.dump(data, f)
        with self.assertRaisesRegex(ValueError, 'differs'):
            fdtd.ReadSetup(self.json)
        self.assertEqual(len(fdtd.GetPorts()), 1)

    def test_look_is_ignored(self):
        _, (orig,) = self._write('Lumped')
        self._edit_xml('<FillColor', '<Attributes RenderMode="edges"/><FillColor')
        self.assertEqual(self._read().GetPort(1).U_filenames, orig.U_filenames)

    def test_unstorable_argument(self):
        fdtd, csx = _make_fdtd()
        fdtd.AddLumpedPort(1, 50, [10, 10, 0], [10, 10, 4], 'z', excite=1, tag=object())
        with self.assertRaises(TypeError):
            fdtd.WriteSetup(self.json)
        with self.assertWarnsRegex(UserWarning, 'setup file not written: Lumped port 1'):
            self.assertTrue(fdtd.Write2XML(self.xml))   # the XML is still written
        self.assertTrue(os.path.exists(self.xml))
        self.assertFalse(os.path.exists(self.json))


if __name__ == '__main__':
    unittest.main()
