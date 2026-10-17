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
# What compare_sim_results() judges a deviation against. Each of these choices
# decides whether two engines are reported as equal, and none of them shows in
# a comparison that passes.

import os
import tempfile
import unittest
import numpy as np
import h5py

from openEMS.compare import compare_sim_results


def write_run(path, probe, time, dump, e_phi=0.0):
    """ A minimal simulation folder: one probe, one field dump with two
    timesteps and a far field result """
    os.mkdir(path)
    np.savetxt(os.path.join(path, 'ut'), np.column_stack([time, probe]), header='t/s\tvoltage', comments='% ')
    with h5py.File(os.path.join(path, 'Et.h5'), 'w') as h:
        h['FieldData/TD/00000000'] = dump
        h['FieldData/TD/00000010'] = dump*1e-6      # a decayed timestep
    with h5py.File(os.path.join(path, 'nf2ff.h5'), 'w') as h:
        h['nf2ff/E_theta/FD/f0'] = np.ones(4)
        h['nf2ff/E_phi/FD/f0'] = np.full(4, e_phi)  # cross polarisation: numerical noise only


class Test_Compare(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.a = os.path.join(self.tmp.name, 'a')
        self.b = os.path.join(self.tmp.name, 'b')
        self.time = np.arange(5)*1e-12
        self.probe = np.array([0, 1, 2, 1, 0], dtype=float)
        self.dump = np.linspace(-1, 1, 24).reshape(3, 2, 2, 2)
        write_run(self.a, self.probe, self.time, self.dump, e_phi=1e-9)

    def tearDown(self):
        self.tmp.cleanup()

    def worst(self):
        return dict(compare_sim_results(self.a, self.b))

    def test_deviations_are_relative_to_the_peak_of_their_kind(self):
        probe = self.probe.copy()
        probe[1] += 2e-3                            # 1e-3 of the probe peak, not of the time values
        dump = self.dump.copy()
        write_run(self.b, probe, self.time, dump, e_phi=3e-9)
        with h5py.File(os.path.join(self.b, 'Et.h5'), 'a') as h:
            h['FieldData/TD/00000010'][0, 0, 0, 0] += 1e-7
        res = self.worst()
        self.assertAlmostEqual(res['ut'], 1e-3, delta=1e-9)
        # 10 % of the decayed timestep itself, but 1e-7 of the whole dump
        self.assertAlmostEqual(res['Et.h5:FieldData/TD/00000010'], 1e-7, delta=1e-12)
        self.assertEqual(res['Et.h5:FieldData/TD/00000000'], 0)
        # three times its own value, but 2e-9 of the co-polarisation
        self.assertAlmostEqual(res['nf2ff.h5:nf2ff/E_phi/FD/f0'], 2e-9, delta=1e-14)

    def test_what_cannot_be_compared_is_infinitely_off(self):
        # a run that stopped at another timestep, and one without the dump
        write_run(self.b, self.probe[:4], self.time[:4], self.dump, e_phi=1e-9)
        os.remove(os.path.join(self.b, 'Et.h5'))
        res = self.worst()
        self.assertEqual(res['ut'], np.inf)
        self.assertEqual(res['Et.h5'], np.inf)
        # the same values at shifted times are no agreement either
        self.tmp.cleanup()
        self.setUp()
        write_run(self.b, self.probe, self.time + 1e-12, self.dump, e_phi=1e-9)
        self.assertEqual(self.worst()['ut'], np.inf)
        # nor is an output that only the second run wrote
        self.tmp.cleanup()
        self.setUp()
        write_run(self.b, self.probe, self.time, self.dump, e_phi=1e-9)
        np.savetxt(os.path.join(self.b, 'it'), np.column_stack([self.time, self.probe]), comments='% ')
        self.assertEqual(self.worst()['it'], np.inf)


if __name__ == '__main__':
    unittest.main()
