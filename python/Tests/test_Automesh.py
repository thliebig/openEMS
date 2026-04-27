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
from openEMS.automesh import mesh_hint_from_primitive, mesh_hint_from_box, mesh_combine


class Test_MeshHintFromBox(unittest.TestCase):
    def setUp(self):
        self.csx = ContinuousStructure()
        self.metal = self.csx.AddMetal('metal')

    def test_basic_hint(self):
        box = self.metal.AddBox([0, 0, 0], [10, 5, 2])
        hint = mesh_hint_from_box(box, 'xyz')
        self.assertIsNotNone(hint)
        self.assertIn(0,  hint[0])
        self.assertIn(10, hint[0])
        self.assertIn(0,  hint[1])
        self.assertIn(5,  hint[1])
        self.assertIn(0,  hint[2])
        self.assertIn(2,  hint[2])

    def test_single_direction(self):
        box = self.metal.AddBox([0, 0, 0], [10, 5, 2])
        hint = mesh_hint_from_box(box, 'x')
        self.assertIsNotNone(hint[0])
        self.assertIsNone(hint[1])
        self.assertIsNone(hint[2])

    def test_two_directions(self):
        box = self.metal.AddBox([0, 0, 0], [10, 5, 2])
        hint = mesh_hint_from_box(box, 'xy')
        self.assertIsNotNone(hint[0])
        self.assertIsNotNone(hint[1])
        self.assertIsNone(hint[2])

    def test_flat_box_single_point(self):
        # z-flat box: start[z]==stop[z], only one hint in z
        box = self.metal.AddBox([0, 0, 3], [10, 5, 3])
        hint = mesh_hint_from_box(box, 'z')
        self.assertEqual(len(hint[2]), 1)
        self.assertAlmostEqual(hint[2][0], 3.0)

    def test_metal_edge_res(self):
        box = self.metal.AddBox([0, 0, 0], [10, 5, 2])
        hint = mesh_hint_from_box(box, 'x', metal_edge_res=1.0)
        # with metal_edge_res we get sub-cell refinement: 4 points per edge
        self.assertGreater(len(hint[0]), 2)

    def test_up_dir_only(self):
        box = self.metal.AddBox([0, 0, 0], [10, 5, 2])
        hint = mesh_hint_from_box(box, 'x', up_dir=True, down_dir=False)
        self.assertIn(10, hint[0])
        self.assertNotIn(0, hint[0])

    def test_down_dir_only(self):
        box = self.metal.AddBox([0, 0, 0], [10, 5, 2])
        hint = mesh_hint_from_box(box, 'x', up_dir=False, down_dir=True)
        self.assertIn(0, hint[0])
        self.assertNotIn(10, hint[0])

    def test_combine_with_existing_mesh(self):
        box = self.metal.AddBox([0, 0, 0], [10, 5, 2])
        existing = [[100.0], None, None]
        hint = mesh_hint_from_box(box, 'x', mesh=existing)
        self.assertIn(100.0, hint[0])
        self.assertIn(0.0,   hint[0])


class Test_MeshHintFromPrimitive(unittest.TestCase):
    def setUp(self):
        self.csx = ContinuousStructure()
        self.metal = self.csx.AddMetal('metal')

    def test_box_primitive(self):
        box = self.metal.AddBox([1, 2, 3], [4, 5, 6])
        hint = mesh_hint_from_primitive(box, 'xyz')
        self.assertIsNotNone(hint)

    def test_unsupported_primitive_returns_none(self):
        sphere = self.metal.AddSphere([0, 0, 0], 5)
        hint = mesh_hint_from_primitive(sphere, 'xyz')
        self.assertIsNone(hint)


class Test_MeshCombine(unittest.TestCase):
    def test_both_none(self):
        result = mesh_combine([None, None, None], [None, None, None])
        self.assertEqual(result, [None, None, None])

    def test_first_none(self):
        result = mesh_combine([None, None, None], [[1, 2], None, None])
        self.assertEqual(result[0], [1, 2])

    def test_second_none(self):
        result = mesh_combine([[1, 2], None, None], [None, None, None])
        self.assertEqual(result[0], [1, 2])

    def test_merge_sorted(self):
        result = mesh_combine([[3, 1], None, None], [[2, 4], None, None])
        self.assertEqual(result[0], [1, 2, 3, 4])


if __name__ == '__main__':
    unittest.main()
