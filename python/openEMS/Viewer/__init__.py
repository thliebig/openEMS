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

"""openEMS.Viewer
Field items for the CSXCAD viewer (:mod:`CSXCAD.Viewer`): near-field slices,
iso-surfaces and far-field patterns of openEMS results.

Everything is imported lazily: ``import openEMS.Viewer`` loads no VTK.

Examples
--------
>>> from CSXCAD.Viewer import Viewer
>>> from openEMS.Viewer import FieldSlices, FarField
>>> v = Viewer(CSX)
>>> v.AddItem(FieldSlices(os.path.join(sim_path, 'Et.h5'), planes=[('z', 1.5)], quantity='E_z', freq=2.4e9))
>>> v.AddItem(FarField(nf2ff_res, quantity='directivity', scale='log'))
>>> v.Show()
"""

__all__ = ['FieldSlices', 'FieldIsoSurface', 'FarField', 'FieldSource']

_MODULES = dict(FieldSlices='Slices', FieldIsoSurface='IsoSurface', FarField='Pattern', FieldSource='Source')


def __getattr__(name):
    if name in _MODULES:
        import importlib
        mod = importlib.import_module('openEMS.Viewer.' + _MODULES[name])
        return getattr(mod, name)
    raise AttributeError("module 'openEMS.Viewer' has no attribute '{}'".format(name))
