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

"""openEMS.Viewer.Colors
Colour maps (by matplotlib name) as VTK lookup tables, and legends.
"""

import numpy as np

from vtkmodules.vtkCommonCore import vtkLookupTable
from vtkmodules.vtkRenderingAnnotation import vtkScalarBarActor

#: colour maps offered in the settings form (any matplotlib name works)
CMAPS = ['viridis', 'magma', 'inferno', 'plasma', 'cividis', 'turbo', 'jet', 'RdBu_r', 'coolwarm',
         'seismic', 'hot', 'gray']


def DefaultCmap(kind, signed):
    """ Diverging for signed values, viridis for magnitudes, magma for SAR. """
    if kind == 'SAR':
        return 'magma'
    return 'RdBu_r' if signed else 'viridis'


def LookupTable(cmap, n=256, lut=None):
    """ A ``vtkLookupTable`` of the matplotlib colour map `cmap` (fills
    `lut` if given). """
    import matplotlib
    try:
        cm = matplotlib.colormaps[cmap]
    except KeyError:
        raise ValueError('unknown colour map "{}"'.format(cmap))
    rgba = cm(np.linspace(0, 1, n))
    if lut is None:
        lut = vtkLookupTable()
    lut.SetNumberOfTableValues(n)
    for i, c in enumerate(rgba):
        lut.SetTableValue(i, *c)
    lut.Build()
    return lut


def ScalarBar(lut, title=''):
    """ A legend for `lut`; the viewer places it (``Viewer.AddLegend``). """
    bar = vtkScalarBarActor()
    bar.SetLookupTable(lut)
    bar.SetTitle(title)
    bar.SetNumberOfLabels(5)
    bar.SetLabelFormat('%-6.3g')
    bar.SetBarRatio(0.25)
    # keep the title clear of the top label
    bar.SetVerticalTitleSeparation(10)
    bar.UnconstrainedFontSizeOn()
    bar.GetTitleTextProperty().SetFontSize(13)
    bar.GetTitleTextProperty().BoldOn()
    bar.GetTitleTextProperty().ItalicOff()
    bar.GetLabelTextProperty().SetFontSize(12)
    bar.GetLabelTextProperty().BoldOff()
    bar.GetLabelTextProperty().ItalicOff()
    return bar
