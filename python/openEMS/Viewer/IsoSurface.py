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

"""openEMS.Viewer.IsoSurface
The FieldIsoSurface item: iso-surfaces of a scalar in a volume dump.
"""

import numpy as np

from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkCommonDataModel import vtkStructuredGrid, vtkPolyData
from vtkmodules.vtkFiltersCore import vtkContourFilter, vtkPolyDataNormals
from vtkmodules.vtkRenderingCore import vtkActor, vtkPolyDataMapper
from vtkmodules.util.numpy_support import numpy_to_vtk, vtk_to_numpy

from openEMS.Viewer import Quantity
from openEMS.Viewer.Item import DumpItem


class FieldIsoSurface(DumpItem):
    """ FieldIsoSurface(source, quantity=None, levels=None, scale='linear', name=None, **settings)

    Iso-surfaces of a scalar quantity of a volume dump. Every keyword
    argument is a setting, stored in the view state and editable in the GUI.

    Parameters
    ----------
    source : str, HDF5Dump or FieldSource
        A volume dump (or computed data).
    quantity : str
        ``'magnitude'`` (default), a component, or ``'value'`` for a scalar
        dump (SAR). A component in linear scale is signed: its
        instantaneous value (``amplitude=False`` in the frequency domain,
        animated over the phase) has surfaces of either sign.
    levels : list of float
        Relative levels: fractions of the maximum in linear scale (default
        ``[0.5]``, ``[-0.5, 0.5]`` when signed), dB below the maximum in log
        scale (default ``[-6]``). In the time domain the maximum is the one
        over all animation steps.
    scale : 'linear' or 'log'
    color : str
        ``'cmap'`` (colour by level) or a colour.
    legend : bool
        Show the colour legend (with ``color='cmap'``).
    domain, freq, time_step, phase, amplitude, norm, unit :
        As for :class:`~openEMS.Viewer.Slices.FieldSlices`; the default is
        the amplitude.
    sample : int
        Spatial sub-sampling of the volume (every n-th line).
    time_sample, n_phase, animate :
        Animation, as for :class:`~openEMS.Viewer.Slices.FieldSlices`.
    name : str
        Item name, default the file stem plus ``' iso'``.
    """

    DEFAULT_GROUP = 'results/fields'
    AMPLITUDE_DEFAULT = True

    HELP = dict(
        DumpItem.HELP,
        levels='Values of the surfaces, relative to the peak: fractions in linear scale (negative ones for\n'
               'negative values of a signed component, e.g. "-0.1, 0.1"), dB below it in log scale\n'
               '(e.g. "-6, -20").',
        color='"cmap" (coloured by level) or one colour for all surfaces.',
    )

    def __init__(self, source, quantity=None, levels=None, scale='linear', dynamic_range=40.0,
                 value_range='auto', cmap=None, color='cmap', domain=None, freq=None, time_step=None, phase=0.0,
                 amplitude=None, sample=1, time_sample=1, n_phase=36, animate=True, norm=1.0, unit=None,
                 legend=True, name=None):
        src, name = self._Source(source, name, ' iso')
        super().__init__(name, src, quantity, domain, freq, time_step, phase, amplitude, sample, time_sample,
                         n_phase, animate, norm, unit, cmap, scale, dynamic_range, value_range, legend)
        self.levels = None if levels is None else Quantity.ParseNumbers(levels)
        self.color = color
        self._grid = None
        self._alpha = None
        self._surfaces = {}
        self.actor = vtkActor()
        m = vtkPolyDataMapper()
        m.SetLookupTable(self._lut)
        m.UseLookupTableScalarRangeOn()
        self.actor.SetMapper(m)
        self.AddProp(self.actor)
        self._BuildGrid()
        self._ResolveTimeStep()

    def Levels(self):
        if self.levels is not None:
            return list(self.levels)
        if self.scale == 'log':
            return [-6.0]
        return [-0.5, 0.5] if self._Signed() else [0.5]

    def _Signed(self):
        """ A component's instantaneous value in linear scale: surfaces of either sign. """
        return (self.quantity not in ('magnitude', 'value') and self.scale == 'linear'
                and not (self.domain == 'frequency' and self.amplitude))

    # ------------------------------------------------------------------
    def _BuildGrid(self):
        src = self.source
        lines = src.SampledLines(None, self.sample, self._delta)
        A, B, C = np.meshgrid(*lines, indexing='ij')
        if src.mesh_type == 'cylindrical':
            pts = np.stack([A * np.cos(B), A * np.sin(B), C], axis=-1)
            self._alpha = B
        else:
            pts = np.stack([A, B, C], axis=-1)
            self._alpha = None
        grid = vtkStructuredGrid()
        grid.SetDimensions(*A.shape)
        points = vtkPoints()
        points.SetData(numpy_to_vtk(pts.transpose(2, 1, 0, 3).reshape(-1, 3).astype(np.float64), deep=1))
        grid.SetPoints(points)
        self._grid = grid
        self._grid_delta = self._delta   # the drawing unit the points are in

    def _ReadTime(self, t):
        """ The quantity in the volume at time index `t` (signed). """
        data = self._Norm(np.asarray(self.source.Read(None, 0, self.sample, t_idx=t), dtype=np.float32))
        return Quantity.Evaluate(data, self.quantity, self.source.mesh_type, self._alpha)

    def _TimePeaks(self):
        """ The peak of every animation step and of the static one, read
        once without keeping the volumes: the levels then stay the same
        over all steps. """
        out = {}
        for t in self._TimeSteps():
            key = ('t', t, self.sample, self.quantity)
            if key not in self._peaks:
                vals = self._raw[key] if key in self._raw else self._ReadTime(t)
                self._peaks[key] = float(np.max(np.abs(vals))) if vals.size else 0.0
            out[t] = self._peaks[key]
        return out

    def _Values(self, step):
        """ The shown (linear) scalar volume at animation `step`. """
        src = self.source
        if self.domain == 'frequency':
            key = ('f', self.freq, self.sample)
            if key not in self._raw:
                self._raw[key] = self._Norm(src.ReadAtFrequency(self.freq, None, 0, self.sample))
            data = self._raw[key]
            vals = Quantity.Evaluate(data, self.quantity, src.mesh_type, self._alpha, self._Phase(step),
                                     self.amplitude)
            if key not in self._peaks:
                self._peaks[key] = Quantity.PeakAbs(data, self.quantity, src.mesh_type, self._alpha)
        else:
            key = ('t', self._TimeIndex(step), self.sample, self.quantity)
            if key not in self._raw:
                self._raw[key] = self._ReadTime(key[1])
            vals = self._raw[key]
        return vals if self._Signed() else np.abs(vals)

    def _Peak(self):
        if self.domain == 'time':
            return max(self._TimePeaks().values(), default=0.0)
        self._Values(self._step)
        return self._peaks[('f', self.freq, self.sample)]

    def _AbsLevels(self):
        """ Linear values of the levels. """
        peak = self._Peak()
        if self.scale == 'log':
            power = Quantity.IsPower(self.source.kind)
            return [peak * 10 ** (l / (10.0 if power else 20.0)) for l in self.Levels()]
        return [peak * l for l in self.Levels()]

    def _Surface(self, step):
        vals = self._Values(step)
        levels = self._AbsLevels()
        key = (step, self.quantity, self.amplitude, self.scale, self.phase, self.n_phase, self.domain, self.freq,
               self.time_step, self.time_sample, self.sample, tuple(levels))
        if key in self._surfaces:
            return self._surfaces[key]
        arr = numpy_to_vtk(np.ascontiguousarray(vals.transpose(2, 1, 0).ravel(), dtype=np.float32), deep=1)
        arr.SetName('value')
        self._grid.GetPointData().SetScalars(arr)
        cf = vtkContourFilter()
        cf.SetInputData(self._grid)
        for i, l in enumerate(levels):
            cf.SetValue(i, l)
        nrm = vtkPolyDataNormals()
        nrm.SetInputConnection(cf.GetOutputPort())
        nrm.Update()
        pd = vtkPolyData()
        pd.DeepCopy(nrm.GetOutput())
        if self.scale == 'log':
            # colour in dB
            sc = pd.GetPointData().GetScalars()
            if sc is not None and sc.GetNumberOfTuples():
                db = Quantity.ToDB(vtk_to_numpy(sc), Quantity.IsPower(self.source.kind), -1e9)
                a = numpy_to_vtk(db, deep=1)
                a.SetName('value')
                pd.GetPointData().SetScalars(a)
        self._surfaces[key] = pd
        return pd

    # ------------------------------------------------------------------
    def Update(self):
        if self.GetViewer() is None:
            return
        if self._grid is None or self._grid_delta != self._delta:   # e.g. the structure's unit once added
            self._BuildGrid()
            self._surfaces = {}
        self._RegisterAnimation()
        self._ShowStep(self._step)

    def _Invalidate(self):
        super()._Invalidate()
        self._surfaces = {}
        self._grid = None    # the sampling may have changed

    def _ShowStep(self, step):
        self._step = step
        if self._grid is None:
            self._BuildGrid()
        pd = self._Surface(step)
        m = self.actor.GetMapper()
        m.SetInputData(pd)
        levels = self._AbsLevels()
        power = Quantity.IsPower(self.source.kind)
        if self.scale == 'log':
            shown = [float(Quantity.ToDB(np.array([l]), power, -1e9)[0]) for l in levels]
        else:
            shown = levels
        signed = self._Signed()
        if self.value_range != 'auto':
            vrange = tuple(self.value_range)
        elif signed:
            top = max(abs(l) for l in shown) or 1.0
            vrange = (-top, top)
        else:
            lo, hi = min(shown), max(shown)
            if hi <= lo:
                lo, hi = lo - 1.0, hi + 1.0
            vrange = (lo, hi)
        self._ApplyLut(signed, vrange)
        if self.color == 'cmap':
            m.ScalarVisibilityOn()
        else:
            from CSXCAD.Viewer.ViewState import ParseColor
            m.ScalarVisibilityOff()
            self.actor.GetProperty().SetColor(*ParseColor(self.color)[0])
        self._ShowLegend(self.color == 'cmap')
        self._bar.SetTitle('{} [{}]'.format(self._Label(), self._Unit()))

    def OnPrepareAnimationStep(self, step):
        self._Values(step)

    def _Label(self):
        kind = self.source.label if self.quantity == 'value' else self.source.kind
        return Quantity.Label(kind, self.quantity, not self._Signed() and self.quantity != 'magnitude')

    def GetDescription(self):
        lv = ', '.join('{:g}{}'.format(l, ' dB' if self.scale == 'log' else '') for l in self.Levels())
        return ', '.join(p for p in ('{} at {}'.format(self._Label(), lv), self._FixedWhenText()) if p)

    def GetInfo(self):
        src = self.source
        n = int(np.prod([len(l) for l in src.SampledLines(None, self.sample)]))
        used = sum(a.nbytes for a in self._raw.values())
        info = dict(source=src.path or src.name, field=src.label, domain=self.domain,
                    levels=', '.join('{:.4g}'.format(l) for l in self._AbsLevels()), points=n,
                    memory='{:.1f} MB'.format(used / 1e6))
        if self.domain == 'time':
            steps = -(-len(src.times) // self.time_sample)
            info['memory'] += ' of {:.1f} MB for {} steps'.format(n * 4 * steps / 1e6, steps)
        return info

    def GetValueText(self, world_pos, prop, cell_id):
        if prop is not self.actor:
            return None
        pd = self.actor.GetMapper().GetInput()
        sc = pd.GetPointData().GetScalars() if pd is not None else None
        if sc is None or cell_id < 0 or cell_id >= pd.GetNumberOfCells():
            return None
        val = sc.GetValue(pd.GetCell(cell_id).GetPointId(0))
        return '{} = {:.4g} {} at {}'.format(self._Label(), val, self._Unit(), self._LengthText(world_pos))

    # ------------------------------------------------------------------
    def GetParameters(self):
        p = self._FieldParameters()
        p.append(dict(name='levels', type='str', value=self.Levels(), section='Field'))
        p += self._CommonParameters(self._Signed())
        p.append(dict(name='color', label='colour', type='str', value=self.color, default='cmap', section='Colours'))
        return p + self._DataParameters() + self._AnimationParameters()

    def _SetOwn(self, name, value):
        if name == 'levels':
            self.levels = Quantity.ParseNumbers(value)
        elif name == 'color':
            self.color = value
        else:
            return False
        return True
