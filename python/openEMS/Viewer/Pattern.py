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

"""openEMS.Viewer.Pattern
The FarField item: a 3D far-field pattern from an NF2FF result.
"""

import os

import numpy as np

from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkCommonDataModel import vtkStructuredGrid, vtkPolyData, vtkCellArray
from vtkmodules.vtkFiltersGeometry import vtkStructuredGridGeometryFilter
from vtkmodules.vtkFiltersCore import vtkPolyDataNormals
from vtkmodules.vtkRenderingCore import vtkActor, vtkPolyDataMapper
from vtkmodules.util.numpy_support import numpy_to_vtk

from CSXCAD.Viewer.Scene import UnionBounds

from openEMS.Viewer import Quantity
from openEMS.Viewer.Item import FieldItem

QUANTITIES = ['directivity', 'gain', 'E_norm', 'E_theta', 'E_phi', 'RHCP', 'LHCP']
_LABELS = dict(directivity='Directivity', gain='Gain', E_norm='|E|', E_theta='|E_theta|', E_phi='|E_phi|',
               RHCP='|E_RHCP|', LHCP='|E_LHCP|')
_SHORT = dict(directivity='D', gain='G', E_norm='|E|', E_theta='|E_θ|', E_phi='|E_φ|', RHCP='|E_R|',
              LHCP='|E_L|')


def _Lines(points, polylines):
    pd = vtkPolyData()
    pts = vtkPoints()
    pts.SetData(numpy_to_vtk(np.asarray(points, dtype=np.float64), deep=1))
    cells = vtkCellArray()
    for line in polylines:
        cells.InsertNextCell(len(line), [int(i) for i in line])
    pd.SetPoints(pts)
    pd.SetLines(cells)
    return pd


def _Actor(pd=None, port=None):
    m = vtkPolyDataMapper()
    if port is not None:
        m.SetInputConnection(port)
    else:
        m.SetInputData(pd if pd is not None else vtkPolyData())
    a = vtkActor()
    a.SetMapper(m)
    return a


class FarField(FieldItem):
    """ FarField(result, quantity='directivity', scale='log', freq=None, origin='auto', size='100%', name=None, **settings)

    A 3D far-field pattern: the radius and the colour show the quantity.
    Every keyword argument is a setting (see :meth:`GetParameters`), stored
    in the view state and editable in the GUI.

    Parameters
    ----------
    result : nf2ff_results or str
        The result of ``nf2ff.CalcNF2FF`` or its result file, computed on a
        full theta x phi grid, e.g. ``theta=np.arange(0, 181, 2)``,
        ``phi=np.arange(0, 361, 4)``.
    quantity : str
        ``'directivity'`` (dBi), ``'gain'`` (needs `accepted_power`),
        ``'E_norm'``, ``'E_theta'``, ``'E_phi'``, ``'RHCP'`` or ``'LHCP'``.
    scale : 'linear' or 'log'
        Linear: radius proportional to the value; log: to the dB value above
        ``max - dynamic_range``.
    dynamic_range : float
        dB range shown in log scale.
    value_range : 'auto' or (vmin, vmax)
        Colour range.
    cmap : str
        A matplotlib colour map name.
    legend : bool
        Show the colour legend.
    freq : float
        Frequency in Hz; default the first.
    sweep : bool
        Animate over all frequencies.
    origin : 'auto' or sequence of 3 float
        Where the pattern is drawn, in drawing units. ``'auto'`` starts at
        the centre of the shown structure and items and moves along the main
        lobe until the pattern clears them (also items added before, e.g. a
        field plane); it is then reported as the resolved numbers.
    size : str or float
        Largest radius: ``'100%'`` of the largest extent of the visible
        structure, or a length in drawing units.
    cuts : bool
        The principal-plane cuts (phi=0, phi=90, theta=90 degree) as lines
        on the pattern.
    accepted_power : float or array
        Accepted power per frequency in W, for the gain.
    name : str
        Item name, default the result file stem or ``'far field'``.
    """

    DEFAULT_GROUP = 'results/farfield'

    def __init__(self, result, quantity='directivity', scale='log', dynamic_range=40.0, value_range='auto',
                 cmap=None, freq=None, sweep=False, origin='auto', size='100%',
                 cuts=False, accepted_power=None, legend=True, name=None):
        if isinstance(result, (str, os.PathLike)):
            from openEMS.nf2ff import nf2ff_results
            result = nf2ff_results(os.fspath(result))
        self.result = result
        th, ph = np.atleast_1d(result.theta), np.atleast_1d(result.phi)
        if len(th) < 3 or len(ph) < 3:
            raise ValueError('the far-field result has {} theta and {} phi values, a 3D pattern needs a full '
                             'grid: call CalcNF2FF with e.g. theta=np.arange(0, 181, 2), '
                             'phi=np.arange(0, 361, 4)'.format(len(th), len(ph)))
        if name is None:
            fn = getattr(result, 'fn', None)
            name = os.path.splitext(os.path.basename(fn))[0] if fn else 'far field'
        super().__init__(name, cmap, scale, dynamic_range, value_range, legend)
        self._freqs = np.atleast_1d(np.asarray(result.freq, dtype=float))
        self.accepted_power = None if accepted_power is None else \
            np.broadcast_to(np.asarray(accepted_power, dtype=float), self._freqs.shape).copy()
        self.quantity = self._CheckQuantity(quantity)
        self.freq = self._SnapFreq(freq if freq is not None else self._freqs[0])
        self.sweep = bool(sweep)
        self.origin = self._Vector(origin)
        self.size = size
        self.cuts = bool(cuts)
        self._step = None
        self._anim = None
        self._polys = {}     # (freq index, settings) -> pattern polydata
        self._grid_vals = None
        self.pattern = vtkActor()
        self.pattern.SetMapper(vtkPolyDataMapper())
        m = self.pattern.GetMapper()
        m.SetLookupTable(self._lut)
        m.UseLookupTableScalarRangeOn()
        m.SetScalarModeToUsePointFieldData()
        m.SelectColorArray('value')
        self.pattern.GetProperty().SetAmbient(0.3)
        self.pattern.GetProperty().SetDiffuse(0.7)
        self.cut_actor = _Actor()
        self.cut_actor.GetProperty().SetColor(0.1, 0.1, 0.1)
        self.cut_actor.GetProperty().SetLineWidth(3.0)
        self.cut_actor.GetMapper().SetRelativeCoincidentTopologyLineOffsetParameters(0, -6)
        for a in (self.pattern, self.cut_actor):
            self.AddProp(a)

    def _CheckQuantity(self, q):
        if q not in QUANTITIES:
            raise ValueError('invalid far-field quantity "{}", use one of {}'.format(q, QUANTITIES))
        if q == 'gain' and self.accepted_power is None:
            raise ValueError('the gain needs `accepted_power` (per frequency, in W)')
        return q

    def _SnapFreq(self, f):
        return float(self._freqs[np.argmin(np.abs(self._freqs - float(f)))])

    def _Kind(self):
        return self.quantity

    def _FreqIndex(self, step=None):
        if self.sweep and step is not None:
            return int(step) % len(self._freqs)
        return int(np.argmin(np.abs(self._freqs - self.freq)))

    # ------------------------------------------------------------------
    # values
    def Values(self, fi):
        """ The quantity on the (theta, phi) grid at frequency index `fi`
        (linear: dimensionless for directivity/gain, V/m for fields). """
        res = self.result
        q = self.quantity
        if q in ('directivity', 'gain'):
            p = np.real(np.asarray(res.P_rad[fi]))
            ref = res.Prad[fi] if q == 'directivity' else self.accepted_power[fi]
            return 4 * np.pi * np.asarray(res.r) ** 2 * p / ref
        arr = dict(E_norm=res.E_norm, E_theta=res.E_theta, E_phi=res.E_phi, RHCP=res.E_cprh,
                   LHCP=res.E_cplh)[q]
        return np.abs(np.asarray(arr[fi]))

    def _Shown(self, vals):
        """ The colour value (dB in log scale) and the radius in 0..1. """
        power = Quantity.IsPower(self.quantity)
        peak = float(np.max(vals))
        if self.scale == 'log':
            top = Quantity.MaxDB(peak, power)
            shown = Quantity.ToDB(vals, power, top - self.dynamic_range)
            radius = (shown - (top - self.dynamic_range)) / self.dynamic_range
        else:
            shown = vals.astype(np.float32)
            radius = vals / peak if peak > 0 else vals
        return shown, radius

    def _Unit(self):
        q = self.quantity
        if self.scale == 'log':
            return 'dBi' if q in ('directivity', 'gain') else 'dBV/m'
        return '' if q in ('directivity', 'gain') else 'V/m'

    @staticmethod
    def _Dir(theta, phi):
        return np.stack([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)], axis=-1)

    def MainLobe(self, fi=None):
        """ ``(theta, phi, value)`` of the maximum, angles in radians. """
        fi = self._FreqIndex(self._step) if fi is None else fi
        vals = self.Values(fi)
        i, j = np.unravel_index(np.argmax(vals), vals.shape)
        return float(self.result.theta[i]), float(self.result.phi[j]), float(vals[i, j])

    # ------------------------------------------------------------------
    # geometry
    def _Radius(self):
        v = self.GetViewer()
        b = v.GetStructureBounds() if v is not None else None
        ref = max(b[1] - b[0], b[3] - b[2], b[5] - b[4]) if b is not None else 1.0
        if ref <= 0:
            ref = 1.0
        return Quantity.ParseSize(self.size, ref)

    def _UnitPattern(self, fi):
        """ Pattern points (radius up to 1) relative to the origin, colour values. """
        th, ph = np.asarray(self.result.theta), np.asarray(self.result.phi)
        TH, PH = np.meshgrid(th, ph, indexing='ij')
        vals = self.Values(fi)
        shown, radius = self._Shown(vals)
        pts = radius[..., None] * self._Dir(TH, PH)
        return pts, shown, vals

    def _ResolveOrigin(self):
        """ 'auto': from the centre of the structure and the other shown
        items along the main lobe, until the pattern clears them. """
        if not isinstance(self.origin, str):
            return
        v = self.GetViewer()
        b = None
        if v is not None:
            b = UnionBounds([v.GetStructureBounds()] + [it.GetBounds() for it in v.GetItems().values()
                                                        if it is not self and it.GetVisible()
                                                        and it.IsIncludedInBounds()])
        if b is None:
            self.origin = [0.0, 0.0, 0.0]
            return
        th, ph, _ = self.MainLobe(self._FreqIndex())
        u = self._Dir(th, ph)
        pts, _, _ = self._UnitPattern(self._FreqIndex())
        R = self._Radius()
        corners = np.array([[x, y, z] for x in b[:2] for y in b[2:4] for z in b[4:]])
        start = corners.mean(axis=0)
        reach = float(np.max((corners - start) @ u))
        back = float(np.min(pts.reshape(-1, 3) @ u)) * R
        # a clear gap, also in perspective views
        d = max(reach - back, 0.0) + 0.2 * R
        self.origin = [float(c) for c in np.round(start + u * d, 9)]

    def _Origin(self):
        return np.array([0.0, 0.0, 0.0] if isinstance(self.origin, str) else self.origin)

    def _BuildPattern(self, fi):
        key = (fi, self.quantity, self.scale, self.dynamic_range, self._Radius(), tuple(self._Origin()))
        if key in self._polys:
            return self._polys[key]
        pts, shown, vals = self._UnitPattern(fi)
        n1, n2 = shown.shape
        xyz = self._Origin() + pts * self._Radius()
        grid = vtkStructuredGrid()
        grid.SetDimensions(n1, n2, 1)
        points = vtkPoints()
        points.SetData(numpy_to_vtk(xyz.transpose(1, 0, 2).reshape(-1, 3).astype(np.float64), deep=1))
        grid.SetPoints(points)
        arr = numpy_to_vtk(np.ascontiguousarray(shown.T.ravel(), dtype=np.float32), deep=1)
        arr.SetName('value')
        grid.GetPointData().AddArray(arr)
        geom = vtkStructuredGridGeometryFilter()
        geom.SetInputData(grid)
        nrm = vtkPolyDataNormals()
        nrm.SetInputConnection(geom.GetOutputPort())
        nrm.SplittingOff()
        nrm.ConsistencyOn()
        nrm.Update()
        pd = vtkPolyData()
        pd.DeepCopy(nrm.GetOutput())
        self._polys[key] = (pd, xyz, shown)
        return self._polys[key]

    def _Cuts(self, xyz):
        th, ph = np.asarray(self.result.theta), np.asarray(self.result.phi)
        n2 = len(ph)
        pts, lines = [], []

        def add(seq):
            lines.append(list(range(len(pts), len(pts) + len(seq))))
            pts.extend(seq)
        for p0 in (0.0, np.pi / 2):
            j0, j1 = self._CutColumns(p0)
            add(list(xyz[:, j0]) + list(xyz[::-1, j1]))
        i = int(np.argmin(np.abs(th - np.pi / 2)))
        add([xyz[i, j] for j in range(n2)] + [xyz[i, 0]])
        return _Lines(pts, lines)

    def Update(self):
        if self.GetViewer() is None:
            return
        self._ResolveOrigin()
        self._RegisterAnimation()
        self._Show(self._step or 0)

    def _RegisterAnimation(self):
        want = self.sweep and len(self._freqs) > 1
        if want and self._anim is None:
            self._anim = self.AddAnimation(len(self._freqs))
        elif not want and self._anim is not None:
            self._anim.Remove()
            self._anim = None
            self._step = None

    def _Show(self, step):
        self._step = step if self.sweep else None
        fi = self._FreqIndex(self._step)
        pd, xyz, shown = self._BuildPattern(fi)
        self._grid_vals = (xyz, shown)
        self.pattern.GetMapper().SetInputData(pd)
        if self.value_range == 'auto':
            vrange = (float(np.min(shown)), float(np.max(shown)))
            if self.scale == 'linear':
                vrange = (0.0, vrange[1])
        else:
            vrange = tuple(self.value_range)
        self._ApplyLut(False, vrange)
        vis = self.GetVisible()
        self.cut_actor.GetMapper().SetInputData(self._Cuts(xyz) if self.cuts else vtkPolyData())
        self.cut_actor.SetVisibility(self.cuts and vis)
        self._bar.SetTitle('{} [{}]\nf = {}'.format(_LABELS[self.quantity], self._Unit(),
                                               Quantity.FormatFreq(self._freqs[fi])))

    def OnVisibilityChanged(self, visible):
        if visible:
            self._Show(self._step or 0)

    def OnAnimationStep(self, step):
        self._Show(step)

    def OnPrepareAnimationStep(self, step):
        self._BuildPattern(self._FreqIndex(step))

    # ------------------------------------------------------------------
    def _MaxText(self, fi):
        val = self.MainLobe(fi)[2]
        if self.quantity in ('directivity', 'gain'):
            return '{:.1f} dBi'.format(10 * np.log10(val)) if val > 0 else '0'
        return '{:.3g} V/m'.format(val)

    def GetDescription(self):
        fi = self._FreqIndex(self._step)
        return '{}, {}, max {}'.format(_LABELS[self.quantity], Quantity.FormatFreq(self._freqs[fi]),
                                       self._MaxText(fi))

    def GetInfo(self):
        fi = self._FreqIndex(self._step)
        res = self.result
        th, ph, _ = self.MainLobe(fi)
        info = {'frequency': Quantity.FormatFreq(self._freqs[fi]),
                'Dmax': '{:.2f} dBi'.format(10 * np.log10(np.atleast_1d(res.Dmax)[fi])),
                'Prad': '{:.4g} W'.format(np.atleast_1d(res.Prad)[fi]),
                'max': '{} at θ={:.1f}°, φ={:.1f}°'.format(self._MaxText(fi), np.degrees(th), np.degrees(ph))}
        for label, width in self._BeamWidths(fi).items():
            info['-3 dB beam width ({})'.format(label)] = '{:.1f}°'.format(width)
        info['origin'] = self.origin if isinstance(self.origin, str) else ', '.join('{:.4g}'.format(c)
                                                                                     for c in self.origin)
        return info

    def _CutColumns(self, p0):
        """ The phi columns nearest `p0` and `p0` + 180 degree: together one
        cut through the pattern. """
        ph = np.asarray(self.result.phi)
        return tuple(int(np.argmin(np.abs(np.angle(np.exp(1j * (ph - a)))))) for a in (p0, p0 + np.pi))

    def _BeamWidths(self, fi):
        """ -3 dB widths around the maximum in the principal planes, degrees. """
        th = np.asarray(self.result.theta)
        vals = self.Values(fi)
        out = {}
        for label, p0 in (('φ=0°', 0.0), ('φ=90°', np.pi / 2)):
            j0, j1 = self._CutColumns(p0)
            # one cut over -180..180 degree of theta
            ang = np.concatenate([-th[::-1], th[1:]])
            cut = np.concatenate([vals[::-1, j1], vals[1:, j0]])
            if not Quantity.IsPower(self.quantity):
                cut = cut ** 2
            k = int(np.argmax(cut))
            ok = cut >= cut[k] / 2
            lo = k
            while lo > 0 and ok[lo - 1]:
                lo -= 1
            hi = k
            while hi < len(cut) - 1 and ok[hi + 1]:
                hi += 1
            if lo > 0 and hi < len(cut) - 1:
                out[label] = float(np.degrees(ang[hi] - ang[lo]))
        return out

    def GetValueText(self, world_pos, prop, cell_id):
        if prop not in (self.pattern, self.cut_actor) or self._grid_vals is None:
            return None
        xyz, shown = self._grid_vals
        d = np.linalg.norm(xyz - np.asarray(world_pos), axis=-1)
        i, j = np.unravel_index(np.argmin(d), d.shape)
        return '{} = {:.3g} {} at θ={:.0f}°, φ={:.0f}°'.format(
            _SHORT[self.quantity], shown[i, j], self._Unit(), np.degrees(self.result.theta[i]),
            np.degrees(self.result.phi[j]))

    # ------------------------------------------------------------------
    HELP = dict(
        FieldItem.HELP,
        quantity='Directivity, gain (needs the accepted power) or a component of the far field.',
        freq='Frequency of the pattern.',
        sweep='Animate over all calculated frequencies.',
        size='Largest radius of the pattern: a percentage of the largest size of the shown structure,\n'
             'or a length in drawing units.',
        origin='Where the pattern is drawn, in drawing units, e.g. "0, 0, 100". "auto" moves it from the centre\n'
               'of the structure along the main lobe until it clears the structure and the items added before.',
        cuts='Draw the principal cuts (the phi = 0 and 90 degree planes and theta = 90 degree) as lines.',
    )

    def GetParameters(self):
        p = [dict(name='quantity', type='choice', value=self.quantity,
                  choices=[q for q in QUANTITIES if q != 'gain' or self.accepted_power is not None],
                  default='directivity', section='Pattern'),
             dict(name='freq', label='frequency', type='choice', value=self.freq,
                  choices=[float(f) for f in self._freqs],
                  labels=[Quantity.FormatFreq(f) for f in self._freqs], section='Pattern')]
        if len(self._freqs) > 1:
            p.append(dict(name='sweep', label='sweep frequencies', type='bool', value=self.sweep, default=False,
                          section='Pattern'))
        p += self._CommonParameters(False, scale_default='log')
        p += [
            dict(name='size', type='str', value=self.size, default='100%', section='Placement'),
            dict(name='origin', type='str', value=self.origin if isinstance(self.origin, str)
                 else list(self.origin), section='Placement'),
            dict(name='cuts', label='principal cuts', type='bool', value=self.cuts, default=False,
                 section='Extras'),
        ]
        return p

    @staticmethod
    def _Vector(value):
        if isinstance(value, str) and value.strip() == 'auto':
            return 'auto'
        vals = Quantity.ParseNumbers(value)
        if len(vals) != 3:
            raise ValueError('expected three numbers')
        return vals

    def SetParameter(self, name, value):
        if name == 'quantity':
            self.quantity = self._CheckQuantity(value)
        elif name == 'freq':
            self.freq = self._SnapFreq(value)
        elif name in ('sweep', 'cuts'):
            setattr(self, name, bool(value))
        elif name == 'size':
            Quantity.ParseSize(value, 1.0)
            self.size = value
        elif name == 'origin':
            self.origin = self._Vector(value)
        elif not self._SetCommon(name, value):
            super().SetParameter(name, value)
        if self.GetViewer() is not None:
            self.Update()
