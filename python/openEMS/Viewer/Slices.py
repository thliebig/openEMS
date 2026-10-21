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

"""openEMS.Viewer.Slices
The FieldSlices item: cutting planes through a field dump.
"""

import numpy as np

from vtkmodules.vtkCommonCore import vtkPoints, vtkLookupTable, reference
from vtkmodules.vtkCommonDataModel import vtkStructuredGrid, vtkPolyData
from vtkmodules.vtkFiltersGeometry import vtkStructuredGridGeometryFilter
from vtkmodules.vtkFiltersGeneral import vtkWarpScalar
from vtkmodules.vtkFiltersCore import vtkContourFilter, vtkGlyph3D
from vtkmodules.vtkFiltersSources import vtkArrowSource, vtkConeSource
from vtkmodules.vtkRenderingCore import vtkActor, vtkPolyDataMapper
from vtkmodules.util.numpy_support import numpy_to_vtk

from openEMS.Viewer import Quantity
from openEMS.Viewer.Item import DumpItem


def _Array(values, name, components=1):
    arr = numpy_to_vtk(np.ascontiguousarray(values, dtype=np.float32).reshape(-1, components)
                       if components > 1 else np.ascontiguousarray(values, dtype=np.float32), deep=1)
    arr.SetName(name)
    return arr


class _Plane:
    """ The VTK pipeline of one plane: grid -> warp -> colour, contours;
    sub-sampled points -> glyphs. """

    def __init__(self, item, ny, idx):
        self.ny, self.idx = ny, idx
        self.grid = vtkStructuredGrid()
        self.geom = vtkStructuredGridGeometryFilter()
        self.geom.SetInputData(self.grid)
        self.warp = vtkWarpScalar()
        self.warp.SetInputConnection(self.geom.GetOutputPort())
        self.warp.SetInputArrayToProcess(0, 0, 0, 0, 'height')
        self.color = self._Actor(self.warp, item._lut)
        self.contour = vtkContourFilter()
        self.contour.SetInputConnection(self.warp.GetOutputPort())
        self.contour.SetInputArrayToProcess(0, 0, 0, 0, 'value')
        self.lines = self._Actor(self.contour, item._lut)
        self.lines.GetProperty().SetLineWidth(2.0)
        self.lines.GetMapper().SetRelativeCoincidentTopologyLineOffsetParameters(0, -6)
        self.points = vtkPolyData()
        self.glyph = vtkGlyph3D()
        self.glyph.SetInputData(self.points)
        self.glyph.SetInputArrayToProcess(1, 0, 0, 0, 'vector')
        self.glyph.OrientOn()
        self.glyph.SetVectorModeToUseVector()
        self.glyph.SetColorModeToColorByVector()
        # an own table: the glyph mapper sets its range to the magnitudes
        self.glyph_lut = vtkLookupTable()
        self.glyphs = self._Actor(self.glyph, self.glyph_lut)
        self.sample_idx = None

    @staticmethod
    def _Actor(src, lut):
        m = vtkPolyDataMapper()
        m.SetInputConnection(src.GetOutputPort())
        m.SetLookupTable(lut)
        m.UseLookupTableScalarRangeOn()
        m.SetScalarModeToUsePointFieldData()
        m.SelectColorArray('value')
        a = vtkActor()
        a.SetMapper(m)
        a.GetProperty().SetAmbient(0.35)
        a.GetProperty().SetDiffuse(0.65)
        return a

    def Actors(self):
        return [self.color, self.lines, self.glyphs]


class FieldSlices(DumpItem):
    """ FieldSlices(source, planes=None, quantity=None, domain=None, freq=None, name=None, **settings)

    One or more cutting planes through a field dump, coloured by a field
    quantity, optionally warped, with contour lines and vector glyphs.
    Every keyword argument is a setting (see :meth:`GetParameters`), stored
    in the view state and editable in the GUI.

    Parameters
    ----------
    source : str, HDF5Dump or FieldSource
        The dump (an openEMS HDF5 field or SAR dump) or computed data.
    planes : list of (dir, pos)
        The planes: normal direction (``'x'``, ``'y'``, ``'z'``, or
        ``'rho'``, ``'alpha'``, ``'z'`` for a cylindrical dump) and position
        in drawing units (``alpha`` in radians), snapped to the nearest dump
        line; None for the middle. Default: the dump's own plane, or the
        middle z plane.
    quantity : str
        A component (``'x'``, ``'E_z'``, ``'rho'``, ...), ``'magnitude'``,
        or ``'value'`` for a scalar dump (SAR). Default: the magnitude.
    domain : 'frequency' or 'time'
        A time-domain dump can be shown at a frequency too (DFT). Default:
        the dump's own domain.
    freq : float
        Frequency in Hz (default: the first stored one).
    time_step : int
        Time-step index shown in the time domain; default the step with the
        largest value. While animated it is the current frame: the clock
        starts there (snapped to every `time_sample`-th step) and moving
        the clock changes it.
    phase : float
        Phase in degrees of the instantaneous frequency-domain value.
    amplitude : bool
        Frequency domain: show the amplitude instead of the instantaneous
        value (default for the magnitude).
    scale : 'linear' or 'log'
        Log is in dB (``10*log10`` for SAR, else ``20*log10``).
    dynamic_range : float
        dB range below the maximum shown in log scale.
    value_range : 'auto' or (vmin, vmax)
        Colour range; ``'auto'`` is the range over all planes and animation
        steps (symmetric for signed values).
    cmap : str
        A matplotlib colour map name; default by quantity.
    legend : bool
        Show the colour legend.
    color, warp, contours, glyphs : bool
        Layers: the coloured plane, its height warp along the plane normal,
        contour lines and vector arrows (the whole field for the magnitude,
        else the shown component).
    warp_scale : str or float
        Height of the largest value: ``'20%'`` of the largest plane extent,
        or a factor.
    contour_levels : int or list
        Number of levels across the colour range, or the levels.
    contour_color : str
        ``'cmap'`` or a colour.
    glyph : 'arrow' or 'cone'
    glyph_density : int
        Arrows per side of the plane.
    glyph_scale : 'magnitude' or 'uniform'
        Arrow length by value (clamped at the peak) or all the same.
    glyph_size : float
        Length factor: 1 makes the longest arrow 1.5 arrow spacings long.
    sample : int
        Spatial sub-sampling (every n-th line).
    time_sample : int
        Time-domain animation: every n-th time step.
    n_phase : int
        Frequency-domain animation: steps per period.
    animate : bool
        Register a phase (time) animation on the viewer clock.
    norm : float or complex
        The values are divided by it, e.g. the incident port voltage at the
        frequency (a complex factor also removes its phase; the time domain
        uses its magnitude) or the square root of the input power.
    unit : str
        Unit of the normalised values, e.g. ``'V/m/V'``; default the dump's.
    name : str
        Item name, default the file stem.

    Examples
    --------
    >>> v.AddItem(FieldSlices('tmp/Et.h5', planes=[('z', 1.5)], quantity='E_z', domain='frequency', freq=2.4e9))
    """

    DEFAULT_GROUP = 'results/fields'

    HELP = dict(
        DumpItem.HELP,
        planes='The planes as [direction, position] pairs.',
        color='Colour the planes by the shown value.',
        warp='Lift the planes along their normal by the value: a height map.',
        warp_scale='Height of the largest value: a percentage of the largest plane size (e.g. 20%),\n'
                   'or a factor in drawing units per value unit.',
        contours='Draw lines of equal value.',
        contour_levels='Number of contour lines across the colour range, or their values (e.g. "-1, 0, 1").',
        contour_color='"cmap" (coloured by value) or one colour for the contour lines.',
        glyphs='Draw arrows: the field vector for the magnitude, else the shown component.',
        glyph='Arrow or cone shape.',
        glyph_density='About this many arrows along each side of a plane; they sit on mesh lines,\n'
                      'at most one per line.',
        glyph_scale='Arrow length by the value (clamped at the peak), or all of the same length.',
        glyph_size='Length factor: 1 makes the longest arrow 1.5 arrow spacings long.',
    )

    def __init__(self, source, planes=None, quantity=None, domain=None, freq=None, time_step=None, phase=0.0,
                 amplitude=None, scale='linear', dynamic_range=40.0, value_range='auto', cmap=None,
                 color=True, warp=False, warp_scale='20%', contours=False, contour_levels=10,
                 contour_color='black', glyphs=False, glyph='arrow', glyph_density=30,
                 glyph_scale='magnitude', glyph_size=1.0, sample=1, time_sample=1, n_phase=36, animate=True,
                 norm=1.0, unit=None, legend=True, name=None):
        src, name = self._Source(source, name)
        super().__init__(name, src, quantity, domain, freq, time_step, phase, amplitude, sample, time_sample,
                         n_phase, animate, norm, unit, cmap, scale, dynamic_range, value_range, legend)
        self.color, self.warp, self.contours, self.glyphs = bool(color), bool(warp), bool(contours), bool(glyphs)
        self.warp_scale = warp_scale
        self.contour_levels = contour_levels
        self.contour_color = contour_color
        self.glyph = glyph
        self.glyph_density = int(glyph_density)
        self.glyph_scale = glyph_scale
        self.glyph_size = float(glyph_size)
        self._planes_in = self._DefaultPlanes() if planes is None else [tuple(p) for p in planes]
        self._planes = []        # [_Plane]
        self._plane_key = None   # what the planes were built for
        self._hidden_planes = set()   # indices of planes hidden in the tree
        self._frames = {}        # step -> [(value, height, vector) vtk arrays per plane]
        self._frame_key = None
        self._shown = None       # [(values, vectors)] of the shown frame, for hover
        self.ResolvePlanes()
        self._ResolveTimeStep()

    def _DefaultPlanes(self):
        shape = self.source.Shape()
        names = self.source.DirNames()
        for n in (2, 1, 0):
            if shape[n] == 1:
                return [(names[n], None)]
        return [(names[2], None)]

    # ------------------------------------------------------------------
    # planes
    def _Lines(self):
        return self.source.Lines(self._delta)

    def ResolvePlanes(self):
        """ Snap the planes to dump lines: ``[(dir name, position)]``. """
        lines = self._Lines()
        names = self.source.DirNames()
        out = []
        for p in self._planes_in:
            if len(p) != 2:
                raise ValueError('a plane is (dir, pos), got {!r}'.format(p))
            ny = self.source.Dir(p[0])
            l = lines[ny]
            idx = len(l) // 2 if p[1] is None else int(np.argmin(np.abs(l - float(p[1]))))
            out.append((names[ny], float(l[idx]), ny, idx))
        self._resolved = out
        return [(d, pos) for d, pos, _, _ in out]

    def GetPlanes(self):
        """ The planes as ``[(dir, pos)]``, positions snapped. """
        return [(d, pos) for d, pos, _, _ in self._resolved]

    def SetPlanes(self, planes):
        self._planes_in = [tuple(p) for p in planes]
        self.ResolvePlanes()
        self._hidden_planes = {n for n in self._hidden_planes if n < len(self._resolved)}
        self._Invalidate()

    def AddPlane(self, plane=None):
        """ Add a plane, by default a copy of the last one. """
        planes = self.GetPlanes()
        self.SetPlanes(planes + [plane if plane is not None else (planes[-1] if planes else ('z', None))])

    def RemovePlane(self, n=-1):
        """ Remove plane `n` (default the last one); the planes after it
        keep whether they are hidden. """
        planes = self.GetPlanes()
        if len(planes) <= 1:
            raise ValueError('the last plane cannot be removed')
        n = range(len(planes))[n]
        del planes[n]
        self._hidden_planes = {k - (k > n) for k in self._hidden_planes if k != n}
        self.SetPlanes(planes)

    # ------------------------------------------------------------------
    # data
    DATA_SETTINGS = DumpItem.DATA_SETTINGS + ('planes',)

    def _Invalidate(self):
        super()._Invalidate()
        self._frames = {}

    def _Raw(self, pi, step=0):
        """ The raw data of plane `pi`: complex at the frequency, or real at
        the time index of animation `step`. """
        _, _, ny, idx = self._resolved[pi]
        src = self.source
        if self.domain == 'frequency':
            key = ('f', ny, idx, self.sample, self.freq)
            if key not in self._raw:
                self._raw[key] = self._Wrap(pi, self._Norm(src.ReadAtFrequency(self.freq, ny, idx, self.sample)))
        else:
            t = self._TimeIndex(step)
            key = ('t', ny, idx, self.sample, t)
            if key not in self._raw:
                self._raw[key] = self._Wrap(pi, self._Norm(np.asarray(src.Read(ny, idx, self.sample, t_idx=t),
                                                                      dtype=np.float32)))
        data = self._raw[key]
        if key not in self._peaks:
            self._peaks[key] = self._Peaks(pi, data)
        return data

    def _Peaks(self, pi, data):
        """ (peak of the quantity, peak of the arrow length) of `data`. """
        src = self.source
        peak = Quantity.PeakAbs(data, self.quantity, src.mesh_type, self._Alpha(pi))
        # the arrows show the quantity as a vector
        return peak, (peak if src.is_vector else None)

    def _TimePeaks(self):
        """ The peaks of every animation step and of the static one, read
        once without keeping the data: the colour range and the warp height
        then stay the same over all steps. """
        src = self.source
        out = {}
        for pi, (_, _, ny, idx) in enumerate(self._resolved):
            for t in self._TimeSteps():
                key = ('t', ny, idx, self.sample, t)
                if key not in self._peaks:
                    data = np.asarray(src.Read(ny, idx, self.sample, t_idx=t), dtype=np.float32)
                    self._peaks[key] = self._Peaks(pi, self._Wrap(pi, self._Norm(data)))
                out[t] = max(out.get(t, 0.0), self._peaks[key][0])
        return out

    def _Alpha(self, pi):
        """ The angle at each point of plane `pi` of a cylindrical dump. """
        if self.source.mesh_type != 'cylindrical':
            return None
        return self._PlaneCoords(pi)[1]

    def _WrapAxis(self, pi):
        """ The in-plane axis (0/1) of a cylindrical plane whose alpha lines
        go almost around: the plane is closed by repeating the first line. """
        _, _, ny, _ = self._resolved[pi]
        if self.source.mesh_type != 'cylindrical' or ny == 1:
            return None
        al = self.source.SampledLines(ny, self.sample)[0 if ny == 0 else 1]
        if len(al) < 3:
            return None
        span = al[-1] - al[0]
        step = span / (len(al) - 1)
        if span < 2 * np.pi - 1e-9 and span + 1.5 * step >= 2 * np.pi:
            return 0 if ny == 0 else 1
        return None

    def _Wrap(self, pi, data):
        """ `data` (``(..., N1, N2)``) closed around, see :meth:`_WrapAxis`. """
        ax = self._WrapAxis(pi)
        if ax is None:
            return data
        ax = data.ndim - 2 + ax
        first = np.take(data, [0], axis=ax)
        return np.concatenate([data, first], axis=ax)

    def _PlaneCoords(self, pi):
        """ The three mesh coordinates of every point of plane `pi`, each
        ``(N1, N2)``. """
        _, _, ny, idx = self._resolved[pi]
        lines = self._Lines()
        sl = self.source.SampledLines(ny, self.sample, self._delta)
        ax = self._WrapAxis(pi)
        if ax is not None:
            sl[ax] = np.append(sl[ax], sl[ax][0] + 2 * np.pi)
        u, v = [n for n in range(3) if n != ny]
        U, V = np.meshgrid(sl[0], sl[1], indexing='ij')
        c = [None] * 3
        c[u], c[v] = U, V
        c[ny] = np.full(U.shape, lines[ny][idx])
        return c

    def _Points(self, pi):
        """ Cartesian points (N1, N2, 3) of plane `pi` and their normals. """
        c = self._PlaneCoords(pi)
        _, _, ny, _ = self._resolved[pi]
        if self.source.mesh_type == 'cylindrical':
            rho, a, z = c
            pts = np.stack([rho * np.cos(a), rho * np.sin(a), z], axis=-1)
            if ny == 0:
                nrm = np.stack([np.cos(a), np.sin(a), np.zeros_like(a)], axis=-1)
            elif ny == 1:
                nrm = np.stack([-np.sin(a), np.cos(a), np.zeros_like(a)], axis=-1)
            else:
                nrm = np.stack([np.zeros_like(a), np.zeros_like(a), np.ones_like(a)], axis=-1)
        else:
            pts = np.stack(c, axis=-1)
            nrm = np.zeros_like(pts)
            nrm[..., ny] = 1.0
        return pts, nrm

    def _GlobalPeaks(self):
        """ Peak of the quantity and of the arrow length over all planes
        and all animation steps. """
        if self.domain == 'time':
            self._TimePeaks()
        keys = [k for k in self._peaks if (k[0] == 'f') == (self.domain == 'frequency')]
        if self.domain == 'frequency':
            keys = [k for k in keys if k[3] == self.sample and k[4] == self.freq]
        q = max([self._peaks[k][0] for k in keys], default=0.0)
        v = max([self._peaks[k][1] or 0.0 for k in keys], default=0.0)
        return q, v

    def _PhaseAnimated(self):
        # arrows show the instantaneous vector: they turn with the phase
        return not self.amplitude or self.glyphs

    def _Signed(self):
        return self.quantity not in ('magnitude', 'value') and not (self.amplitude and self.domain == 'frequency') \
            and self.scale == 'linear'

    # ------------------------------------------------------------------
    # frames
    def _FrameKey(self):
        # dB values depend on the (growing) global peak
        peak = self._GlobalPeaks()[0] if self.scale == 'log' else None
        return (peak, self.quantity, self.amplitude, self.scale, self.dynamic_range, self.phase, self.n_phase,
                self.domain, self.freq, self.time_step, self.time_sample, self.sample, self.glyphs,
                self.warp, tuple(self.GetPlanes()))

    def _Frame(self, step):
        """ VTK arrays (value, height, vector) per plane for `step`. """
        if self._frame_key != self._FrameKey():
            self._frames = {}
            self._frame_key = self._FrameKey()
        if step in self._frames:
            return self._frames[step]
        src = self.source
        fd = self.domain == 'frequency'
        phase = self._Phase(step) if fd else 0.0
        raws = [self._Raw(pi, step) for pi in range(len(self._resolved))]
        peak, _ = self._GlobalPeaks()
        power = Quantity.IsPower(src.kind)
        floor = Quantity.MaxDB(peak, power) - self.dynamic_range
        out = []
        for pi, data in enumerate(raws):
            alpha = self._Alpha(pi)
            vals = Quantity.Evaluate(data, self.quantity, src.mesh_type, alpha, phase, self.amplitude and fd)
            if self.scale == 'log':
                vals = Quantity.ToDB(vals, power, floor)
                height = vals - floor
            else:
                height = vals
            vec = None
            if src.is_vector and self.glyphs:
                inst = Quantity.Instantaneous(data, phase)
                vec = Quantity.ArrowVector(inst, self.quantity, src.mesh_type, alpha)
                vec = np.moveaxis(vec, 0, -1).astype(np.float32)
            # VTK point order: first in-plane direction fastest
            arrays = (_Array(vals.ravel(order='F'), 'value'), _Array(height.ravel(order='F'), 'height'),
                      None if vec is None else vec.transpose(1, 0, 2).reshape(-1, 3))
            out.append(arrays)
        # the shown values for hover come from the frame
        self._frames[step] = out
        return out

    # ------------------------------------------------------------------
    # build and update
    def Update(self):
        """ Rebuild the pipelines and show the current settings. """
        if self.GetViewer() is None:
            return
        self.ResolvePlanes()
        self._BuildPlanes()
        self._RegisterAnimation()
        self._ShowStep(self._step)
        self.PartsChanged()   # e.g. a plane added from the context menu

    def _BuildPlanes(self):
        # the arrow positions depend on the density
        key = [(ny, idx, self.sample, self.glyph_density) for _, _, ny, idx in self._resolved]
        if self._plane_key == key and self._planes:
            return
        for pl in self._planes:
            for a in pl.Actors():
                self.RemoveProp(a)
        self._planes = []
        for pi, (_, _, ny, idx) in enumerate(self._resolved):
            pl = _Plane(self, ny, idx)
            pts, nrm = self._Points(pi)
            n1, n2 = pts.shape[:2]
            dims = [1, 1, 1]
            dims[0], dims[1] = n1, n2
            pl.grid.SetDimensions(*dims)
            points = vtkPoints()
            points.SetData(numpy_to_vtk(pts.transpose(1, 0, 2).reshape(-1, 3).astype(np.float64), deep=1))
            pl.grid.SetPoints(points)
            pl.grid.GetPointData().SetNormals(_Array(nrm.transpose(1, 0, 2).reshape(-1, 3), 'normals', 3))
            # glyph positions: a sub-grid with about glyph_density points per side
            s1 = max(1, int(round(n1 / max(self.glyph_density, 1))))
            s2 = max(1, int(round(n2 / max(self.glyph_density, 1))))
            I, J = np.meshgrid(np.arange(0, n1, s1), np.arange(0, n2, s2), indexing='ij')
            pl.sample_idx = (I.ravel(order='F') + n1 * J.ravel(order='F'))
            gp = vtkPoints()
            gp.SetData(numpy_to_vtk(pts.transpose(1, 0, 2).reshape(-1, 3)[pl.sample_idx].astype(np.float64), deep=1))
            pl.points.SetPoints(gp)
            ext = np.ptp(pts.reshape(-1, 3), axis=0)
            pl.extent = float(max(ext.max(), 1e-12))
            pl.spacing = pl.extent / max(self.glyph_density, 1)
            self._planes.append(pl)
            for a in pl.Actors():
                self.AddProp(a)
        self._plane_key = key

    def _ShowStep(self, step):
        self._step = step
        frame = self._Frame(step)
        peak, vpeak = self._GlobalPeaks()
        signed = self._Signed()
        vrange = self._Range(peak, signed)
        self._ApplyLut(signed, vrange)
        self._shown = []
        for pi, (pl, (val, height, vec)) in enumerate(zip(self._planes, frame)):
            shown = self.GetVisible() and pi not in self._hidden_planes
            pd = pl.grid.GetPointData()
            pd.AddArray(val)
            pd.SetActiveScalars('value')
            pd.AddArray(height)
            pl.grid.Modified()
            # warp: the largest value gets warp_scale of the plane extent
            top = max(abs(vrange[0]), abs(vrange[1])) if self.scale == 'linear' else self.dynamic_range
            ws = self.warp_scale
            if not self.warp:
                factor = 0.0
            elif isinstance(ws, str) and ws.strip().endswith('%'):
                factor = Quantity.ParseSize(ws, pl.extent) / top if top > 0 else 0.0
            else:
                factor = float(ws)
            pl.warp.SetScaleFactor(factor)
            pl.color.SetVisibility(self.color and shown)
            # contours
            pl.contour.SetNumberOfContours(0)
            levels = self.contour_levels
            if isinstance(levels, (int, np.integer)) or (isinstance(levels, str) and levels.strip().isdigit()):
                n = int(levels)
                levels = np.linspace(vrange[0], vrange[1], n + 2)[1:-1]
            for i, lv in enumerate(levels):
                pl.contour.SetValue(i, float(lv))
            lm = pl.lines.GetMapper()
            if self.contour_color == 'cmap':
                lm.ScalarVisibilityOn()
            else:
                from CSXCAD.Viewer.ViewState import ParseColor
                lm.ScalarVisibilityOff()
                pl.lines.GetProperty().SetColor(*ParseColor(self.contour_color)[0])
            pl.lines.SetVisibility(self.contours and shown)
            # glyphs
            show_glyphs = self.glyphs and vec is not None
            if show_glyphs:
                sub = vec[pl.sample_idx]
                pl.points.GetPointData().SetVectors(_Array(sub, 'vector', 3))
                pl.points.GetPointData().SetScalars(_Array(np.linalg.norm(sub, axis=1), 'value'))
                pl.points.Modified()
                src = vtkConeSource() if self.glyph == 'cone' else vtkArrowSource()
                if self.glyph == 'cone':
                    src.SetResolution(12)
                pl.glyph.SetSourceConnection(src.GetOutputPort())
                if self.glyph_scale == 'uniform':
                    pl.glyph.SetScaleModeToDataScalingOff()
                else:
                    pl.glyph.SetScaleModeToScaleByVector()
                    pl.glyph.ClampingOn()
                    pl.glyph.SetRange(0.0, vpeak if vpeak > 0 else 1.0)
                pl.glyph.SetScaleFactor(1.5 * pl.spacing * self.glyph_size)
                pl.glyph_lut.DeepCopy(self._lut)
                gm = pl.glyphs.GetMapper()
                gm.UseLookupTableScalarRangeOff()
                gm.SetScalarModeToUsePointData()
                gm.SetScalarRange(0.0, vpeak if vpeak > 0 else 1.0)
            pl.glyphs.SetVisibility(show_glyphs and shown)
            self._shown.append((val, vec))
        self._bar.SetTitle(self._LegendTitle())

    def OnPrepareAnimationStep(self, step):
        for pi in range(len(self._resolved)):
            self._Raw(pi, step)

    # ------------------------------------------------------------------
    # texts
    def _QuantityLabel(self):
        fd = self.domain == 'frequency'
        return Quantity.Label(self.source.label if self.quantity == 'value' else self.source.kind,
                              self.quantity, (self.amplitude and fd) or self.scale == 'log')

    def _LegendTitle(self):
        # 'f = ...' so the frequency is not read as a value of the bar;
        # VTK's font has no Greek letters
        when = self._WhenText(phase='phase ').replace(', ', '\n')
        if self.domain == 'frequency':
            when = 'f = ' + when
        return '{} [{}]\n{}'.format(self._QuantityLabel(), self._Unit(), when.replace('t=', 't = '))

    def GetDescription(self):
        parts = [self._QuantityLabel(), self._FixedWhenText()]
        if self.scale == 'log':
            parts.append('log')
        return ', '.join(p for p in parts if p)

    def GetValueText(self, world_pos, prop, cell_id):
        for pl, shown in zip(self._planes, self._shown or []):
            if prop not in pl.Actors():
                continue
            pd = pl.warp.GetOutput()
            if prop is pl.color and 0 <= cell_id < pd.GetNumberOfCells():
                cell = pd.GetCell(cell_id)
                n = cell.GetNumberOfPoints()
                closest, pc, w = [0.0] * 3, [0.0] * 3, [0.0] * n
                sub, d2 = reference(0), reference(0.0)
                cell.EvaluatePosition(list(world_pos), closest, sub, pc, d2, w)
                ids = [cell.GetPointId(i) for i in range(n)]
            else:
                # on a contour line or an arrow: the nearest plane point
                pid = pd.FindPoint(list(world_pos))
                if pid < 0:
                    return None
                n, w, ids = 1, [1.0], [pid]
            vals = pd.GetPointData().GetArray('value')
            val = sum(w[i] * vals.GetValue(ids[i]) for i in range(n))
            text = '{} = {:.4g} {} at {}'.format(self._QuantityLabel(), val, self._Unit(),
                                                 self._LengthText(world_pos))
            vec = shown[1]
            # the arrows of a component repeat the value
            if vec is not None and self.quantity == 'magnitude':
                vv = sum(w[i] * vec[ids[i]] for i in range(n))
                text += '   {} = ({})'.format(self.source.kind, ', '.join('{:.3g}'.format(c) for c in vv))
            return text
        return None

    def GetInfo(self):
        src = self.source
        info = dict(source=src.path or src.name, field=src.label, mesh=src.mesh_type,
                    domain=self.domain, shown=self._FixedWhenText(),
                    planes=', '.join('{}={:.4g}'.format(d, p) for d, p in self.GetPlanes()))
        if not info['shown']:
            del info['shown']
        peak, _ = self._GlobalPeaks()
        signed = self._Signed()
        info['value range'] = '{:.4g} … {:.4g} {}'.format(*self._Range(peak, signed), self._Unit())
        if self.domain == 'time':
            pts = sum(np.prod([len(l) for l in self.source.SampledLines(ny, self.sample)])
                      for _, _, ny, _ in self._resolved)
            per_step = pts * (3 if src.is_vector else 1) * 4
            n = -(-self._NumTimes() // self.time_sample)
            used = sum(a.nbytes for k, a in self._raw.items() if k[0] == 't')
            info['memory'] = '{:.1f} MB of {:.1f} MB for {} steps'.format(used / 1e6, per_step * n / 1e6, n)
        return info

    # ------------------------------------------------------------------
    # settings
    def GetParameters(self):
        src = self.source
        lines = self._Lines()
        p = self._FieldParameters()
        p.append(dict(name='planes', type='str', value=[list(x) for x in self.GetPlanes()], form=False))
        for n, (d, pos) in enumerate(self.GetPlanes()):
            ny = src.Dir(d)
            p.append(dict(name='plane{}_dir'.format(n), label='plane {} normal'.format(n + 1), type='choice',
                          value=d, choices=list(src.DirNames()), state=False, section='Planes',
                          help='Normal direction of this plane (more planes: the item\'s context menu).'))
            p.append(dict(name='plane{}_pos'.format(n), label='plane {} position'.format(n + 1), type='slider',
                          value=pos, values=[float(x) for x in lines[ny]], state=False, section='Planes',
                          help='Position of this plane, on the mesh lines of the dump.'))
        p += self._CommonParameters(self._Signed())
        p += [
            dict(name='color', label='coloured plane', type='bool', value=self.color, default=True,
                 section='Layers'),
            dict(name='warp', type='bool', value=self.warp, default=False, section='Layers'),
            dict(name='warp_scale', label='warp height', type='str', value=self.warp_scale, default='20%',
                 section='Layers'),
            dict(name='contours', type='bool', value=self.contours, default=False, section='Layers'),
            dict(name='contour_levels', label='contour levels', type='int' if isinstance(
                self.contour_levels, (int, np.integer)) else 'str', value=self.contour_levels, min=1, max=100,
                 default=10, section='Layers'),
            dict(name='contour_color', label='contour colour', type='str', value=self.contour_color,
                 default='black', section='Layers'),
        ]
        if src.is_vector:
            p += [
                dict(name='glyphs', label='arrows', type='bool', value=self.glyphs, default=False,
                     section='Layers'),
                dict(name='glyph', label='arrow shape', type='choice', value=self.glyph,
                     choices=['arrow', 'cone'], default='arrow', section='Layers'),
                dict(name='glyph_density', label='arrows per side', type='int', value=self.glyph_density,
                     min=2, max=200, default=30, section='Layers'),
                dict(name='glyph_scale', label='arrow length', type='choice', value=self.glyph_scale,
                     choices=['magnitude', 'uniform'], default='magnitude', section='Layers'),
                dict(name='glyph_size', label='arrow size', type='float', value=self.glyph_size, min=0.05,
                     max=50.0, step=0.25, decimals=2, default=1.0, section='Layers'),
            ]
        return p + self._DataParameters() + self._AnimationParameters()

    def _SetOwn(self, name, value):
        if name.startswith('plane') and name.endswith(('_dir', '_pos')):
            n = int(name[5:-4])
            planes = self.GetPlanes()
            d, pos = planes[n]
            if name.endswith('_dir'):
                planes[n] = (value, None)
            else:
                planes[n] = (d, float(value))
            self.SetPlanes(planes)
        elif name == 'planes':
            self.SetPlanes(value)
        elif name == 'glyph_density':
            self.glyph_density = max(int(value), 1)
        elif name == 'glyph_size':
            self.glyph_size = float(value)
        elif name in ('color', 'warp', 'contours', 'glyphs'):
            setattr(self, name, bool(value))
        elif name in ('warp_scale', 'contour_color', 'glyph', 'glyph_scale'):
            setattr(self, name, value)
        elif name == 'contour_levels':
            if isinstance(value, str):
                vals = Quantity.ParseNumbers(value)
                value = int(vals[0]) if len(vals) == 1 and vals[0].is_integer() else vals
            self.contour_levels = value
        else:
            return False
        return True

    # parts: the planes, each shown or hidden in the tree
    def GetParts(self):
        planes = self.GetPlanes()
        if len(planes) < 2:
            return []
        return [('plane {}'.format(n + 1), 'plane {}: {} = {:.4g}'.format(n + 1, d, pos))
                for n, (d, pos) in enumerate(planes)]

    def _PartIndex(self, part):
        n = int(part.split()[-1]) - 1 if part.startswith('plane ') and part.split()[-1].isdigit() else -1
        if not 0 <= n < len(self.GetPlanes()):
            raise KeyError('no part "{}"'.format(part))
        return n

    def SetPartVisible(self, part, val):
        n = self._PartIndex(part)
        if val:
            self._hidden_planes.discard(n)
        else:
            self._hidden_planes.add(n)
        if self.GetViewer() is not None and self._planes:
            self._ShowStep(self._step)

    def GetPartVisible(self, part):
        return self._PartIndex(part) not in self._hidden_planes

    def GetActions(self):
        return [('Add plane', lambda: (self.AddPlane(), self.Update(), self.Render())),
                ('Remove plane', lambda: (self.RemovePlane(), self.Update(), self.Render()))]
