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
Results of a simulation: what to compute from it and how to plot it.

A :class:`Results` object knows the simulation (the openEMS object with its
ports and NF2FF boxes, and the simulation path) and hands out *quantities*:
small objects naming one result, e.g. ``res.S(2, 1)`` or
``res.FarFieldCut('nf2ff', f0, phi=0)``. A quantity computes its data on
first use (``GetXY()``, ``GetMap()``, ``GetValue()``) and caches it.
Figures collect quantities as curves or maps; :meth:`Results.Show` draws
them with matplotlib. The whole description is stored in the setup file
next to the XML (:meth:`Results.Write`), so that it can be shown again
without the script.

>>> res = Results(FDTD, Sim_Path)
>>> fig = res.AddFigure('S-Parameter', ylabel='S (dB)')
>>> fig.AddCurve(res.S(1, 1), fmt='dB')
>>> fig.AddCurve(res.S(2, 1), fmt='dB')
>>> f, s21 = res.S(2, 1).GetXY()
>>> res.Show()
"""

import hashlib
import json
import os

import numpy as np

from openEMS import setup_file
from openEMS.physical_constants import Z0


# ----------------------------------------------------------------------
# formats

FORMATS = ('abs', 'amplitude', 'dB', 'dBi', 'real', 'imag', 'phase', 'phase_unwrap')


def Format(values, fmt, power=False):
    """ Format(values, fmt, power=False)

    Apply the curve format `fmt` (see FORMATS) to `values`. `power` selects
    10 log10 instead of 20 log10 for the dB formats; ``'amplitude'`` is the
    square root of a power (e.g. the field pattern from the directivity),
    else the magnitude. """
    values = np.asarray(values)
    if fmt is None:
        return np.real(values) if np.isrealobj(values) or np.all(np.imag(values) == 0) else np.abs(values)
    if fmt == 'abs':
        return np.abs(values)
    if fmt == 'amplitude':
        return np.sqrt(np.abs(values)) if power else np.abs(values)
    if fmt in ('dB', 'dBi'):
        with np.errstate(divide='ignore'):
            return (10 if power else 20) * np.log10(np.abs(values))
    if fmt == 'real':
        return np.real(values)
    if fmt == 'imag':
        return np.imag(values)
    if fmt == 'phase':
        return np.rad2deg(np.angle(values))
    if fmt == 'phase_unwrap':
        return np.rad2deg(np.unwrap(np.angle(values)))
    raise ValueError('unknown format "{}", use one of {}'.format(fmt, FORMATS))


def _EngUnit(values, unit):
    """ Scale factor and prefixed unit for an axis of `values` in `unit`. """
    peak = np.nanmax(np.abs(values)) if np.size(values) else 0
    for fac, pre in ((1e12, 'T'), (1e9, 'G'), (1e6, 'M'), (1e3, 'k'), (1, ''), (1e-3, 'm'),
                     (1e-6, 'µ'), (1e-9, 'n'), (1e-12, 'p')):
        if peak >= fac:
            return fac, pre + unit
    return 1, unit


# ----------------------------------------------------------------------
# quantities

QUANTITIES = {}


def _register(cls):
    QUANTITIES[cls.NAME] = cls
    return cls


class Quantity:
    """ One result of a simulation, see the module documentation. Created
    through the methods of :class:`Results`. """
    NAME = None
    KIND = 'curve'      # 'curve', 'map' or 'value'
    POWER = False       # dB as 10 log10
    UNIT = ''
    YNAME = ''          # for the axis label, e.g. 'S-parameter'
    X = 'f'             # 'f', 't' or 'angle' (curves)
    DEFAULT_FMT = None

    def __init__(self, res, **params):
        self.res = res
        self.params = params

    def __repr__(self):
        return '{}({})'.format(self.NAME, ', '.join('{}={!r}'.format(k, v) for k, v in self.params.items()))

    def ToDict(self):
        return dict(q=self.NAME, **setup_file._to_json(self.params, self.NAME))

    def Key(self):
        return json.dumps(self.ToDict(), sort_keys=True)

    def Label(self):
        return self.NAME

    def XName(self):
        """ Name of the x axis of a curve. """
        if self.X == 'angle':
            if 'phi' in self.params:
                return 'theta' if self.params.get('phi') is not None else 'phi'
            return 'angle'
        return {'f': 'frequency', 't': 'time'}.get(self.X, self.X or '')

    def _Cached(self):
        key = self.Key()
        if key not in self.res._cache:
            self.res._cache[key] = self._Compute()
        return self.res._cache[key]

    def GetXY(self):
        """ GetXY()

        The data of a curve: x (frequency in Hz, time in s or angle in
        degree) and the (possibly complex) values. """
        if self.KIND != 'curve':
            raise TypeError('{} is a {}, not a curve'.format(self.NAME, self.KIND))
        return self._Cached()

    def GetMap(self):
        """ GetMap()

        The data of a map: the lines x and y (drawing units) and the values
        of shape (len(x), len(y)). """
        if self.KIND != 'map':
            raise TypeError('{} is a {}, not a map'.format(self.NAME, self.KIND))
        return self._Cached()

    def GetValue(self):
        """ GetValue()

        The value of a single-valued quantity. """
        if self.KIND != 'value':
            raise TypeError('{} is a {}, not a value'.format(self.NAME, self.KIND))
        return self._Cached()

    def At(self, freq):
        """ The value at frequency `freq`: a value as is, a curve over
        frequency interpolated (complex values component-wise). """
        if self.KIND == 'value':
            return self.GetValue()
        if self.KIND == 'curve' and self.X == 'f':
            x, y = self.GetXY()
            y = np.asarray(y)
            if np.iscomplexobj(y):
                return np.interp(freq, x, y.real) + 1j * np.interp(freq, x, y.imag)
            return np.interp(freq, x, y)
        raise TypeError('{} has no value at a frequency'.format(self.NAME))

    def _Compute(self):
        raise NotImplementedError


class _PortQuantity(Quantity):
    LABEL = None   # the label before the port number, default NAME

    def Label(self):
        return '{}{}'.format(self.LABEL or self.NAME, self.params['port'])

    def _p(self, nr):
        return self.res._Port(nr)


@_register
class S(_PortQuantity):
    NAME = 'S'
    YNAME = 'S-parameter'
    DEFAULT_FMT = 'dB'

    def Label(self):
        return 'S{}{}'.format(self.params['i'], self.params['j'])

    def _Compute(self):
        pi, pj = self._p(self.params['i']), self._p(self.params['j'])
        return self.res.f, pi.uf_ref / pj.uf_inc


@_register
class Zin(_PortQuantity):
    NAME = 'Zin'
    YNAME = 'impedance'
    UNIT = 'Ohm'

    def _Compute(self):
        p = self._p(self.params['port'])
        return self.res.f, p.uf_tot / p.if_tot


@_register
class Yin(_PortQuantity):
    NAME = 'Yin'
    UNIT = 'S'
    YNAME = 'admittance'

    def _Compute(self):
        p = self._p(self.params['port'])
        return self.res.f, p.if_tot / p.uf_tot


class _Power(_PortQuantity):
    POWER = True
    YNAME = 'power'
    UNIT = 'W'
    ATTR = None

    def _Compute(self):
        return self.res.f, getattr(self._p(self.params['port']), self.ATTR)


@_register
class P_inc(_Power):
    NAME, ATTR = 'P_inc', 'P_inc'


@_register
class P_ref(_Power):
    NAME, ATTR = 'P_ref', 'P_ref'


@_register
class P_acc(_Power):
    NAME, ATTR = 'P_acc', 'P_acc'


class _UI(_PortQuantity):
    LETTER = None

    def __init__(self, res, **params):
        super().__init__(res, **params)
        if params.get('part', 'tot') not in ('tot', 'inc', 'ref'):
            raise ValueError('part must be "tot", "inc" or "ref"')
        if params.get('domain', 'f') not in ('f', 't'):
            raise ValueError('domain must be "f" or "t"')
        self.X = params.get('domain', 'f')

    def Label(self):
        return '{}{}_{}'.format(self.LETTER, self.params['port'], self.params.get('part', 'tot'))

    def _Compute(self):
        p = self._p(self.params['port'])
        part = self.params.get('part', 'tot')
        if self.X == 'f':
            return self.res.f, getattr(p, '{}f_{}'.format(self.LETTER.lower(), part))
        name = '{}t_{}'.format(self.LETTER.lower(), part)
        if not hasattr(p, name):
            raise ValueError('port {}: no time-domain {} (needs a real reference impedance)'.format(
                self.params['port'], part))
        t = p.u_time if self.LETTER == 'U' else p.i_time
        return t, getattr(p, name)


@_register
class U(_UI):
    NAME, LETTER, UNIT, YNAME = 'U', 'U', 'V', 'voltage'


@_register
class I(_UI):
    NAME, LETTER, UNIT, YNAME = 'I', 'I', 'A', 'current'


@_register
class Beta(_PortQuantity):
    NAME = 'Beta'
    LABEL = 'beta'
    YNAME = 'propagation constant'
    UNIT = '1/m'

    def _Compute(self):
        p = self._p(self.params['port'])
        if not hasattr(p, 'beta'):
            raise ValueError('port {} has no propagation constant'.format(self.params['port']))
        return self.res.f, p.beta * np.ones_like(self.res.f)


@_register
class ZL(_PortQuantity):
    NAME = 'ZL'
    YNAME = 'line impedance'
    UNIT = 'Ohm'

    def _Compute(self):
        p = self._p(self.params['port'])
        return self.res.f, getattr(p, 'ZL', p.Z_ref) * np.ones_like(self.res.f)


@_register
class Probe(Quantity):
    NAME = 'Probe'
    YNAME = 'probe'

    def __init__(self, res, **params):
        super().__init__(res, **params)
        self.X = params.get('domain', 'f')

    def Label(self):
        return self.params['name']

    def _Compute(self):
        from openEMS.ports import UI_data
        d = UI_data([self.params['name']], self.res.sim_path, self.res.f)
        if self.X == 'f':
            return self.res.f, np.asarray(d.ui_f_val[0])
        return d.ui_time[0], d.ui_val[0]


def _Angles(angles, default=(-180, 180, 2)):
    if angles is None:
        return np.arange(default[0], default[1] + default[2] / 2, default[2])
    return np.asarray(angles, dtype=float)


@_register
class FarFieldCut(Quantity):
    """ Directivity along a cut with a fixed `phi` (over theta) or a fixed
    `theta` (over phi). `component`: ``'norm'`` (total), ``'theta'``,
    ``'phi'``, ``'cprh'``, ``'cplh'``. """
    NAME = 'FarFieldCut'
    YNAME = 'directivity'
    X = 'angle'
    POWER = True
    DEFAULT_FMT = 'dBi'

    def __init__(self, res, **params):
        super().__init__(res, **params)
        if (params.get('phi') is None) == (params.get('theta') is None):
            raise ValueError('FarFieldCut: give either phi or theta')

    def Label(self):
        fixed = 'phi' if self.params.get('phi') is not None else 'theta'
        return 'D ({}={:g}°, {})'.format(fixed, self.params[fixed], _FormatFreq(self.params['freq']))

    def _Compute(self):
        p = self.params
        angles = _Angles(p.get('angles'))
        if p.get('phi') is not None:
            theta, phi = angles, [p['phi']]
        else:
            theta, phi = [p['theta']], angles
        r = self.res._NF2FF(p['box'], [p['freq']], theta, phi, p.get('center'))
        comp = p.get('component') or 'norm'
        field = {'norm': r.E_norm, 'theta': r.E_theta, 'phi': r.E_phi,
                 'cprh': r.E_cprh, 'cplh': r.E_cplh}[comp][0]
        D = r.Dmax[0] * (np.abs(field) / np.max(r.E_norm[0])) ** 2
        return angles, np.squeeze(D)


class _Sphere(Quantity):
    """ A result of the far field over the given angles (default the full
    sphere in 5° steps): a value for one frequency, a curve for several. """

    def __init__(self, res, **params):
        super().__init__(res, **params)
        self.KIND = 'curve' if np.ndim(params['freq']) else 'value'

    def _Compute(self):
        p = self.params
        theta = np.arange(0, 181, 5.0) if p.get('theta') is None else np.asarray(p['theta'], dtype=float)
        phi = np.arange(0, 360, 5.0) if p.get('phi') is None else np.asarray(p['phi'], dtype=float)
        freq = np.atleast_1d(p['freq'])
        vals = np.asarray(self._Of(self.res._NF2FF(p['box'], freq, theta, phi, p.get('center')), freq))
        return float(vals[0]) if self.KIND == 'value' else (freq, vals)

    def _Of(self, r, freq):
        """ The values per frequency of the NF2FF result `r`. """
        raise NotImplementedError


@_register
class Dmax(_Sphere):
    """ Maximum directivity, see :class:`_Sphere`. """
    NAME = 'Dmax'
    YNAME = 'directivity'
    POWER = True
    DEFAULT_FMT = 'dBi'

    def _Of(self, r, freq):
        return r.Dmax


@_register
class Prad(_Sphere):
    """ Radiated power (W), see :class:`_Sphere`. """
    NAME = 'Prad'
    POWER = True
    UNIT = 'W'
    YNAME = 'radiated power'

    def _Of(self, r, freq):
        return r.Prad


@_register
class Efficiency(_Sphere):
    """ Radiation efficiency Prad / P_acc of `port`, see :class:`_Sphere`. """
    NAME = 'Efficiency'
    YNAME = 'efficiency'

    def _Of(self, r, freq):
        return np.asarray(r.Prad) / np.real(self.res.P_acc(self.params['port']).At(freq))


@_register
class RCS(Quantity):
    """ Radar cross section (m²) of a plane-wave excitation of amplitude
    `E_amp` (V/m), whose incident field is recorded by the probe `probe`
    (``'et'``, the excitation signal): along a cut with fixed `phi` or
    `theta` at one frequency (curve over angle), or in the direction
    (`theta`, `phi`), e.g. the backscatter, over frequency (curve; a value
    for one frequency). """
    NAME = 'RCS'
    YNAME = 'RCS'
    POWER = True
    UNIT = 'm²'

    def __init__(self, res, **params):
        super().__init__(res, **params)
        if params.get('phi') is None and params.get('theta') is None:
            raise ValueError('RCS: give phi, theta or both')
        self.direction = params.get('phi') is not None and params.get('theta') is not None
        if self.direction:
            self.KIND = 'curve' if np.ndim(params['freq']) else 'value'
        else:
            if np.ndim(params['freq']):
                raise ValueError('RCS: a cut over angle needs a single frequency')
            self.X = 'angle'

    def Label(self):
        if self.direction:
            return 'RCS (theta={:g}°, phi={:g}°)'.format(self.params['theta'], self.params['phi'])
        fixed = 'phi' if self.params.get('phi') is not None else 'theta'
        return 'RCS ({}={:g}°, {})'.format(fixed, self.params[fixed], _FormatFreq(self.params['freq']))

    def _Compute(self):
        from openEMS.ports import UI_data
        p = self.params
        freq = np.atleast_1d(p['freq'])
        if self.direction:
            theta, phi = [p['theta']], [p['phi']]
        elif p.get('phi') is not None:
            angles = _Angles(p.get('angles'))
            theta, phi = angles, [p['phi']]
        else:
            angles = _Angles(p.get('angles'))
            theta, phi = [p['theta']], angles
        ef = UI_data([p['probe']], self.res.sim_path, freq)
        P_in = 0.5 * (p.get('E_amp') or 1.0) ** 2 / Z0 * np.abs(np.asarray(ef.ui_f_val[0])) ** 2
        r = self.res._NF2FF(p['box'], freq, theta, phi, p.get('center'))
        if not self.direction:
            return angles, np.squeeze(4 * np.pi / P_in[0] * r.P_rad[0])
        rcs = np.array([4 * np.pi / P_in[n] * r.P_rad[n][0][0] for n in range(len(freq))])
        if self.KIND == 'value':
            return float(rcs[0])
        return freq, rcs


@_register
class DumpSlice(Quantity):
    """ A plane of an HDF5 field or SAR dump: `plane` = (direction,
    position in drawing units), `component` ``'magnitude'``, ``'value'``
    (scalar dumps such as SAR) or a direction (``'x'``, ``'rho'``, ...), at
    `freq` (frequency-domain or DFT of a time-domain dump) or `time_step`
    (index). """
    NAME = 'DumpSlice'
    KIND = 'map'

    def Label(self):
        return self._Cached()[3]['label']

    def _Compute(self):
        from openEMS.Viewer.Source import FieldSource
        from openEMS.Viewer.Quantity import Evaluate
        p = self.params
        fn = p['file'] if os.path.isabs(p['file']) else os.path.join(self.res.sim_path, p['file'])
        src = FieldSource(fn)
        try:
            ny = src.Dir(p['plane'][0])
            lines = src.Lines(self.res._DeltaUnit())
            idx = int(np.argmin(np.abs(lines[ny] - p['plane'][1])))
            if p.get('time_step') is not None:
                data = src.Read(ny, idx, t_idx=p['time_step'])
            elif p.get('freq') is not None:
                data = src.ReadAtFrequency(p['freq'], ny, idx)
            else:
                data = src.Read(ny, idx, f_idx=0) if src.HasFD() else src.Read(ny, idx, t_idx=0)
            comp = p.get('component') or ('magnitude' if src.is_vector else 'value')
            if comp not in ('magnitude', 'value'):
                comp = src.DirNames()[src.Dir(comp)]
            values = Evaluate(data, comp, src.mesh_type, amplitude=True)
            x, y = [l for n, l in enumerate(lines) if n != ny]
            label = src.label if comp == 'value' else '{} {}'.format(src.label, comp)
            meta = dict(power=src.kind == 'SAR', axes=[d for n, d in enumerate(src.DirNames()) if n != ny],
                        label='{} ({})'.format(label, src.unit) if src.unit else label)
        finally:
            src.Close()
        return x, y, np.asarray(values, dtype=float), meta

    def GetMap(self):
        return self._Cached()[:3]

    def IsPower(self):
        """ SAR is a power density, fields are amplitudes (dB scale). """
        return self._Cached()[3]['power']

    def AxisNames(self):
        return self._Cached()[3]['axes']


@_register
class Array(Quantity):
    """ Curve data computed by the script; stored with the results. """
    NAME = 'Array'

    def __init__(self, res, x=None, y=None, label='', domain='f', power=False, **params):
        super().__init__(res, label=label, domain=domain, power=power, **params)
        self.X = domain
        self.POWER = power
        self.data = (np.asarray(x), np.asarray(y))

    def ToDict(self):
        return dict(q=self.NAME, id=self.res._ArrayId(self), label=self.params['label'],
                    domain=self.params['domain'], power=self.POWER)

    def Key(self):
        return 'Array:{}'.format(id(self))

    def Label(self):
        return self.params['label']

    def _Compute(self):
        return self.data


@_register
class ArrayMap(Quantity):
    """ Map data computed by the script; stored with the results. """
    NAME = 'ArrayMap'
    KIND = 'map'

    def __init__(self, res, x=None, y=None, c=None, label='', **params):
        super().__init__(res, label=label, **params)
        self.data = (np.asarray(x), np.asarray(y), np.asarray(c))

    def ToDict(self):
        return dict(q=self.NAME, id=self.res._ArrayId(self), label=self.params['label'])

    def Key(self):
        return 'ArrayMap:{}'.format(id(self))

    def Label(self):
        return self.params['label']

    def _Compute(self):
        return self.data


def _PartLabel(fmt):
    return {'real': 'Re ', 'imag': 'Im ', 'phase': 'arg ', 'phase_unwrap': 'arg '}.get(fmt, '')


def _YLabel(e):
    """ Axis label of a curve entry from its quantity and format, e.g.
    'S-parameter (dB)'; the real and imaginary part share one. """
    name, unit = e['q'].YNAME, e['q'].UNIT
    if e['fmt'] in ('phase', 'phase_unwrap'):
        name, unit = (name + ' phase').strip(), 'deg'
    elif e['fmt'] in ('dB', 'dBi'):
        unit = e['fmt']
    if name and unit:
        return '{} ({})'.format(name, unit)
    return name or unit


def _FormatFreq(f):
    fac, unit = _EngUnit([f], 'Hz')
    return '{:g} {}'.format(f / fac, unit)


# ----------------------------------------------------------------------
# figures

class Figure:
    """ A figure: one or several axes (`grid`) with curves or maps. Created
    with :meth:`Results.AddFigure`. """
    KINDS = ('xy', 'polar', 'map')
    #: options of the figure (title: of the whole figure) and of its axes
    OPTIONS = ('title', 'xlabel', 'ylabel', 'xscale', 'yscale', 'xlim', 'ylim',
               'ylabel_right', 'ylim_right',          # the right y axis (AddCurve(right=True))
               'theta_zero', 'theta_direction')       # polar: e.g. 'N' and -1 (clockwise)

    def __init__(self, res, name, kind='xy', grid=(1, 1), shared_clim=False, sharex=False, **opts):
        if kind not in self.KINDS:
            raise ValueError('unknown figure kind "{}", use one of {}'.format(kind, self.KINDS))
        self.res = res
        self.name = name
        self.kind = kind
        self.grid = tuple(grid)
        self.shared_clim = shared_clim
        self.sharex = sharex
        self.opts = opts          # OPTIONS, for all axes (title: the figure's)
        self.axes_opts = {}       # {axes: OPTIONS} of single axes
        self.lines = []           # reference lines: dicts kind ('h'/'v'), value, axes, kw
        self.entries = []         # dicts: q, fmt, label, axes, style, ...

    def AddCurve(self, q, fmt=None, label=None, axes=0, style=None, right=False, **kw):
        """ AddCurve(q, fmt=None, label=None, axes=0, style=None, right=False, **kw)

        Add the curve of quantity `q` in format `fmt` (see FORMATS, default
        per quantity), to the axes number `axes`; `right`: on a second y
        axis at the right (e.g. a phase next to a magnitude). `style` is a
        matplotlib format string (``'--'``, ``'r-o'``), further keywords
        (``color``, ``linewidth``) go to ``plot``. An axes shows a legend for
        several curves, or for a single one with a `label` of its own. """
        if q.KIND != 'curve':
            raise TypeError('{} is a {}; AddCurve needs a curve'.format(q.NAME, q.KIND))
        self._Axes(axes)
        e = dict(q=q, fmt=fmt or q.DEFAULT_FMT, label=label, axes=axes, style=style, kw=kw)
        if right:
            e['right'] = True
        self.entries.append(e)
        return e

    def AddHLine(self, y, axes=0, **kw):
        """ AddHLine(y, axes=0, **kw)

        A horizontal reference line at `y` (in the plotted unit, e.g. -10
        for -10 dB); keywords go to matplotlib's ``axhline`` (default a thin
        dashed grey line). """
        self._Axes(axes)
        self.lines.append(dict(kind='h', value=y, axes=axes, kw=kw))

    def AddVLine(self, x, axes=0, **kw):
        """ AddVLine(x, axes=0, **kw)

        A vertical reference line at `x` (in the axis' unit as plotted,
        e.g. GHz), see :meth:`AddHLine`. """
        self._Axes(axes)
        self.lines.append(dict(kind='v', value=x, axes=axes, kw=kw))

    def AddMap(self, q, fmt=None, norm=None, label=None, axes=0, cmap=None, clim=None, **kw):
        """ AddMap(q, fmt=None, norm=None, label=None, axes=0, cmap=None, clim=None, **kw)

        Add the map of quantity `q`, divided by the quantity `norm` at the
        map's frequency (e.g. the accepted power ``res.P_acc(1)``). """
        if q.KIND != 'map':
            raise TypeError('{} is a {}; AddMap needs a map'.format(q.NAME, q.KIND))
        self._Axes(axes)
        e = dict(q=q, fmt=fmt, norm=norm, label=label, axes=axes, cmap=cmap, clim=clim, kw=kw)
        self.entries.append(e)
        return e

    def SetAxes(self, axes, **opts):
        """ SetAxes(axes, title=None, xlabel=None, ylabel=None, xscale=None, yscale=None, xlim=None, ylim=None)

        Options of the axes number `axes` (in a grid; they override the
        figure's), e.g. a title per subplot or a fixed range
        (``ylim=(-40, 2)``). """
        self._Axes(axes)
        for opt in opts:
            if opt not in self.OPTIONS:
                raise ValueError('unknown axes option "{}", use one of {}'.format(opt, self.OPTIONS))
        self.axes_opts.setdefault(axes, {}).update(opts)

    def Size(self):
        """ The default size in inches: 6.4 x 4.8 per axes. """
        return 6.4 * self.grid[1], 4.8 * self.grid[0]

    def _Axes(self, axes):
        if not 0 <= axes < self.grid[0] * self.grid[1]:
            raise ValueError('axes {} outside the grid {}'.format(axes, self.grid))

    # ------------------------------------------------------------------
    def MapValues(self, e):
        """ x, y and the formatted, normalised values of the map entry `e`. """
        q = e['q']
        x, y, c = q.GetMap()
        if e.get('norm') is not None:
            if q.params.get('freq') is None:
                raise ValueError('{}: normalising needs the frequency of the map'.format(q))
            c = c / np.real(e['norm'].At(q.params['freq']))
        if e.get('fmt'):
            c = Format(c, e['fmt'], power=q.IsPower() if hasattr(q, 'IsPower') else False)
        return x, y, c

    def Draw(self, fig=None):
        """ Draw(fig=None)

        Draw into the matplotlib figure `fig` (a new one by default) and
        return it. """
        import matplotlib.pyplot as plt
        if fig is None:
            fig = plt.figure(figsize=self.Size())
        n_axes = self.grid[0] * self.grid[1]
        axs = []
        for n in range(n_axes):
            share = axs[0] if self.sharex and self.kind == 'xy' and axs else None
            axs.append(fig.add_subplot(self.grid[0], self.grid[1], n + 1, sharex=share,
                                       projection='polar' if self.kind == 'polar' else None))
        # a shared x axis is labelled at the bottom row only
        inner = set(range((self.grid[0] - 1) * self.grid[1])) if self.sharex and self.kind == 'xy' else set()
        twins = {}   # the right y axes

        def target(e):
            if not e.get('right'):
                return axs[e['axes']]
            if e['axes'] not in twins:
                twins[e['axes']] = axs[e['axes']].twinx()
            return twins[e['axes']]

        clim = None
        if self.shared_clim:
            vals = [self.MapValues(e)[2] for e in self.entries if e['q'].KIND == 'map']
            if vals:
                clim = (min(np.nanmin(v) for v in vals), max(np.nanmax(v) for v in vals))
        xunits = {}
        for e in self.entries:
            q = e['q']
            label = e['label'] or _PartLabel(e['fmt']) + q.Label()
            if q.KIND == 'map':
                ax = axs[e['axes']]
                x, y, c = self.MapValues(e)
                lim = e.get('clim') or clim or (None, None)
                kw = dict(dict(shading='auto'), **e['kw'])
                m = ax.pcolormesh(x, y, c if np.ndim(x) == 2 else c.T, cmap=e.get('cmap'),
                                  vmin=lim[0], vmax=lim[1], **kw)
                if e['label'] is None and e.get('norm') is not None:
                    label = '{} / {}'.format(label, e['norm'].Label())
                fig.colorbar(m, ax=ax, label=label)
                ax.set_aspect('equal')
                if hasattr(q, 'AxisNames'):
                    ax.set_xlabel(q.AxisNames()[0])
                    ax.set_ylabel(q.AxisNames()[1])
                continue
            x, y = q.GetXY()
            y = Format(y, e['fmt'], q.POWER)
            if q.X == 'angle' and self.kind == 'polar':
                x = np.deg2rad(x)
            elif q.X == 'angle':
                xunits[e['axes']] = (q.XName(), 'deg')
            elif q.X in ('f', 't'):
                fac, unit = _EngUnit(x, 'Hz' if q.X == 'f' else 's')
                xunits[e['axes']] = (q.XName(), unit)
                x = np.asarray(x) / fac
            args = (x, y) if e['style'] is None else (x, y, e['style'])
            target(e).plot(*args, label=label, **e['kw'])
        for line in self.lines:
            ax = axs[line['axes']]
            kw = dict(dict(color='gray', linestyle='--', linewidth=1), **line['kw'])
            (ax.axhline if line['kind'] == 'h' else ax.axvline)(line['value'], **kw)
        for n, ax in enumerate(axs):
            curves = [e for e in self.entries if e['axes'] == n and e['q'].KIND == 'curve']
            left = [e for e in curves if not e.get('right')]
            right = [e for e in curves if e.get('right')]
            if self.kind != 'map' and curves:
                ax.grid(True)
                if self.kind == 'xy':
                    ax.set_xmargin(0)   # the curves span the axis
                # a legend for several curves or a label of your own; one for both y axes
                if len(curves) > 1 or any(e['label'] for e in curves):
                    h, l = ax.get_legend_handles_labels()
                    if n in twins:
                        h2, l2 = twins[n].get_legend_handles_labels()
                        h, l = h + h2, l + l2
                    ax.legend(h, l)
            if n in inner:
                ax.tick_params(labelbottom=False)
            elif self.sharex and self.kind == 'xy' and xunits:
                ax.set_xlabel('{} ({})'.format(*xunits.get(n, next(iter(xunits.values())))))
            elif n in xunits and self.kind == 'xy':
                ax.set_xlabel('{} ({})'.format(*xunits[n]))
            o = {k: v for k, v in self.opts.items() if k != 'title'}
            o.update(self.axes_opts.get(n, {}))
            for opt in ('title', 'xlabel', 'ylabel', 'xscale', 'yscale'):
                if opt in o and not (opt == 'xlabel' and n in inner):
                    getattr(ax, 'set_' + opt)(o[opt])
            if self.kind == 'xy':
                for side, es, a in (('', left, ax), ('_right', right, twins.get(n))):
                    if a is None or 'ylabel' + side in o:
                        continue
                    labels = {_YLabel(e) for e in es}
                    if len(labels) == 1:
                        a.set_ylabel(labels.pop())
            if n in twins:
                tw = twins[n]
                if 'ylabel_right' in o:
                    tw.set_ylabel(o['ylabel_right'])
                if 'ylim_right' in o:
                    tw.set_ylim(o['ylim_right'])
                colors = {l.get_color() for l in tw.get_lines()}
                if len(colors) == 1:   # the axis in the colour of its curves
                    c = colors.pop()
                    tw.yaxis.label.set_color(c)
                    tw.tick_params(axis='y', colors=c)
            if self.kind == 'polar':
                if 'theta_zero' in o:
                    ax.set_theta_zero_location(o['theta_zero'])
                if 'theta_direction' in o:
                    ax.set_theta_direction(o['theta_direction'])
                dB = [e for e in curves if e['fmt'] in ('dB', 'dBi')]
                if dB:
                    if 'ylim' in o:
                        rmin = o['ylim'][0]
                    else:
                        top = max(np.nanmax(Format(e['q'].GetXY()[1], e['fmt'], e['q'].POWER)) for e in dB)
                        rmin = top - 40
                        ax.set_ylim(rmin, top + 0.05 * 40)
                    for l in ax.get_lines():   # below the centre would wrap around
                        l.set_ydata(np.maximum(l.get_ydata(), rmin))
            for opt in ('xlim', 'ylim'):
                if opt in o:
                    getattr(ax, 'set_' + opt)(o[opt])
        fig.suptitle(self.opts.get('title', self.name))
        return fig

    # ------------------------------------------------------------------
    def ToDict(self):
        out = dict(name=self.name, kind=self.kind, grid=list(self.grid),
                   opts=setup_file._to_json(self.opts, self.name), entries=[])
        if self.axes_opts:
            out['axes'] = {str(n): setup_file._to_json(o, self.name) for n, o in self.axes_opts.items()}
        if self.lines:
            out['lines'] = setup_file._to_json(self.lines, self.name)
        if self.shared_clim:
            out['shared_clim'] = True
        if self.sharex:
            out['sharex'] = True
        for e in self.entries:
            d = {k: v for k, v in e.items() if v is not None and k not in ('q', 'norm', 'kw')}
            d['q'] = e['q'].ToDict()
            if e.get('norm') is not None:
                d['norm'] = e['norm'].ToDict()
            if e['kw']:
                d['kw'] = e['kw']
            out['entries'].append(d)
        return out


# ----------------------------------------------------------------------

class Results:
    """ Results(FDTD, sim_path, f=None, num=401)

    The results of the simulation of `FDTD` (an :class:`openEMS.openEMS`
    object with its ports and NF2FF boxes, also one recreated with
    `ReadFromXML` + `ReadSetup`), run in `sim_path`.

    :param f: frequencies (Hz) of the frequency-domain quantities; default
        `num` points over the excitation band
    """
    def __init__(self, FDTD, sim_path, f=None, num=401):
        self.FDTD = FDTD
        self.sim_path = sim_path
        self.f = np.asarray(f, dtype=float) if f is not None else self._DefaultFreq(num)
        self.figures = []
        self._port_opts = {}
        self._cache = {}
        self._ports = {}
        self._arrays = []

    def _DefaultFreq(self, num):
        exc = self.FDTD.GetExcitation()
        if exc is None:
            raise ValueError('Results: no excitation set, give the frequencies f')
        if exc['type'] == 'gauss':
            f_max = exc['f0'] + exc['fc']
            return np.linspace(max(exc['f0'] - exc['fc'], f_max / num), f_max, num)
        if exc['type'] == 'sinus':
            return np.array([exc['f0']])
        return np.linspace(exc['f_max'] / num, exc['f_max'], num)

    def SetFrequency(self, f):
        """ SetFrequency(f)

        The frequencies (Hz) of all frequency-domain quantities. """
        self.f = np.asarray(f, dtype=float)
        self._cache.clear()
        self._ports.clear()

    def SetPort(self, port_nr, **opts):
        """ SetPort(port_nr, ref_impedance=None, ref_plane_shift=None, signal_type='pulse')

        Options for `CalcPort` of port `port_nr`, used by all its
        quantities. """
        allowed = ('ref_impedance', 'ref_plane_shift', 'signal_type')
        for k in opts:
            if k not in allowed:
                raise ValueError('unknown port option "{}", use one of {}'.format(k, allowed))
        self._port_opts[port_nr] = opts
        self._cache.clear()
        self._ports.pop(port_nr, None)

    # ------------------------------------------------------------------
    # quantities
    def S(self, i, j):
        """ S-parameter S_ij: wave out of port `i` over wave into port `j`. """
        return S(self, i=i, j=j)

    def Zin(self, port):
        """ Input impedance of `port`. """
        return Zin(self, port=port)

    def Yin(self, port):
        """ Input admittance of `port`. """
        return Yin(self, port=port)

    def P_inc(self, port):
        return P_inc(self, port=port)

    def P_ref(self, port):
        return P_ref(self, port=port)

    def P_acc(self, port):
        """ Accepted power of `port` (incident minus reflected). """
        return P_acc(self, port=port)

    def U(self, port, part='tot', domain='f'):
        """ Voltage of `port`: `part` ``'tot'``, ``'inc'`` or ``'ref'``,
        `domain` ``'f'`` (frequency) or ``'t'`` (time). """
        return U(self, port=port, part=part, domain=domain)

    def I(self, port, part='tot', domain='f'):
        """ Current of `port`, see :meth:`U`. """
        return I(self, port=port, part=part, domain=domain)

    def Beta(self, port):
        """ Propagation constant of a transmission line or waveguide port. """
        return Beta(self, port=port)

    def ZL(self, port):
        """ Line or wave impedance of a transmission line or waveguide port. """
        return ZL(self, port=port)

    def Probe(self, name, domain='f'):
        """ The probe `name` (its file in the simulation path). """
        return Probe(self, name=name, domain=domain)

    def FarFieldCut(self, box, freq, phi=None, theta=None, angles=None, component='norm', center=None):
        """ Directivity of NF2FF box `box` (its name) at `freq` along a cut
        with fixed `phi` or `theta` (degree) over `angles` (default -180 to
        180 in 2° steps). """
        return FarFieldCut(self, box=box, freq=freq, phi=phi, theta=theta, angles=angles,
                           component=component, center=center)

    def Dmax(self, box, freq, theta=None, phi=None, center=None):
        """ Maximum directivity of NF2FF box `box` at `freq` (a value), or at
        each of an array of frequencies (a curve). """
        return Dmax(self, box=box, freq=freq, theta=theta, phi=phi, center=center)

    def Prad(self, box, freq, theta=None, phi=None, center=None):
        """ Radiated power of NF2FF box `box` at `freq` (a value), or at each
        of an array of frequencies (a curve), see :meth:`Dmax`. """
        return Prad(self, box=box, freq=freq, theta=theta, phi=phi, center=center)

    def Efficiency(self, box, freq, port, theta=None, phi=None, center=None):
        """ Radiation efficiency: radiated power of `box` over the power
        accepted by `port`. """
        return Efficiency(self, box=box, freq=freq, port=port, theta=theta, phi=phi, center=center)

    def RCS(self, box, freq, probe, phi=None, theta=None, angles=None, E_amp=1.0, center=None):
        """ Radar cross section of a plane-wave excitation: a cut over angle
        (`phi` or `theta`) or a direction (both), see :class:`RCS`. """
        return RCS(self, box=box, freq=freq, probe=probe, phi=phi, theta=theta, angles=angles,
                   E_amp=E_amp, center=center)

    def DumpSlice(self, file, plane, component=None, freq=None, time_step=None):
        """ A plane of a field or SAR dump, see :class:`DumpSlice`. """
        return DumpSlice(self, file=file, plane=list(plane), component=component, freq=freq,
                         time_step=time_step)

    def Array(self, x, y, label='', domain='f', power=False):
        """ A curve of data computed by the script: x over `domain` ``'f'``
        (Hz), ``'t'`` (s), ``'angle'`` (degree) or None (as is, label the
        axis with ``xlabel``); `power`: a power quantity (10 log10 in dB). """
        return Array(self, x=x, y=y, label=label, domain=domain, power=power)

    def ArrayMap(self, x, y, c, label=''):
        """ A map of data computed by the script: lines `x` and `y` with `c`
        of shape (len(x), len(y)), or 2D coordinate arrays `x`, `y` of the
        shape of `c` (a curvilinear grid, e.g. a cylindrical mesh). """
        return ArrayMap(self, x=x, y=y, c=c, label=label)

    def Quantity(self, d):
        """ The quantity described by the dict `d` (from ToDict()). """
        d = setup_file._from_json(dict(d))
        cls = QUANTITIES.get(d.pop('q'))
        if cls is None:
            raise ValueError('unknown quantity {}'.format(d))
        if cls in (Array, ArrayMap):
            return self._arrays[d['id']]
        return cls(self, **d)

    # ------------------------------------------------------------------
    # data
    def _Port(self, nr):
        if nr not in self._ports:
            import copy
            port = self.FDTD.GetPort(nr)
            port.CalcPort(self.sim_path, self.f, **self._port_opts.get(nr, {}))
            # a copy: the script may call CalcPort on the port again (other
            # frequencies) before a figure reads it
            self._ports[nr] = copy.copy(port)
        return self._ports[nr]

    def _DeltaUnit(self):
        return self.FDTD.GetCSX().GetGrid().GetDeltaUnit()

    def _NF2FF(self, name, freq, theta, phi, center=None, radius=1):
        boxes = [b for b in self.FDTD.GetNF2FFBoxes() if b.name == name]
        if not boxes:
            raise ValueError('no NF2FF box "{}"'.format(name))
        center = [0, 0, 0] if center is None else list(center)
        key = json.dumps([list(np.atleast_1d(freq).astype(float)), list(np.asarray(theta, float)),
                          list(np.asarray(phi, float)), center, radius])
        if key not in self._cache:
            # a file per parameter set, so that a cached file always matches
            outfile = '{}_{}.h5'.format(name, hashlib.sha1(key.encode()).hexdigest()[:10])
            self._cache[key] = boxes[0].CalcNF2FF(self.sim_path, freq, theta, phi, radius=radius,
                                                  center=center, outfile=outfile)
        return self._cache[key]

    def _ArrayId(self, q):
        for n, a in enumerate(self._arrays):
            if a is q:
                return n
        self._arrays.append(q)
        return len(self._arrays) - 1

    # ------------------------------------------------------------------
    # figures
    def AddFigure(self, name, kind='xy', grid=(1, 1), shared_clim=False, sharex=False, **opts):
        """ AddFigure(name, kind='xy', grid=(1, 1), shared_clim=False, sharex=False, title=None, xlabel=None, ylabel=None, xscale=None, yscale=None, xlim=None, ylim=None)

        Add a figure: `kind` ``'xy'`` (curves), ``'polar'`` (curves over
        angle) or ``'map'``; `grid` (rows, columns) of axes; `shared_clim`
        one colour range for all maps; `sharex` one x axis for all axes of
        an ``'xy'`` figure (e.g. S-parameters above the impedance, zoomed
        together), labelled at the bottom only. The axis options apply to all axes,
        :meth:`Figure.SetAxes` sets them per axes. Returns the
        :class:`Figure`. """
        for opt in opts:
            if opt not in Figure.OPTIONS:
                raise ValueError('unknown figure option "{}", use one of {}'.format(opt, Figure.OPTIONS))
        fig = Figure(self, name, kind, grid, shared_clim, sharex, **opts)
        self.figures.append(fig)
        return fig

    def GetFigure(self, name):
        for fig in self.figures:
            if fig.name == name:
                return fig
        raise KeyError('no figure "{}"'.format(name))

    def Show(self, block=True):
        """ Show(block=True)

        Draw all figures with matplotlib and show them. """
        import matplotlib.pyplot as plt
        figs = [fig.Draw() for fig in self.figures]
        plt.show(block=block)
        return figs

    # ------------------------------------------------------------------
    # persistence
    def ToDict(self):
        f = self.f
        if len(f) > 2 and np.allclose(np.diff(f), (f[-1] - f[0]) / (len(f) - 1), rtol=1e-9, atol=0):
            freq = dict(start=float(f[0]), stop=float(f[-1]), num=len(f))
        else:
            freq = [float(v) for v in f]
        self._arrays = [a for a in self._arrays if any(
            e['q'] is a for fig in self.figures for e in fig.entries)]
        return dict(f=freq, ports={str(k): v for k, v in self._port_opts.items()},
                    figures=[fig.ToDict() for fig in self.figures], items=[])

    def Write(self, file):
        """ Write(file)

        Store the results description under ``"results"`` in the setup file
        `file` (written by :meth:`openEMS.openEMS.WriteSetup`; created if
        missing). Array data goes to the same name ending in ``.npz``. """
        data = {'format': setup_file.FORMAT, 'version': setup_file.VERSION}
        if os.path.exists(file):
            with open(file) as fh:
                data = json.load(fh)
        data['results'] = self.ToDict()
        text = json.dumps(data, indent=1, allow_nan=False)
        npz = os.path.splitext(file)[0] + '.npz'
        if self._arrays:
            arrays = {}
            for n, a in enumerate(self._arrays):
                for k, v in zip('xyc', a.data):
                    arrays['{}_{}'.format(n, k)] = v
            np.savez(npz, **arrays)
        elif os.path.exists(npz):
            os.remove(npz)
        with open(file, 'w') as fh:
            fh.write(text)

    def Read(self, file):
        """ Read(file)

        Restore frequencies, port options and figures from the setup file
        `file`, replacing the current ones. """
        with open(file) as fh:
            data = json.load(fh).get('results')
        if data is None:
            raise ValueError('{}: no results'.format(file))
        f = data['f']
        f = np.linspace(f['start'], f['stop'], f['num']) if isinstance(f, dict) else np.asarray(f, float)
        self.SetFrequency(f)
        self._port_opts = {int(k): setup_file._from_json(v) for k, v in data.get('ports', {}).items()}
        self._arrays = []
        npz = os.path.splitext(file)[0] + '.npz'
        if os.path.exists(npz):
            with np.load(npz) as z:
                n = 0
                while '{}_x'.format(n) in z:
                    if '{}_c'.format(n) in z:
                        self._arrays.append(ArrayMap(self, z['{}_x'.format(n)], z['{}_y'.format(n)],
                                                     z['{}_c'.format(n)]))
                    else:
                        self._arrays.append(Array(self, z['{}_x'.format(n)], z['{}_y'.format(n)]))
                    n += 1
        self.figures = []
        for fd in data.get('figures', []):
            fig = self.AddFigure(fd['name'], fd.get('kind', 'xy'), fd.get('grid', (1, 1)),
                                 fd.get('shared_clim', False), fd.get('sharex', False), **fd.get('opts', {}))
            for n, o in fd.get('axes', {}).items():
                fig.SetAxes(int(n), **o)
            for line in fd.get('lines', []):
                (fig.AddHLine if line['kind'] == 'h' else fig.AddVLine)(line['value'], line['axes'], **line['kw'])
            for e in fd.get('entries', []):
                q = self.Quantity(e['q'])
                if q.NAME in ('Array', 'ArrayMap'):   # label, domain and power from the file
                    q.params.update({k: v for k, v in e['q'].items() if k in ('label', 'domain', 'power')})
                    q.X = q.params.get('domain', q.X)
                    q.POWER = q.params.get('power', False)
                opts = {k: v for k, v in e.items() if k not in ('q', 'norm', 'kw')}
                if q.KIND == 'map':
                    norm = self.Quantity(e['norm']) if 'norm' in e else None
                    fig.AddMap(q, norm=norm, **opts, **e.get('kw', {}))
                else:
                    fig.AddCurve(q, **opts, **e.get('kw', {}))

    @classmethod
    def Load(cls, FDTD, sim_path, file):
        """ Load(FDTD, sim_path, file)

        A :class:`Results` object restored from the setup file `file`, for
        `FDTD` with its ports and boxes recreated by
        :meth:`openEMS.openEMS.ReadSetup`. """
        res = cls(FDTD, sim_path, f=[1.0])
        res.Read(file)
        return res
