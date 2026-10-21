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

"""openEMS.Viewer.Quantity
What is shown of a field: component or magnitude, instantaneous or
amplitude, linear or in dB, and the value range.
"""

import numpy as np

#: quantities whose dB value is 10*log10 (power-like); others use 20*log10
POWER_KINDS = ('SAR', 'P_rad', 'directivity', 'gain')


def Components(mesh_type, is_vector):
    """ The quantities of a source: the components and ``'magnitude'``, or
    ``'value'`` for a scalar. A cylindrical source offers the Cartesian
    components too. """
    if not is_vector:
        return ['value']
    if mesh_type == 'cylindrical':
        return ['rho', 'alpha', 'z', 'x', 'y', 'magnitude']
    return ['x', 'y', 'z', 'magnitude']


def Normalize(quantity, kind, choices):
    """ Accept ``'E_z'`` as well as ``'z'``; None picks the default. """
    if quantity is None:
        return 'magnitude' if 'magnitude' in choices else choices[0]
    q = str(quantity)
    for prefix in (kind + '_', kind.lower() + '_'):
        if q.startswith(prefix):
            q = q[len(prefix):]
    if q in ('abs', 'norm', 'mag'):
        q = 'magnitude'
    if q not in choices:
        raise ValueError('invalid quantity "{}", use one of {}'.format(quantity, choices))
    return q


def Label(kind, quantity, amplitude=False):
    """ Short label such as ``E_z``, ``|E_z|`` or ``|E|``. """
    if quantity == 'value':
        return kind
    if quantity == 'magnitude':
        return '|{}|'.format(kind)
    name = '{}_{}'.format(kind, quantity)
    return '|{}|'.format(name) if amplitude else name


def CartesianVector(data, mesh_type, alpha=None):
    """ The vector `data` ``(3, ...)`` in Cartesian components; a
    cylindrical one is rotated by the angle `alpha` (same shape as one
    component). """
    if mesh_type != 'cylindrical':
        return data
    c, s = np.cos(alpha), np.sin(alpha)
    return np.stack([data[0] * c - data[1] * s, data[0] * s + data[1] * c, data[2]])


def ArrowVector(data, quantity, mesh_type='cartesian', alpha=None):
    """ Cartesian vectors ``(3, ...)`` of real field `data` for arrows: the
    whole field for the magnitude, else only the component `quantity`
    along its own direction. """
    if quantity == 'magnitude':
        return CartesianVector(data, mesh_type, alpha)
    if mesh_type == 'cylindrical' and quantity in ('rho', 'alpha'):
        only = np.zeros_like(data)
        n = ('rho', 'alpha').index(quantity)
        only[n] = data[n]
        return CartesianVector(only, mesh_type, alpha)
    vec = CartesianVector(data, mesh_type, alpha)
    out = np.zeros_like(vec)
    n = ('x', 'y', 'z').index(quantity)
    out[n] = vec[n]
    return out


def Instantaneous(data, phase_deg):
    """ ``Re(F e^{j phase})`` of complex data, real data unchanged. """
    if not np.iscomplexobj(data):
        return data
    return np.real(data * np.exp(1j * np.deg2rad(phase_deg)))


def Evaluate(data, quantity, mesh_type='cartesian', alpha=None, phase=0.0, amplitude=False):
    """ Evaluate(data, quantity, mesh_type='cartesian', alpha=None, phase=0.0, amplitude=False)

    The shown real value of field `data` (complex for frequency domain).

    Parameters
    ----------
    data : array
        ``(3, ...)`` vector or ``(...)`` scalar field.
    quantity : str
        A component (``'x'``, ``'rho'``, ...), ``'magnitude'`` or ``'value'``.
    alpha : array, optional
        The angle at each point, for Cartesian components of a cylindrical
        field.
    phase : float
        Phase in degrees for the instantaneous value of complex data.
    amplitude : bool
        Complex data: the amplitude ``|F|`` (peak amplitude
        ``sqrt(sum |F_i|^2)`` for the magnitude) instead of the
        instantaneous value.
    """
    if quantity == 'value':
        v = np.abs(data) if np.iscomplexobj(data) else data
        return np.asarray(v, dtype=np.float32)
    if quantity == 'magnitude':
        if amplitude or not np.iscomplexobj(data):
            if np.iscomplexobj(data):
                return np.sqrt(np.sum(np.abs(data) ** 2, axis=0)).astype(np.float32)
            return np.sqrt(np.sum(data ** 2, axis=0)).astype(np.float32)
        inst = Instantaneous(data, phase)
        return np.sqrt(np.sum(inst ** 2, axis=0)).astype(np.float32)
    names = ['rho', 'alpha', 'z'] if mesh_type == 'cylindrical' else ['x', 'y', 'z']
    if quantity in names:
        comp = data[names.index(quantity)]
    elif quantity in ('x', 'y') and mesh_type == 'cylindrical':
        c, s = np.cos(alpha), np.sin(alpha)
        comp = data[0] * c - data[1] * s if quantity == 'x' else data[0] * s + data[1] * c
    else:
        raise ValueError('invalid quantity "{}"'.format(quantity))
    if np.iscomplexobj(comp):
        comp = np.abs(comp) if amplitude else Instantaneous(comp, phase)
    return np.asarray(comp, dtype=np.float32)


def PeakAbs(data, quantity, mesh_type='cartesian', alpha=None):
    """ The largest absolute value the quantity reaches over all phases. """
    if np.iscomplexobj(data):
        v = Evaluate(data, quantity, mesh_type, alpha, amplitude=True)
    else:
        v = np.abs(Evaluate(data, quantity, mesh_type, alpha))
    return float(np.nanmax(np.abs(v))) if v.size else 0.0


def IsPower(kind):
    return kind in POWER_KINDS


def ToDB(values, power, floor_db):
    """ ``|values|`` in dB (``10*log10`` for power-like, else ``20*log10``),
    clipped below at `floor_db`. """
    factor = 10.0 if power else 20.0
    with np.errstate(divide='ignore'):
        db = factor * np.log10(np.abs(values))
    return np.maximum(db, floor_db).astype(np.float32)


def MaxDB(peak, power):
    if peak <= 0:
        return 0.0
    return (10.0 if power else 20.0) * np.log10(peak)


def AutoRange(peak, signed, scale='linear', power=False, dynamic_range=40.0):
    """ The colour range for the largest absolute value `peak`: symmetric
    around 0 for signed linear values, ``max - dynamic_range ... max`` in
    dB. """
    if scale == 'log':
        top = MaxDB(peak, power)
        return (top - float(dynamic_range), top)
    if peak <= 0:
        peak = 1.0
    return (-peak, peak) if signed else (0.0, peak)


def ParseNumbers(value):
    """ A list of floats from a sequence or from its text in the settings
    form, e.g. ``'-1, 2'``, ``'[0, 0, 100]'`` or ``'(1 2)'``. """
    if isinstance(value, str):
        for c in '()[],':
            value = value.replace(c, ' ')
        value = value.split()
    return [float(v) for v in np.atleast_1d(value)]


def ParseRange(value_range):
    """ ``'auto'`` or a ``(vmin, vmax)`` pair, also from a string. """
    if value_range is None or (isinstance(value_range, str) and value_range.strip().lower() in ('', 'auto')):
        return 'auto'
    vals = ParseNumbers(value_range)
    if len(vals) != 2 or not vals[0] < vals[1]:
        raise ValueError('value_range must be "auto" or (vmin, vmax) with vmin < vmax')
    return vals


def ParseSize(size, reference):
    """ ``'20%'`` of `reference`, or an absolute number. """
    if isinstance(size, str) and size.strip().endswith('%'):
        return float(size.strip()[:-1]) / 100.0 * reference
    return float(size)


def FormatFreq(f):
    for div, unit in ((1e12, 'THz'), (1e9, 'GHz'), (1e6, 'MHz'), (1e3, 'kHz')):
        if abs(f) >= div:
            return '{:.4g} {}'.format(f / div, unit)
    return '{:.4g} Hz'.format(f)


def FormatTime(t):
    if t == 0:
        return '0 s'
    for div, unit in ((1, 's'), (1e-3, 'ms'), (1e-6, 'µs'), (1e-9, 'ns'), (1e-12, 'ps')):
        if abs(t) >= div:
            return '{:.4g} {}'.format(t / div, unit)
    return '{:.4g} fs'.format(t / 1e-15)


def LengthUnit(delta_unit):
    """ Name of a drawing unit in m, or '' if it has none. """
    for val, name in ((1.0, 'm'), (1e-2, 'cm'), (1e-3, 'mm'), (1e-6, 'µm'), (1e-9, 'nm')):
        if np.isclose(delta_unit, val, rtol=1e-9):
            return name
    return ''
