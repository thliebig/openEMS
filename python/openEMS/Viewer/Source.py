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

"""openEMS.Viewer.Source
Field data for the viewer items: an openEMS HDF5 dump or numpy arrays.
"""

import os

import numpy as np

#: dump_type % 10 of a field dump -> field kind
_KINDS = {0: 'E', 1: 'H', 2: 'J', 3: 'rotH', 4: 'D', 5: 'B'}
#: SAR dump types -> label
_SAR = {20: 'SAR', 21: 'SAR 1g', 22: 'SAR 10g', 29: 'SAR'}
MESH_TYPES = ('cartesian', 'cylindrical')
_DIR_NAMES = {'cartesian': ('x', 'y', 'z'), 'cylindrical': ('rho', 'alpha', 'z')}


def AsSource(src):
    """ A :class:`FieldSource` for a path, an ``HDF5Dump`` or a FieldSource. """
    return src if isinstance(src, FieldSource) else FieldSource(src)


class FieldSource:
    """ FieldSource(dump)

    Field data shown by the viewer items, read plane by plane (or as a
    volume) so that only what is shown is read from disk.

    Parameters
    ----------
    dump : str or openEMS.utilities.HDF5Dump
        An openEMS HDF5 field or SAR dump file, or an open dump.

    See Also
    --------
    FromArrays : a source for fields computed in Python
    """

    def __init__(self, dump=None):
        self._dump = None
        self._arrays = None
        self._owns = False
        self.path = None
        if dump is None:
            return
        from openEMS.utilities import HDF5Dump
        self._owns = not isinstance(dump, HDF5Dump)
        if self._owns:
            self.path = os.fspath(dump)
            dump = HDF5Dump(self.path)
        else:
            self.path = dump.file.filename
        self._dump = dump
        mesh = dump.GetMesh()
        if mesh['type'] not in (0, 1):
            raise ValueError('"{}": only Cartesian and cylindrical dumps can be shown'.format(self.path))
        #: 'cartesian' or 'cylindrical'
        self.mesh_type = MESH_TYPES[mesh['type']]
        self._lines = [np.asarray(l, dtype=float) for l in mesh['lines']]
        #: the length unit of the simulation mesh, in m
        self.delta_unit = mesh['scaling']
        self.is_vector = dump.IsVector()
        dt = dump.dump_type
        if dt is not None and dt in _SAR:
            self.kind, self.label = 'SAR', _SAR[dt]
        elif dt is None and not self.is_vector and 'mass' in dump.file.attrs:
            # a result of SAR_Calculation: no dump type, but the total mass
            self.kind, self.label = 'SAR', 'SAR'
        else:
            self.kind = _KINDS.get(dt % 10 if dt is not None else 0, 'E')
            self.label = self.kind
        self.unit = UNITS.get(self.kind, '')
        self.freqs = np.asarray(dump.frequencies, dtype=float)
        self.times = np.asarray(dump.GetTimes(), dtype=float)
        self.name = os.path.splitext(os.path.basename(self.path))[0]

    @classmethod
    def FromArrays(cls, lines, field, domain='frequency', freqs=None, times=None, mesh_type='cartesian',
                   kind='E', unit=None, delta_unit=1.0, name=None, label=None):
        """ FromArrays(lines, field, domain='frequency', freqs=None, times=None, mesh_type='cartesian', kind='E', unit=None, delta_unit=1.0, name=None)

        A source for a field computed in Python.

        Parameters
        ----------
        lines : sequence of 3 arrays
            The mesh lines in drawing units (``alpha`` in radians).
        field : array
            ``(n, 3, Nx, Ny, Nz)`` for a vector field, ``(n, Nx, Ny, Nz)``
            for a scalar one, with `n` frequencies (complex) or time steps
            (real).
        domain : 'frequency' or 'time'
        freqs, times : array
            The `n` frequencies (Hz) or times (s).
        mesh_type : 'cartesian' or 'cylindrical'
        kind : str
            ``'E'``, ``'H'``, ``'J'``, ``'SAR'``, ... (sets label and unit).
        unit : str, optional
            Unit of the values, default from `kind`.
        delta_unit : float
            The drawing unit in m.
        """
        self = cls()
        if mesh_type not in MESH_TYPES:
            raise ValueError('mesh_type must be one of {}'.format(MESH_TYPES))
        self.mesh_type = mesh_type
        self.delta_unit = float(delta_unit)
        self._lines = [np.atleast_1d(np.asarray(l, dtype=float)) for l in lines]
        for n in (0, 2) if mesh_type == 'cylindrical' else (0, 1, 2):
            self._lines[n] = self._lines[n] * self.delta_unit
        field = np.asarray(field)
        shape = tuple(len(l) for l in self._lines)
        if field.shape[-3:] != shape or field.ndim not in (4, 5):
            raise ValueError('field shape {} does not match (n, [3,] {}, {}, {})'.format(field.shape, *shape))
        self.is_vector = field.ndim == 5
        if domain == 'frequency':
            self.freqs = np.atleast_1d(np.asarray(freqs if freqs is not None else np.zeros(len(field)), float))
            self.times = np.array([])
        elif domain == 'time':
            self.times = np.atleast_1d(np.asarray(times if times is not None else np.arange(len(field)), float))
            self.freqs = np.array([])
        else:
            raise ValueError('domain must be "frequency" or "time"')
        if len(field) != max(len(self.freqs), len(self.times)):
            raise ValueError('field has {} samples, but {} frequencies/times are given'.format(
                len(field), max(len(self.freqs), len(self.times))))
        self._arrays = field
        self.kind = kind
        self.label = label or kind
        self.unit = UNITS.get(kind, '') if unit is None else unit
        self.name = name
        return self

    # ------------------------------------------------------------------
    def __repr__(self):
        return '<FieldSource {}: {} {}, {}, {} frequencies, {} time steps>'.format(
            self.name, self.label, 'vector' if self.is_vector else 'scalar', self.mesh_type,
            len(self.freqs), len(self.times))

    def HasFD(self):
        """ True if frequency-domain data is stored. """
        return len(self.freqs) > 0

    def HasTD(self):
        """ True if time-domain data is stored. """
        return len(self.times) > 0

    def Domains(self):
        """ The domains the data can be shown in; a time-domain source can be
        shown at a frequency too (by a DFT). """
        out = []
        if self.HasFD() or self.HasTD():
            out.append('frequency')
        if self.HasTD():
            out.append('time')
        return out

    def DirNames(self):
        return _DIR_NAMES[self.mesh_type]

    def Dir(self, d):
        """ Direction index 0/1/2 of a name (``'x'``, ``'rho'``, ...) or index. """
        names = self.DirNames()
        if isinstance(d, str):
            if d.lower() not in names:
                raise ValueError('invalid direction "{}", use one of {}'.format(d, names))
            return names.index(d.lower())
        d = int(d)
        if d not in (0, 1, 2):
            raise ValueError('invalid direction {}'.format(d))
        return d

    def Shape(self):
        return tuple(len(l) for l in self._lines)

    def Lines(self, delta_unit=None):
        """ The mesh lines in drawing units (`delta_unit` in m, default the
        dump's own unit); ``alpha`` stays in radians. """
        delta = self.delta_unit if delta_unit is None else delta_unit
        out = []
        for n, l in enumerate(self._lines):
            out.append(l.copy() if (self.mesh_type == 'cylindrical' and n == 1) else l / delta)
        return out

    # ------------------------------------------------------------------
    # reading
    def _Select(self, ny, idx, sample):
        """ Index tuple into (Nx, Ny, Nz) for a plane (ny, idx) or the volume. """
        sel = []
        for n in range(3):
            if ny is not None and n == ny:
                sel.append(int(idx))
            else:
                sel.append(slice(None, None, max(int(sample), 1)))
        return tuple(sel)

    def SampledLines(self, ny=None, sample=1, delta_unit=None):
        """ The lines matching a read of plane normal `ny` (or the volume)
        with `sample`, in drawing units: two arrays for a plane, three for
        the volume. """
        lines = self.Lines(delta_unit)
        sel = self._Select(ny, 0, sample)
        return [l[s] for n, (l, s) in enumerate(zip(lines, sel)) if n != ny]

    def Read(self, ny=None, idx=0, sample=1, f_idx=None, t_idx=None):
        """ Read(ny=None, idx=0, sample=1, f_idx=None, t_idx=None)

        Read the plane with normal `ny` at line `idx` (the volume for
        ``ny=None``), every `sample`-th line, at frequency index `f_idx`
        (complex) or time step `t_idx` (real).

        Returns
        -------
        array
            ``(3, N1, N2)`` / ``(N1, N2)`` for a plane, ``(3, Nx, Ny, Nz)`` /
            ``(Nx, Ny, Nz)`` for the volume.
        """
        sel = self._Select(ny, idx, sample)
        if self._arrays is not None:
            n = f_idx if f_idx is not None else t_idx
            data = self._arrays[n]
            return np.asarray(data[(slice(None),) + sel] if self.is_vector else data[sel])
        d = self._Region(ny, idx, sample)
        if f_idx is not None:
            return d.GetFieldAtIndex(f_idx=int(f_idx))
        return d.GetFieldAtIndex(t_idx=int(t_idx))

    def ReadAtFrequency(self, freq, ny=None, idx=0, sample=1):
        """ The field at frequency `freq` (Hz): a stored frequency, or the
        DFT of the time-domain data (reads every time step of the plane). """
        if self.HasFD():
            n = int(np.argmin(np.abs(self.freqs - freq)))
            if np.isclose(self.freqs[n], freq, rtol=1e-6) or not self.HasTD():
                return self.Read(ny, idx, sample, f_idx=n)
        if self._arrays is not None:
            dt = self.times[1] - self.times[0]
            acc = 0
            for n, t in enumerate(self.times):
                acc = acc + self.Read(ny, idx, sample, t_idx=n) * np.exp(-2j * np.pi * freq * t)
            return 2 * acc * dt
        return self._Region(ny, idx, sample).GetFieldAtFrequency(freq)

    def _Region(self, ny, idx, sample):
        """ The dump, set to read the plane (ny, idx) or the volume with `sample`. """
        if self._dump is None:   # closed: open it again
            from openEMS.utilities import HDF5Dump
            self._dump = HDF5Dump(self.path)
        d = self._dump
        d.ResetRegion()
        d.SetSampling(*[1 if n == ny else max(int(sample), 1) for n in range(3)])
        if ny is not None:
            d.SetPlane(ny, idx=int(idx))
        return d

    def Close(self):
        """ Close the dump file this source opened; the next read opens it
        again. A dump given as ``HDF5Dump`` stays open. """
        if self._dump is not None and self._owns:
            self._dump.Close()
            self._dump = None


#: units of the field kinds
UNITS = dict(E='V/m', H='A/m', J='A/m²', rotH='A/m²', D='As/m²', B='Vs/m²', SAR='W/kg')
