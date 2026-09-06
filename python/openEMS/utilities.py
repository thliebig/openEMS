# -*- coding: utf-8 -*-
#
# Copyright (C) 2015,20016 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

import numpy as np
import h5py

from CSXCAD.Utilities import CheckNyDir

def DFT_time2freq( t, val, freq, signal_type='pulse'):
    assert len(t)==len(val)
    assert len(freq)>0
    f_val = np.zeros(len(freq))*1j
    for n_f in range(len(freq)):
        f_val[n_f] = np.sum( val*np.exp( -1j*2*np.pi*freq[n_f] * t ) )

    if signal_type == 'pulse':
        f_val *= t[1]-t[0]
    elif signal_type == 'periodic':
        f_val /= len(t)
    else:
        raise Exception('Unknown signal type: "{}"'.format(signal_type))

    return 2*f_val  # single-sided spectrum

def Check_Array_Equal(a,b, tol, relative=False):
    a = np.array(a)
    b = np.array(b)
    if a.shape!=b.shape:
        return False
    if tol==0:
        return (a==b).all()
    if relative:
        d = np.abs((a-b)/a)
    else:
        d = np.abs((a-b))
    return np.max(d)<tol

def check_mode_purity(label, signal, purity, threshold=0.99, sig_frac=0.01):
    """Assert mode purity > threshold where the signal exceeds sig_frac * peak.

    Parameters
    ----------
    label : str
        Descriptive name used in the assertion message.
    signal : array
        Time-domain signal amplitude (column 1 of probe file).
    purity : array or None
        Mode purity time series (column 2 of probe file), or None if unavailable.
    threshold : float
        Minimum acceptable mode purity (default 0.99 = 99 %).
    sig_frac : float
        Ignore time steps where ``|signal| < sig_frac * max(|signal|)``.

    Notes
    -----
    Purity can be negative when the wave travels in the opposite direction
    (e.g. the receive port seeing the transmitted wave), so abs(purity) is used.
    """
    if purity is None:
        return
    mask = np.abs(signal) >= sig_frac * np.max(np.abs(signal))
    if not np.any(mask):
        return
    min_purity = np.min(np.abs(purity[mask]))
    print('{}: min mode purity = {:.1f}% ({:.1f}% of samples considered)'.format(
        label, 100*min_purity, 100*np.sum(mask)/len(signal)))
    assert min_purity >= threshold, \
        '{}: mode purity {:.1f}% below {:.0f}% threshold'.format(
            label, 100*min_purity, 100*threshold)


class HDF5Dump:
    """Reader for an openEMS HDF5 field dump file.

    Opens the file, reads the (cheap) metadata up front and keeps the file
    open for subsequent field reads.  All metadata -- the dump type, whether
    the file holds time or frequency domain data, the number of samples and
    the mesh size -- is available *before* any field data is read.

    The region of interest is configured on the object and then applies to
    every field access.  `SetPlane`, `SetRange` and `SetSampling` restrict
    the read to a plane, a sub-range or a sub-sampled grid; the selection is
    passed down to HDF5 so that only the requested data is read from disk.
    Field data is always returned in the natural, x-inner order
    ``(3, Nx, Ny, Nz)`` regardless of how the file is stored on disk.

    Parameters
    ----------
    filename : str or h5py.File
        Path to the dump file, or an already open file.  An already open file
        is *not* closed by this class.

    Examples
    --------
    >>> with HDF5Dump('Ef.h5') as dump:
    ...     print(dump)
    ...     dump.SetPlane('z', pos=10e-3)          # nearest z-line to 10 mm
    ...     dump.SetSampling(2, 2, 1)              # every other x and y line
    ...     for freq, field in dump.IterFD():
    ...         pass

    A time domain dump can be evaluated at a frequency directly; the DFT is
    then done on the fly, holding only one timestep in memory at a time:

    >>> with HDF5Dump('Et.h5') as dump:
    ...     E = dump.GetFieldAtFrequency(2.4e9)

    See Also
    --------
    openEMS.sar_utils.readSAR : convenience reader for SAR result files
    """

    #: Relative tolerance when matching a frequency requested from
    #: `GetFieldAtFrequency` against the frequencies stored in the file.
    FREQ_RTOL = 1e-6

    # names of the /Mesh datasets, keyed by the mesh_type attribute;
    # mirrors the openEMS HDF5 writer (tools/hdf5_file_writer.cpp)
    _MESH_NAMES = {0: ('x', 'y', 'z'),            # Cartesian
                   1: ('rho', 'alpha', 'z'),      # cylindrical
                   2: ('r', 'theta', 'phi')}      # spherical (NF2FF)

    # dump_type attribute --> human readable name
    _DUMP_TYPE_NAMES = {  0: 'E-field (TD)',      10: 'E-field (FD)',
                          1: 'H-field (TD)',      11: 'H-field (FD)',
                          2: 'J-field (TD)',      12: 'J-field (FD)',
                          3: 'rot(H)-field (TD)', 13: 'rot(H)-field (FD)',
                          4: 'D-field (TD)',      14: 'D-field (FD)',
                          5: 'B-field (TD)',      15: 'B-field (FD)',
                         20: 'local SAR',         21: '1g averaged SAR',
                         22: '10g averaged SAR',  29: 'raw SAR data'}

    def __init__(self, filename):
        self._owns_file = not isinstance(filename, h5py.File)
        self._h5 = h5py.File(filename, 'r') if self._owns_file else filename
        try:
            self._Init()
        except Exception:
            if self._owns_file:
                self._h5.close()
            raise

    ###########################################################################
    # setup / teardown
    ###########################################################################

    def _Init(self):
        h5 = self._h5
        self._root_attrs = dict(h5.attrs)

        if 'Mesh' not in h5:
            raise KeyError('"{}" does not contain a /Mesh group, this does not '
                           'look like an openEMS HDF5 dump'.format(h5.filename))
        self._mesh = self._ReadMesh()

        self._td_names = sorted(h5['FieldData/TD'].keys()) if 'FieldData/TD' in h5 else []
        self._fd_indices = self._ScanFDIndices()
        if not self._td_names and not self._fd_indices:
            raise KeyError('"{}" does not contain any /FieldData/TD or '
                           '/FieldData/FD data'.format(h5.filename))

        self._frequencies = self._ReadFrequencies()
        self._legacy, self._shape, self._is_vector = self._ProbeLayout()

        self.ResetRegion()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.Close()
        return False

    def Close(self):
        """Close the file, unless it was handed in already open."""
        if self._owns_file and self._h5:
            self._h5.close()
        self._h5 = None

    def __repr__(self):
        if self._h5 is None:
            return '<HDF5Dump (closed)>'
        domain = []
        if self._td_names:
            domain.append('{} timesteps'.format(len(self._td_names)))
        if self._fd_indices:
            domain.append('{} frequencies'.format(len(self._fd_indices)))
        region = 'full' if self._region == [None]*3 else str(self._region)
        if self._sampling != [1, 1, 1]:
            region += ', sampling {}'.format(tuple(self._sampling))
        return '<HDF5Dump {!r}: {}, {}, shape {}, region {}>'.format(
            self._h5.filename, self.DumpTypeName, ', '.join(domain),
            self._shape, region)

    ###########################################################################
    # metadata gathered at open time
    ###########################################################################

    def _ReadMesh(self):
        grp = self._h5['Mesh']
        attrs = dict(grp.attrs)
        m_type = int(attrs.get('mesh_type', 0))
        names = self._MESH_NAMES.get(m_type, self._MESH_NAMES[0])
        if not all(n in grp for n in names):
            # no (or a wrong) mesh_type attribute: identify by the datasets present
            for t, cand in self._MESH_NAMES.items():
                if all(n in grp for n in cand):
                    m_type, names = t, cand
                    break
            else:
                raise KeyError('"{}" does not contain a valid /Mesh group'.format(self._h5.filename))
        return {'lines': [np.array(grp[n]) for n in names],
                'names': list(names), 'type': m_type,
                'scaling': float(attrs.get('mesh_scaling', 1.0))}

    def _ScanFDIndices(self):
        """Collect the numeric indices of the f<n> datasets present."""
        if 'FieldData/FD' not in self._h5:
            return []
        found = set()
        for key in self._h5['FieldData/FD'].keys():
            if not key.startswith('f'):
                continue
            stem = key[1:]
            for suffix in ('_real', '_imag'):
                if stem.endswith(suffix):
                    stem = stem[:-len(suffix)]
                    break
            if stem.isdigit():
                found.add(int(stem))
        return sorted(found)

    def _ReadFrequencies(self):
        if not self._fd_indices:
            return np.array([])
        grp = self._h5['FieldData/FD']
        if 'frequency' in grp.attrs:
            return np.atleast_1d(np.array(grp.attrs['frequency'], dtype=float))
        # fall back to the per-dataset attribute
        freqs = []
        for n in self._fd_indices:
            ds = self._FDDataset(grp, n)
            freqs.append(float(np.asarray(ds.attrs['frequency']).flatten()[0])
                         if 'frequency' in ds.attrs else np.nan)
        return np.array(freqs)

    def _ProbeLayout(self):
        """Determine axis order, logical grid shape and vector/scalar layout."""
        if self._fd_indices:
            ds = self._FDDataset(self._h5['FieldData/FD'], self._fd_indices[0])
        else:
            ds = self._h5['FieldData/TD'][self._td_names[0]]
        # the d_order attribute is authoritative where present, otherwise the
        # file version decides (mirrors the C++ reader): no version or <= 0.2
        # is always legacy, above that the legacy_fmt attribute says so
        if 'd_order' in ds.attrs:
            legacy = self._DecodeStr(ds.attrs['d_order']).upper() in ('NZYX', 'ZYX')
        elif float(self._root_attrs.get('openEMS_HDF5_version', 0.0)) <= 0.2:
            legacy = True
        else:
            legacy = bool(self._root_attrs.get('legacy_fmt', False))

        # len(shape) instead of Dataset.ndim: the latter is missing in the
        # ancient h5py of the oldest supported distributions
        is_vector = len(ds.shape) == 4
        shape = tuple(ds.shape[1:] if is_vector else ds.shape)
        if legacy:                    # stored (Nz,Ny,Nx) --> logical (Nx,Ny,Nz)
            shape = shape[::-1]
        return legacy, shape, is_vector

    @staticmethod
    def _DecodeStr(val):
        """Decode an HDF5 string attribute (str, bytes or 1-element array)."""
        val = np.asarray(val).flatten()
        if len(val) == 0:
            return ''
        val = val[0]
        if isinstance(val, bytes):
            return val.decode('utf-8', 'replace')
        return str(val)

    @staticmethod
    def _FDDataset(grp, n):
        """Return the dataset holding frequency index *n* (compound or legacy)."""
        for name in ('f{}'.format(n), 'f{}_real'.format(n)):
            if name in grp:
                return grp[name]
        raise KeyError('no dataset for frequency index {}'.format(n))

    ###########################################################################
    # public metadata
    ###########################################################################

    @property
    def File(self):
        """The open `h5py.File`, for direct access to anything not wrapped here
        (e.g. the ``/CellData`` and ``/CellWidth`` groups of a raw SAR dump)."""
        return self._h5

    @property
    def DumpType(self):
        """The openEMS dump type as an integer, or None if not stored."""
        if 'dump_type' not in self._root_attrs:
            return None
        return int(np.asarray(self._root_attrs['dump_type']).flatten()[0])

    @property
    def DumpTypeName(self):
        """Human readable name of the dump type."""
        return self._DUMP_TYPE_NAMES.get(self.DumpType, 'unknown dump type')

    @property
    def IsTD(self):
        """True if the file holds time domain data."""
        return len(self._td_names) > 0

    @property
    def IsFD(self):
        """True if the file holds frequency domain data."""
        return len(self._fd_indices) > 0

    @property
    def IsVector(self):
        """True for a vector field dump, False for a scalar one (e.g. SAR)."""
        return self._is_vector

    @property
    def Shape(self):
        """Grid size (Nx, Ny, Nz) of the full dump, in logical axis order."""
        return self._shape

    @property
    def NumTimesteps(self):
        """Number of recorded timesteps (0 for an FD-only dump)."""
        return len(self._td_names)

    @property
    def NumFrequencies(self):
        """Number of recorded frequencies (0 for a TD-only dump)."""
        return len(self._fd_indices)

    @property
    def Frequencies(self):
        """The frequencies stored in the file, in Hz (empty for a TD-only dump)."""
        return self._frequencies

    @property
    def Times(self):
        """The simulation times of the recorded timesteps, in s."""
        if not self.IsTD:
            return np.array([])
        grp = self._h5['FieldData/TD']
        return np.array([self._TimeOf(grp, n) for n in self._td_names])

    @staticmethod
    def _TimeOf(grp, name):
        ds = grp[name]
        if 'time' not in ds.attrs:
            return np.nan
        return float(np.asarray(ds.attrs['time']).flatten()[0])

    ###########################################################################
    # region of interest
    ###########################################################################

    def _CheckDir(self, ny):
        """Resolve a direction to 0/1/2, accepting the mesh coordinate names."""
        if isinstance(ny, str):
            names = [n.lower() for n in self._mesh['names']]
            if ny.lower() in names:
                return names.index(ny.lower())
        try:
            return CheckNyDir(ny)
        except (ValueError, AssertionError):
            raise ValueError('Invalid direction `{}`, valid directions are '
                             '0/1/2 or {}'.format(ny, tuple(self._mesh['names'])))

    def ResetRegion(self, ny=None):
        """Drop plane/range/sampling settings and read the full dump again.

        Parameters
        ----------
        ny : int or str, optional
            Reset only this direction, given as 0/1/2 or a coordinate name.
            By default all three directions are reset.
        """
        # per direction the region is one of:
        #   None           -- the full extent
        #   int            -- collapsed to a single line (a plane or line)
        #   (start, stop)  -- a range of lines, half open
        if ny is None:
            self._region   = [None, None, None]
            self._sampling = [1, 1, 1]
        else:
            ny = self._CheckDir(ny)
            self._region[ny]   = None
            self._sampling[ny] = 1

    def _CheckNoConflict(self, ny, wanted):
        """Refuse to silently replace a plane by a range or vice versa."""
        current = self._region[ny]
        if current is None:
            return
        has = 'range' if isinstance(current, tuple) else 'plane'
        if has == wanted:
            return          # replacing a setting of the same kind is fine
        raise ValueError(
            'direction "{}" already has a {} set; call ResetRegion("{}") first '
            'if you want a {} instead'.format(self._mesh['names'][ny], has,
                                              self._mesh['names'][ny], wanted))

    def NearestIndex(self, ny, coord):
        """Index of the mesh line closest to *coord* along direction *ny*.

        Parameters
        ----------
        ny : int or str
            Direction, as 0/1/2 or a coordinate name ('x', 'rho', ...).
        coord : float
            Coordinate in SI units, i.e. metres for lengths and radians for
            angles -- the same units as ``GetMesh()['lines']``.
        """
        ny = self._CheckDir(ny)
        return int(np.argmin(np.abs(self._mesh['lines'][ny] - coord)))

    def SetPlane(self, ny, pos=None, idx=None):
        """Restrict all field reads to a single plane normal to *ny*.

        The plane axis collapses, so a 3D vector dump then yields data of
        shape (3, N1, N2) instead of (3, Nx, Ny, Nz).

        Parameters
        ----------
        ny : int or str
            Normal direction, as 0/1/2 or a coordinate name ('x', 'rho', ...).
        pos : float, optional
            Position in SI units; the nearest mesh line is used.
        idx : int, optional
            Mesh line index.  Exactly one of `pos` and `idx` must be given.

        There is only ever one plane: setting a plane for a different
        direction replaces the previous one, so switching the slice
        orientation needs no reset.  Ranges and sampling on the other
        directions are kept.  It raises if *this* direction already has a
        range set; use ``ResetRegion(ny)`` to clear it first.
        """
        ny = self._CheckDir(ny)
        self._CheckNoConflict(ny, 'plane')
        if (pos is None) == (idx is None):
            raise ValueError('SetPlane: give exactly one of `pos` or `idx`')
        if pos is not None:
            idx = self.NearestIndex(ny, pos)
        n_max = self._shape[ny]
        if not -n_max <= idx < n_max:
            raise IndexError('SetPlane: index {} out of range for direction {} '
                             'with {} lines'.format(idx, ny, n_max))
        # there is only ever one plane: drop a plane set on another direction,
        # ranges and sampling on the other directions are kept
        for n in range(3):
            if n != ny and isinstance(self._region[n], int):
                self._region[n] = None
        self._region[ny] = int(idx) % n_max

    def SetLine(self, ny, pos=None, idx=None):
        """Restrict all field reads to a single line along direction *ny*.

        The two directions perpendicular to *ny* collapse, so a 3D vector dump
        yields data of shape (3, N) instead of (3, Nx, Ny, Nz).

        Parameters
        ----------
        ny : int or str
            Direction the line runs along, as 0/1/2 or a coordinate name.
        pos : sequence of three floats, optional
            Position in SI units for each direction; the nearest mesh line is
            used.  The entry for the line direction *ny* must be None.
        idx : sequence of three ints, optional
            The same as mesh line indices.  Exactly one of `pos` and `idx`
            must be given.

        The entry for each direction is given by its **position in the
        sequence**, so there is no ambiguity about the order:
        ``SetLine('y', idx=(3, None, 5))`` is a line along y at x-index 3 and
        z-index 5.

        The two perpendicular directions always collapse, replacing whatever
        was set for them.  For the line direction itself, a setting that would
        leave fewer than two lines -- a plane normal, or a previous line -- is
        reset to the full extent, while an existing range of two or more lines
        is kept.  Sampling is never changed.

        Examples
        --------
        >>> dump.SetLine('x', idx=(None, 3, 5))       # along x at y=3, z=5
        >>> dump.SetLine('z', pos=(0.0, 1e-3, None))  # along z at x=0, y=1 mm
        """
        ny = self._CheckDir(ny)
        if (pos is None) == (idx is None):
            raise ValueError('SetLine: give exactly one of `pos` or `idx`')
        vals = list(pos if pos is not None else idx)
        names = self._mesh['names']
        if len(vals) != 3:
            raise ValueError('SetLine: `{}` must have three entries, one per '
                             'direction {}, with None for the line direction '
                             '"{}"'.format('pos' if pos is not None else 'idx',
                                           tuple(names), names[ny]))
        if vals[ny] is not None:
            raise ValueError('SetLine: entry {} is the line direction "{}" and '
                             'must be None, got {!r}'.format(ny, names[ny],
                                                             vals[ny]))
        for n in range(3):
            if n == ny:
                continue
            if vals[n] is None:
                raise ValueError('SetLine: no position given for direction '
                                 '"{}"'.format(names[n]))
            i = self.NearestIndex(n, vals[n]) if pos is not None else vals[n]
            n_max = self._shape[n]
            if not -n_max <= i < n_max:
                raise IndexError('SetLine: index {} out of range for direction '
                                 '"{}" with {} lines'.format(i, names[n], n_max))
            self._region[n] = int(i) % n_max
        # the line direction has to extend: drop a setting that would leave a
        # single line (a former plane normal or line), keep a real range
        if self._SelectionLength(ny) < 2:
            self._region[ny] = None

    def SetRange(self, ny, start=None, stop=None, idx_start=None, idx_stop=None):
        """Restrict all field reads to a sub-range along direction *ny*.

        Parameters
        ----------
        ny : int or str
            Direction, as 0/1/2 or a coordinate name ('x', 'rho', ...).
        start, stop : float, optional
            Range limits in SI units.  The nearest mesh lines are used and
            both limits are *inclusive*.
        idx_start, idx_stop : int, optional
            Range limits as mesh line indices, half open as usual in Python
            (`idx_stop` is not included).  Cannot be combined with
            `start`/`stop`.

        Calling `SetRange` again for the same direction replaces the range.  It
        raises if that direction already has a plane set; use
        ``ResetRegion(ny)`` to clear it first.
        """
        ny = self._CheckDir(ny)
        self._CheckNoConflict(ny, 'range')
        by_coord = start is not None or stop is not None
        by_index = idx_start is not None or idx_stop is not None
        if by_coord and by_index:
            raise ValueError('SetRange: give either `start`/`stop` or '
                             '`idx_start`/`idx_stop`, not both')
        if not by_coord and not by_index:
            raise ValueError('SetRange: no range given')
        if by_coord:
            i0 = 0 if start is None else self.NearestIndex(ny, start)
            i1 = self._shape[ny] if stop is None else self.NearestIndex(ny, stop) + 1
        else:
            i0, i1 = idx_start, idx_stop
        self._region[ny] = (i0, i1)

    def SetSampling(self, *factors):
        """Sub-sample the field data by the given step in each direction.

        Accepts either a single factor applied to all three directions, or
        one factor per direction, e.g. ``SetSampling(2)`` or
        ``SetSampling(2, 2, 1)``.  Sub-sampling is applied by HDF5 while
        reading, it does not interpolate.
        """
        if len(factors) == 1:
            factors = factors * 3
        if len(factors) != 3:
            raise ValueError('SetSampling: give one factor, or one per direction')
        if any(int(f) < 1 for f in factors):
            raise ValueError('SetSampling: factors must be >= 1')
        self._sampling = [int(f) for f in factors]

    def GetMesh(self, region=False):
        """Return the mesh of the dump.

        Parameters
        ----------
        region : bool, optional
            If True, return only the mesh lines covered by the currently
            configured region and sampling, i.e. the lines matching the field
            data returned by the getters.  A plane direction then yields a
            single-element array.

        Returns
        -------
        mesh : dict
            * `lines`   : list of the three coordinate vectors, in SI units,
              i.e. metres for lengths and radians for angles
            * `names`   : the corresponding coordinate names, e.g. ``['x','y','z']``
            * `type`    : 0 --> Cartesian, 1 --> cylindrical, 2 --> spherical
            * `scaling` : the simulation length unit in metres (e.g. 1e-3 for a
              mm mesh); divide `lines` by it to get the drawing units
        """
        mesh = dict(self._mesh)
        mesh['lines'] = list(self._mesh['lines'])
        if region:
            sel = self._SpatialSelection()
            mesh['lines'] = [np.atleast_1d(l[s]) for l, s in zip(mesh['lines'], sel)]
        return mesh

    ###########################################################################
    # reading field data
    ###########################################################################

    def _SpatialSelection(self):
        """The configured region as a tuple of three slices/ints, x-y-z order."""
        sel = []
        for n in range(3):
            r, step = self._region[n], self._sampling[n]
            if isinstance(r, tuple):
                sel.append(slice(r[0], r[1], step))
            elif r is None:
                sel.append(slice(None, None, step))
            else:
                sel.append(r)
        return tuple(sel)

    def _SelectionLength(self, ny):
        """Number of lines direction *ny* currently contributes."""
        sel = self._SpatialSelection()[ny]
        if not isinstance(sel, slice):
            return 1
        return len(range(*sel.indices(self._shape[ny])))

    def _BuildIndex(self, component):
        """Map the configured region onto the on-disk axis order."""
        sel = list(self._SpatialSelection())
        if self._legacy:                        # on disk as (..,Nz,Ny,Nx)
            sel = sel[::-1]
        if self._is_vector:
            sel = [slice(None) if component is None else component] + sel
        return tuple(sel)

    def _FixOrder(self, data, index):
        """Reverse the spatial axes of a legacy array back to x-inner order."""
        if not self._legacy or data.ndim < 2:
            return data
        # the component axis survives only if it was not indexed with an int
        offset = 1 if (self._is_vector and isinstance(index[0], slice)) else 0
        perm = list(range(offset)) + list(range(data.ndim - 1, offset - 1, -1))
        return np.transpose(data, perm)

    def _ReadDataset(self, grp, name, component):
        index = self._BuildIndex(component)
        if name in grp:
            # h5py maps the compound {r,i} type of a complex FD dump onto a
            # native complex array by itself, see h5py's complex_names config
            data = grp[name][index]
        elif name + '_real' in grp:
            # legacy format: real and imaginary part in two separate datasets.
            # Combined by hand rather than via `real + 1j*imag`, which would
            # promote float32 data to complex128.
            real = grp[name + '_real'][index]
            imag = grp[name + '_imag'][index]
            data = np.empty(real.shape, dtype=np.complex128
                            if real.dtype == np.float64 else np.complex64)
            data.real = real
            data.imag = imag
        else:
            raise KeyError('"{}" does not contain the dataset {}/{}'.format(
                self._h5.filename, grp.name, name))
        return self._FixOrder(data, index)

    def _TDName(self, t_idx):
        if isinstance(t_idx, str):
            return t_idx
        if not -len(self._td_names) <= t_idx < len(self._td_names):
            raise IndexError('timestep index {} out of range, "{}" holds {} '
                             'timesteps'.format(t_idx, self._h5.filename,
                                                len(self._td_names)))
        return self._td_names[t_idx]

    def _FDName(self, f_idx):
        if isinstance(f_idx, str):
            return f_idx
        if not -len(self._fd_indices) <= f_idx < len(self._fd_indices):
            raise IndexError('frequency index {} out of range, "{}" holds {} '
                             'frequencies'.format(f_idx, self._h5.filename,
                                                  len(self._fd_indices)))
        return 'f{}'.format(self._fd_indices[f_idx])

    def _Group(self, domain):
        path = 'FieldData/' + domain
        if path not in self._h5:
            raise KeyError('"{}" does not contain any /{} data'.format(
                self._h5.filename, path))
        return self._h5[path]

    def GetFieldAtIndex(self, f_idx=None, t_idx=None, component=None):
        """Read one field sample, addressed by its index in the file.

        Parameters
        ----------
        f_idx : int or str, optional
            Frequency index, or a dataset name such as ``'f0'``.
        t_idx : int or str, optional
            Position of the timestep in the file, or a dataset name such as
            ``'000100'``.  Exactly one of `f_idx` and `t_idx` must be given.
        component : int or str, optional
            Read only this vector component; by default all three are read.

        Returns
        -------
        data : ndarray
            Shape (3, Nx, Ny, Nz) for a vector dump, (Nx, Ny, Nz) for a scalar
            one, reduced by the configured region, sampling and `component`.
            Frequency domain data is complex, time domain data is real.
        """
        if (f_idx is None) == (t_idx is None):
            raise ValueError('give exactly one of `f_idx` or `t_idx`')
        if component is not None:
            component = self._CheckDir(component)
        if f_idx is not None:
            grp = self._Group('FD')
            return self._ReadDataset(grp, self._FDName(f_idx), component)
        grp = self._Group('TD')
        return self._ReadDataset(grp, self._TDName(t_idx), component)

    def GetFieldAtFrequency(self, freq, component=None):
        """Read the field at a given frequency, in Hz.

        If the file holds frequency domain data, the matching dataset is read
        directly; the frequency must match one stored in the file to within
        `FREQ_RTOL`.  Otherwise the frequency is computed from the time domain
        data by an on-the-fly DFT, which reads every timestep but holds only
        one of them in memory at a time.

        Parameters
        ----------
        freq : float
            Frequency in Hz.
        component : int or str, optional
            Read only this vector component; by default all three are read.
        """
        if self.IsFD:
            f_idx = self._MatchFrequency(freq)
            if f_idx is not None:
                return self.GetFieldAtIndex(f_idx=f_idx, component=component)
            if not self.IsTD:
                raise ValueError(
                    '{:g} Hz is not stored in "{}", available frequencies are '
                    '{}'.format(freq, self._h5.filename, self._frequencies))
        return self._DFT(freq, component)

    def _MatchFrequency(self, freq):
        """Index of the stored frequency matching *freq*, or None."""
        if not self.NumFrequencies:
            return None
        n = int(np.argmin(np.abs(self._frequencies - freq)))
        if np.isclose(self._frequencies[n], freq, rtol=self.FREQ_RTOL, atol=0.0):
            return n
        return None

    def _DFT(self, freq, component):
        """Single frequency DFT of the time domain data, one timestep at a time."""
        if not self.IsTD:
            raise ValueError('"{}" contains no time domain data to transform'.format(
                self._h5.filename))
        grp = self._Group('TD')
        times = self.Times
        if len(times) < 2:
            raise ValueError('need at least two timesteps for a DFT, "{}" has '
                             '{}'.format(self._h5.filename, len(times)))
        dt = times[1] - times[0]
        accum = None
        for name, t in zip(self._td_names, times):
            data = self._ReadDataset(grp, name, component)
            if accum is None:
                c_type = np.complex64 if data.dtype == np.float32 else np.complex128
                accum = np.zeros(data.shape, dtype=c_type)
            accum += data * np.exp(-2j*np.pi*freq*t)
        return 2 * accum * dt     # single-sided spectrum

    def IterFD(self, component=None):
        """Iterate over all frequency domain samples.

        Yields
        ------
        (freq, data) : tuple of float and ndarray
            The frequency in Hz and the field data, see `GetFieldAtIndex`.
        """
        for n in range(self.NumFrequencies):
            yield float(self._frequencies[n]), \
                  self.GetFieldAtIndex(f_idx=n, component=component)

    def IterTD(self, component=None):
        """Iterate over all recorded timesteps.

        Yields
        ------
        (time, data) : tuple of float and ndarray
            The simulation time in s and the field data, see `GetFieldAtIndex`.
        """
        grp = self._Group('TD')
        if component is not None:
            component = self._CheckDir(component)
        for name in self._td_names:
            yield self._TimeOf(grp, name), self._ReadDataset(grp, name, component)

    def GetAttributes(self, f_idx=None, t_idx=None):
        """Collect the attributes applying to one sample, without reading it.

        Attributes are taken from three levels in order of increasing
        precedence: the file root, the ``/FieldData/{FD,TD}`` group and the
        dataset itself.  Later levels overwrite earlier ones when the same key
        appears at several levels; in particular the `frequency` array of the
        FD group is overwritten by the scalar `frequency` of the dataset.
        """
        if (f_idx is None) == (t_idx is None):
            raise ValueError('give exactly one of `f_idx` or `t_idx`')
        domain = 'FD' if f_idx is not None else 'TD'
        grp = self._Group(domain)
        name = self._FDName(f_idx) if f_idx is not None else self._TDName(t_idx)
        if name not in grp:
            name = name + '_real'
        if name not in grp:
            raise KeyError('"{}" does not contain the dataset /FieldData/{}/{}'.format(
                self._h5.filename, domain, name))
        attrs = dict(self._root_attrs)
        attrs.update(grp.attrs)
        attrs.update(grp[name].attrs)
        return attrs


if __name__=="__main__":
    import pylab as plt

    t = np.linspace(0,2,201)

    s = np.sin(2*np.pi*2*t)
    plt.plot(t,s)

    f = np.linspace(0,3,101)
    sf = DFT_time2freq(t, s, f, 'periodic')

    plt.figure()
    plt.plot(f, np.abs(sf))

    plt.show()


