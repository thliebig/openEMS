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
from __future__ import absolute_import

import os
import numpy as np
from CSXCAD.Utilities import CheckNyDir
from CSXCAD.CSRectGrid import CoordinateSystem
from openEMS import utilities

from openEMS.physical_constants import *

def _load_ui_file(filepath):
    """Read an openEMS probe output file in a single pass.

    Returns (data, col_names) where data is a float64 array of shape (N, ncols)
    and col_names is the list of column-name strings from the last header line
    (e.g. ['t/s', 'voltage', 'mode_purity']), or None if no header was found.
    """
    comments = []
    rows = []
    with open(filepath) as f:
        for line in f:
            if line.startswith('%'):
                comments.append(line[1:].strip())
            else:
                s = line.strip()
                if s:
                    rows.append(s.split())
    col_names = comments[-1].split() if comments else None
    if col_names is not None and col_names[0] != 't/s':
        raise ValueError('{}: unexpected first column "{}", expected "t/s"'.format(filepath, col_names[0]))
    return np.array(rows, dtype=np.double), col_names


class UI_data:
    def __init__(self, fns, path, freq, signal_type='pulse', **kw):
        self.path = path
        if type(fns)==str:
            fns = [fns]
        self.fns  = fns

        if np.isscalar(freq):
            freq = [freq]
        self.freq = freq

        self.ui_time        = []
        self.ui_val         = []
        self.ui_f_val       = []
        self.col_names      = []  # column names per file, from the file header
        self.ui_mode_purity = []  # mode-purity time series or None

        for fn in fns:
            data, col_names = _load_ui_file(os.path.join(path, fn))
            self.col_names.append(col_names)
            self.ui_time.append(data[:, 0])
            self.ui_val.append(data[:, 1])
            self.ui_f_val.append(utilities.DFT_time2freq(data[:, 0], data[:, 1], freq, signal_type=signal_type))
            has_purity = (col_names is not None and len(col_names) > 1
                          and col_names[-1] == 'mode_purity')
            self.ui_mode_purity.append(data[:, -1] if has_purity else None)

# Port Base-Class
class Port(object):
    """
    The port base class.

    :param CSX: Continuous Structure
    :param port_nr: int -- port number
    :param R: float -- port reference impedance, e.g. 50 (Ohms)
    :param start, stop: (3,) array -- Start/Stop box coordinates
    :param p_dir: int -- port direction
    :param excite: float -- port excitation amplitude
    :param priority: int -- priority of all contained primtives
    :param PortNamePrefix: str -- a prefix for all ports-names
    :param delay: float -- a positive delay value to e.g. emulate a phase shift
    """
    def __init__(self, CSX, port_nr, start, stop, excite, **kw):
        self.CSX      = CSX
        self.number   = port_nr
        self.excite   = excite
        self.start    = np.array(start, np.double)
        self.stop     = np.array(stop, np.double)
        self.Z_ref    = None
        self.U_filenames = kw.get('U_filenames', [])
        self.I_filenames = kw.get('I_filenames', [])
        self.port_props = []

        self.priority = 0
        if 'priority' in kw:
            self.priority = kw['priority']

        self.prefix = ''
        if 'PortNamePrefix' in kw:
            self.prefix = kw['PortNamePrefix']
        self.delay = 0

        if 'delay' in kw:
            self.delay = kw['delay']

        self.lbl_temp = self.prefix + 'port_{}' +  '_{}'.format(self.number)

    def SetEnabled(self, val):
        from CSXCAD.CSProperties import CSPropExcitation
        found_any = False
        for prop in self.port_props:
            if type(prop) == CSPropExcitation:
                prop.SetEnabled(val)
                found_any = True
        if not found_any and val:
            # if we attempt to activate this port and it does not have any excitation set, raise an exception!
            raise Exception('Unable to enable port! No excitation found!')

    def ReadUIData(self, sim_path, freq, signal_type ='pulse'):
        self.u_data = UI_data(self.U_filenames, sim_path, freq, signal_type )
        self.uf_tot = 0
        self.ut_tot = 0
        for n in range(len(self.u_data.fns)):
            self.uf_tot += self.u_data.ui_f_val[n]
            self.ut_tot += self.u_data.ui_val[n]
        self.u_time = self.u_data.ui_time[0]

        self.i_data = UI_data(self.I_filenames, sim_path, freq, signal_type )
        self.if_tot = 0
        self.it_tot = 0
        for n in range(len(self.i_data.fns)):
            self.if_tot += self.i_data.ui_f_val[n]
            self.it_tot += self.i_data.ui_val[n]
        self.i_time = self.i_data.ui_time[0]

        # mode purity: extra column written by ProcessModeMatch (index 0 = purity)
        self.u_mode_purity = self.u_data.ui_mode_purity
        self.i_mode_purity = self.i_data.ui_mode_purity


    def CalcPort(self, sim_path, freq, ref_impedance=None, ref_plane_shift=None, signal_type='pulse'):
        self.ReadUIData(sim_path, freq, signal_type)

        if ref_impedance is not None:
            self.Z_ref = ref_impedance
        if self.Z_ref is None:
            raise Exception('Port Z_ref should not be None!')

        if ref_plane_shift is not None:
            if not hasattr(self, 'beta'):
                raise Exception('Port has no beta attribute!')
            shift = ref_plane_shift
            if self.measplane_shift:
                shift -= self.measplane_shift
            shift *= self.CSX.GetGrid().GetDeltaUnit()
            phase = np.real(self.beta)*shift
            uf_tot = self.uf_tot * np.cos(-phase) + 1j * self.if_tot * self.Z_ref * np.sin(-phase)
            if_tot = self.if_tot * np.cos(-phase) + 1j * self.uf_tot / self.Z_ref * np.sin(-phase)
            self.uf_tot = uf_tot
            self.if_tot = if_tot

        self.uf_inc = 0.5 * ( self.uf_tot + self.if_tot * self.Z_ref )
        self.if_inc = 0.5 * ( self.if_tot + self.uf_tot / self.Z_ref )
        self.uf_ref = self.uf_tot - self.uf_inc
        self.if_ref = self.if_inc - self.if_tot

        if type(self.Z_ref) in [int, float]:
            self.ut_inc = 0.5 * ( self.ut_tot + self.it_tot * self.Z_ref )
            self.it_inc = 0.5 * ( self.it_tot + self.ut_tot / self.Z_ref )
            self.ut_ref = self.ut_tot - self.ut_inc
            self.it_ref = self.it_inc - self.it_tot

        # calc some more port parameter
        # incoming power
        self.P_inc = 0.5*np.real(self.uf_inc*np.conj(self.if_inc))
        # reflected power
        self.P_ref = 0.5*np.real(self.uf_ref*np.conj(self.if_ref))
        # accepted power (incoming - reflected)
        self.P_acc = 0.5*np.real(self.uf_tot*np.conj(self.if_tot))

class LumpedPort(Port):
    """
    The lumped port.

    See Also
    --------
    Port
    """
    def __init__(self, CSX,  port_nr, R, start, stop, exc_dir, excite=0, **kw):
        super(LumpedPort, self).__init__(CSX, port_nr=port_nr, start=start, stop=stop, excite=excite, **kw)
        self.R = R
        self.exc_ny  = CheckNyDir(exc_dir)

        self.direction = np.sign(self.stop[self.exc_ny]-self.start[self.exc_ny])
        if not self.start[self.exc_ny]!=self.stop[self.exc_ny]:
            raise Exception('LumpedPort: start and stop may not be identical in excitation direction')

        if self.R > 0:
            lumped_R = CSX.AddLumpedElement(self.lbl_temp.format('resist'), ny=self.exc_ny, caps=True, R=self.R)
        elif self.R==0:
            lumped_R = CSX.AddMetal(self.lbl_temp.format('resist'))
        self.port_props.append(lumped_R)
        lumped_R.AddBox(self.start, self.stop, priority=self.priority)

        if excite!=0:
            exc_vec = np.zeros(3)
            exc_vec[self.exc_ny] = -1*self.direction*excite
            exc = CSX.AddExcitation(self.lbl_temp.format('excite'), exc_type=0, exc_val=exc_vec, delay=self.delay)
            exc.AddBox(self.start, self.stop, priority=self.priority)
            self.port_props.append(exc)

        self.U_filenames = [self.lbl_temp.format('ut'), ]
        u_start = 0.5*(self.start+self.stop)
        u_start[self.exc_ny] = self.start[self.exc_ny]
        u_stop  = 0.5*(self.start+self.stop)
        u_stop[self.exc_ny]  = self.stop[self.exc_ny]
        u_probe = CSX.AddProbe(self.U_filenames[0], p_type=0, weight=-1)
        u_probe.AddBox(u_start, u_stop)
        self.port_props.append(u_probe)

        self.I_filenames = [self.lbl_temp.format('it'), ]
        i_start = np.array(self.start)
        i_start[self.exc_ny] = 0.5*(self.start[self.exc_ny]+self.stop[self.exc_ny])
        i_stop  = np.array(self.stop)
        i_stop[self.exc_ny]  = 0.5*(self.start[self.exc_ny]+self.stop[self.exc_ny])
        i_probe = CSX.AddProbe(self.I_filenames[0], p_type=1, weight=self.direction, norm_dir=self.exc_ny)
        i_probe.AddBox(i_start, i_stop)
        self.port_props.append(i_probe)

    def CalcPort(self, sim_path, freq, ref_impedance=None, ref_plane_shift=None, signal_type='pulse'):
        if ref_impedance is None:
            self.Z_ref = self.R
        if ref_plane_shift is not None:
            Warning('A lumped port does not support a reference plane shift! Ignoring...')
        super(LumpedPort, self).CalcPort(sim_path, freq, ref_impedance, ref_plane_shift, signal_type)

class MSLPort(Port):
    """
    The microstrip transmission line port.

    :param prop_dir: int/str -- direction of propagation

    See Also
    --------
    Port
    """
    def __init__(self, CSX, port_nr, metal_prop, start, stop, prop_dir, exc_dir, excite=0, **kw):
        super(MSLPort, self).__init__(CSX, port_nr=port_nr, start=start, stop=stop, excite=excite, **kw)
        self.exc_ny  = CheckNyDir(exc_dir)
        self.prop_ny = CheckNyDir(prop_dir)
        self.direction   = np.sign(stop[self.prop_ny]-start[self.prop_ny])
        self.upside_down = np.sign(stop[self.exc_ny]  -start[self.exc_ny])
        if not (self.start!=self.stop).all():
            raise Exception('Start coordinate must not be equal to stop coordinate')
#        assert stop[self.prop_ny]!=start[self.prop_ny], 'port length in propergation direction may not be zero!'
#        assert stop[self.exc_ny] !=start[self.exc_ny], 'port length in propergation direction may not be zero!'
        if not self.exc_ny!=self.prop_ny:
            raise Exception('Excitation direction must not be equal to propagation direction')

        self.feed_shift = 0
        if 'FeedShift' in kw:
            self.feed_shift = kw['FeedShift']
        self.measplane_shift = 0.5*np.abs(self.start[self.prop_ny]-self.stop[self.prop_ny])
        if 'MeasPlaneShift' in kw:
            self.measplane_shift =  kw['MeasPlaneShift']
        self.measplane_pos = self.start[self.prop_ny] + self.measplane_shift*self.direction
        self.feed_R = np.inf
        if 'Feed_R' in kw:
            self.feed_R = kw['Feed_R']

        # add metal msl-plane
        MSL_start = np.array(self.start)
        MSL_stop  = np.array(self.stop)
        MSL_stop[self.exc_ny] = MSL_start[self.exc_ny]
        metal_prop.AddBox(MSL_start, MSL_stop, priority=self.priority )

        mesh = CSX.GetGrid()
        prop_lines = mesh.GetLines(self.prop_ny)
        if not len(prop_lines)>5:
            raise Exception('At least 5 lines in propagation direction required!')
        meas_pos_idx = np.argmin(np.abs(prop_lines-self.measplane_pos))
        if meas_pos_idx==0:
            meas_pos_idx=1
        if meas_pos_idx>=len(prop_lines)-1:
            meas_pos_idx=len(prop_lines)-2
        self.measplane_shift = np.abs(self.start[self.prop_ny]-prop_lines[meas_pos_idx])
        prope_idx = np.array([meas_pos_idx-1, meas_pos_idx, meas_pos_idx+1], int)
        if self.direction<0:
            prope_idx = np.flipud(prope_idx)
        u_prope_pos = prop_lines[prope_idx]
        self.U_filenames = []
        self.U_delta = np.diff(u_prope_pos)
        suffix = ['A', 'B', 'C']
        for n in range(len(prope_idx)):
            u_start = 0.5*(self.start+self.stop)
            u_stop  = 0.5*(self.start+self.stop)
            u_start[self.prop_ny] = u_prope_pos[n]
            u_stop[self.prop_ny]  = u_prope_pos[n]
            u_start[self.exc_ny]  = self.start[self.exc_ny]
            u_stop[self.exc_ny]   = self.stop [self.exc_ny]
            u_name = self.lbl_temp.format('ut') + suffix[n]
            self.U_filenames.append(u_name)
            u_probe = CSX.AddProbe(u_name, p_type=0)
            u_probe.AddBox(u_start, u_stop)
            self.port_props.append(u_probe)

        i_prope_pos = u_prope_pos[0:2] + np.diff(u_prope_pos)/2.0
        self.I_filenames = []
        self.I_delta = np.diff(i_prope_pos)
        i_start = np.array(self.start)
        i_stop  = np.array(self.stop)
        i_stop[self.exc_ny] = self.start[self.exc_ny]
        for n in range(len(i_prope_pos)):
            i_start[self.prop_ny] = i_prope_pos[n]
            i_stop[self.prop_ny]  = i_prope_pos[n]
            i_name = self.lbl_temp.format('it') + suffix[n]
            self.I_filenames.append(i_name)
            i_probe = CSX.AddProbe(i_name, p_type=1, weight=self.direction, norm_dir=self.prop_ny)
            i_probe.AddBox(i_start, i_stop)
            self.port_props.append(i_probe)

        if excite!=0:
            excide_pos_idx = np.argmin(np.abs(prop_lines-(self.start[self.prop_ny] + self.feed_shift*self.direction)))
            exc_start = np.array(self.start)
            exc_stop  = np.array(self.stop)
            exc_start[self.prop_ny] = prop_lines[excide_pos_idx]
            exc_stop [self.prop_ny] = prop_lines[excide_pos_idx]
            exc_vec = np.zeros(3)
            exc_vec[self.exc_ny] = -1*self.upside_down*excite
            exc = CSX.AddExcitation(self.lbl_temp.format('excite'), exc_type=0, exc_val=exc_vec, delay=self.delay)
            exc.AddBox(exc_start, exc_stop, priority=self.priority)
            self.port_props.append(exc)

        if self.feed_R>=0 and not np.isinf(self.feed_R):
            R_start = np.array(self.start)
            R_stop  = np.array(self.stop)
            R_stop [self.prop_ny] = R_start[self.prop_ny]
            if self.feed_R==0:
                metal_prop.AddBox(R_start, R_stop)
            else:
                lumped_R = CSX.AddLumpedElement(self.lbl_temp.format('resist'), ny=self.exc_ny, caps=True, R=self.feed_R)
                lumped_R.AddBox(R_start, R_stop)
                self.port_props.append(lumped_R)

    def ReadUIData(self, sim_path, freq, signal_type ='pulse'):
        self.u_data = UI_data(self.U_filenames, sim_path, freq, signal_type )
        self.uf_tot = self.u_data.ui_f_val[1]
        self.ut_tot = self.u_data.ui_val[1]
        self.u_time = self.u_data.ui_time[0]

        self.i_data = UI_data(self.I_filenames, sim_path, freq, signal_type )
        self.if_tot = 0.5*(self.i_data.ui_f_val[0]+self.i_data.ui_f_val[1])
        self.it_tot = 0.5*(self.i_data.ui_val[0]+self.i_data.ui_val[1])
        self.i_time = self.i_data.ui_time[0]

        unit = self.CSX.GetGrid().GetDeltaUnit()
        Et = self.u_data.ui_f_val[1]
        dEt = (self.u_data.ui_f_val[2] - self.u_data.ui_f_val[0]) / (np.sum(np.abs(self.U_delta)) * unit)
        Ht = self.if_tot # space averaging: Ht is now defined at the same pos as Et
        dHt = (self.i_data.ui_f_val[1] - self.i_data.ui_f_val[0]) / (np.abs(self.I_delta[0]) * unit)

        beta = np.sqrt( - dEt * dHt / (Ht * Et) )
        beta[np.real(beta) < 0] *= -1 # determine correct sign (unlike the paper)
        self.beta = beta

        # determine ZL
        self.Z_ref = np.sqrt(Et * dEt / (Ht * dHt))

class WaveguidePort(Port):
    """
    Base class for any waveguide port.

    The mode shape can be supplied either as analytic expressions (E_WG_func /
    H_WG_func) or as an HDF5 mode file (E_WG_file / H_WG_file).  Exactly one
    of the two sources must be provided.

    Parameters
    ----------
    exc_dir : int or str
        Propagation direction of the waveguide (0/'x', 1/'y', 2/'z').
    E_WG_func : list of str or None
        Electric field mode profile as a list of three fparser expressions,
        one per Cartesian component.  Use ``None`` when supplying a mode file.
    H_WG_func : list of str or None
        Magnetic field mode profile as a list of three fparser expressions.
        Use ``None`` when supplying a mode file.
    kc : float
        Cut-off wavenumber of the mode in drawing units (e.g. pi/a for TE10).
        Used by :meth:`CalcPort` to compute the propagation constant beta and
        the analytic waveguide impedance.
    E_WG_file : str or None
        Path to an HDF5 file containing the electric field mode profile.
        Required when E_WG_func / H_WG_func are ``None``.
    H_WG_file : str or None
        Path to an HDF5 file containing the magnetic field mode profile.
        Required when E_WG_func / H_WG_func are ``None``.
    local_origin : array-like, 'corner', 'center', or None
        Coordinate origin used when evaluating mode functions or looking up
        mode file data, given as a **Cartesian** point for any mesh type.
        Coordinates are converted to Cartesian and shifted by this amount, and
        only then are ``x``/``y``/``z``/``rho``/``a``/``r``/``t`` derived, so
        that ``rho`` and ``a`` are measured from this point.  Both the
        excitation and the two mode-match probes receive the same origin, which
        is what keeps the excited and the probed mode identical.  Shorthands:
        ``'corner'`` resolves to ``min(start, stop)`` per axis; ``'center'``
        resolves to the midpoint.  Because they derive from ``start``/``stop``,
        which are given in mesh coordinates, the shorthands are Cartesian-mesh
        only -- on a cylindrical mesh pass an explicit vector (``[0, 0, z]``
        keeps the mode centered on the mesh axis).  Default ``None`` means no
        shift (global coordinates are used as-is).

    See Also
    --------
    Port, RectWGPort

    """
    def __init__(self, CSX, port_nr, start, stop, exc_dir, E_WG_func, H_WG_func, kc, excite = 0, excite_type = 0, E_WG_file = None, H_WG_file = None, local_origin = None, **kw):
        super(WaveguidePort, self).__init__(CSX, port_nr=port_nr, start=start, stop=stop, excite=excite, excite_type=excite_type, **kw)
        self.exc_ny  = CheckNyDir(exc_dir)
        self.ny_P  = (self.exc_ny+1)%3
        self.ny_PP = (self.exc_ny+2)%3
        self.direction = np.sign(stop[self.exc_ny]-start[self.exc_ny])
        self.ref_index = 1

        if (self.excite!=0 and stop[self.exc_ny]==start[self.exc_ny]):
            raise Exception('Port length in excitation direction may not be zero if port is excited!')

        self.kc = kc
        self.E_func = E_WG_func
        self.H_func = H_WG_func
        self.E_file = E_WG_file
        self.H_file = H_WG_file

        # Resolve local_origin shorthand into a Cartesian coordinate array (or
        # keep None). start/stop are in mesh coordinates, so on a cylindrical
        # mesh their corner/midpoint is a (rho,a,z) triple rather than the
        # Cartesian point the origin has to be -- reject the shorthands there
        # instead of silently handing cylindrical components to the solver.
        if isinstance(local_origin, str):
            if CSX.GetGrid().GetMeshType() == CoordinateSystem.CYLINDRICAL:
                raise Exception("local_origin '{}' is not supported on a cylindrical mesh; "
                                'pass an explicit Cartesian [x,y,z] vector '
                                '(use [0,0,z] to stay on the mesh axis)'.format(local_origin))
            if local_origin == 'corner':
                local_origin = np.minimum(start, stop)
            elif local_origin == 'center':
                local_origin = 0.5 * (np.array(start) + np.array(stop))
            else:
                raise Exception("unknown local_origin value '{}'; use 'corner', "
                                "'center' or an [x,y,z] vector".format(local_origin))
        elif local_origin is not None:
            local_origin = np.array(local_origin, dtype=float)

        use_function_expr = (self.E_func is not None) and (self.H_func is not None)

        if excite != 0:
            e_start = np.array(start)
            e_stop  = np.array(stop)
            e_stop[self.exc_ny] = e_start[self.exc_ny]
            e_vec = excite*np.ones(3)
            e_vec[self.exc_ny] = 0
            exc = CSX.AddExcitation(self.lbl_temp.format('excite'), exc_type=excite_type, exc_val=e_vec, delay=self.delay)


            if not use_function_expr:
                if not ((type(self.E_file) is str) and (type(self.H_file) is str)):
                    raise Exception('Both E_WG_file and H_WG_file must be strings (HDF5 file paths)')

                # E-field (TE case)
                if excite_type in [0,1]:
                    exc.SetWeightFile(self.E_file)
                # H-field (TM case)
                elif excite_type in [2,3]:
                    exc.SetWeightFile(self.H_file)
                else:
                    raise Exception('Unsupported excitation type. Only 0 or 2 for WaveguidePort')
            else:
                if not (type(self.E_func) is list):
                    raise Exception('Unsupported input type for "E_WG_func" or "H_WG_func". Expected a list of strings')

                if excite_type in [0,1]:
                    exc.SetWeightFunction([str(x) for x in self.E_func])
                elif excite_type in [2,3]:
                    exc.SetWeightFunction([str(x) for x in self.H_func])
                else:
                    raise Exception('Unsupported excitation type. Only 0 or 2 for WaveguidePort')

            # For the mode file to be used correctly, the direction of
            # propagation has to be explicitly set.
            if not use_function_expr:
                dirVect = [0,0,0]
                dirVect[self.exc_ny] = 1
                exc.SetPropagationDir(dirVect)

            if local_origin is not None:
                exc.SetWeightOrigin(local_origin)

            # Finally, add the box
            exc.AddBox(e_start, e_stop, priority=self.priority)
            self.port_props.append(exc)


        # voltage/current planes
        m_start = np.array(start)
        m_stop  = np.array(stop)
        m_start[self.exc_ny] = m_stop[self.exc_ny]
        self.measplane_shift = np.abs(stop[self.exc_ny] - start[self.exc_ny])

        self.U_filenames = [self.lbl_temp.format('ut'), ]
        u_probe_kw = {'mode_function': self.E_func} if use_function_expr else {}
        u_probe = CSX.AddProbe(self.U_filenames[0], p_type=10, **u_probe_kw)
        if not use_function_expr:
            u_probe.SetModeFile(self.E_file)
        if local_origin is not None:
            u_probe.SetModeOrigin(local_origin)
        u_probe.AddBox(m_start, m_stop)
        self.port_props.append(u_probe)

        self.I_filenames = [self.lbl_temp.format('it'), ]
        i_probe_kw = {'mode_function': self.H_func} if use_function_expr else {}
        i_probe = CSX.AddProbe(self.I_filenames[0], p_type=11, weight=self.direction, **i_probe_kw)
        if not use_function_expr:
            i_probe.SetModeFile(self.H_file)
        if local_origin is not None:
            i_probe.SetModeOrigin(local_origin)
        i_probe.AddBox(m_start, m_stop)
        self.port_props.append(i_probe)

    def CalcPort(self, sim_path, freq, ref_impedance=None, ref_plane_shift=None, signal_type='pulse', ZL = -1):
        k = 2.0*np.pi*freq/C0*self.ref_index
        self.beta = np.sqrt(k**2 - self.kc**2)
        if ZL <= 0:
            self.ZL = k * Z0 / self.beta    #analytic waveguide impedance
        else:
            self.ZL = ZL
        if ref_impedance is None:
            self.Z_ref = self.ZL
        super(WaveguidePort, self).CalcPort(sim_path, freq, ref_impedance, ref_plane_shift, signal_type)

class RectWGPort(WaveguidePort):
    """
    Rectangular waveguide port with automatic TE mode profile generation.

    Constructs the analytic TE_mn mode functions (Pozar, Microwave Engineering)
    and forwards them to :class:`WaveguidePort`.  Only TE modes are currently
    supported.

    Parameters
    ----------
    a : float
        Waveguide width in metres (dimension along the first transverse axis).
    b : float
        Waveguide height in metres (dimension along the second transverse axis).
    mode_name : str
        Four-character mode string, e.g. ``'TE10'`` or ``'TE11'``.
    local_origin : array-like, 'corner', 'center', or None
        Defaults to ``'corner'`` so that the generated mode functions, which
        are defined relative to the lower-left corner of the port, are
        evaluated correctly regardless of where the port is placed in the mesh.
        See :class:`WaveguidePort` for the full description.

    See Also
    --------
    Port, WaveguidePort

    """
    def __init__(self, CSX, port_nr, start, stop, exc_dir, a, b, mode_name, excite=0, local_origin='corner', **kw):
        Port.__init__(self, CSX, port_nr, start, stop, excite=0, **kw)
        self.exc_ny  = CheckNyDir(exc_dir)
        self.ny_P  = (self.exc_ny+1)%3
        self.ny_PP = (self.exc_ny+2)%3
        self.WG_size = [a, b]

        self.WG_mode = mode_name
        if not len(self.WG_mode)==4:
            raise Exception('Invalid mode definition')
        self.unit = self.CSX.GetGrid().GetDeltaUnit()
        if self.WG_mode.startswith('TE'):
            self.TE = True
            self.TM = False
        else:
            self.TE = False
            self.TM = True
        self.M = float(self.WG_mode[2])
        self.N = float(self.WG_mode[3])

        if not self.TE:
            raise Exception('Currently only TE-modes are supported! Mode found: {}'.format(self.WG_mode))

        # values by David M. Pozar, Microwave Engineering, third edition
        a = self.WG_size[0]
        b = self.WG_size[1]

        xyz = 'xyz'
        name_P  = xyz[self.ny_P]
        name_PP = xyz[self.ny_PP]

        kc = np.sqrt((self.M*np.pi/a)**2 + (self.N*np.pi/b)**2)

        a /= self.unit
        b /= self.unit
        E_func = [0,0,0]
        H_func = [0,0,0]
        if self.N>0:
            E_func[self.ny_P]  = '{}*cos({}*{})*sin({}*{})'.format(self.N/b   , self.M*np.pi/a, name_P, self.N*np.pi/b, name_PP)
        if self.M>0:
            E_func[self.ny_PP] = '{}*sin({}*{})*cos({}*{})'.format(-1*self.M/a, self.M*np.pi/a, name_P, self.N*np.pi/b, name_PP)

        if self.M>0:
            H_func[self.ny_P]  = '{}*sin({}*{})*cos({}*{})'.format(self.M/a, self.M*np.pi/a, name_P, self.N*np.pi/b, name_PP)
        if self.N>0:
            H_func[self.ny_PP] = '{}*cos({}*{})*sin({}*{})'.format(self.N/b, self.M*np.pi/a, name_P, self.N*np.pi/b, name_PP)

        super(RectWGPort, self).__init__(CSX, port_nr=port_nr, start=start, stop=stop, exc_dir=exc_dir, E_WG_func=E_func, H_WG_func=H_func, kc=kc, excite=excite, local_origin=local_origin, **kw)

class CircWGPort(WaveguidePort):
    """Circular waveguide port with analytic TE mode profile.

    Generates Bessel-function mode functions (Pozar 3rd ed.) and forwards
    them to :class:`WaveguidePort`.  Only TE modes with the listed (n, m)
    indices are supported.  The fparser variables ``rho`` and ``a`` describe
    cylindrical coordinates in the xy-plane, so **the propagation axis must
    be z** (``exc_dir=2``) for the mode functions to be evaluated correctly
    in Cartesian meshes.

    Parameters
    ----------
    exc_dir : int or str
        Propagation direction: 0/'x', 1/'y', or 2/'z'.  The mode functions
        are expressed in the two transverse Cartesian coordinates, so all
        three directions are supported.
    radius : float
        Waveguide radius in metres.
    mode_name : str
        Four-character TE mode string, e.g. ``'TE11'`` or ``'TE21'``.
    pol_ang : float, optional
        Polarisation angle in radians (0 = horizontal, pi/2 = vertical).
    local_origin : array-like, 'corner', or 'center', optional
        Forwarded to :class:`WaveguidePort`.  Defaults to ``'center'`` so
        that the mode functions are evaluated relative to the port midpoint.

    Notes
    -----
    Supported modes and their Bessel zeros p'_nm (zeros of J_n'):

      TE01 3.832  TE11 1.841  TE21 3.054
      TE02 7.016  TE12 5.331  TE22 6.706
      TE03 10.174 TE13 8.536  TE23 9.970

    For each propagation direction the two transverse Cartesian coordinates
    (u, v) are determined by ``ny_P = (exc_ny+1)%3`` and
    ``ny_PP = (exc_ny+2)%3``.  The transverse radius and angle are then
    ``rho_t = sqrt(u²+v²)`` and ``a_t = atan2(v, u)``, which replace the
    xy-plane built-ins ``rho`` and ``a`` used in earlier versions.

    See Also
    --------
    Port, WaveguidePort, RectWGPort
    """

    _pnm = {
        (0, 1): 3.832,  (1, 1): 1.841,  (2, 1): 3.054,
        (0, 2): 7.016,  (1, 2): 5.331,  (2, 2): 6.706,
        (0, 3): 10.174, (1, 3): 8.536,  (2, 3): 9.970,
    }

    def __init__(self, CSX, port_nr, start, stop, exc_dir, radius, mode_name,
                 pol_ang=0, excite=0, local_origin='center', **kw):
        Port.__init__(self, CSX, port_nr, start, stop, excite=0, **kw)

        if not mode_name.upper().startswith('TE'):
            raise Exception('CircWGPort: only TE modes are supported')
        n = int(mode_name[2])
        m = int(mode_name[3])
        pnm = self._pnm.get((n, m))
        if pnm is None:
            raise Exception('CircWGPort: unsupported TE_nm mode: {}'.format(mode_name))

        unit    = CSX.GetGrid().GetDeltaUnit()
        kc      = pnm / radius       # cut-off wavenumber, 1/m
        kc_draw = kc * unit          # 1/drawing-unit (for fparser)

        exc_ny = CheckNyDir(exc_dir)
        ny_P   = (exc_ny + 1) % 3
        ny_PP  = (exc_ny + 2) % 3

        # Transverse cylindrical coordinates as fparser expressions for any propagation axis
        coords = 'xyz'
        u     = coords[ny_P]
        v     = coords[ny_PP]
        rho_t = 'sqrt({0}*{0}+{1}*{1})'.format(u, v)
        a_t   = 'atan2({},{})'.format(v, u)
        ang   = '({})-{:.15g}'.format(a_t, pol_ang)

        # Cylindrical E and H components (Pozar 3rd ed., TE_nm pattern n=1)
        c_r = -1.0 / kc_draw**2
        c_a =  1.0 / kc_draw
        Er = '{C:.15g}/({r})*cos({a})*j1({k:.15g}*({r}))'.format(C=c_r,  r=rho_t, a=ang, k=kc_draw)
        Ea = '{C:.15g}*sin({a})*0.5*(j0({k:.15g}*({r}))-jn(2,{k:.15g}*({r})))'.format(C=c_a,  r=rho_t, a=ang, k=kc_draw)
        Hr = '{C:.15g}*sin({a})*0.5*(j0({k:.15g}*({r}))-jn(2,{k:.15g}*({r})))'.format(C=-c_a, r=rho_t, a=ang, k=kc_draw)
        Ha = '{C:.15g}/({r})*cos({a})*j1({k:.15g}*({r}))'.format(C=c_r,  r=rho_t, a=ang, k=kc_draw)

        # Cartesian conversion: cos(a_t)=u/rho_t, sin(a_t)=v/rho_t; map to transverse axes
        r_draw = radius / unit
        mask  = '(({r})<{d:.15g})'.format(r=rho_t, d=r_draw)
        cos_a = '{}/({})'.format(u, rho_t)
        sin_a = '{}/({})'.format(v, rho_t)
        E_func = ['0', '0', '0']
        H_func = ['0', '0', '0']
        E_func[ny_P]  = '(({Er})*({ca})-({Ea})*({sa}))*{m}'.format(Er=Er, Ea=Ea, ca=cos_a, sa=sin_a, m=mask)
        E_func[ny_PP] = '(({Er})*({sa})+({Ea})*({ca}))*{m}'.format(Er=Er, Ea=Ea, ca=cos_a, sa=sin_a, m=mask)
        H_func[ny_P]  = '(({Hr})*({ca})-({Ha})*({sa}))*{m}'.format(Hr=Hr, Ha=Ha, ca=cos_a, sa=sin_a, m=mask)
        H_func[ny_PP] = '(({Hr})*({sa})+({Ha})*({ca}))*{m}'.format(Hr=Hr, Ha=Ha, ca=cos_a, sa=sin_a, m=mask)

        super(CircWGPort, self).__init__(
            CSX, port_nr=port_nr, start=start, stop=stop,
            exc_dir=exc_dir, E_WG_func=E_func, H_WG_func=H_func,
            kc=kc, excite=excite, local_origin=local_origin, **kw)


class CoaxialPort(Port):
    """Coaxial transmission line port.

    Creates the coaxial geometry (inner conductor, outer shell, optional
    dielectric fill), voltage/current probes for the differential TL method,
    and an optional weighted radial-field excitation.

    Parameters
    ----------
    pec_prop : CSProperties
        Metal property for the inner conductor and outer shell.
    mat_prop : CSProperties or None
        Dielectric property for the filling between the conductors.  Pass
        ``None`` for an air-filled coaxial line.
    prop_dir : int or str
        Direction of wave propagation (0/1/2 or 'x'/'y'/'z').
    r_i : float
        Inner conductor radius in drawing units.
    r_o : float
        Outer conductor inner radius in drawing units.
    r_os : float
        Outer conductor outer radius in drawing units.
    excite : float, optional
        Excitation amplitude of the transverse E-field profile.  Set to 0
        (default) for a passive port.
    FeedShift : float, optional
        Distance in drawing units to shift the excitation from ``start``.
    Feed_R : float, optional
        Lumped port resistance.  Default is ``numpy.inf`` (open).  Only 0
        (short/metal termination) and ``inf`` (open) are currently supported.
    MeasPlaneShift : float, optional
        Distance in drawing units from ``start`` to the measurement plane.
        Default is the midpoint of the port.

    See Also
    --------
    Port, MSLPort
    """

    def __init__(self, CSX, port_nr, pec_prop, mat_prop, start, stop,
                 prop_dir, r_i, r_o, r_os, excite=0, **kw):
        super(CoaxialPort, self).__init__(CSX, port_nr=port_nr, start=start, stop=stop, excite=excite, **kw)

        self.prop_ny = CheckNyDir(prop_dir)
        self.ny_P    = (self.prop_ny + 1) % 3
        self.ny_PP   = (self.prop_ny + 2) % 3
        self.direction = np.sign(stop[self.prop_ny] - start[self.prop_ny])
        if self.direction == 0:
            raise Exception('CoaxialPort: start and stop must differ in propagation direction')

        self.r_i = r_i
        self.r_o = r_o

        feed_shift = kw.get('FeedShift', 0)
        feed_R     = kw.get('Feed_R', np.inf)

        # Default measurement plane at midpoint
        measplane_pos = 0.5 * (start[self.prop_ny] + stop[self.prop_ny])
        if 'MeasPlaneShift' in kw:
            measplane_pos = start[self.prop_ny] + self.direction * kw['MeasPlaneShift']

        # Coaxial geometry
        pec_prop.AddCylinder(start, stop, r_i, priority=self.priority)
        pec_prop.AddCylindricalShell(start, stop, 0.5*(r_o+r_os), r_os-r_o, priority=self.priority)
        if mat_prop is not None:
            mat_prop.AddCylindricalShell(start, stop, 0.5*(r_i+r_o), r_o-r_i, priority=self.priority-1)
        self.port_props.append(pec_prop)

        # Snap measurement plane to mesh and extract 3 adjacent lines
        mesh       = CSX.GetGrid()
        prop_lines = mesh.GetLines(self.prop_ny)
        idx = np.argmin(np.abs(prop_lines - measplane_pos))
        idx = max(1, min(idx, len(prop_lines) - 2))
        meshlines = prop_lines[idx-1:idx+2]
        if self.direction < 0:
            meshlines = meshlines[::-1]

        self.measplane_shift = abs(meshlines[1] - start[self.prop_ny])
        self.U_delta = np.diff(meshlines)
        i_pos = meshlines[:2] + np.diff(meshlines) / 2.0
        self.I_delta = np.diff(i_pos)

        # Voltage probes: radial line from r_i to r_o at three prop positions
        suffix = ['A', 'B', 'C']
        self.U_filenames = []
        for n in range(3):
            v_start = np.array(start, dtype=float)
            v_start[self.prop_ny] = meshlines[n]
            v_start[self.ny_P]    = start[self.ny_P] + r_i
            v_start[self.ny_PP]   = start[self.ny_PP]
            v_stop = v_start.copy()
            v_stop[self.ny_P]     = start[self.ny_P] + r_o

            u_name = self.lbl_temp.format('ut') + suffix[n]
            self.U_filenames.append(u_name)
            u_probe = CSX.AddProbe(u_name, p_type=0, weight=1)
            u_probe.AddBox(v_start, v_stop)
            self.port_props.append(u_probe)

        # Current probes: square loop encircling the inner conductor
        margin = 0.1 * (r_o - r_i)
        self.I_filenames = []
        for n in range(2):
            i_start = np.zeros(3)
            i_start[self.prop_ny] = 0.5 * (meshlines[n] + meshlines[n+1])
            i_start[self.ny_P]    = start[self.ny_P]  - r_i - margin
            i_start[self.ny_PP]   = start[self.ny_PP] - r_i - margin
            i_stop = i_start.copy()
            i_stop[self.ny_P]     = start[self.ny_P]  + r_i + margin
            i_stop[self.ny_PP]    = start[self.ny_PP] + r_i + margin

            i_name = self.lbl_temp.format('it') + suffix[n]
            self.I_filenames.append(i_name)
            i_probe = CSX.AddProbe(i_name, p_type=1, weight=self.direction, norm_dir=self.prop_ny)
            i_probe.AddBox(i_start, i_stop)
            self.port_props.append(i_probe)

        # Excitation: thin cylindrical shell with radial E-field weighting
        if excite != 0:
            prop_feed_idx = np.argmin(np.abs(prop_lines - (start[self.prop_ny] + feed_shift*self.direction)))
            min_cell = np.min(np.diff(prop_lines))
            ex_start = np.array(start, dtype=float)
            ex_start[self.prop_ny] = prop_lines[prop_feed_idx] - 0.01*min_cell
            ex_stop = ex_start.copy()
            ex_stop[self.prop_ny]  = prop_lines[prop_feed_idx] + 0.01*min_cell

            xyz = 'xyz'
            nP_name  = xyz[self.ny_P]
            nPP_name = xyz[self.ny_PP]
            cx = start[self.ny_P]
            cy = start[self.ny_PP]
            dX  = '({}-{:.15g})'.format(nP_name,  cx)
            dY  = '({}-{:.15g})'.format(nPP_name, cy)
            r2  = '({}*{}+{}*{})'.format(dX, dX, dY, dY)
            mask = '*(sqrt({r2})<{r_o:.15g})*(sqrt({r2})>{r_i:.15g})'.format(r2=r2, r_o=r_o, r_i=r_i)
            func_E = ['0', '0', '0']
            func_E[self.ny_P]  = '{}/{}{}' .format(dX, r2, mask)
            func_E[self.ny_PP] = '{}/{}{}' .format(dY, r2, mask)

            exc_val = excite*np.ones(3)
            exc_val[self.prop_ny] = 0
            exc = CSX.AddExcitation(self.lbl_temp.format('excite'), exc_type=0,
                                    exc_val=exc_val, delay=self.delay)
            exc.SetWeightFunction(func_E)
            exc.AddCylindricalShell(ex_start, ex_stop, 0.5*(r_i+r_o), r_o-r_i, priority=0)
            self.port_props.append(exc)

        # Termination at start of line
        term_start = np.array(start, dtype=float)
        term_stop  = np.array(stop,  dtype=float)
        term_stop[self.prop_ny] = term_start[self.prop_ny]
        if feed_R == 0:
            pec_prop.AddBox(term_start, term_stop, priority=self.priority)
            pec_prop.AddCylindricalShell(term_start, term_stop, 0.5*(r_i+r_o), r_o-r_i, priority=self.priority)
        elif np.isinf(feed_R):
            pass  # open termination
        elif feed_R > 0:
            raise NotImplementedError('CoaxialPort: finite Feed_R > 0 is not yet supported; use Feed_R=0 or Feed_R=inf')
        else:
            raise Exception('CoaxialPort: Feed_R <= 0 is not allowed')

    def ReadUIData(self, sim_path, freq, signal_type='pulse'):
        self.u_data = UI_data(self.U_filenames, sim_path, freq, signal_type)
        self.uf_tot = self.u_data.ui_f_val[1]
        self.ut_tot = self.u_data.ui_val[1]

        self.i_data = UI_data(self.I_filenames, sim_path, freq, signal_type)
        self.if_tot = 0.5 * (self.i_data.ui_f_val[0] + self.i_data.ui_f_val[1])
        self.it_tot = 0.5 * (self.i_data.ui_val[0]   + self.i_data.ui_val[1])

        unit = self.CSX.GetGrid().GetDeltaUnit()
        Et   = self.u_data.ui_f_val[1]
        dEt  = (self.u_data.ui_f_val[2] - self.u_data.ui_f_val[0]) / (np.sum(np.abs(self.U_delta)) * unit)
        Ht   = self.if_tot
        dHt  = (self.i_data.ui_f_val[1] - self.i_data.ui_f_val[0]) / (np.abs(self.I_delta[0]) * unit)

        beta = np.sqrt(-dEt * dHt / (Ht * Et))
        beta[np.real(beta) < 0] *= -1
        self.beta  = beta
        self.Z_ref = np.sqrt(Et * dEt / (Ht * dHt))


class StripLinePort(Port):
    """Stripline transmission line port.

    Creates the stripline metal, symmetric voltage probes (above and below the
    strip toward the upper and lower ground planes), current probes spanning the
    full strip width, and an optional symmetric excitation.

    Parameters
    ----------
    metal_prop : CSProperties
        Metal property for the stripline conductor.
    prop_dir : int or str
        Direction of wave propagation (0/1/2 or 'x'/'y'/'z').
    exc_dir : int or str or (3,) array
        E-field direction (must have exactly one non-zero component).
    height : float
        Distance from the strip to each ground plane, in drawing units.
    excite : bool or float, optional
        Set to True (or non-zero) to enable the port excitation.
    FeedShift : float, optional
        Excitation shift from ``start`` in drawing units.
    Feed_R : float, optional
        Lumped resistance in Ohms.  Default ``inf`` (open).
    MeasPlaneShift : float, optional
        Measurement plane distance from ``start`` in drawing units.

    See Also
    --------
    Port, MSLPort, CPWPort
    """

    def __init__(self, CSX, port_nr, metal_prop, start, stop, prop_dir, exc_dir,
                 height, excite=0, **kw):
        super(StripLinePort, self).__init__(CSX, port_nr=port_nr, start=start, stop=stop, excite=excite, **kw)

        self.prop_ny = CheckNyDir(prop_dir)

        # Determine height (E-field) and width directions from exc_dir
        exc_vec = np.zeros(3)
        exc_vec[CheckNyDir(exc_dir)] = 1.0
        self.height_ny = int(np.argmax(np.abs(exc_vec)))

        prop_vec = np.zeros(3)
        prop_vec[self.prop_ny] = 1.0
        self.width_ny = int(np.argmax(np.abs(np.cross(prop_vec, exc_vec))))

        if start[self.height_ny] != stop[self.height_ny]:
            raise Exception('StripLinePort: start/stop in height direction must be equal')

        self.direction = np.sign(stop[self.prop_ny] - start[self.prop_ny])
        if self.direction == 0:
            raise Exception('StripLinePort: start/stop in prop direction must differ')

        feed_shift    = kw.get('FeedShift', 0)
        feed_R        = kw.get('Feed_R', np.inf)

        nstart = np.minimum(start, stop)
        nstop  = np.maximum(start, stop)

        # Stripline metal layer
        metal_prop.AddBox(np.array(start), np.array(stop), priority=self.priority)
        self.port_props.append(metal_prop)

        # Measurement plane position
        measplane_pos = 0.5 * (nstart[self.prop_ny] + nstop[self.prop_ny])
        if 'MeasPlaneShift' in kw:
            measplane_pos = start[self.prop_ny] + self.direction * kw['MeasPlaneShift']

        mesh       = CSX.GetGrid()
        prop_lines = mesh.GetLines(self.prop_ny)
        idx = np.argmin(np.abs(prop_lines - measplane_pos))
        idx = max(1, min(idx, len(prop_lines) - 2))
        meshlines = prop_lines[idx-1:idx+2]
        if self.direction < 0:
            meshlines = meshlines[::-1]

        self.measplane_shift = abs(meshlines[1] - start[self.prop_ny])
        self.U_delta = np.diff(meshlines)
        i_pos = meshlines[:2] + np.diff(meshlines) / 2.0
        self.I_delta = np.diff(i_pos)

        # Width center (nearest mesh line to strip center)
        width_lines = mesh.GetLines(self.width_ny)
        w_center = 0.5 * (nstart[self.width_ny] + nstop[self.width_ny])
        w_center_idx = np.argmin(np.abs(width_lines - w_center))
        SL_w2 = width_lines[w_center_idx]

        # Height direction vector
        height_vec = np.zeros(3)
        height_vec[self.height_ny] = height

        # Voltage probes: pairs (upper + lower) at three prop positions
        suffix_pairs = [('A1', 'A2'), ('B1', 'B2'), ('C1', 'C2')]
        self.U_filenames = []
        for n, (s1, s2) in enumerate(suffix_pairs):
            v_pt = np.zeros(3)
            v_pt[self.prop_ny]   = meshlines[n]
            v_pt[self.width_ny]  = SL_w2
            v_pt[self.height_ny] = start[self.height_ny]

            for s, sign in [(s1, +1), (s2, -1)]:
                u_name = self.lbl_temp.format('ut') + s
                self.U_filenames.append(u_name)
                u_probe = CSX.AddProbe(u_name, p_type=0, weight=0.5)
                u_probe.AddBox(v_pt, v_pt + sign * height_vec, priority=self.priority)
                self.port_props.append(u_probe)

        # Current probe boundaries: span full strip width, one cell above/below
        height_lines = mesh.GetLines(self.height_ny)
        h_idx = np.argmin(np.abs(height_lines - start[self.height_ny]))
        h_idx = max(2, min(h_idx, len(height_lines) - 3))

        w_idx_start = np.argmin(np.abs(width_lines - nstart[self.width_ny]))
        w_idx_start = max(1, min(w_idx_start, len(width_lines) - 1))
        w_idx_stop  = np.argmin(np.abs(width_lines - nstop[self.width_ny]))
        w_idx_stop  = max(0, min(w_idx_stop,  len(width_lines) - 2))

        i_base_start = np.zeros(3)
        i_base_start[self.width_ny]  = 0.5 * (width_lines[w_idx_start-1] + width_lines[w_idx_start])
        i_base_start[self.height_ny] = 0.5 * (height_lines[h_idx-2] + height_lines[h_idx-1])
        i_base_stop = np.zeros(3)
        i_base_stop[self.width_ny]   = 0.5 * (width_lines[w_idx_stop]   + width_lines[w_idx_stop+1])
        i_base_stop[self.height_ny]  = 0.5 * (height_lines[h_idx+1] + height_lines[h_idx+2])

        self.I_filenames = []
        for n, s in enumerate(['A', 'B']):
            i_start = i_base_start.copy()
            i_stop  = i_base_stop.copy()
            i_start[self.prop_ny] = 0.5 * (meshlines[n]   + meshlines[n+1])
            i_stop[self.prop_ny]  = i_start[self.prop_ny]

            i_name = self.lbl_temp.format('it') + s
            self.I_filenames.append(i_name)
            i_probe = CSX.AddProbe(i_name, p_type=1, weight=self.direction, norm_dir=self.prop_ny)
            i_probe.AddBox(i_start, i_stop)
            self.port_props.append(i_probe)

        # Excitation: two symmetric boxes (above/below strip) at feed position
        if excite != 0:
            feed_idx = np.argmin(np.abs(prop_lines - (start[self.prop_ny] + feed_shift * self.direction)))
            ex_start = np.zeros(3)
            ex_start[self.prop_ny]   = prop_lines[feed_idx]
            ex_start[self.width_ny]  = nstart[self.width_ny]
            ex_start[self.height_ny] = nstart[self.height_ny]
            ex_stop = ex_start.copy()
            ex_stop[self.prop_ny]    = prop_lines[feed_idx]
            ex_stop[self.width_ny]   = nstop[self.width_ny]
            ex_stop[self.height_ny]  = nstop[self.height_ny]

            exc_val = np.zeros(3)
            exc_val[self.height_ny] = excite
            for lbl_s, sign in [('excite_1', +1), ('excite_2', -1)]:
                exc = CSX.AddExcitation(self.lbl_temp.format(lbl_s), exc_type=0,
                                        exc_val=sign * exc_val, delay=self.delay)
                exc.AddBox(ex_start, ex_stop + sign * height_vec, priority=self.priority)
                self.port_props.append(exc)

        # Termination resistance at start
        r_start = np.array(start, dtype=float)
        r_stop  = np.array(stop,  dtype=float)
        r_stop[self.prop_ny] = r_start[self.prop_ny]
        if feed_R > 0 and not np.isinf(feed_R):
            lumped = CSX.AddLumpedElement(self.lbl_temp.format('resist'),
                                          ny=self.height_ny, R=2*feed_R)
            lumped.AddBox(r_start, r_stop + height_vec,  priority=self.priority)
            lumped.AddBox(r_start, r_stop - height_vec,  priority=self.priority)
            self.port_props.append(lumped)
        elif np.isinf(feed_R):
            pass
        elif feed_R == 0:
            metal_prop.AddBox(r_start, r_stop + height_vec, priority=self.priority)
            metal_prop.AddBox(r_start, r_stop - height_vec, priority=self.priority)
        else:
            raise Exception('StripLinePort: Feed_R must be >= 0')

    def ReadUIData(self, sim_path, freq, signal_type='pulse'):
        all_u = UI_data(self.U_filenames, sim_path, freq, signal_type)

        # Sum paired (upper+lower) probes at each of the three positions
        uf_A = all_u.ui_f_val[0] + all_u.ui_f_val[1]
        uf_B = all_u.ui_f_val[2] + all_u.ui_f_val[3]
        uf_C = all_u.ui_f_val[4] + all_u.ui_f_val[5]
        ut_B = all_u.ui_val[2]   + all_u.ui_val[3]

        self.uf_tot = uf_B
        self.ut_tot = ut_B

        self.i_data = UI_data(self.I_filenames, sim_path, freq, signal_type)
        self.if_tot = 0.5 * (self.i_data.ui_f_val[0] + self.i_data.ui_f_val[1])
        self.it_tot = 0.5 * (self.i_data.ui_val[0]   + self.i_data.ui_val[1])

        unit = self.CSX.GetGrid().GetDeltaUnit()
        Et   = uf_B
        dEt  = (uf_C - uf_A) / (np.sum(np.abs(self.U_delta)) * unit)
        Ht   = self.if_tot
        dHt  = (self.i_data.ui_f_val[1] - self.i_data.ui_f_val[0]) / (np.abs(self.I_delta[0]) * unit)

        beta = np.sqrt(-dEt * dHt / (Ht * Et))
        beta[np.real(beta) < 0] *= -1
        self.beta  = beta
        self.Z_ref = np.sqrt(Et * dEt / (Ht * dHt))


class CPWPort(Port):
    """Coplanar waveguide (CPW) port.

    Creates the CPW metal, gap-spanning voltage probes (left and right gaps),
    current probes, and an optional symmetric excitation across the two gaps.

    Parameters
    ----------
    metal_prop : CSProperties
        Metal property for the CPW conductor.
    prop_dir : int or str
        Direction of wave propagation (0/1/2 or 'x'/'y'/'z').
    exc_dir : int or str or (3,) array
        E-field direction across the gaps (one non-zero component).
    gap_width : float
        Width of each CPW gap in drawing units.
    excite : bool or float, optional
        Enable port excitation.
    FeedShift : float, optional
        Excitation shift from ``start``.
    Feed_R : float, optional
        Lumped resistance in Ohms (applied to each gap as 2*R).
    MeasPlaneShift : float, optional
        Measurement plane distance from ``start``.

    See Also
    --------
    Port, MSLPort, StripLinePort
    """

    def __init__(self, CSX, port_nr, metal_prop, start, stop, prop_dir, exc_dir,
                 gap_width, excite=0, **kw):
        super(CPWPort, self).__init__(CSX, port_nr=port_nr, start=start, stop=stop, excite=excite, **kw)

        self.prop_ny = CheckNyDir(prop_dir)

        # Height direction = E-field direction; width direction = cross product
        exc_vec = np.zeros(3)
        exc_vec[CheckNyDir(exc_dir)] = 1.0
        self.height_ny = int(np.argmax(np.abs(exc_vec)))

        prop_vec = np.zeros(3)
        prop_vec[self.prop_ny] = 1.0
        self.width_ny = int(np.argmax(np.abs(np.cross(prop_vec, exc_vec))))

        if start[self.height_ny] != stop[self.height_ny]:
            raise Exception('CPWPort: start/stop in height direction must be equal')

        self.direction = np.sign(stop[self.prop_ny] - start[self.prop_ny])
        if self.direction == 0:
            raise Exception('CPWPort: start/stop in prop direction must differ')

        feed_shift = kw.get('FeedShift', 0)
        feed_R     = kw.get('Feed_R', np.inf)

        nstart = np.minimum(start, stop)
        nstop  = np.maximum(start, stop)

        # CPW metal layer
        metal_prop.AddBox(np.array(start), np.array(stop), priority=self.priority)
        self.port_props.append(metal_prop)

        # Measurement plane
        measplane_pos = 0.5 * (nstart[self.prop_ny] + nstop[self.prop_ny])
        if 'MeasPlaneShift' in kw:
            measplane_pos = start[self.prop_ny] + self.direction * kw['MeasPlaneShift']

        mesh       = CSX.GetGrid()
        prop_lines = mesh.GetLines(self.prop_ny)
        idx = np.argmin(np.abs(prop_lines - measplane_pos))
        idx = max(1, min(idx, len(prop_lines) - 2))
        meshlines = prop_lines[idx-1:idx+2]
        if self.direction < 0:
            meshlines = meshlines[::-1]

        self.measplane_shift = abs(meshlines[1] - start[self.prop_ny])
        self.U_delta = np.diff(meshlines)
        i_pos = meshlines[:2] + np.diff(meshlines) / 2.0
        self.I_delta = np.diff(i_pos)

        # Half-width and gap offsets in width direction
        w_center = 0.5 * (nstart[self.width_ny] + nstop[self.width_ny])
        half_w   = 0.5 * (nstop[self.width_ny] - nstart[self.width_ny])

        w_add_start = np.zeros(3)
        w_add_stop  = np.zeros(3)
        w_add_start[self.width_ny] = half_w
        w_add_stop[self.width_ny]  = half_w + gap_width

        # Voltage probes: pairs (left/right gap) at three prop positions
        suffix_pairs = [('A1', 'A2'), ('B1', 'B2'), ('C1', 'C2')]
        self.U_filenames = []
        for n, (s1, s2) in enumerate(suffix_pairs):
            v_pt = np.zeros(3)
            v_pt[self.prop_ny]   = meshlines[n]
            v_pt[self.width_ny]  = w_center
            v_pt[self.height_ny] = start[self.height_ny]

            for s, sign in [(s1, -1), (s2, +1)]:
                u_name = self.lbl_temp.format('ut') + s
                self.U_filenames.append(u_name)
                u_probe = CSX.AddProbe(u_name, p_type=0, weight=0.5)
                u_probe.AddBox(v_pt + sign*w_add_start, v_pt + sign*w_add_stop,
                               priority=self.priority)
                self.port_props.append(u_probe)

        # Current probes: span CPW + gaps width, one cell in height
        height_lines = mesh.GetLines(self.height_ny)
        width_lines  = mesh.GetLines(self.width_ny)
        h_idx = np.argmin(np.abs(height_lines - start[self.height_ny]))
        h_idx = max(2, min(h_idx, len(height_lines) - 3))

        w_idx_start = np.argmin(np.abs(width_lines - nstart[self.width_ny]))
        w_idx_start = max(1, min(w_idx_start, len(width_lines) - 1))
        w_idx_stop  = np.argmin(np.abs(width_lines - nstop[self.width_ny]))
        w_idx_stop  = max(0, min(w_idx_stop,  len(width_lines) - 2))

        i_base_start = np.zeros(3)
        i_base_start[self.width_ny]  = 0.5 * (width_lines[w_idx_start-1] + width_lines[w_idx_start])
        i_base_start[self.height_ny] = 0.5 * (height_lines[h_idx-2] + height_lines[h_idx-1])
        i_base_stop = np.zeros(3)
        i_base_stop[self.width_ny]   = 0.5 * (width_lines[w_idx_stop]   + width_lines[w_idx_stop+1])
        i_base_stop[self.height_ny]  = 0.5 * (height_lines[h_idx+1] + height_lines[h_idx+2])

        self.I_filenames = []
        for n, s in enumerate(['A', 'B']):
            i_start = i_base_start.copy()
            i_stop  = i_base_stop.copy()
            i_start[self.prop_ny] = 0.5 * (meshlines[n]   + meshlines[n+1])
            i_stop[self.prop_ny]  = i_start[self.prop_ny]

            i_name = self.lbl_temp.format('it') + s
            self.I_filenames.append(i_name)
            i_probe = CSX.AddProbe(i_name, p_type=1, weight=self.direction, norm_dir=self.prop_ny)
            i_probe.AddBox(i_start, i_stop)
            self.port_props.append(i_probe)

        # Excitation: two gap boxes, E-field in the width (gap) direction
        if excite != 0:
            feed_idx = np.argmin(np.abs(prop_lines - (start[self.prop_ny] + feed_shift * self.direction)))
            ex_pt = np.zeros(3)
            ex_pt[self.prop_ny]   = prop_lines[feed_idx]
            ex_pt[self.width_ny]  = w_center
            ex_pt[self.height_ny] = start[self.height_ny]

            for lbl_s, sign in [('excite_1', -1), ('excite_2', +1)]:
                exc_val = np.zeros(3)
                exc_val[self.width_ny] = sign*excite
                exc = CSX.AddExcitation(self.lbl_temp.format(lbl_s), exc_type=0,
                                        exc_val=exc_val, delay=self.delay)
                exc.AddBox(ex_pt + sign*w_add_start, ex_pt + sign*w_add_stop,
                           priority=self.priority)
                self.port_props.append(exc)

        # Termination resistance at start — centred on w_center, one box per gap
        r_pt = np.zeros(3)
        r_pt[self.prop_ny]   = start[self.prop_ny]
        r_pt[self.width_ny]  = w_center
        r_pt[self.height_ny] = start[self.height_ny]
        if feed_R > 0 and not np.isinf(feed_R):
            lumped = CSX.AddLumpedElement(self.lbl_temp.format('resist'),
                                          ny=self.width_ny, R=2*feed_R)
            for sign in [-1, +1]:
                lumped.AddBox(r_pt + sign*w_add_start, r_pt + sign*w_add_stop, priority=self.priority)
            self.port_props.append(lumped)
        elif np.isinf(feed_R):
            pass
        elif feed_R == 0:
            for sign in [-1, +1]:
                metal_prop.AddBox(r_pt + sign*w_add_start, r_pt + sign*w_add_stop, priority=self.priority)
        else:
            raise Exception('CPWPort: Feed_R must be >= 0')

    def ReadUIData(self, sim_path, freq, signal_type='pulse'):
        all_u = UI_data(self.U_filenames, sim_path, freq, signal_type)

        # Sum paired (left+right gap) probes at each of the three positions
        uf_A = all_u.ui_f_val[0] + all_u.ui_f_val[1]
        uf_B = all_u.ui_f_val[2] + all_u.ui_f_val[3]
        uf_C = all_u.ui_f_val[4] + all_u.ui_f_val[5]
        ut_B = all_u.ui_val[2]   + all_u.ui_val[3]

        self.uf_tot = uf_B
        self.ut_tot = ut_B

        self.i_data = UI_data(self.I_filenames, sim_path, freq, signal_type)
        self.if_tot = 0.5 * (self.i_data.ui_f_val[0] + self.i_data.ui_f_val[1])
        self.it_tot = 0.5 * (self.i_data.ui_val[0]   + self.i_data.ui_val[1])

        unit = self.CSX.GetGrid().GetDeltaUnit()
        Et   = uf_B
        dEt  = (uf_C - uf_A) / (np.sum(np.abs(self.U_delta)) * unit)
        Ht   = self.if_tot
        dHt  = (self.i_data.ui_f_val[1] - self.i_data.ui_f_val[0]) / (np.abs(self.I_delta[0]) * unit)

        beta = np.sqrt(-dEt * dHt / (Ht * Et))
        beta[np.real(beta) < 0] *= -1
        self.beta  = beta
        self.Z_ref = np.sqrt(Et * dEt / (Ht * dHt))


class CurvePort(Port):
    """One-dimensional curve-based lumped port.

    Creates a single-cell lumped port aligned with the nearest mesh edge.
    When ``start`` and ``stop`` span more than one mesh cell, PEC curves are
    added to connect the user coordinates to the one-cell port location at the
    midpoint.

    Parameters
    ----------
    R : float
        Port reference impedance in Ohms.
    excite : bool or float, optional
        Enable port excitation.

    See Also
    --------
    Port, LumpedPort
    """

    def __init__(self, CSX, port_nr, R, start, stop, excite=0, **kw):
        super(CurvePort, self).__init__(CSX, port_nr=port_nr, start=start, stop=stop, excite=excite, **kw)
        self.R     = R
        self.Z_ref = R

        mesh  = CSX.GetGrid()
        unit  = mesh.GetDeltaUnit()
        lines = [mesh.GetLines(n) for n in range(3)]

        # Find dominant direction
        delta = np.abs(self.stop - self.start)
        self.port_dir = int(np.argmax(delta))
        dir1 = (self.port_dir + 1) % 3
        dir2 = (self.port_dir + 2) % 3

        # Normalise so start < stop in port direction
        if self.start[self.port_dir] <= self.stop[self.port_dir]:
            nstart, nstop = np.array(self.start), np.array(self.stop)
        else:
            nstart, nstop = np.array(self.stop),  np.array(self.start)

        # Snap all coordinates to nearest mesh lines
        snap = lambda v, ln: (np.argmin(np.abs(ln - v)), ln[np.argmin(np.abs(ln - v))])
        s_idx = [snap(nstart[d], lines[d])[0] for d in range(3)]
        e_idx = [snap(nstop[d],  lines[d])[0] for d in range(3)]

        if abs(s_idx[self.port_dir] - e_idx[self.port_dir]) != 1:
            # Port spans more than one cell; find one-cell edge at midpoint
            mid = 0.5 * (nstart[self.port_dir] + nstop[self.port_dir])
            p_idx  = np.argmin(np.abs(lines[self.port_dir] - mid))
            p1_idx = np.argmin(np.abs(lines[dir1] - 0.5*(nstart[dir1]+nstop[dir1])))
            p2_idx = np.argmin(np.abs(lines[dir2] - 0.5*(nstart[dir2]+nstop[dir2])))

            port_start_idx = [0, 0, 0]
            port_start_idx[self.port_dir] = p_idx
            port_start_idx[dir1]          = p1_idx
            port_start_idx[dir2]          = p2_idx
            port_stop_idx = port_start_idx[:]
            port_stop_idx[self.port_dir]  = p_idx + 1

            # PEC curves connecting user coordinates to port edge
            edge_start = np.array([lines[d][port_start_idx[d]] for d in range(3)])
            edge_stop  = np.array([lines[d][port_stop_idx[d]]  for d in range(3)])

            metal_name = self.lbl_temp.format('PEC')
            metal_prop = CSX.AddMetal(metal_name)
            self.port_props.append(metal_prop)
            for user_pt in [nstart, nstop]:
                near = edge_start if np.linalg.norm(user_pt - edge_start) <= np.linalg.norm(user_pt - edge_stop) else edge_stop
                metal_prop.AddCurve(
                    [[user_pt[0], near[0]], [user_pt[1], near[1]], [user_pt[2], near[2]]],
                    priority=self.priority)
        else:
            port_start_idx = s_idx
            port_stop_idx  = e_idx
            edge_start = np.array([lines[d][port_start_idx[d]] for d in range(3)])
            edge_stop  = np.array([lines[d][port_stop_idx[d]]  for d in range(3)])

        # Cell-size margins for the current probe box in the transverse plane
        delta1_n = lines[dir1][port_start_idx[dir1]] - lines[dir1][port_start_idx[dir1]-1]
        delta1_p = lines[dir1][port_start_idx[dir1]+1] - lines[dir1][port_start_idx[dir1]]
        delta2_n = lines[dir2][port_start_idx[dir2]] - lines[dir2][port_start_idx[dir2]-1]
        delta2_p = lines[dir2][port_start_idx[dir2]+1] - lines[dir2][port_start_idx[dir2]]

        i_start = edge_start.copy()
        i_stop  = edge_stop.copy()
        i_start[dir1] -= delta1_n / 2
        i_stop[dir1]  += delta1_p / 2
        i_start[dir2] -= delta2_n / 2
        i_stop[dir2]  += delta2_p / 2
        # Current probe is in transverse plane; set prop coordinate to midpoint
        mid_prop = 0.5 * (edge_start[self.port_dir] + edge_stop[self.port_dir])
        i_start[self.port_dir] = mid_prop
        i_stop[self.port_dir]  = mid_prop

        # Voltage probe (with weight=-1 to get correct sign) and current probe
        self.U_filenames = [self.lbl_temp.format('ut')]
        u_probe = CSX.AddProbe(self.U_filenames[0], p_type=0, weight=-1)
        u_probe.AddBox(edge_start, edge_stop)
        self.port_props.append(u_probe)

        self.I_filenames = [self.lbl_temp.format('it')]
        i_probe = CSX.AddProbe(self.I_filenames[0], p_type=1, weight=1)
        i_probe.AddBox(i_start, i_stop)
        self.port_props.append(i_probe)

        # Lumped resistance or short
        if R > 0 and not np.isinf(R):
            lumped = CSX.AddLumpedElement(self.lbl_temp.format('resist'),
                                          ny=self.port_dir, R=R)
            lumped.AddBox(edge_start, edge_stop, priority=self.priority)
            self.port_props.append(lumped)
        elif R == 0:
            metal_name2 = self.lbl_temp.format('short')
            short_prop  = CSX.AddMetal(metal_name2)
            short_prop.AddBox(edge_start, edge_stop, priority=self.priority)
            self.port_props.append(short_prop)

        # Excitation
        if excite:
            exc_dir_vec = excite*(np.array(port_stop_idx) != np.array(port_start_idx)).astype(float)
            exc = CSX.AddExcitation(self.lbl_temp.format('excite'), exc_type=0,
                                    exc_val=exc_dir_vec, delay=self.delay)
            exc.AddBox(edge_start, edge_stop, priority=self.priority)
            self.port_props.append(exc)

    def CalcPort(self, sim_path, freq, ref_impedance=None, ref_plane_shift=None, signal_type='pulse'):
        if ref_impedance is None:
            self.Z_ref = self.R
        if ref_plane_shift is not None:
            Warning('CurvePort does not support a reference plane shift; ignoring')
        super(CurvePort, self).CalcPort(sim_path, freq, ref_impedance, None, signal_type)

