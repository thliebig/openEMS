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
from openEMS import utilities

from openEMS.physical_constants import *

class UI_data:
    def __init__(self, fns, path, freq, signal_type='pulse', **kw):
        self.path = path
        if type(fns)==str:
            fns = [fns]
        self.fns  = fns

        if np.isscalar(freq):
            freq = [freq]
        self.freq = freq

        self.ui_time = []
        self.ui_val  = []
        self.ui_f_val = []

        for fn in fns:
            tmp = np.loadtxt(os.path.join(path, fn),comments='%')
            self.ui_time.append(tmp[:,0])
            self.ui_val.append(tmp[:,1])
            self.ui_f_val.append(utilities.DFT_time2freq(tmp[:,0], tmp[:,1], freq, signal_type=signal_type))

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
    :param delay: float -- a positiv delay value to e.g. emulate a phase shift
    """
    def __init__(self, CSX, port_nr, start, stop, excite, **kw):
        self.CSX      = CSX
        self.number   = port_nr
        self.excite   = excite
        self.start    = np.array(start, np.float)
        self.stop     = np.array(stop, np.float)
        self.Z_ref    = None
        self.U_filenames = kw.get('U_filenames', [])
        self.I_filenames = kw.get('I_filenames', [])

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

    def ReadUIData(self, sim_path, freq, signal_type ='pulse'):
        self.u_data = UI_data(self.U_filenames, sim_path, freq, signal_type )
        self.uf_tot = 0
        self.ut_tot = 0
        for n in range(len(self.U_filenames)):
            self.uf_tot += self.u_data.ui_f_val[n]
            self.ut_tot += self.u_data.ui_val[n]

        self.i_data = UI_data(self.I_filenames, sim_path, freq, signal_type )
        self.if_tot = 0
        self.it_tot = 0
        for n in range(len(self.U_filenames)):
            self.if_tot += self.i_data.ui_f_val[n]
            self.it_tot += self.i_data.ui_val[n]


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

        if type(self.Z_ref) == float:
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

        lumped_R.AddBox(self.start, self.stop, priority=self.priority)

        if excite!=0:
            exc_vec = np.zeros(3)
            exc_vec[self.exc_ny] = -1*self.direction*excite
            exc = CSX.AddExcitation(self.lbl_temp.format('excite'), exc_type=0, exc_val=exc_vec, delay=self.delay)
            exc.AddBox(self.start, self.stop, priority=self.priority)

        self.U_filenames = [self.lbl_temp.format('ut'), ]
        u_start = 0.5*(self.start+self.stop)
        u_start[self.exc_ny] = self.start[self.exc_ny]
        u_stop  = 0.5*(self.start+self.stop)
        u_stop[self.exc_ny]  = self.stop[self.exc_ny]
        u_probe = CSX.AddProbe(self.U_filenames[0], p_type=0, weight=-1)
        u_probe.AddBox(u_start, u_stop)

        self.I_filenames = [self.lbl_temp.format('it'), ]
        i_start = np.array(self.start)
        i_start[self.exc_ny] = 0.5*(self.start[self.exc_ny]+self.stop[self.exc_ny])
        i_stop  = np.array(self.stop)
        i_stop[self.exc_ny]  = 0.5*(self.start[self.exc_ny]+self.stop[self.exc_ny])
        i_probe = CSX.AddProbe(self.I_filenames[0], p_type=1, weight=self.direction, norm_dir=self.exc_ny)
        i_probe.AddBox(i_start, i_stop)

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
        prope_idx = np.array([meas_pos_idx-1, meas_pos_idx, meas_pos_idx+1], np.int)
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

        if self.feed_R>=0 and not np.isinf(self.feed_R):
            R_start = np.array(self.start)
            R_stop  = np.array(self.stop)
            R_stop [self.prop_ny] = R_start[self.prop_ny]
            if self.feed_R==0:
                metal_prop.AddBox(R_start, R_stop)
            else:
                lumped_R = CSX.AddLumpedElement(self.lbl_temp.format('resist'), ny=self.exc_ny, caps=True, R=self.feed_R)
                lumped_R.AddBox(R_start, R_stop)

    def ReadUIData(self, sim_path, freq, signal_type ='pulse'):
        self.u_data = UI_data(self.U_filenames, sim_path, freq, signal_type )
        self.uf_tot = self.u_data.ui_f_val[1]

        self.i_data = UI_data(self.I_filenames, sim_path, freq, signal_type )
        self.if_tot = 0.5*(self.i_data.ui_f_val[0]+self.i_data.ui_f_val[1])

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

    See Also
    --------
    Port, RectWGPort

    """
    def __init__(self, CSX, port_nr, start, stop, exc_dir, E_WG_func, H_WG_func, kc, excite=0, **kw):
        super(WaveguidePort, self).__init__(CSX, port_nr=port_nr, start=start, stop=stop, excite=excite, **kw)
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

        if excite!=0:
            e_start = np.array(start)
            e_stop  = np.array(stop)
            e_stop[self.exc_ny] = e_start[self.exc_ny]
            e_vec = np.ones(3)
            e_vec[self.exc_ny]=0
            exc = CSX.AddExcitation(self.lbl_temp.format('excite'), exc_type=0, exc_val=e_vec, delay=self.delay)
            exc.SetWeightFunction([str(x) for x in self.E_func])
            exc.AddBox(e_start, e_stop, priority=self.priority)

        # voltage/current planes
        m_start = np.array(start)
        m_stop  = np.array(stop)
        m_start[self.exc_ny] = m_stop[self.exc_ny]
        self.measplane_shift = np.abs(stop[self.exc_ny] - start[self.exc_ny])

        self.U_filenames = [self.lbl_temp.format('ut'), ]

        u_probe = CSX.AddProbe(self.U_filenames[0], p_type=10, mode_function=self.E_func)
        u_probe.AddBox(m_start, m_stop)

        self.I_filenames = [self.lbl_temp.format('it'), ]
        i_probe = CSX.AddProbe(self.I_filenames[0], p_type=11, weight=self.direction, mode_function=self.H_func)
        i_probe.AddBox(m_start, m_stop)


    def CalcPort(self, sim_path, freq, ref_impedance=None, ref_plane_shift=None, signal_type='pulse'):
        k = 2.0*np.pi*freq/C0*self.ref_index
        self.beta = np.sqrt(k**2 - self.kc**2)
        self.ZL = k * Z0 / self.beta    #analytic waveguide impedance
        if ref_impedance is None:
            self.Z_ref = self.ZL
        super(WaveguidePort, self).CalcPort(sim_path, freq, ref_impedance, ref_plane_shift, signal_type)

class RectWGPort(WaveguidePort):
    """
    Rectangular waveguide port.

    :param a,b: float -- Width/Height of rectangular waveguide port

    See Also
    --------
    Port, WaveguidePort

    """
    def __init__(self, CSX, port_nr, start, stop, exc_dir, a, b, mode_name, excite=0, **kw):
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
        if self.start[self.ny_P]!=0:
            name_P = '({}-{})'.format(xyz[self.ny_P], self.start[self.ny_P])
        else:
            name_P = xyz[self.ny_P]
        if self.start[self.ny_PP]!=0:
            name_PP = '({}-{})'.format(xyz[self.ny_P], self.start[self.ny_P])
        else:
            name_PP = xyz[self.ny_P]

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

        super(RectWGPort, self).__init__(CSX, port_nr=port_nr, start=start, stop=stop, exc_dir=exc_dir, E_WG_func=E_func, H_WG_func=H_func, kc=kc, excite=excite, **kw)

