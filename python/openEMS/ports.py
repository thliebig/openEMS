# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 22:53:39 2015

@author: thorsten
"""

import os
import numpy as np
from CSXCAD.Utilities import CheckNyDir
from openEMS import utilities

from openEMS.physical_constants import *

class UI_data:
    def __init__(self, fns, path, freq, signal_type='pulse', **kw):
        self.path = path
        self.fns  = fns
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
class Port:
    def __init__(self, CSX, port_nr, start, stop, excite, **kw):
        self.CSX      = CSX
        self.number   = port_nr
        self.excite   = excite
        self.start    = np.array(start, np.float)
        self.stop     = np.array(stop, np.float)
        self.Z_ref    = None
        self.U_filenames = []
        self.I_filenames = []

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
        assert self.Z_ref  is not None

        if ref_plane_shift is not None:
            assert hasattr(self, 'beta')
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

# Lumped-Port
class LumpedPort(Port):
    def __init__(self, CSX,  port_nr, R, start, stop, exc_dir, excite=0, **kw):
        super(LumpedPort, self).__init__(CSX, port_nr=port_nr, start=start, stop=stop, excite=excite, **kw)
        self.R = R
        self.exc_ny  = CheckNyDir(exc_dir)

        self.direction = np.sign(self.stop[self.exc_ny]-self.start[self.exc_ny])
        assert self.start[self.exc_ny]!=self.stop[self.exc_ny]

        if self.R > 0:
            lumped_R = CSX.AddLumpedElement(self.lbl_temp.format('resist'), ny=self.exc_ny, caps=True, R=self.R)
        elif self.R==0:
            lumped_R = CSX.AddMetal(self.lbl_temp.format('resist'))

        CSX.AddBox(lumped_R, self.start, self.stop, priority=self.priority)

        if excite!=0:
            exc_vec = np.zeros(3)
            exc_vec[self.exc_ny] = -1*self.direction*excite
            exc = CSX.AddExcitation(self.lbl_temp.format('excite'), exc_type=0, exc_val=exc_vec, delay=self.delay)
            CSX.AddBox(exc, self.start, self.stop, priority=self.priority)

        self.U_filenames = [self.lbl_temp.format('ut'), ]
        u_start = 0.5*(self.start+self.stop)
        u_start[self.exc_ny] = self.start[self.exc_ny]
        u_stop  = 0.5*(self.start+self.stop)
        u_stop[self.exc_ny]  = self.stop[self.exc_ny]
        u_probe = CSX.AddProbe(self.U_filenames[0], p_type=0, weight=-1*self.direction)
        CSX.AddBox(u_probe, u_start, u_stop)

        self.I_filenames = [self.lbl_temp.format('it'), ]
        i_start = np.array(self.start)
        i_start[self.exc_ny] = 0.5*(self.start[self.exc_ny]+self.stop[self.exc_ny])
        i_stop  = np.array(self.stop)
        i_stop[self.exc_ny]  = 0.5*(self.start[self.exc_ny]+self.stop[self.exc_ny])
        i_probe = CSX.AddProbe(self.I_filenames[0], p_type=1, weight=self.direction, norm_dir=self.exc_ny)
        CSX.AddBox(i_probe, i_start, i_stop)

    def CalcPort(self, sim_path, freq, ref_impedance=None, ref_plane_shift=None, signal_type='pulse'):
        if ref_impedance is None:
            self.Z_ref = self.R
        if ref_plane_shift is not None:
            Warning('A lumped port does not support a reference plane shift! Ignoring...')
        super(LumpedPort, self).CalcPort(sim_path, freq, ref_impedance, ref_plane_shift, signal_type)

class MSLPort(Port):
    def __init__(self, CSX, port_nr, metal_prop, start, stop, prop_dir, exc_dir, excite=0, **kw):
        super(MSLPort, self).__init__(CSX, port_nr=port_nr, start=start, stop=stop, excite=excite, **kw)
        self.exc_ny  = CheckNyDir(exc_dir)
        self.prop_ny = CheckNyDir(prop_dir)
        self.direction   = np.sign(stop[self.prop_ny]-start[self.prop_ny])
        self.upside_down = np.sign(stop[self.exc_ny]  -start[self.exc_ny])
        assert (self.start!=self.stop).all()
#        assert stop[self.prop_ny]!=start[self.prop_ny], 'port length in propergation direction may not be zero!'
#        assert stop[self.exc_ny] !=start[self.exc_ny], 'port length in propergation direction may not be zero!'
        assert self.exc_ny!=self.prop_ny

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
        CSX.AddBox( metal_prop, MSL_start, MSL_stop, priority=self.priority )

        mesh = CSX.GetGrid()
        prop_lines = mesh.GetLines(self.prop_ny)
        assert len(prop_lines)>5, 'At least 5 lines in propagation direction required!'
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
            u_probe = CSX.AddProbe(u_name, p_type=0, weight=self.upside_down)
            CSX.AddBox(u_probe, u_start, u_stop)

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
            CSX.AddBox(i_probe, i_start, i_stop)

        if excite!=0:
            excide_pos_idx = np.argmin(np.abs(prop_lines-(self.start[self.prop_ny] + self.feed_shift*self.direction)))
            exc_start = np.array(self.start)
            exc_stop  = np.array(self.stop)
            exc_start[self.prop_ny] = prop_lines[excide_pos_idx]
            exc_stop [self.prop_ny] = prop_lines[excide_pos_idx]
            exc_vec = np.zeros(3)
            exc_vec[self.exc_ny] = -1*self.upside_down*excite
            exc = CSX.AddExcitation(self.lbl_temp.format('excite'), exc_type=0, exc_val=exc_vec, delay=self.delay)
            CSX.AddBox(exc, exc_start, exc_stop, priority=self.priority)

        if self.feed_R>=0 and not np.isinf(self.feed_R):
            R_start = np.array(self.start)
            R_stop  = np.array(self.stop)
            R_stop [self.prop_ny] = R_start[self.prop_ny]
            if self.feed_R==0:
                CSX.AddBox(metal_prop, R_start, R_stop)
            else:
                lumped_R = CSX.AddLumpedElement(self.lbl_temp.format('resist'), ny=self.exc_ny, caps=True, R=self.feed_R)
                CSX.AddBox(lumped_R, R_start, R_stop)

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
    def __init__(self, CSX, port_nr, start, stop, exc_dir, mode_name, excite=0, **kw):
        super(WaveguidePort, self).__init__(CSX, port_nr=port_nr, start=start, stop=stop, excite=excite, **kw)
        self.exc_ny  = CheckNyDir(exc_dir)
        self.ny_P  = (self.exc_ny+1)%3
        self.ny_PP = (self.exc_ny+2)%3
        self.direction = np.sign(stop[self.exc_ny]-start[self.exc_ny])
        self.ref_index = 1

        assert not (self.excite!=0 and stop[self.exc_ny]==start[self.exc_ny]), 'port length in excitation direction may not be zero if port is excited!'

        self.InitMode(mode_name)

        if excite!=0:
            e_start = np.array(start)
            e_stop  = np.array(stop)
            e_stop[self.exc_ny] = e_start[self.exc_ny]
            e_vec = np.ones(3)
            e_vec[self.exc_ny]=0
            exc = CSX.AddExcitation(self.lbl_temp.format('excite'), exc_type=0, exc_val=e_vec, delay=self.delay)
            CSX.AddBox(exc, e_start, e_stop, priority=self.priority)

        # voltage/current planes
        m_start = np.array(start)
        m_stop  = np.array(stop)
        m_start[self.exc_ny] = m_stop[self.exc_ny]
        self.measplane_shift = np.abs(stop[self.exc_ny] - start[self.exc_ny])

        self.U_filenames = [self.lbl_temp.format('ut'), ]

        u_probe = CSX.AddProbe(self.U_filenames[0], p_type=10, mode_function=self.E_func)
        CSX.AddBox(u_probe, m_start, m_stop)

        self.I_filenames = [self.lbl_temp.format('it'), ]
        i_probe = CSX.AddProbe(self.I_filenames[0], p_type=11, weight=self.direction, mode_function=self.H_func)
        CSX.AddBox(i_probe, m_start, m_stop)

    def InitMode(self, wg_mode):
        self.WG_mode = wg_mode
        assert len(self.WG_mode)==4, 'Invalid mode definition'
        self.unit = self.CSX.GetGrid().GetDeltaUnit()
        if self.WG_mode.startswith('TE'):
            self.TE = True
            self.TM = False
        else:
            self.TE = False
            self.TM = True
        self.M = float(self.WG_mode[2])
        self.N = float(self.WG_mode[3])
        self.kc = None
        self.E_func = [0,0,0]
        self.H_func = [0,0,0]

    def CalcPort(self, sim_path, freq, ref_impedance=None, ref_plane_shift=None, signal_type='pulse'):
        k = 2.0*np.pi*freq/C0*self.ref_index
        self.beta = np.sqrt(k**2 - self.kc**2)
        self.ZL = k * Z0 / self.beta    #analytic waveguide impedance
        if ref_impedance is None:
            self.Z_ref = self.ZL
        super(WaveguidePort, self).CalcPort(sim_path, freq, ref_impedance, ref_plane_shift, signal_type)

class RectWGPort(WaveguidePort):
    def __init__(self, CSX, port_nr, start, stop, exc_dir, a, b, mode_name, excite=0, **kw):
        self.WG_size = [a, b]
        super(RectWGPort, self).__init__(CSX, port_nr=port_nr, start=start, stop=stop, exc_dir=exc_dir, mode_name=mode_name, excite=excite, **kw)

    def InitMode(self, wg_mode):
        super(RectWGPort, self).InitMode(wg_mode)
        assert self.TE, 'Currently only TE-modes are supported! Mode found: {}'.format(self.WG_mode)

        # values by David M. Pozar, Microwave Engineering, third edition
        a = self.WG_size[0]
        b = self.WG_size[1]
        self.kc = np.sqrt((self.M*np.pi/a)**2 + (self.N*np.pi/b)**2)

        xyz = 'xyz'
        if self.start[self.ny_P]!=0:
            name_P = '({}-{})'.format(xyz[self.ny_P], self.start[self.ny_P])
        else:
            name_P = xyz[self.ny_P]
        if self.start[self.ny_PP]!=0:
            name_PP = '({}-{})'.format(xyz[self.ny_P], self.start[self.ny_P])
        else:
            name_PP = xyz[self.ny_P]

        a /= self.unit
        b /= self.unit
        if self.N>0:
            self.E_func[self.ny_P]  = '{}*cos({}*{})*sin({}*{})'.format(self.N/b   , self.M*np.pi/a, name_P, self.N*np.pi/b, name_PP)
        if self.M>0:
            self.E_func[self.ny_PP] = '{}*sin({}*{})*cos({}*{})'.format(-1*self.M/a, self.M*np.pi/a, name_P, self.N*np.pi/b, name_PP)

        if self.M>0:
            self.H_func[self.ny_P]  = '{}*sin({}*{})*cos({}*{})'.format(self.M/a, self.M*np.pi/a, name_P, self.N*np.pi/b, name_PP)
        if self.N>0:
            self.H_func[self.ny_PP] = '{}*cos({}*{})*sin({}*{})'.format(self.N/b, self.M*np.pi/a, name_P, self.N*np.pi/b, name_PP)
