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
import h5py
from openEMS import _nf2ff
from openEMS import utilities

class nf2ff:
    """
    Create an nf2ff recording box. The nf2ff can either record in time-domain
    or frequency-domain. Further more certain directions and boundary condition
    mirroring can be enabled/disabled.

    :param name: str -- Name for this recording box.
    :param start/stop: (3,) array -- Box start/stop coordinates.
    :param directions: (6,) bool array -- Enable/Disables directions.
    :param mirror: (6,) int array -- 0 (Off), 1 (PEC) or 2 (PMC) boundary mirroring
    :param frequency: array like -- List of frequencies (FD-domain recording)
    """
    def __init__(self, CSX, name, start, stop, **kw):
        self.CSX   = CSX
        self.name  = name
        self.start = start
        self.stop  = stop

        self.freq  = None
        self.theta = None
        self.phi   = None
        self.center = None

        self.directions = [True]*6 # all directions by default
        if 'directions' in kw:
            self.directions = kw['directions']
            del kw['directions']
            assert len(self.directions)==6

        self.mirror = [0]*6
        if 'mirror' in kw:
            self.mirror = kw['mirror']
            del kw['mirror']
            assert len(self.mirror)==6

        self.dump_type = 0  # default Et/Ht
        self.dump_mode = 1  # default cell interpolated

        self.freq = None   # broadband recording by defualt
        if 'frequency' in kw:
            self.freq = kw['frequency']
            del kw['frequency']
            self.dump_type = 10  # Ef/Hf

            if np.isscalar(self.freq):
                self.freq = [self.freq]

        self.e_file = '{}_E'.format(self.name)
        self.h_file = '{}_H'.format(self.name)

        self.e_dump = CSX.AddDump(self.e_file, dump_type=self.dump_type  , dump_mode=self.dump_mode, file_type=1, **kw)
        self.h_dump = CSX.AddDump(self.h_file, dump_type=self.dump_type+1, dump_mode=self.dump_mode, file_type=1, **kw)
        if self.freq is not None:
            self.e_dump.SetFrequency(self.freq)
            self.h_dump.SetFrequency(self.freq)

#        print(self.directions)
        for ny in range(3):
            pos = 2*ny
            if self.directions[pos]:
                l_start = np.array(start)
                l_stop  = np.array(stop)
                l_stop[ny] = l_start[ny]
                self.e_dump.AddBox(l_start, l_stop)
                self.h_dump.AddBox(l_start, l_stop)
            if self.directions[pos+1]:
                l_start = np.array(start)
                l_stop  = np.array(stop)
                l_start[ny] = l_stop[ny]
                self.e_dump.AddBox(l_start, l_stop)
                self.h_dump.AddBox(l_start, l_stop)

    def CalcNF2FF(self, sim_path, freq, theta, phi, radius=1, center=[0,0,0], outfile=None, read_cached=False, verbose=0):
        """ CalcNF2FF(sim_path, freq, theta, phi, center=[0,0,0], outfile=None, read_cached=True, verbose=0):

        Calculate the far-field after the simulation is done.

        :param sim_path: str -- Simulation path
        :param freq: array like -- list of frequency for transformation
        :param theta/phi: array like -- Theta/Phi angles to calculate the far-field
        :param radius: float -- Radius to calculate the far-field (default is 1m)
        :param center: (3,) array -- phase center, must be inside the recording box
        :param outfile: str -- File to save results in. (defaults to recording name)
        :param read_cached: bool -- enable/disable read already existing results (default off)
        :param verbose: int -- set verbose level (default 0)

        :returns: nf2ff_results class instance
        """
        if np.isscalar(freq):
            freq = [freq]
        self.freq  = freq
        if np.isscalar(theta):
            theta = [theta]
        self.theta = theta
        if np.isscalar(phi):
            phi = [phi]
        self.phi   = phi
        self.center = center

        if outfile is None:
            fn = os.path.join(sim_path, self.name + '.h5')
        else:
            fn = os.path.join(sim_path,  outfile)
        if  not read_cached or not os.path.exists(fn):
            nfc = _nf2ff._nf2ff(self.freq, np.deg2rad(theta), np.deg2rad(phi), center, verbose=verbose)

            for ny in range(3):
                nfc.SetMirror(self.mirror[2*ny]  , ny, self.start[ny])
                nfc.SetMirror(self.mirror[2*ny+1], ny, self.stop[ny])

            nfc.SetRadius(radius)

            for n in range(6):
                fn_e = os.path.join(sim_path, self.e_file + '_{}.h5'.format(n))
                fn_h = os.path.join(sim_path, self.h_file + '_{}.h5'.format(n))
                if os.path.exists(fn_e) and os.path.exists(fn_h):
                    if not nfc.AnalyseFile(fn_e, fn_h):
                        raise Exception('CalcNF2FF:: Unable to analyse files!')

            nfc.Write2HDF5(fn)

        result = nf2ff_results(fn)
        if result.phi is not None:
            if not np.abs((result.r-radius)/radius)<1e-6:
                raise Exception('Radius does not match. Did you read an invalid chached result? Try "read_cached=False"')
            if not utilities.Check_Array_Equal(np.rad2deg(result.theta), self.theta, 1e-4):
                raise Exception('Theta array does not match. Did you read an invalid chached result? Try "read_cached=False"')
            if not utilities.Check_Array_Equal(np.rad2deg(result.phi), self.phi, 1e-4):
                raise Exception('Phi array does not match. Did you read an invalid chached result? Try "read_cached=False"')
            if not utilities.Check_Array_Equal(result.freq, self.freq, 1e-6, relative=True):
                raise Exception('Frequency array does not match. Did you read an invalid chached result? Try "read_cached=False"')
        return result

class nf2ff_results:
    """
    nf2ff result class containing all results obtained by the nf2ff calculation.
    Usueally returned from nf2ff.CalcNF2FF

    Available attributes:

    * `fn`   : file name
    * `theta`: theta angles
    * `phi`  : phi angles
    * `r`    : radius
    * `freq` : frequencies
    * `Dmax` : directivity over frequency
    * `Prad` : total radiated power over frequency

    * `E_theta` : theta component of electric field over frequency/theta/phi
    * `E_phi`   : phi   component of electric field over frequency/theta/phi
    * `E_norm`  : abs   component of electric field over frequency/theta/phi
    * `E_cprh`  : theta component of electric field over frequency/theta/phi
    * `E_cplh`  : theta component of electric field over frequency/theta/phi
    * `P_rad`   : radiated power (S) over frequency/theta/phi
    """
    def __init__(self, fn):
        self.fn  = fn
        h5_file  = h5py.File(fn, 'r')
        mesh_grp = h5_file['Mesh']
        self.phi   = np.array(mesh_grp['phi'])
        self.theta = np.array(mesh_grp['theta'])
        self.r     = np.array(mesh_grp['r'])

        data  = h5_file['nf2ff']
        self.freq = np.array(data.attrs['Frequency'])

        self.Dmax = np.array(data.attrs['Dmax'])
        self.Prad = np.array(data.attrs['Prad'])

        THETA, PHI = np.meshgrid(self.theta, self.phi, indexing='ij')
        cos_phi = np.cos(PHI)
        sin_phi = np.sin(PHI)

        self.E_theta = []
        self.E_phi   = []
        self.P_rad   = []
        self.E_norm  = []
        self.E_cprh  = []
        self.E_cplh  = []
        for n in range(len(self.freq)):
            E_theta = np.array(h5_file['/nf2ff/E_theta/FD/f{}_real'.format(n)]) + 1j*np.array(h5_file['/nf2ff/E_theta/FD/f{}_imag'.format(n)])
            E_theta = np.swapaxes(E_theta, 0, 1)
            E_phi   = np.array(h5_file['/nf2ff/E_phi/FD/f{}_real'.format(n)])   + 1j*np.array(h5_file['/nf2ff/E_phi/FD/f{}_imag'.format(n)])
            E_phi   = np.swapaxes(E_phi, 0, 1)
            self.P_rad  .append(np.swapaxes(np.array(h5_file['/nf2ff/P_rad/FD/f{}'.format(n)]), 0, 1))

            self.E_theta.append(E_theta)
            self.E_phi  .append(E_phi)
            self.E_norm .append(np.sqrt(np.abs(E_theta)**2 + np.abs(E_phi)**2))
            self.E_cprh .append((cos_phi+1j*sin_phi) * (E_theta+1j*E_phi)/np.sqrt(2.0))
            self.E_cplh .append((cos_phi-1j*sin_phi) * (E_theta-1j*E_phi)/np.sqrt(2.0))
