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

cimport _nf2ff
import numpy as np
import os
from CSXCAD.Utilities import CheckNyDir

cdef class _nf2ff:
    def __cinit__(self, freq, theta, phi, center, numThreads=0, **kw):
        if type(freq) in [float, int]:
            freq = list(float(freq))
        if type(theta) in [float, int]:
            theta = list(float(theta))
        if type(phi) in [float, int]:
            phi = list(float(phi))
        self.thisptr = new cpp_nf2ff(freq, theta, phi, center, numThreads)

        if 'verbose' in kw:
            self.SetVerboseLevel(kw['verbose'])
            del kw['verbose']

        assert len(kw)==0, 'Unknown keyword(s): {}'.format(kw)

    def AnalyseFile(self, e_file, h_file):
        if not os.path.exists(e_file):
            raise Exception('File "{}" does not exist'.format(e_file))
        if not os.path.exists(h_file):
            raise Exception('File "{}" does not exist'.format(e_file))
        cdef string e_fn = e_file.encode('UTF-8')
        cdef string h_fn = h_file.encode('UTF-8')
        with nogil:
            ok = self.thisptr.AnalyseFile(e_fn, h_fn)
        return ok

    def SetMirror(self, mirr_type, ny, pos):
        if mirr_type<=0:
            return
        if not mirr_type<3:
            raise Exception('SetMirror: invalid mirror type!')
        ny = CheckNyDir(ny)
        self.thisptr.SetMirror(mirr_type, ny, pos)

    def SetRadius(self, radius):
        self.thisptr.SetRadius(radius)

    def Write2HDF5(self, filename):
        return self.thisptr.Write2HDF5(filename.encode('UTF-8'))

    def SetVerboseLevel(self, level):
        self.thisptr.SetVerboseLevel(level)
