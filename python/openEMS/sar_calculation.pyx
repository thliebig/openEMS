# -*- coding: utf-8 -*-
#
# Copyright (C) 2025 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

cimport openEMS.sar_calculation
import os

cdef class SAR_Calculation:
    def __cinit__(self):
        self.thisptr = new _SAR_Calculation()

    def SetDebugLevel(self, level):
        self.thisptr.SetDebugLevel(level)

    def SetAveragingMass(self, mass):
        ### Set the averaging mass in g
        # convert from g to kg
        self.thisptr.SetAveragingMass(float(mass)/1000)

    def SetAveragingMethod(self, method, silent=True):
         #cdef string c_method = method.encode('UTF-8')
         return self.thisptr.SetAveragingMethod(method.encode('UTF-8'), silent)

    def CalcFromHDF5(self, h5_fn, out_name):
        if not os.path.exists(h5_fn):
            raise Exception('File "{}" does not exist'.format(h5_fn))
        cdef string in_fn = h5_fn.encode('UTF-8')
        cdef string out_fn = out_name.encode('UTF-8')
        with nogil:
            ok = self.thisptr.CalcFromHDF5(in_fn, out_fn)
        return ok

