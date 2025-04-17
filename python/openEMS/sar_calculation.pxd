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

from libcpp.string cimport string
from libcpp cimport bool

cdef extern from "openEMS/sar_calculation.h":
    cdef cppclass _SAR_Calculation "SAR_Calculation":
        _SAR_Calculation() nogil except +

        void SetDebugLevel(int level)
        void SetAveragingMass(float mass)
        bool SetAveragingMethod(string method, bool silent)
        bool CalcFromHDF5(string h5_fn, string out_name) nogil


cdef class SAR_Calculation:
    cdef _SAR_Calculation *thisptr
