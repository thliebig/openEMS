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

from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.complex cimport complex
from libcpp cimport bool
cimport cython.numeric

cdef extern from "openEMS/nf2ff.h":
    cdef cppclass cpp_nf2ff "nf2ff":
        cpp_nf2ff(vector[float] freq, vector[float] theta, vector[float] phi, vector[float] center, unsigned int numThreads) nogil except +

        bool AnalyseFile(string E_Field_file, string H_Field_file) nogil

        void SetRadius(float radius)
        void SetPermittivity(vector[float] permittivity);
        void SetPermeability(vector[float] permeability);

        void SetMirror(int _type, int _dir, float pos);

        double GetTotalRadPower(size_t f_idx)
        double GetMaxDirectivity(size_t f_idx)

        complex[double]** GetETheta(size_t f_idx)
        complex[double]** GetEPhi(size_t f_idx)
        double** GetRadPower(size_t f_idx)

        bool Write2HDF5(string filename) nogil

        void SetVerboseLevel(int level)

cdef class _nf2ff:
    cdef  cpp_nf2ff *thisptr
