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
from libcpp cimport bool

from CSXCAD.CSXCAD cimport _ContinuousStructure, ContinuousStructure

cdef extern from "openEMS/openems.h":
    cdef cppclass _openEMS "openEMS":
        _openEMS() nogil except +
        void Reset()

        void SetNumberOfTimeSteps(unsigned int val)
        void SetCSX(_ContinuousStructure* csx)

        void SetEndCriteria(double val)
        void SetOverSampling(int val)
        void SetCellConstantMaterial(bool val)

        void SetCylinderCoords(bool val)
        void SetupCylinderMultiGrid(string val)

        void SetTimeStepMethod(int val)
        void SetTimeStep(double val)
        void SetTimeStepFactor(double val)
        void SetMaxTime(double val)

        void SetLibraryArguments(vector[string] allOptions) except +

        void Set_BC_Type(int idx, int _type)
        int Get_BC_Type(int idx)
        void Set_BC_PML(int idx, unsigned int size)
        int Get_PML_Size(int idx)
        void Set_Mur_PhaseVel(int idx, double val)

        void SetGaussExcite(double f0, double fc)
        void SetSinusExcite(double f0)
        void SetDiracExcite(double f_max)
        void SetStepExcite(double f_max)
        void SetCustomExcite(string _str, double f0, double fmax)

        void SetAbort(bool val)

        int SetupFDTD() nogil
        void RunFDTD()  nogil

        @staticmethod
        void WelcomeScreen()

cdef class openEMS:
    cdef  _openEMS *thisptr
    cdef readonly ContinuousStructure __CSX      # hold a C++ instance which we're wrapping
