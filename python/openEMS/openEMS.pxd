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
from libcpp cimport bool

from CSXCAD.CSXCAD cimport _ContinuousStructure, ContinuousStructure

cdef extern from "openEMS/openems.h":
    cdef cppclass _openEMS "openEMS":
        _openEMS() nogil except +
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

        void SetNumberOfThreads(int val)

        void Set_BC_Type(int idx, int _type)
        int Get_BC_Type(int idx)
        void Set_BC_PML(int idx, unsigned int size)
        int Get_PML_Size(int idx)
        void Set_Mur_PhaseVel(int idx, double val)

        void SetGaussExcite(double f0, double fc)

        void SetAbort(bool val)

        void SetVerboseLevel(int level)
        void DebugPEC()      nogil
        void DebugMaterial() nogil
        void DebugCSX()      nogil

        int SetupFDTD() nogil
        void RunFDTD()  nogil

        @staticmethod
        void WelcomeScreen()

cdef class openEMS:
    cdef  _openEMS *thisptr
    cdef readonly ContinuousStructure __CSX      # hold a C++ instance which we're wrapping
