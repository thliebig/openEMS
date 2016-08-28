# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 00:12:43 2015

@author: thorsten
"""

from libcpp.string cimport string
from libcpp cimport bool

from CSXCAD.CSXCAD cimport _ContinuousStructure, ContinuousStructure

cdef extern from "openEMS/openems.h":
    cdef cppclass _openEMS "openEMS":
        _openEMS() except +
        void SetNumberOfTimeSteps(unsigned int val)
        void SetCSX(_ContinuousStructure* csx)

        void SetEndCriteria(double val)
        void SetOverSampling(int val)
        void SetCellConstantMaterial(bool val)

        void SetCylinderCoords(bool val)
        #void SetupCylinderMultiGrid(std::vector<double> val)

        void SetTimeStepMethod(int val)
        void SetTimeStep(double val)
        void SetTimeStepFactor(double val)
        void SetMaxTime(double val)

        void Set_BC_Type(int idx, int _type)
        int Get_BC_Type(int idx)
        void Set_BC_PML(int idx, unsigned int size)
        int Get_PML_Size(int idx)
        void Set_Mur_PhaseVel(int idx, double val)

        void SetGaussExcite(double f0, double fc)

        void SetVerboseLevel(int level)

        int SetupFDTD()
        void RunFDTD()

        @staticmethod
        void WelcomeScreen()

cdef class openEMS:
    cdef  _openEMS *thisptr
    cdef readonly ContinuousStructure CSX      # hold a C++ instance which we're wrapping
