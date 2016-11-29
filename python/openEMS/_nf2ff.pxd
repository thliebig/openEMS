# -*- coding: utf-8 -*-
"""
Created on Sun Dec 20 22:43:35 2015

@author: thorsten
"""

from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.complex cimport complex
from libcpp cimport bool
cimport cython.numeric

cdef extern from "openEMS/nf2ff.h":
    cdef cppclass cpp_nf2ff "nf2ff":
        cpp_nf2ff(vector[float] freq, vector[float] theta, vector[float] phi, vector[float] center, unsigned int numThreads) except +

        bool AnalyseFile(string E_Field_file, string H_Field_file)

        void SetRadius(float radius)
        void SetPermittivity(vector[float] permittivity);
        void SetPermeability(vector[float] permeability);

        void SetMirror(int _type, int _dir, float pos);

        double GetTotalRadPower(size_t f_idx)
        double GetMaxDirectivity(size_t f_idx)

        complex[double]** GetETheta(size_t f_idx)
        complex[double]** GetEPhi(size_t f_idx)
        double** GetRadPower(size_t f_idx)

        bool Write2HDF5(string filename)

        void SetVerboseLevel(int level)

cdef class _nf2ff:
    cdef  cpp_nf2ff *thisptr
