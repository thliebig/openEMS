# -*- coding: utf-8 -*-
"""
Created on Sun Dec 20 22:42:19 2015

@author: thorsten
"""

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
        assert os.path.exists(e_file)
        assert os.path.exists(h_file)
        return self.thisptr.AnalyseFile(e_file.encode('UTF-8'), h_file.encode('UTF-8'))

    def SetMirror(self, mirr_type, ny, pos):
        if mirr_type<=0:
            return
        assert mirr_type<3
        ny = CheckNyDir(ny)
        self.thisptr.SetMirror(mirr_type, ny, pos)

    def Write2HDF5(self, filename):
        return self.thisptr.Write2HDF5(filename.encode('UTF-8'))

    def SetVerboseLevel(self, level):
        self.thisptr.SetVerboseLevel(level)
