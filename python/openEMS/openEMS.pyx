# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 23:50:24 2015

@author: thorsten
"""

import os, sys, shutil
import numpy as np
cimport openEMS
from . import ports, nf2ff

cdef class openEMS:
    @staticmethod
    def WelcomeScreen():
        _openEMS.WelcomeScreen()

    """ __cinit__
    NrTS:           max. number of timesteps to simulate (e.g. default=1e9)
    EndCriteria:    end criteria, e.g. 1e-5, simulations stops if energy has decayed by this value (<1e-4 is recommended, default=1e-5)
    MaxTime:        max. real time in seconds to simulate
    OverSampling:   nyquist oversampling of time domain dumps
    CoordSystem:    choose coordinate system (0 Cartesian, 1 Cylindrical)
    #MultiGrid:      define a cylindrical sub-grid radius
    TimeStep:       force to use a given timestep (dangerous!)
    TimeStepFactor: reduce the timestep by a given factor (>0 to <=1)
    TimeStepMethod: 1 or 3 chose timestep method (1=CFL, 3=Rennigs (default))
    CellConstantMaterial: set to 1 to assume a material is constant inside a cell (material probing in cell center)
    """
    def __cinit__(self, *args, **kw):
        self.thisptr = new _openEMS()
        self.CSX = None

        if 'NrTS' in kw:
            self.SetNumberOfTimeSteps(kw['NrTS'])
        else:
            self.SetNumberOfTimeSteps(1e9)
        if 'EndCriteria' in kw:
            self.SetEndCriteria(kw['EndCriteria'])
        if 'MaxTime' in kw:
            self.SetMaxTime(kw['MaxTime'])
        if 'OverSampling' in kw:
            self.SetOverSampling(kw['OverSampling'])
        if 'CoordSystem' in kw:
            self.SetCoordSystem(kw['CoordSystem'])
        if 'TimeStep' in kw:
            self.SetTimeStep(kw['TimeStep'])
        if 'TimeStepFactor' in kw:
            self.SetTimeStepFactor(kw['TimeStepFactor'])
        if 'TimeStepMethod' in kw:
            self.SetTimeStepMethod(kw['TimeStepMethod'])
        if 'CellConstantMaterial' in kw:
            self.SetCellConstantMaterial(kw['CellConstantMaterial'])

    def __dealloc__(self):
        del self.thisptr
        if self.CSX is not None:
            self.CSX.thisptr = NULL

    def SetNumberOfTimeSteps(self, val):
        self.thisptr.SetNumberOfTimeSteps(val)

    def SetEndCriteria(self, val):
        self.thisptr.SetEndCriteria(val)

    def SetOverSampling(self, val):
        self.thisptr.SetOverSampling(val)

    def SetCellConstantMaterial(self, val):
        self.thisptr.SetCellConstantMaterial(val)

    def SetCoordSystem(self, val):
        assert (val==0 or val==1)
        if val==0:
            pass
        elif val==1:
            self.SetCylinderCoords()

    def SetCylinderCoords(self):
        self.thisptr.SetCylinderCoords(True)

    def SetTimeStepMethod(self, val):
        self.thisptr.SetTimeStepMethod(val)
    def SetTimeStep(self, val):
        self.thisptr.SetTimeStep(val)
    def SetTimeStepFactor(self, val):
        self.thisptr.SetTimeStepFactor(val)
    def SetMaxTime(self, val):
        self.thisptr.SetMaxTime(val)

    def SetGaussExcite(self, f0, fc):
        self.thisptr.SetGaussExcite(f0, fc)


    def SetBoundaryCond(self, BC):
        assert len(BC)==6
        for n in range(len(BC)):
            if type(BC[n])==int:
                self.thisptr.Set_BC_Type(n, BC[n])
                continue
            if BC[n] in ['PEC', 'PMC', 'MUR']:
                self.thisptr.Set_BC_Type(n, ['PEC', 'PMC', 'MUR'].index(BC[n]))
                continue
            if BC[n].startswith('PML_'):
                size = int(BC[n].strip('PML_'))
                self.thisptr.Set_BC_PML(n, size)
                continue
            raise Exception('Unknown boundary condition')

    def AddLumpedPort(self, port_nr, R, start, stop, p_dir, excite=0, **kw):
        assert self.CSX is not None
        return ports.LumpedPort(self.CSX, port_nr, R, start, stop, p_dir, excite, **kw)

    def AddRectWaveGuidePort(self, port_nr, start, stop, p_dir, a, b, mode_name, excite=0, **kw):
        assert self.CSX is not None
        return ports.RectWGPort(self.CSX, port_nr, start, stop, p_dir, a, b, mode_name, excite, **kw)

    def AddMSLPort(self, port_nr, metal_prop, start, stop, prop_dir, exc_dir, excite=0, **kw):
        assert self.CSX is not None
        return ports.MSLPort(self.CSX, port_nr, metal_prop, start, stop, prop_dir, exc_dir, excite, **kw)

    def CreateNF2FFBox(self, name='nf2ff', start=None, stop=None, **kw):
        assert self.CSX is not None
        directions = [True]*6
        mirror     = [0]*6
        BC_size = [0]*6
        BC_type = [0]*6
        for n in range(6):
            BC_type[n] = self.thisptr.Get_BC_Type(n)
            if BC_type[n]==0:
                directions[n]= False
                mirror[n]    = 1  # PEC mirror
            elif BC_type[n]==1:
                directions[n]= False
                mirror[n]    = 2  # PMC mirror
            elif BC_type[n]==2:
                BC_size[n] = 2
            elif BC_type[n]==3:
                BC_size[n] = self.thisptr.Get_PML_Size(n)+1

        if start is None or stop is None:
            grid = self.CSX.GetGrid()
            assert grid.IsValid(), 'Error::CreateNF2FFBox: Grid is invalid'
            start = np.zeros(3)
            stop  = np.zeros(3)
            for n in range(3):
                l = grid.GetLines(n)
                BC_type = self.thisptr.Get_BC_Type(2*n)
                assert len(l)>(BC_size[2*n]+BC_size[2*n+1]), 'Error::CreateNF2FFBox: not enough lines in some direction'
                start[n] = l[BC_size[2*n]]
                stop[n]  = l[-1*BC_size[2*n+1]-1]
        return nf2ff.nf2ff(self.CSX, name, start, stop, directions=directions, mirror=mirror, **kw)

    def SetCSX(self, ContinuousStructure CSX):
        self.CSX = CSX
        self.thisptr.SetCSX(CSX.thisptr)

    def Run(self, sim_path, cleanup=False, setup_only=False, verbose=None):
        if cleanup and os.path.exists(sim_path):
            shutil.rmtree(sim_path)
            os.mkdir(sim_path)
        if not os.path.exists(sim_path):
            os.mkdir(sim_path)
        cwd = os.getcwd()
        os.chdir(sim_path)
        if verbose is not None:
            self.thisptr.SetVerboseLevel(verbose)
        assert os.getcwd() == sim_path
        _openEMS.WelcomeScreen()
        cdef int EC
        EC = self.thisptr.SetupFDTD()
        if EC!=0:
            print('Run: Setup failed, error code: {}'.format(EC))
        if setup_only or EC!=0:
            return EC
        self.thisptr.RunFDTD()
