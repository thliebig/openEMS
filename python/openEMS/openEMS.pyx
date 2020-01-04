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

import os, sys, shutil
import numpy as np
cimport openEMS
from . import ports, nf2ff, automesh

from CSXCAD.Utilities import GetMultiDirs

cdef class openEMS:
    """ openEMS

    This class is the main control class for the FDTD options and setup and
    to run the final simulation.

    Examples
    --------

    >>> CSX = CSXCAD.ContinuousStructure()
    >>>
    >>> grid = CSX.GetGrid()
    >>> grid.SetLines('x', np.arange(-50,50,1))
    >>> grid.SetLines('y', np.arange(-50,50,1))
    >>> grid.SetLines('z', np.arange(-2,2.1,1))
    >>> grid.SetDeltaUnit(1e-3)
    >>>
    >>> FDTD = openEMS(NrTS=1e4, EndCriteria=1e-4)
    >>>
    >>> FDTD.SetCSX(CSX)
    >>> FDTD.SetBoundaryCond(['PML_8', 'PML_8', 'PML_8', 'PML_8', 'PEC', 'PEC'])
    >>> FDTD.SetGaussExcite(0, 10e9)
    >>>
    >>> FDTD.AddLumpedPort(port_nr=1, R=50, start=[10, 0, -2], stop=[10, 0, 2], p_dir='z', excite=1)
    >>>
    >>> FDTD.Run(sim_path='/tmp/test')

    :param NrTS:           max. number of timesteps to simulate (e.g. default=1e9)
    :param EndCriteria:    end criteria, e.g. 1e-5, simulations stops if energy has decayed by this value (<1e-4 is recommended, default=1e-5)
    :param MaxTime:        max. real time in seconds to simulate
    :param OverSampling:   nyquist oversampling of time domain dumps
    :param CoordSystem:    choose coordinate system (0 Cartesian, 1 Cylindrical)
    :param MultiGrid:      define a cylindrical sub-grid radius
    :param TimeStep:       force to use a given timestep (dangerous!)
    :param TimeStepFactor: reduce the timestep by a given factor (>0 to <=1)
    :param TimeStepMethod: 1 or 3 chose timestep method (1=CFL, 3=Rennigs (default))
    :param CellConstantMaterial: set to 1 to assume a material is constant inside a cell (material probing in cell center)
    """
    @staticmethod
    def WelcomeScreen():
        """
        Show the openEMS welcome screen.
        """
        _openEMS.WelcomeScreen()

    def __cinit__(self, *args, **kw):
        self.thisptr = new _openEMS()
        self.__CSX = None

        if 'NrTS' in kw:
            self.SetNumberOfTimeSteps(kw['NrTS'])
            del kw['NrTS']
        else:
            self.SetNumberOfTimeSteps(1e9)
        if 'EndCriteria' in kw:
            self.SetEndCriteria(kw['EndCriteria'])
            del kw['EndCriteria']
        if 'MaxTime' in kw:
            self.SetMaxTime(kw['MaxTime'])
            del kw['MaxTime']
        if 'OverSampling' in kw:
            self.SetOverSampling(kw['OverSampling'])
            del kw['OverSampling']
        if 'CoordSystem' in kw:
            self.SetCoordSystem(kw['CoordSystem'])
            del kw['CoordSystem']
        if 'TimeStep' in kw:
            self.SetTimeStep(kw['TimeStep'])
            del kw['TimeStep']
        if 'TimeStepFactor' in kw:
            self.SetTimeStepFactor(kw['TimeStepFactor'])
            del kw['TimeStepFactor']
        if 'TimeStepMethod' in kw:
            self.SetTimeStepMethod(kw['TimeStepMethod'])
            del kw['TimeStepMethod']
        if 'CellConstantMaterial' in kw:
            self.SetCellConstantMaterial(kw['CellConstantMaterial'])
            del kw['CellConstantMaterial']
        if 'MultiGrid' in kw:
            self.SetMultiGrid(kw['MultiGrid'])
            del kw['MultiGrid']

        assert len(kw)==0, 'Unknown keyword arguments: "{}"'.format(kw)

    def __dealloc__(self):
        del self.thisptr
        if self.__CSX is not None:
            self.__CSX.thisptr = NULL

    def SetNumberOfTimeSteps(self, val):
        """ SetNumberOfTimeSteps(val)

        Set the number of timesteps. E.g. 5e4 (default is 1e9)
        """
        self.thisptr.SetNumberOfTimeSteps(val)

    def SetEndCriteria(self, val):
        """ SetEndCriteria(val)

        Set the end critera value. E.g. 1e-6 for -60dB
        """
        self.thisptr.SetEndCriteria(val)

    def SetOverSampling(self, val):
        """ SetOverSampling(val)

        Set the time domain signal oversampling as multiple of the Nyquist-rate.
        """
        self.thisptr.SetOverSampling(val)

    def SetCellConstantMaterial(self, val):
        """ SetCellConstantMaterial(val)

        Set cell material averaging to assume constant material inside each primary cell. (Advanced option)

        :param val: bool -- Enable or Disable (default disabled)
        """
        self.thisptr.SetCellConstantMaterial(val)

    def SetCoordSystem(self, val):
        """ SetCoordSystem(val)

        Set the coordinate system. 0 --> Cartesian (default), 1 --> cylindrical
        """
        if not (val==0 or val==1):
            raise Exception('SetCoordSystem: Invalid coordinate system')
        if val==0:
            self.thisptr.SetCylinderCoords(False)
        elif val==1:
            self.thisptr.SetCylinderCoords(True)

    def SetMultiGrid(self, radii):
        """ SetMultiGrid(radii)

        Define radii at which a cylindrical multi grid should be defined.

        :param radii: array like, multigrid radii

        See Also
        --------
        openEMS.SetCylinderCoords
        """
        if not len(radii)>0:
            raise Exception('SetMultiGrid: invalid multi grid definition')

        grid_str = ','.join(['{}'.format(x) for x in radii])
        self.thisptr.SetupCylinderMultiGrid(grid_str.encode('UTF-8'))

    def SetCylinderCoords(self):
        """ SetCylinderCoords()

        Enable use of cylindircal coordinates.

        See Also
        --------
        openEMS.SetMultiGrid
        """
        self.thisptr.SetCylinderCoords(True)

    def SetTimeStepMethod(self, val):
        """ SetTimeStepMethod(val)

        Set the time step calculation method. (Advanced option)

        Options:

        * 1: CFL criteria
        * 3: Advanced Rennings criteria (default)

        :param val: int -- 1 or 3 (See above)
        """
        self.thisptr.SetTimeStepMethod(val)

    def SetTimeStep(self, val):
        """ SetTimeStep(val)

        Set/force the timestep. (Advanced option)

        It is highly recommended to not use this method! You may use the
        SetTimeStepFactor instead to reduce the time step if necessary!
        """
        self.thisptr.SetTimeStep(val)

    def SetTimeStepFactor(self, val):
        """ SetTimeStepFactor(val)

        Set a time step factor (>0..1) to increase FDTD stability.

        :param val: float -- >0..1
        """
        self.thisptr.SetTimeStepFactor(val)

    def SetMaxTime(self, val):
        """ SetMaxTime(val)

        Set max simulation time for a max. number of timesteps.
        """
        self.thisptr.SetMaxTime(val)

    def SetGaussExcite(self, f0, fc):
        """ SetGaussExcite(f0, fc)

        Set a Gaussian pulse as excitation signal.

        :param f0: float -- Center frequency in Hz.
        :param fc: float -- -20dB bandwidth in Hz.
        """
        self.thisptr.SetGaussExcite(f0, fc)


    def SetBoundaryCond(self, BC):
        """ SetBoundaryCond(BC)

        Set the boundary conditions for all six FDTD directions.

        Options:

        * 0 or 'PEC' : perfect electric conductor (default)
        * 1 or 'PMC' : perfect magnetic conductor, useful for symmetries
        * 2 or 'MUR' : simple MUR absorbing boundary conditions
        * 3 or 'PML-8' : PML absorbing boundary conditions

        :param BC: (8,) array or list -- see options above
        """
        if not len(BC)==6:
            raise Exception('Invalid boundary condition size!')
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
        """ AddLumpedPort(port_nr, R, start, stop, p_dir, excite=0, **kw)

        Add a lumped port wit the given values and location.

        See Also
        --------
        openEMS.ports.LumpedPort
        """
        if self.__CSX is None:
            raise Exception('AddLumpedPort: CSX is not set!')
        port = ports.LumpedPort(self.__CSX, port_nr, R, start, stop, p_dir, excite, **kw)
        edges2grid = kw.get('edges2grid', None)
        if edges2grid is not None:
            grid = self.__CSX.GetGrid()
            for n in GetMultiDirs(edges2grid):
                grid.AddLine(n, start[n])
                if start[n] != stop[n]:
                    grid.AddLine(n, stop[n])
        return port

    def AddWaveGuidePort(self, port_nr, start, stop, p_dir, E_func, H_func, kc, excite=0, **kw):
        """ AddWaveGuidePort(self, port_nr, start, stop, p_dir, E_func, H_func, kc, excite=0, **kw)

        Add a arbitrary waveguide port.

        See Also
        --------
        openEMS.ports.WaveguidePort
        """
        if self.__CSX is None:
            raise Exception('AddWaveGuidePort: CSX is not set!')
        return ports.WaveguidePort(self.__CSX, port_nr, start, stop, p_dir, E_func, H_func, kc, excite, **kw)

    def AddRectWaveGuidePort(self, port_nr, start, stop, p_dir, a, b, mode_name, excite=0, **kw):
        """ AddRectWaveGuidePort(port_nr, start, stop, p_dir, a, b, mode_name, excite=0, **kw)

        Add a rectilinear waveguide port.

        See Also
        --------
        openEMS.ports.RectWGPort
        """
        if self.__CSX is None:
            raise Exception('AddRectWaveGuidePort: CSX is not set!')
        return ports.RectWGPort(self.__CSX, port_nr, start, stop, p_dir, a, b, mode_name, excite, **kw)

    def AddMSLPort(self, port_nr, metal_prop, start, stop, prop_dir, exc_dir, excite=0, **kw):
        """ AddMSLPort(port_nr, metal_prop, start, stop, prop_dir, exc_dir, excite=0, **kw)

        Add a microstrip transmission line port.

        See Also
        --------
        openEMS.ports.MSLPort
        """
        if self.__CSX is None:
            raise Exception('AddMSLPort: CSX is not set!')
        return ports.MSLPort(self.__CSX, port_nr, metal_prop, start, stop, prop_dir, exc_dir, excite, **kw)

    def CreateNF2FFBox(self, name='nf2ff', start=None, stop=None, **kw):
        """ CreateNF2FFBox(name='nf2ff', start=None, stop=None, **kw)

        Create a near-field to far-field box.

        This method will automatically adept the recording box to the current
        FDTD grid and boundary conditions.

        Notes
        -----
        * Make sure the mesh grid and all boundary conditions are finially defined.

        See Also
        --------
        openEMS.nf2ff.nf2ff
        """
        if self.__CSX is None:
            raise Exception('CreateNF2FFBox: CSX is not set!')
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
            grid = self.__CSX.GetGrid()
            if not grid.IsValid():
                raise Exception('Error::CreateNF2FFBox: Grid is invalid')
            start = np.zeros(3)
            stop  = np.zeros(3)
            for n in range(3):
                l = grid.GetLines(n)
                BC_type = self.thisptr.Get_BC_Type(2*n)
                if not len(l)>(BC_size[2*n]+BC_size[2*n+1]):
                    raise Exception('Error::CreateNF2FFBox: not enough lines in some direction')
                start[n] = l[BC_size[2*n]]
                stop[n]  = l[-1*BC_size[2*n+1]-1]
        return nf2ff.nf2ff(self.__CSX, name, start, stop, directions=directions, mirror=mirror, **kw)

    def SetCSX(self, ContinuousStructure CSX):
        """ SetCSX(CSX)

        Set the CSXCAD Continuous Structure for CAD data handling.

        See Also
        --------
        CSXCAD.ContinuousStructure
        """
        self.__CSX = CSX
        self.thisptr.SetCSX(CSX.thisptr)

    def GetCSX(self):
        return self.__CSX

    def AddEdges2Grid(self, dirs, primitives=None, properties=None, **kw):
        """ AddEdges2Grid(primitives, dirs, **kw)

        Add the edges of the given primitives to the FDTD grid.

        :param dirs: primitives -- one or more primitives
        :param dirs: str -- 'x','y','z' or 'xy', 'yz' or 'xyz' or 'all'
        """
        csx = self.GetCSX()
        if csx is None:
            raise Exception('AddEdges2Grid: Unable to access CSX!')
        prim_list = []
        if primitives is not None and  type(primitives) is not list:
            prim_list.append(primitives)
        elif primitives is not None:
            prim_list += primitives

        if properties is not None and  type(properties) is not list:
            prim_list += properties.GetAllPrimitives()
        elif properties is not None:
            for prop in properties:
                prim_list += prop.GetAllPrimitives()

        grid = csx.GetGrid()
        for prim in prim_list:
            hint = automesh.mesh_hint_from_primitive(prim, dirs, **kw)
            if hint is None:
                continue
            for n in range(3):
                if hint[n] is None:
                    continue
                grid.AddLine(n, hint[n])

    def Run(self, sim_path, cleanup=False, setup_only=False, debug_pec=False, verbose=None, **kw):
        """ Run(sim_path, cleanup=False, setup_only=False, verbose=None)

        Run the openEMS FDTD simulation.

        :param sim_path: str -- path to run in and create result data
        :param cleanup: bool -- remove exisiting sim_path to cleanup old results
        :param setup_only: bool -- only perform FDTD setup, do not run simulation
        :param verbose: int -- set the openEMS verbosity level 0..3

        Additional keyword parameter:
        :param numThreads: int -- set the number of threads (default 0 --> max)
        """
        if cleanup and os.path.exists(sim_path):
            shutil.rmtree(sim_path)
            os.mkdir(sim_path)
        if not os.path.exists(sim_path):
            os.mkdir(sim_path)
        os.chdir(sim_path)
        if verbose is not None:
            self.thisptr.SetVerboseLevel(verbose)
        if debug_pec:
            with nogil:
                self.thisptr.DebugPEC()
        if 'numThreads' in kw:
            self.thisptr.SetNumberOfThreads(int(kw['numThreads']))
        assert os.getcwd() == sim_path
        _openEMS.WelcomeScreen()
        cdef int EC
        with nogil:
            EC = self.thisptr.SetupFDTD()
        if EC!=0:
            print('Run: Setup failed, error code: {}'.format(EC))
        if setup_only or EC!=0:
            return EC
        with nogil:
            self.thisptr.RunFDTD()

    def SetAbort(self, val):
        self.thisptr.SetAbort(val)
