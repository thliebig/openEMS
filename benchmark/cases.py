# -*- coding: utf-8 -*-
#
# Copyright (C) 2026 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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
"""The benchmark cases: each one isolates one cost of the solver.

A case is written to an openEMS XML file at a given size, which sets the cell
count. The number of timesteps is practically unlimited: bench.py stops each
run after a fixed time, so every engine and machine runs equally long. Changing
a case or a size changes the results: bump `CASES_VERSION`, results of
different versions are not compared.

All cases use a sinusoidal excitation: its energy never decays, so no run
stops early on the end criteria (an end criteria of 0 in the XML means the
default of 1e-6), and the probes carry a signal over the whole run. A sinus
only stores two periods, whatever the number of timesteps. With a sinus
openEMS adds its steady-state detection to every case.
"""

import numpy as np

from CSXCAD import ContinuousStructure
from CSXCAD.CSProperties import CSPropLorentzMaterial, CSPropDebyeMaterial
from openEMS import openEMS

CASES_VERSION = 2

SIZES = ('XS', 'S', 'M', 'L')
DEFAULT_SIZES = ('S', 'M')

# stopped by bench.py long before
NRTS = 1e9

# mesh lines per direction of the Cartesian cases: ~50 k, 1 M, 8 M, 30 M cells
LINES = {'XS': 37, 'S': 100, 'M': 200, 'L': 311}

# r, alpha and z cells of the cylindrical case, about the same cell counts
CYL_CELLS = {'XS': (40, 32, 40), 'S': (125, 64, 125), 'M': (250, 128, 250), 'L': (390, 192, 390)}

unit = 1e-3     # drawing unit: mm, the uniform meshes have 1 mm cells
f0   = 5e9      # excitation frequency, 60 mm wavelength


def _setup(bc, **kw):
    FDTD = openEMS(NrTS=NRTS, EndCriteria=1e-30, **kw)
    FDTD.SetSinusExcite(f0)
    FDTD.SetBoundaryCond(bc)
    CSX = ContinuousStructure(**kw)
    FDTD.SetCSX(CSX)
    CSX.GetGrid().SetDeltaUnit(unit)
    return FDTD, CSX


def _uniform(CSX, size):
    """ a cube of 1 mm cells centred at the origin, returns its half width """
    n = LINES[size]
    lines = np.arange(n) - (n - 1) / 2
    mesh = CSX.GetGrid()
    for d in 'xyz':
        mesh.AddLine(d, lines)
    return lines[-1]


def _dipole(CSX):
    """ a short z-directed soft source at the origin with field probes two cells away """
    CSX.AddExcitation('dipole', exc_type=0, exc_val=[0, 0, 1]).AddBox([0, 0, -1], [0, 0, 1])
    CSX.AddProbe('et', p_type=2).AddPoint([2, 1, 0])
    CSX.AddProbe('ht', p_type=3).AddPoint([-1, 2, 1])


def _vacuum(size, bc):
    FDTD, CSX = _setup(bc)
    _uniform(CSX, size)
    _dipole(CSX)
    return FDTD


def ideal(size):
    """ uniform vacuum, PEC and PMC walls: the bare update kernel. Both walls are
        only zeroed operator coefficients, no work in the time-stepping """
    return _vacuum(size, ['PEC', 'PMC', 'PMC', 'PEC', 'PEC', 'PMC'])


def pml(size):
    """ uniform vacuum, PML on all sides: the UPML extension """
    return _vacuum(size, ['PML_8'] * 6)


def _graded(n):
    """ n mesh lines, symmetric to 0: 0.25 mm in the centre, graded by 1.4 to 1 mm """
    half = [0.0]
    d = 0.25
    while len(half) < (n + 1) // 2:
        if len(half) > 12 and d < 1:
            d = min(1.0, d * 1.4)
        half.append(half[-1] + d)
    half = np.array(half)
    return np.concatenate([-half[:0:-1], half])


def inhomogeneous(size):
    """ graded mesh, dielectric, lossy and metal objects, a microstrip line between two
        lumped ports, PML: many distinct coefficients and the usual extensions """
    FDTD, CSX = _setup(['PML_8'] * 6)
    lines = _graded(LINES[size])
    mesh = CSX.GetGrid()
    for d in 'xyz':
        mesh.AddLine(d, lines)
    E = lines[-1]

    # microstrip in the fine centre: substrate, ground plane and a strip between two
    # lumped ports, which sit on a mesh line (a zero-width box between two is not used)
    h = 1.5
    px = lines[np.argmin(np.abs(lines - 0.5*E))]
    CSX.AddMaterial('substrate', epsilon=3.66).AddBox([-0.6*E, -0.6*E, -h], [0.6*E, 0.6*E, 0])
    metal = CSX.AddMetal('metal')
    metal.AddBox([-0.6*E, -0.6*E, -h], [0.6*E, 0.6*E, -h], priority=10)
    metal.AddBox([-px, -1, 0], [px, 1, 0], priority=10)
    FDTD.AddLumpedPort(1, 50, [-px, -1, -h], [-px, 1, 0], 'z', excite=1, priority=5)
    FDTD.AddLumpedPort(2, 50, [ px, -1, -h], [ px, 1, 0], 'z', priority=5)

    CSX.AddMaterial('eps2_2', epsilon=2.2).AddBox([0.2*E, 0.2*E, 0.2*E], [0.5*E, 0.6*E, 0.7*E])
    CSX.AddMaterial('eps10', epsilon=10).AddBox([-0.6*E, 0.2*E, 0.1*E], [-0.2*E, 0.6*E, 0.5*E])
    CSX.AddMaterial('eps4_5', epsilon=4.5).AddBox([-0.3*E, -0.7*E, -0.6*E], [0.3*E, -0.3*E, -0.2*E])
    CSX.AddMaterial('lossy', epsilon=6, kappa=0.5).AddBox([0.2*E, -0.6*E, 0.2*E], [0.6*E, -0.2*E, 0.6*E])
    metal.AddBox([-0.7*E, -0.2*E, -0.5*E], [-0.6*E, 0.2*E, -0.2*E], priority=10)

    CSX.AddProbe('et', p_type=2).AddPoint([0, 0, -h/2])
    return FDTD


def dispersive(size):
    """ uniform mesh, a Lorentz and a Debye block, Mur ABC: the dispersive material extensions """
    FDTD, CSX = _setup(['MUR'] * 6)
    E = _uniform(CSX, size)

    lorentz = CSPropLorentzMaterial(CSX.GetParameterSet(), order=1)
    lorentz.SetName('lorentz')
    lorentz.SetDispersiveMaterialProperty(0, eps_plasma=5e9, eps_pole_freq=3e9, eps_relax=1e-9)
    CSX.AddProperty(lorentz)
    lorentz.AddBox([-0.8*E, -0.8*E, -0.8*E], [-1, 0.8*E, 0.8*E], priority=10)

    debye = CSPropDebyeMaterial(CSX.GetParameterSet(), order=1, epsilon=4)
    debye.SetName('debye')
    debye.SetDispersiveMaterialProperty(0, eps_delta=1, eps_relax=1/(2*np.pi*4e9))
    CSX.AddProperty(debye)
    debye.AddBox([3, -0.8*E, -0.8*E], [0.8*E, 0.8*E, 0.8*E], priority=10)

    _dipole(CSX)
    return FDTD


def cylindrical(size):
    """ closed cylindrical mesh including r=0, a dielectric ring, PEC/Mur: the cylinder
        operator and extension """
    FDTD, CSX = _setup(['PEC', 'PEC', 'PEC', 'PEC', 'MUR', 'MUR'], CoordSystem=1)
    nr, na, nz = CYL_CELLS[size]
    mesh = CSX.GetGrid()
    mesh.AddLine('r', np.arange(nr))
    mesh.AddLine('a', np.linspace(-np.pi, np.pi, na + 1))
    mesh.AddLine('z', np.arange(nz) - (nz - 1) / 2)
    R = nr - 1
    Z = (nz - 1) / 2

    CSX.AddMaterial('ring', epsilon=4).AddBox([0.4*R, -np.pi, -0.6*Z], [0.6*R, np.pi, 0.6*Z])
    CSX.AddExcitation('line', exc_type=0, exc_val=[0, 0, 1]).AddBox([3, 0, -1], [3, 0, 1])
    CSX.AddProbe('et', p_type=2).AddPoint([5, 0.2, 0])
    CSX.AddProbe('ht', p_type=3).AddPoint([2, -0.3, 1])
    return FDTD


CASES = {
    'ideal':         ideal,
    'pml':           pml,
    'inhomogeneous': inhomogeneous,
    'dispersive':    dispersive,
    'cylindrical':   cylindrical,
}

# the engine choice of the CPU engines has no effect on a cylindrical mesh
CARTESIAN_ONLY_ENGINES = ('basic', 'sse', 'sse-compressed')


def write_case(name, size, xml_file):
    """ Write case `name` at `size` to `xml_file` """
    CASES[name](size).Write2XML(xml_file)
