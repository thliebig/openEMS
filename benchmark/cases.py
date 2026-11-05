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

Each case also runs with a field dump, see `DUMPS`: the dumps are part of the
case, a change of them bumps `CASES_VERSION` too.
"""

import numpy as np

from CSXCAD import ContinuousStructure
from CSXCAD.CSProperties import CSPropLorentzMaterial, CSPropDebyeMaterial
from openEMS import openEMS

CASES_VERSION = 3

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
    """ a cube of 1 mm cells around the origin, returns its half width. The lines
        are integers, so the dipole and the probes sit on the mesh: with an even
        number of lines the cube is half a cell off centre """
    n = LINES[size]
    lines = np.arange(n) - n // 2
    mesh = CSX.GetGrid()
    for d in 'xyz':
        mesh.AddLine(d, lines)
    return lines[-1]


def _dipole(CSX):
    """ a short z-directed soft source at the origin with field probes two cells away """
    CSX.AddExcitation('dipole', exc_type=0, exc_val=[0, 0, 1]).AddBox([0, 0, -1], [0, 0, 1])
    CSX.AddProbe('E_probe', p_type=2).AddPoint([2, 1, 0])
    CSX.AddProbe('H_probe', p_type=3).AddPoint([-1, 2, 1])


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

    CSX.AddProbe('E_probe', p_type=2).AddPoint([0, 0, -h/2])
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
    """ closed cylindrical mesh including r=0, a dielectric tube, PEC/Mur: the cylinder
        operator and extension. The tube runs through both Mur ends: a first-order Mur
        ABC loads an evanescent TM field with a negative resistance, and the ends of a
        shorter tube reflect that field back to it, so the fields grew exponentially """
    FDTD, CSX = _setup(['PEC', 'PEC', 'PEC', 'PEC', 'MUR', 'MUR'], CoordSystem=1)
    nr, na, nz = CYL_CELLS[size]
    mesh = CSX.GetGrid()
    mesh.AddLine('r', np.arange(nr))
    mesh.AddLine('a', np.linspace(-np.pi, np.pi, na + 1))
    mesh.AddLine('z', np.arange(nz) - (nz - 1) / 2)
    R = nr - 1
    Z = (nz - 1) / 2

    CSX.AddMaterial('tube', epsilon=4).AddBox([0.4*R, -np.pi, -Z], [0.6*R, np.pi, Z])
    CSX.AddExcitation('line', exc_type=0, exc_val=[0, 0, 1]).AddBox([3, 0, -1], [3, 0, 1])
    CSX.AddProbe('E_probe', p_type=2).AddPoint([5, 0.2, 0])
    CSX.AddProbe('H_probe', p_type=3).AddPoint([2, -0.3, 1])
    return FDTD


CASES = {
    'ideal':         ideal,
    'pml':           pml,
    'inhomogeneous': inhomogeneous,
    'dispersive':    dispersive,
    'cylindrical':   cylindrical,
}


# --- field dumps --------------------------------------------------------------
#
# Written during the time-stepping, they cost the engine the field transfer,
# the interpolation and the file output. All are time-domain HDF5 files; a dump
# of the full volume is not among them, already at size M it writes several GB/s.

def _lines(CSX):
    """ the mesh lines of each direction """
    grid = CSX.GetGrid()
    return [grid.GetLines(n) for n in range(3)]


def _plane(FDTD, CSX):
    """ the E-field in the plane z=0 across the whole mesh, cell interpolated, at the
        default rate of a time-domain dump, Nyquist/4 """
    l = _lines(CSX)
    CSX.AddDump('Et_plane', dump_type=0, dump_mode=2, file_type=1).AddBox(
        [l[0][0], l[1][0], 0], [l[0][-1], l[1][-1], 0])


def _partial(FDTD, CSX):
    """ the E-field in the central box of half the extent in each direction (an eighth of
        the volume), not interpolated, at the default rate of a time-domain dump """
    start, stop = [], []
    for l in _lines(CSX):
        start.append(l[0] + 0.25 * (l[-1] - l[0]))
        stop.append(l[0] + 0.75 * (l[-1] - l[0]))
    CSX.AddDump('Et_partial', dump_type=0, dump_mode=0, file_type=1).AddBox(start, stop)


def _nf2ff(FDTD, CSX):
    """ a near-field to far-field box as CreateNF2FFBox() records it by default: the E- and
        H-fields on the faces of a box 10 mesh lines inside the boundary, node interpolated,
        in the time domain at the default rate. Also on the closed cases, where it has no
        far-field meaning; on the cylindrical mesh the outer cylinder and both ends. """
    l = _lines(CSX)
    start = [l[n][10] for n in range(3)]
    stop = [l[n][-11] for n in range(3)]
    directions = [True] * 6
    if CSX.GetGrid().GetMeshType() == 1:
        start[0], start[1], stop[1] = l[0][0], l[1][0], l[1][-1]
        directions = [False, True, False, False, True, True]
    FDTD.CreateNF2FFBox(start=start, stop=stop, directions=directions, mirror=[0] * 6)


DUMPS = {
    'none':    None,
    'plane':   _plane,
    'partial': _partial,
    'nf2ff':   _nf2ff,
}


def write_case(name, size, xml_file, dumps='none', nrts=None):
    """ Write case `name` at `size` with the field dumps `dumps` to `xml_file`; `nrts`
        sets a number of timesteps instead of a practically unlimited one """
    FDTD = CASES[name](size)
    if DUMPS[dumps]:
        DUMPS[dumps](FDTD, FDTD.GetCSX())
    if nrts:
        FDTD.SetNumberOfTimeSteps(nrts)
    FDTD.Write2XML(xml_file)
