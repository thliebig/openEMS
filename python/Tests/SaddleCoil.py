# -*- coding: utf-8 -*-
"""
NMR saddle coil — inductance and the effect of the tuning capacitor.

The geometry of python/Tutorials/SaddleCoil.py without its field dumps: a one
turn saddle coil of copper wire inside a conducting bore shield, fed by a
lumped port through a 1 kOhm resistor. The coil is small against the wavelength, so
its input impedance is that of an inductor with a small loss resistance.

Test cases
----------
1. Bare coil                 — the impedance is inductive over the whole band
2. With the tuning capacitor — the capacitor sits across the coil, so the pair
                               approaches a parallel resonance above the band
                               and the apparent inductance rises with it

Pass criteria
-------------
Bare coil        : reactance positive across the band, inductance 15 to 45 nH
                   and flat within 20 % (well below any resonance),
                   loss resistance under 5 Ohm
Tuning capacitor : apparent inductance above the bare one at every frequency,
                   rising with frequency (the parallel resonance is above the
                   band, as 1/(2*pi*sqrt(L*C)) says for this L and C), and no
                   resonance inside the band

The second case is the interesting one: it checks that the lumped element is
actually in the circuit and in parallel, not that a particular value came out.
"""

import os
import tempfile
import numpy as np

from CSXCAD import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import *

Sim_Path = os.path.join(tempfile.gettempdir(), 'SaddleCoil_Test')

unit = 1e-3
f0, fc = 500.13e6, 200e6

# the cells have to resolve the feed gap of the port (wire_space_c below): with
# coarser cells its two ends land on one mesh line, the port has no extent and
# nothing is excited at all
coil_height   = 20
coil_radius   = 2
coil_angle    = np.pi / 2
wire_radius   = 0.25
wire_space    = 0.25
wire_space_c  = 2 * wire_radius + wire_space
coil_rotation = coil_angle / 2
port_resist   = 1000
bore_radius   = 10
bore_shield   = 0.5
capaT_value   = 1.2e-12

mesh_step = wire_radius               # mm, the cell size around the coil
Airbox    = C0 / (f0 - fc) / unit / 25


def coil_points():
    """the wire path of the coil, one point per row, in mm"""
    ang_space_c = wire_space_c / coil_radius
    pts = []

    def at(angle, radius, z):
        a = angle - coil_rotation
        return [np.cos(a) * radius, np.sin(a) * radius, z]

    a = (np.pi + coil_angle - ang_space_c) / 2 - coil_rotation
    r = coil_radius + wire_space_c
    pts.append([np.cos(a) * r, np.sin(a) * r + wire_space_c, 0])
    for a in np.linspace((np.pi + coil_angle - ang_space_c) / 2, coil_angle - ang_space_c, 9):
        pts.append(at(a, coil_radius + wire_space_c, 0))
    for a in np.linspace(coil_angle - ang_space_c, 0, 9):
        pts.append(at(a, coil_radius, 0))
    for a in np.linspace(0, coil_angle, 9):
        pts.append(at(a, coil_radius, coil_height))
    for a in np.linspace(coil_angle, np.pi + coil_angle, 13):
        pts.append(at(a, coil_radius, 0))
    for a in np.linspace(np.pi + coil_angle, np.pi, 9):
        pts.append(at(a, coil_radius, coil_height))
    pts.append(at(np.pi, coil_radius, wire_space_c))
    pts.append(at(np.pi, coil_radius + wire_space_c, wire_space_c))
    for a in np.linspace(np.pi, (np.pi + coil_angle + ang_space_c) / 2, 9):
        pts.append(at(a, coil_radius + wire_space_c, 0))
    pts.append([pts[-1][0], pts[-1][1] + wire_space_c, pts[-1][2]])
    return np.round(np.array(pts), 3)


def merge(values, tol):
    """sort and drop values closer together than tol, which would make tiny cells"""
    out = []
    for v in sorted(values):
        if not out or v - out[-1] > tol:
            out.append(v)
    return out


def run(with_capacitor, sim_path):
    """simulate the coil, give the frequencies and the input impedance"""
    FDTD = openEMS(NrTS=200000, EndCriteria=1e-4)
    FDTD.SetGaussExcite(f0, fc)
    FDTD.SetBoundaryCond(['MUR'] * 4 + ['PML_8'] * 2)

    CSX = ContinuousStructure()
    FDTD.SetCSX(CSX)
    mesh = CSX.GetGrid()
    mesh.SetDeltaUnit(unit)

    points = coil_points()
    coil = CSX.AddMaterial('saddle', kappa=56e6)
    coil.AddWire(points.T, wire_radius, priority=10)

    start, stop = points[0], points[-1]
    if with_capacitor:
        capaT = CSX.AddLumpedElement('CapaT', ny='x', caps=True, C=capaT_value)
        capaT.AddBox([points[1][0], points[1][1] - 0.25, points[1][2] - 0.25],
                     [points[-2][0], points[-2][1] + 0.25, points[-2][2] + 0.25])

    port = FDTD.AddLumpedPort(1, port_resist, start, stop, 'x', 1.0, priority=100)

    shield = CSX.AddMetal('BoreShield')
    shield.AddCylindricalShell([0, 0, -wire_radius - Airbox],
                               [0, 0, coil_height + wire_radius + Airbox],
                               bore_radius, bore_shield, priority=10)

    edges = {'x': [-bore_radius - bore_shield, bore_radius + bore_shield],
             'y': [-bore_radius - bore_shield, bore_radius + bore_shield],
             'z': [-wire_radius, coil_height + wire_radius]}
    for n, ny in enumerate('xyz'):
        lines = merge(np.concatenate([points[:, n], edges[ny]]), mesh_step / 2)
        refined = list(lines)
        for a, b in zip(lines[:-1], lines[1:]):
            N = int(round((b - a) / mesh_step))
            if N >= 2:
                refined.extend(np.linspace(a, b, N))
        refined = merge(refined, mesh_step / 2)
        mesh.AddLine(ny, refined)
        mesh.AddLine(ny, [refined[0] - Airbox, refined[-1] + Airbox])
        mesh.SmoothMeshLines(ny, C0 / (f0 + fc) / unit / 20, 1.4)

    z = mesh.GetLines('z')
    mesh.AddLine('z', [z[0] - (z[1] - z[0]) * n for n in range(1, 9)])
    mesh.AddLine('z', [z[-1] + (z[-1] - z[-2]) * n for n in range(1, 9)])

    FDTD.Run(sim_path, cleanup=True, verbose=0,
             engine=os.environ.get('OPENEMS_TEST_ENGINE', 'multithreaded'))

    f = np.linspace(f0 - fc, f0 + fc, 201)
    port.CalcPort(sim_path, f)
    return f, port.uf_tot / port.if_tot


def inductance(f, Zin):
    return np.imag(Zin) / (2 * np.pi * f)


print('--- bare coil')
f, Z_bare = run(False, Sim_Path + '_bare')
L_bare = inductance(f, Z_bare)
R_bare = np.real(Z_bare)
print('    inductance %.1f nH at %.1f MHz, %.1f to %.1f nH over the band'
      % (np.interp(f0, f, L_bare) * 1e9, f0 * 1e-6, L_bare.min() * 1e9, L_bare.max() * 1e9))
print('    loss resistance %.2f Ohm at %.1f MHz' % (np.interp(f0, f, R_bare), f0 * 1e-6))

assert np.all(np.imag(Z_bare) > 0), 'the bare coil has to be inductive over the band'
L0 = np.interp(f0, f, L_bare)
assert 15e-9 < L0 < 45e-9, 'inductance %.1f nH is not that of this coil' % (L0 * 1e9)
assert L_bare.max() / L_bare.min() < 1.2, \
    'the inductance of the bare coil varies by %.0f %% over the band, so it is not ' \
    'well below its resonance' % (100 * (L_bare.max() / L_bare.min() - 1))
assert np.interp(f0, f, R_bare) < 5, 'the loss resistance of a copper coil is small'

print('--- with the tuning capacitor (%.1f pF)' % (capaT_value * 1e12))
f, Z_capa = run(True, Sim_Path + '_capa')
L_capa = inductance(f, Z_capa)
print('    apparent inductance %.1f nH at %.1f MHz, %.1f to %.1f nH over the band'
      % (np.interp(f0, f, L_capa) * 1e9, f0 * 1e-6, L_capa.min() * 1e9, L_capa.max() * 1e9))
print('    parallel resonance of %.1f nH and %.1f pF: %.0f MHz'
      % (L0 * 1e9, capaT_value * 1e12, 1 / (2 * np.pi * np.sqrt(L0 * capaT_value)) / 1e6))

assert np.all(np.imag(Z_capa) > 0), \
    'the resonance of this L and C is above the band, so the pair stays inductive in it'
assert np.all(L_capa > L_bare), \
    'the capacitor is across the coil, so it has to raise the apparent inductance'
assert L_capa[-1] / L_capa[0] > L_bare[-1] / L_bare[0], \
    'the apparent inductance has to rise towards the resonance'

print('PASS')
