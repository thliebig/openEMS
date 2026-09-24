# -*- coding: utf-8 -*-
"""
Re-entrant coaxial cavity resonator — a long, high-Q run.

A quarter-wave coaxial line, shorted at one end and capacitively loaded at the
other, inside a closed PEC can. The resonance is narrow, so the run is long: the
field has to ring down before the spectrum resolves the peak, and a small error
in the update coefficients shows up as a frequency shift after a few hundred
cycles that a short run would hide. That is what this test is for; the other
tests in this directory stop after 10k to 40k timesteps, this one needs a few
hundred thousand.

Geometry
--------
A can of inner radius b, shorted at z=0, with a post of radius a from z=0 to
z=L. The post does not reach the lid: the gap between the post and the lid is
where the fields concentrate, and where a lumped capacitor can be placed across
it. Looking down into the gap, the shorted line of length L has the input
reactance of a transmission line, so the structure resonates where

    Z0 * tan(2*pi*f*L/c) = 1 / (2*pi*f*C)      Z0 = 60*ln(b/a) for an air line

with C the total capacitance across the gap.

Test cases
----------
1. Bare gap          — calibration: the measured resonance gives the gap
                       capacitance C_gap, which fringes and cannot be had from
                       the plate formula alone
2. Gap + 1 pF        — C_gap + 1 pF predicts the resonance
3. Gap + 2 pF        — C_gap + 2 pF predicts it again, further down

Pass criteria
-------------
- the fitted C_gap is within a factor of three of the parallel-plate estimate,
  which catches a geometry or mesh that is not the intended one
- cases 2 and 3 land within 4 % of the frequency the fitted C_gap predicts, so
  the line physics and the lumped element agree over a 2:1 change of load
- the resonance moves down monotonically as the load grows
- the peak is narrower than f/300. The walls are PEC, so nothing but the run
  length limits the width: this is the check that the run was long enough for
  the spectrum to resolve the resonance at all
"""

import os
import tempfile
import numpy as np

from CSXCAD import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import *

unit = 1e-3                  # mm

post_r    = 3.0              # post radius a
can_r     = 7.0              # can inner radius b
can_t     = 1.0              # wall thickness
post_len  = 60.0             # shorted line length L
gap       = 1.0              # post to lid
port_R    = 10e3             # weak coupling: the port must not damp the cavity

Z0 = 60 * np.log(can_r / post_r)          # air line, Ohm
C_gap_plate = EPS0 * np.pi * (post_r * unit)**2 / (gap * unit)

f0, fc = 1.0e9, 0.6e9        # excitation: covers 0.4 to 1.6 GHz
freq   = np.linspace(0.4e9, 1.6e9, 2001)


def resonance(C):
    """ The f solving Z0*tan(2*pi*f*L/c) = 1/(2*pi*f*C), by bisection. """
    L = post_len * unit
    def residual(f):
        return Z0 * np.tan(2 * np.pi * f * L / C0) - 1.0 / (2 * np.pi * f * C)
    lo, hi = 0.05e9, 0.999 * C0 / (4 * L)   # below the quarter-wave pole
    assert residual(lo) < 0 < residual(hi)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if residual(mid) < 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def fit_gap_capacitance(f_res):
    """ The C that puts the resonance at f_res, i.e. resonance() inverted. """
    L = post_len * unit
    return 1.0 / (2 * np.pi * f_res * Z0 * np.tan(2 * np.pi * f_res * L / C0))


def run(Sim_Path, C_load=None):
    """ Build and run the cavity, with an optional lumped capacitor across the
        gap, and return the port impedance over freq. """
    FDTD = openEMS(NrTS=400000, EndCriteria=1e-5)
    FDTD.SetGaussExcite(f0, fc)
    FDTD.SetBoundaryCond(['PEC'] * 6)      # the can is closed, nothing leaves

    CSX = ContinuousStructure()
    FDTD.SetCSX(CSX)
    mesh = CSX.GetGrid()
    mesh.SetDeltaUnit(unit)

    lid = post_len + gap
    metal = CSX.AddMetal('metal')
    # the can: wall, floor and lid
    metal.AddCylindricalShell([0, 0, 0], [0, 0, lid], can_r + can_t/2, can_t)
    metal.AddCylinder([0, 0, -can_t], [0, 0, 0],   can_r + can_t)
    metal.AddCylinder([0, 0, lid],    [0, 0, lid + can_t], can_r + can_t)
    # the post, shorted to the floor
    metal.AddCylinder([0, 0, 0], [0, 0, post_len], post_r)

    # the gap: the port, and the capacitor across it
    start = [0, 0, post_len]
    stop  = [0, 0, lid]
    if C_load is not None:
        capa = CSX.AddLumpedElement('C_load', ny='z', caps=True, C=C_load)
        capa.AddBox([-post_r, -post_r, post_len], [post_r, post_r, lid], priority=5)
    port = FDTD.AddLumpedPort(1, port_R, start, stop, 'z', 1.0, priority=10)

    # transverse: the radii on mesh lines, then smoothed outwards
    for ny in ('x', 'y'):
        mesh.AddLine(ny, [-can_r - can_t, -can_r, -post_r, 0, post_r, can_r, can_r + can_t])
        mesh.SmoothMeshLines(ny, 0.5, 1.4)
    # axial: the gap resolved, the line coarser
    mesh.AddLine('z', [-can_t, 0, post_len, lid, lid + can_t])
    mesh.AddLine('z', np.linspace(post_len, lid, 3))
    mesh.SmoothMeshLines('z', 1.0, 1.4)

    FDTD.Run(Sim_Path, cleanup=True, exact_endcriteria=True)

    port.CalcPort(Sim_Path, freq)
    return port.uf_tot / port.if_tot


def peak(Z):
    """ The resonance: the frequency of the |Z| maximum, and f/width at -3 dB.
        With PEC walls the width is set by the length of the run, not by a loss. """
    mag = np.abs(Z)
    i   = np.argmax(mag)
    f_r = freq[i]
    half = mag[i] / np.sqrt(2)
    lo = np.where(mag[:i] < half)[0]
    hi = np.where(mag[i:] < half)[0]
    if len(lo) == 0 or len(hi) == 0:
        return f_r, np.inf
    f_lo = np.interp(half, [mag[lo[-1]], mag[lo[-1] + 1]], [freq[lo[-1]], freq[lo[-1] + 1]])
    f_hi = np.interp(-half, [-mag[i + hi[0] - 1], -mag[i + hi[0]]], [freq[i + hi[0] - 1], freq[i + hi[0]]])
    return f_r, f_r / (f_hi - f_lo)


print('Coaxial resonator: Z0 = %.1f Ohm, quarter-wave at %.3f GHz'
      % (Z0, C0 / (4 * post_len * unit) / 1e9))

print('Running the bare cavity (calibration) ...')
f_bare, Q_bare = peak(run(os.path.join(tempfile.gettempdir(), 'CoaxRes_bare')))
C_gap = fit_gap_capacitance(f_bare)
print('  resonance %.4f GHz, f/width %.0f, fitted gap capacitance %.3f pF'
      % (f_bare / 1e9, Q_bare, C_gap * 1e12))

# the plate value is a lower bound: the fringing field at the rim only adds to it
assert C_gap_plate < C_gap < 5 * C_gap_plate, \
    'fitted gap capacitance %.3f pF against the %.3f pF of the plate formula' \
    % (C_gap * 1e12, C_gap_plate * 1e12)
assert Q_bare > 300, 'the peak is %.0f wide: the cavity is not resonating' % Q_bare

results = []
for C_load in (3e-12, 6e-12):
    print('Running with a %.0f pF capacitor across the gap ...' % (C_load * 1e12))
    f_meas, Q = peak(run(os.path.join(tempfile.gettempdir(),
                                      'CoaxRes_%.0fpF' % (C_load * 1e12)), C_load))
    f_pred = resonance(C_gap + C_load)
    err = abs(f_meas - f_pred) / f_pred
    print('  resonance %.4f GHz, predicted %.4f GHz (%.1f %%), f/width %.0f'
          % (f_meas / 1e9, f_pred / 1e9, err * 100, Q))
    assert err < 0.04, \
        '%.0f pF: resonance %.4f GHz against the predicted %.4f GHz (%.1f %% > 4 %%)' \
        % (C_load * 1e12, f_meas / 1e9, f_pred / 1e9, err * 100)
    results.append(f_meas)

assert f_bare > results[0] > results[1], \
    'the resonance has to move down as the gap is loaded: %.4f, %.4f, %.4f GHz' \
    % (f_bare / 1e9, results[0] / 1e9, results[1] / 1e9)

print('PASS')
