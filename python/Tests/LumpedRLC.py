# -*- coding: utf-8 -*-
"""
Lumped RLC element validation — CurvePort loop test.

Geometry: a small rectangular loop in the xz-plane.
  - CurvePort  at x=0,         z ∈ [0, h]  — measures input impedance
  - LumpedRLC  at x=w,         z ∈ [0, h]  — element under test
  - PEC curves at z=0 and z=h, x ∈ [0, w]  — close the circuit

A reference run with a PEC short at the element position yields Z_loop(f).
Subtracting it isolates the element impedance:
    Z_element(f) = Z_in(f) - Z_loop(f)
This is exact (FDTD is linear), removing the parasitic loop inductance.

Test cases
----------
1. Parallel RLC  (R, L, C all set) — resonance at f0, peak  |Z| ≈ R_par
2. Parallel RL   (C = NaN)         — absent capacitor, RL parallel impedance
3. Series   RLC  (R, L, C all set) — resonance at f0, min   |Z| ≈ R_ser
4. Series   RL   (C = NaN)         — absent capacitor, RL series impedance
5. Parallel RC   (L = NaN)         — absent inductor, handled by base FDTD (not the extension)
6. Series   RC   (L = NaN)         — absent inductor, RC series impedance via extension with dL=0

Pass criteria
-------------
Resonant cases : f_res within 5 % of f0, |Z| at resonance within 15 % of R
Non-resonant   : |Z_measured − Z_analytic| / |Z_analytic| < 10 % at spot freqs

(c) 2026 Thorsten Liebig <thorsten.liebig@gmx.de>
"""

import os, tempfile
import numpy as np

from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.ports import CurvePort

# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------
unit   = 1e-3   # drawing unit: mm
loop_w = 2.0    # mm — port-to-element separation (x)
loop_h = 1.0    # mm — element / port height (z)
small  = 0.05   # mm — half-thickness of element snap box
pad    = 15.0   # mm — airbox clearance to boundary

# ---------------------------------------------------------------------------
# Frequency sweep
# ---------------------------------------------------------------------------
f_start = 0.5e9
f_stop  = 2.0e9
freq    = np.linspace(f_start, f_stop, 401)

# ---------------------------------------------------------------------------
# Component values — target resonance at f0 = 1 GHz
# ---------------------------------------------------------------------------
f_0   = 1.0e9
L_val = 10e-9                                      # 10 nH
C_val = 1.0 / ((2*np.pi*f_0)**2 * L_val)           # ≈ 2.53 pF → f0 = 1 GHz
R_par = 1000.0   # Ω  — parallel resonator peak impedance
R_ser = 10.0     # Ω  — series   resonator minimum impedance

# ---------------------------------------------------------------------------
# Analytic impedance formulas (NaN component = absent)
# ---------------------------------------------------------------------------
def z_par(f, R=np.nan, L=np.nan, C=np.nan):
    """Parallel RLC admittance; absent components contribute zero admittance."""
    f = np.asarray(f, dtype=float)
    Y = np.zeros_like(f, dtype=complex)
    if not np.isnan(R): Y += 1.0 / R
    if not np.isnan(L): Y += 1.0 / (1j * 2*np.pi*f * L)
    if not np.isnan(C): Y += 1j * 2*np.pi*f * C
    return 1.0 / Y

def z_ser(f, R=np.nan, L=np.nan, C=np.nan):
    """Series RLC impedance; absent components contribute zero impedance (wire)."""
    f = np.asarray(f, dtype=float)
    Z = np.zeros_like(f, dtype=complex)
    if not np.isnan(R): Z += R
    if not np.isnan(L): Z += 1j * 2*np.pi*f * L
    if not np.isnan(C): Z += 1.0 / (1j * 2*np.pi*f * C)
    return Z

# ---------------------------------------------------------------------------
# Build + run one loop simulation
# ---------------------------------------------------------------------------
def run_loop(Sim_Path, element_fn):
    """
    Build the rectangular loop, apply element_fn(CSX) at the element position,
    run FDTD and return Z_in(freq) as a complex numpy array.

    element_fn must add its geometry inside the box
        [loop_w ± small, ± small, 0 .. loop_h]
    """
    FDTD = openEMS(EndCriteria=1e-4)
    FDTD.SetGaussExcite(0.5*(f_start + f_stop), 0.5*(f_stop - f_start))
    FDTD.SetBoundaryCond(['MUR'] * 6)

    CSX = ContinuousStructure()
    FDTD.SetCSX(CSX)
    mesh = CSX.GetGrid()
    mesh.SetDeltaUnit(unit)

    mesh.AddLine('x', [-pad, 0, loop_w, loop_w + pad])
    mesh.AddLine('y', [-pad, 0, pad])
    mesh.AddLine('z', [-pad, 0, loop_h, loop_h + pad])
    mesh.SmoothMeshLines('all', 0.5, ratio=1.4)

    # CurvePort at x=0, spanning z=[0, loop_h]
    port = CurvePort(CSX, 1, R=50, start=[0, 0, 0], stop=[0, 0, loop_h], excite=1)

    # PEC curves closing the loop at z=0 and z=loop_h
    conn = CSX.AddMetal('connections')
    conn.AddCurve([[0, loop_w], [0, 0], [0,      0     ]])   # bottom
    conn.AddCurve([[0, loop_w], [0, 0], [loop_h, loop_h]])   # top

    element_fn(CSX)

    FDTD.Run(Sim_Path, cleanup=True, exact_endcriteria=True)

    port.CalcPort(Sim_Path, freq)
    s11  = port.uf_ref / port.uf_inc
    return port.Z_ref * (1 + s11) / (1 - s11)

# ---------------------------------------------------------------------------
# Element box coordinates (same for element and calibration short)
# ---------------------------------------------------------------------------
elem_start = [loop_w - small, -small, 0       ]
elem_stop  = [loop_w + small,  small, loop_h  ]

# ---------------------------------------------------------------------------
# Calibration: PEC short at element position → Z_loop(f)
# ---------------------------------------------------------------------------
def add_short(CSX):
    CSX.AddMetal('short').AddBox(elem_start, elem_stop, priority=20)

print('Running calibration (PEC short)...')
Z_loop = run_loop(os.path.join(tempfile.gettempdir(), 'LRLC_calib'), add_short)

# ---------------------------------------------------------------------------
# Check helpers
# ---------------------------------------------------------------------------
def check_spot(label, Z_in, Z_ref_fn, f_spot, tol=0.10):
    """Compare calibrated Z at spot frequencies against analytic formula."""
    Z_corr = Z_in - Z_loop
    Z_ref  = Z_ref_fn(f_spot)
    errs   = np.abs(Z_corr[[np.argmin(np.abs(freq - f)) for f in f_spot]] - Z_ref) \
             / np.maximum(np.abs(Z_ref), 1.0)
    for f, err in zip(f_spot, errs):
        assert err < tol, (
            'FAIL [{}] @ {:.2f} GHz: err={:.1%} (limit {:.0%})'.format(label, f/1e9, err, tol)
        )
    print('  PASS [{}]: spot-frequency check OK (max err {:.1%})'.format(label, np.max(errs)))

def check_resonance(label, Z_in, expect_peak, R_ref, f_lo=0.7e9, f_hi=1.4e9):
    """Check resonant frequency and impedance value at resonance."""
    Z_corr = np.abs(Z_in - Z_loop)
    mask   = (freq >= f_lo) & (freq <= f_hi)
    i_res  = np.argmax(Z_corr[mask]) if expect_peak else np.argmin(Z_corr[mask])
    f_res  = freq[mask][i_res]
    Z_res  = Z_corr[mask][i_res]

    f_err = abs(f_res - f_0) / f_0
    Z_err = abs(Z_res - R_ref) / R_ref
    assert f_err < 0.05, (
        'FAIL [{}]: f_res={:.3f} GHz, expected {:.3f} GHz ({:.1%} > 5 %)'.format(
            label, f_res/1e9, f_0/1e9, f_err)
    )
    assert Z_err < 0.15, (
        'FAIL [{}]: |Z|_res={:.1f} Ohm, expected {:.0f} Ohm ({:.1%} > 15 %)'.format(
            label, Z_res, R_ref, Z_err)
    )
    print('  PASS [{}]: f_res={:.3f} GHz, |Z|_res={:.1f} Ohm'.format(label, f_res/1e9, Z_res))

# ---------------------------------------------------------------------------
# Test 1: Parallel RLC  (R, L, C all set)
# ---------------------------------------------------------------------------
def add_par_rlc(CSX):
    e = CSX.AddLumpedElement('par_rlc', ny='z', caps=False,
                             R=R_par, L=L_val, C=C_val, LEtype=0)
    e.AddBox(elem_start, elem_stop, priority=20)

print('Test 1: Parallel RLC ...')
Z_par_rlc = run_loop(os.path.join(tempfile.gettempdir(), 'LRLC_par_rlc'), add_par_rlc)
check_resonance('Parallel RLC', Z_par_rlc, expect_peak=True, R_ref=R_par)

# ---------------------------------------------------------------------------
# Test 2: Parallel RL  (C = NaN — absent)
# ---------------------------------------------------------------------------
def add_par_rl(CSX):
    e = CSX.AddLumpedElement('par_rl', ny='z', caps=False,
                             R=R_par, L=L_val, LEtype=0)   # C not set → NaN
    e.AddBox(elem_start, elem_stop, priority=20)

f_spot = np.array([0.7e9, 1.0e9, 1.5e9])

print('Test 2: Parallel RL (NaN C) ...')
Z_par_rl = run_loop(os.path.join(tempfile.gettempdir(), 'LRLC_par_rl'), add_par_rl)
check_spot('Parallel RL', Z_par_rl,
           lambda f: z_par(f, R=R_par, L=L_val),
           f_spot)

# ---------------------------------------------------------------------------
# Test 3: Series RLC  (R, L, C all set)
# ---------------------------------------------------------------------------
def add_ser_rlc(CSX):
    e = CSX.AddLumpedElement('ser_rlc', ny='z', caps=False,
                             R=R_ser, L=L_val, C=C_val, LEtype=1)
    e.AddBox(elem_start, elem_stop, priority=20)

print('Test 3: Series RLC ...')
Z_ser_rlc = run_loop(os.path.join(tempfile.gettempdir(), 'LRLC_ser_rlc'), add_ser_rlc)
check_resonance('Series RLC', Z_ser_rlc, expect_peak=False, R_ref=R_ser)

# ---------------------------------------------------------------------------
# Test 4: Series RL  (C = NaN — absent)
# ---------------------------------------------------------------------------
def add_ser_rl(CSX):
    e = CSX.AddLumpedElement('ser_rl', ny='z', caps=False,
                             R=R_ser, L=L_val, LEtype=1)    # C not set → NaN
    e.AddBox(elem_start, elem_stop, priority=20)

print('Test 4: Series RL (NaN C) ...')
Z_ser_rl = run_loop(os.path.join(tempfile.gettempdir(), 'LRLC_ser_rl'), add_ser_rl)
check_spot('Series RL', Z_ser_rl,
           lambda f: z_ser(f, R=R_ser, L=L_val),
           f_spot)

# ---------------------------------------------------------------------------
# Test 5: Parallel RC  (L = NaN — absent)
#         IsLElumpedRLC returns False for parallel without L, so this element
#         is handled entirely by the base FDTD operator (standard EC_G / EC_C).
# ---------------------------------------------------------------------------
R_rc = 100.0                                # Ω
C_rc = 1.0 / (2*np.pi*f_0 * R_rc)         # → -3 dB corner at f_0 = 1 GHz

def add_par_rc(CSX):
    e = CSX.AddLumpedElement('par_rc', ny='z', caps=False,
                             R=R_rc, C=C_rc, LEtype=0)   # L not set → NaN
    e.AddBox(elem_start, elem_stop, priority=20)

print('Test 5: Parallel RC (NaN L, base FDTD) ...')
Z_par_rc = run_loop(os.path.join(tempfile.gettempdir(), 'LRLC_par_rc'), add_par_rc)
check_spot('Parallel RC', Z_par_rc,
           lambda f: z_par(f, R=R_rc, C=C_rc),
           f_spot)

# ---------------------------------------------------------------------------
# Test 6: Series RC  (L = NaN — absent)
#         Goes through the extension with dL = 0 — exercises NaN-L handling.
# ---------------------------------------------------------------------------
def add_ser_rc(CSX):
    e = CSX.AddLumpedElement('ser_rc', ny='z', caps=False,
                             R=R_ser, C=C_val, LEtype=1)  # L not set → NaN
    e.AddBox(elem_start, elem_stop, priority=20)

print('Test 6: Series RC (NaN L) ...')
Z_ser_rc = run_loop(os.path.join(tempfile.gettempdir(), 'LRLC_ser_rc'), add_ser_rc)
check_spot('Series RC', Z_ser_rc,
           lambda f: z_ser(f, R=R_ser, C=C_val),
           f_spot)

print('\nAll lumped RLC tests PASSED')

if 1:  # set to 1 for debugging plots
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 3, figsize=(16, 8))
    plt.tight_layout()

    for ax, Z, label, Z_analytic in [
        (axes[0, 0], Z_par_rlc, 'Parallel RLC',      z_par(freq, R=R_par, L=L_val, C=C_val)),
        (axes[0, 1], Z_par_rl,  'Parallel RL (NaN C)', z_par(freq, R=R_par, L=L_val)),
        (axes[0, 2], Z_par_rc,  'Parallel RC (NaN L)', z_par(freq, R=R_rc,  C=C_rc)),
        (axes[1, 0], Z_ser_rlc, 'Series RLC',        z_ser(freq, R=R_ser, L=L_val, C=C_val)),
        (axes[1, 1], Z_ser_rl,  'Series RL (NaN C)', z_ser(freq, R=R_ser, L=L_val)),
        (axes[1, 2], Z_ser_rc,  'Series RC (NaN L)', z_ser(freq, R=R_ser, C=C_val)),
    ]:
        Z_corr = Z - Z_loop
        ax.plot(freq/1e9, np.abs(Z_corr),    label='FDTD')
        ax.plot(freq/1e9, np.abs(Z_analytic), label='Analytic', linestyle='--')
        ax.set_title(label)
        ax.set_xlabel('Frequency (GHz)')
        ax.set_ylabel('|Z| (Ω)')
        ax.grid(True)
        ax.legend()

    plt.show()
