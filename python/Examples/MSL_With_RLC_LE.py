# -*- coding: utf-8 -*-
"""
Microstrip line with a lumped series/parallel RLC load.

A 50 Ω lumped port excites the MSL at y=0.
An RLC lumped element terminates it at y=substrate_length.
S11 is de-embedded (naive matched-line assumption) and compared
with the analytic reflection coefficient of the load.
"""

import os, tempfile
import numpy as np
import matplotlib.pyplot as plt

from CSXCAD  import ContinuousStructure
from CSXCAD.CSProperties import LEtype
from openEMS import openEMS
from openEMS.physical_constants import C0

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
post_proc_only    = False   # set True to skip simulation and only re-plot
# show the structure before simulating; set True here, or run with SHOW_STRUCTURE=1
show_structure    = os.environ.get('SHOW_STRUCTURE', '0') == '1'

substrate_epsR      = 4.5
substrate_width     = 8       # mm
substrate_length    = 10      # mm
substrate_thickness = 1       # mm
substrate_cells     = 5
cu_thick            = 0.1     # mm
microstrip_W        = 1.875   # mm
Airbox_Add          = 7.5     # mm — airbox padding to boundary

R = 10      # Ω
L = 2e-9    # H  (np.nan → absent)
C = 1e-12   # F  (np.nan → absent)
Z0 = 50     # Ω — port reference impedance

networkType = LEtype.LE_SERIES   # or LEtype.LE_PARALLEL

f0 = 2e9   # Hz — Gaussian centre
fc = 1e9   # Hz — 20 dB corner

# ---------------------------------------------------------------------------
# FDTD setup
# ---------------------------------------------------------------------------
Sim_Path = os.path.join(tempfile.gettempdir(), 'Simp_MSL_W_RLC')

SimBox = np.array([
    -substrate_width*0.5 - Airbox_Add,
     substrate_width*0.5 + Airbox_Add,
    -Airbox_Add,
     substrate_length + Airbox_Add,
    -cu_thick - Airbox_Add,
     substrate_thickness + cu_thick + Airbox_Add])

FDTD = openEMS(NrTS=160000, EndCriteria=1e-4)
FDTD.SetGaussExcite(f0, fc)
FDTD.SetBoundaryCond(['MUR'] * 6)

CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(1e-3)
# ~λ/150 at the highest excited frequency — fine enough for the thin substrate
mesh_res = C0 / (f0 + fc) / 1e-3 / 150

# Seed the outer boundary lines; SmoothMeshLines fills in the interior
mesh.AddLine('x', SimBox[0:2])
mesh.AddLine('y', SimBox[2:4])
mesh.AddLine('z', SimBox[4:6])

# Top copper trace (microstrip signal line)
line = CSX.AddMaterial('cu_top', kappa=56000000)
line.AddBox([-microstrip_W/2, 0.0, substrate_thickness],
            [ microstrip_W/2, substrate_length, substrate_thickness + cu_thick],
            priority=10)
FDTD.AddEdges2Grid(dirs='xyz', properties=line)

# FR4 dielectric substrate
sub = CSX.AddMaterial('FR4', epsilon=substrate_epsR)
sub.AddBox([-substrate_width/2, 0.0, 0.0],
           [ substrate_width/2,  substrate_length, substrate_thickness],
           priority=2)

# Force uniform z-discretisation through the substrate
mesh.AddLine('z', np.linspace(0, substrate_thickness, substrate_cells + 1))

# Bottom copper ground plane
gnd = CSX.AddMaterial('cu_bot', kappa=56000000)
gnd.AddBox([-substrate_width/2, 0.0, -cu_thick],
           [ substrate_width/2, substrate_length, 0.0],
           priority=10)
FDTD.AddEdges2Grid(dirs='xyz', properties=gnd)

# Lumped 50 Ω excitation port at y=0
port = FDTD.AddLumpedPort(1, Z0,
                          [-microstrip_W*0.5, 0, 0],
                          [ microstrip_W*0.5, 0, substrate_thickness],
                          'z', 1.0, priority=15, edges2grid='xy')

# RLC load at y=substrate_length — the element under test
LE = CSX.AddLumpedElement('RLC_load', ny='z', caps=False,
                          R=R, L=L, C=C, LEtype=networkType)
LE.AddBox([-microstrip_W*0.5, substrate_length, 0],
          [ microstrip_W*0.5, substrate_length, substrate_thickness],
          priority=25)

mesh.SmoothMeshLines('all', mesh_res, 1.4)

### Show the structure (see show_structure above)
if show_structure:
    view = dict(view='iso', render_mode='solid+edges')
    try:
        from CSXCAD.Viewer import Show
    except ImportError:
        print('Showing the structure needs vtk: pip install vtk')
    else:
        Show(CSX, **view)

if not post_proc_only:
    FDTD.Run(Sim_Path, verbose=3, cleanup=True)

# ---------------------------------------------------------------------------
# Post-processing
# ---------------------------------------------------------------------------
f = np.linspace(max(1e9, f0 - fc), f0 + fc, 401)
port.CalcPort(Sim_Path, f)
s11 = port.uf_ref / port.uf_inc

# Analytic load impedance (NaN component = absent)
if networkType == LEtype.LE_SERIES:
    Zref = np.zeros(len(f), dtype=complex)
    if not np.isnan(R): Zref += R
    if not np.isnan(L): Zref += 1j*2*np.pi*f*L
    if not np.isnan(C): Zref += 1.0 / (1j*2*np.pi*f*C)
elif networkType == LEtype.LE_PARALLEL:
    Yref = np.zeros(len(f), dtype=complex)
    if not np.isnan(R): Yref += 1.0 / R
    if not np.isnan(L): Yref += 1.0 / (1j*2*np.pi*f*L)
    if not np.isnan(C): Yref += 1j*2*np.pi*f*C
    Zref = 1.0 / Yref
else:
    raise ValueError('networkType must be LEtype.LE_SERIES or LEtype.LE_PARALLEL')

Gref = (Zref - Z0) / (Zref + Z0)

# Naive matched-line de-embedding: remove the electrical length of the MSL.
# Dk_eff uses the Hammerstad-Jensen quasi-static effective permittivity formula.
Dk_eff = 0.5*(substrate_epsR + 1.0) * (1 + 1/np.sqrt(1 + 12*(substrate_thickness/microstrip_W)))
s11_deemb = s11 / np.exp(-2j * substrate_length*1e-3 * 2*np.pi*f / (C0/np.sqrt(Dk_eff)))

# Recover load impedance from de-embedded S11
Z_deemb = Z0 * (1 + s11_deemb) / (1 - s11_deemb)

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
label = 'Series RLC' if networkType == LEtype.LE_SERIES else 'Parallel RLC'

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
plt.tight_layout()

ax1.plot(f/1e9, 20*np.log10(np.abs(s11_deemb)), label='openEMS (de-embedded)')
ax1.plot(f/1e9, 20*np.log10(np.abs(Gref)), '--', label='Analytic')
ax1.set_title('{} - S11'.format(label))
ax1.set_xlabel('Frequency (GHz)')
ax1.set_ylabel('|S11| (dB)')
ax1.legend()
ax1.grid(True)

ax2.plot(f/1e9, np.abs(Z_deemb), label='openEMS (de-embedded)')
ax2.plot(f/1e9, np.abs(Zref), '--', label='Analytic')
ax2.set_title('{} - Load impedance'.format(label))
ax2.set_xlabel('Frequency (GHz)')
ax2.set_ylabel('|Z| (Ω)')
ax2.legend()
ax2.grid(True)

plt.show()
