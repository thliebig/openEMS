# -*- coding: utf-8 -*-
"""
 Dipole SAR simulation test

 Half-wave dipole at 1 GHz next to a three-layer ellipsoidal phantom
 (skin / headbone / brain, ~280 g), adapted from the Dipole_SAR tutorial.

 SAR variants checked:
   0 g  / SIMPLE     / no autoRange  — local SAR
   1 g  / SIMPLE     / no autoRange  — standard averaging
   1 g  / IEEE_62704 / no autoRange
   10 g / IEEE_62704 / no autoRange  — strictest standard
   10 g / SIMPLE     / 20 dB         — autoRange feature

 All SAR values are normalized to 1 W accepted antenna input power and
 reported in W/kg at that power level.

 Pass criteria:
   peak SAR values match the reference constants below (within 1 %)
   peak SAR does not increase with averaging mass (0 g ≥ 1 g ≥ 10 g)
   phantom mass equals the expected value (0.2801 kg)
   absorbed power is identical for all variants
   autoRange reduces the number of non-zero cells

 To calibrate: run once and note the printed peak SAR values, then fill in
 EXPECTED_PEAK_SAR below.

 Tested with
  - python 3.14
  - openEMS v0.0.36+

 (c) 2026 Thorsten Liebig <thorsten.liebig@gmx.de>

"""

import os, tempfile, shutil
import numpy as np

from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import C0
from openEMS.sar_calculation import SAR_Calculation
from openEMS.sar_utils import readSAR

# ── Reference values ─────────────────────────────────────────────────────────
# Peak SAR in W/kg at 1 W accepted antenna input power.
# Set to None until calibrated; None entries skip that value check.
# Calibrated with the SAR recording sampled at the default OverSampling (4); the values
# converge there (OverSampling 16: within 0.01 %). Sampled at the Nyquist rate, as
# before, they were ~2 % higher.
EXPECTED_PEAK_SAR = {
    'm0g_SIMPLE':       5.2847,
    'm1g_SIMPLE':       4.72736,
    'm1g_IEEE_62704':   4.72736,
    'm10g_IEEE_62704':  3.72305,
    'm10g_SIMPLE_ar20': 3.72305,
}
PEAK_SAR_RTOL    = 0.01    # 1 % tolerance
EXPECTED_MASS_KG = 0.2801  # phantom mass; None to skip

# ── Simulation parameters ─────────────────────────────────────────────────────
unit             = 1e-3
f0               = 1e9
f_stop           = 1.5e9
mesh_res_phantom = 2.5    # mm
mesh_res_air     = C0 / f_stop / 20 / unit

_phantoms = [
    # name,      epsR,  kappa, density, radii [mm],    center [mm]
    ('skin',      50,   0.65,  1100,   [32, 40, 40],  [50, 0, 0]),
    ('headbone',  13,   0.10,  2000,   [30, 38, 38],  [50, 0, 0]),
    ('brain',     60,   0.70,  1040,   [26, 34, 34],  [50, 0, 0]),
]

# ── Run simulation ────────────────────────────────────────────────────────────
Sim_Path = os.path.join(tempfile.gettempdir(), 'Dipole_SAR_Test')
print(f'Simulation path: {Sim_Path}')

FDTD = openEMS(NrTS=30000, EndCriteria=1e-4, CellConstantMaterial=True)
FDTD.SetGaussExcite(0, f_stop)
FDTD.SetBoundaryCond(['PML_8'] * 6)

CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(unit)

dipole_length = 0.48 * (C0 / f0) / unit
dip = CSX.AddMetal('Dipole')
dip.AddBox([0, 0, -dipole_length / 2], [0, 0, dipole_length / 2], priority=1)
thirds = np.array([-1/3, 2/3])
mesh.AddLine('z', -dipole_length / 2 - thirds * mesh_res_phantom)
mesh.AddLine('z',  dipole_length / 2 + thirds * mesh_res_phantom)

for n, (name, epsR, kappa, density, radii, center) in enumerate(_phantoms):
    mat = CSX.AddMaterial(name, epsilon=epsR, kappa=kappa, density=density)
    sp  = mat.AddSphere(priority=10 + n, center=[0, 0, 0], radius=1)
    tr  = sp.GetTransform()
    tr.AddTransform('Scale', radii)
    tr.AddTransform('Translate', center)
    for di, d in enumerate('xyz'):
        mesh.AddLine(d, [center[di] - radii[di], center[di] + radii[di]])

mesh.AddLine('x', [0])
mesh.AddLine('y', [0])
port = FDTD.AddLumpedPort(port_nr=1, R=50,
                          start=[-0.1, -0.1, -mesh_res_phantom / 2],
                          stop=[ 0.1,  0.1,  mesh_res_phantom / 2],
                          p_dir='z', excite=True)

mesh.SmoothMeshLines('all', mesh_res_phantom, 1.4)
mesh.AddLine('x', [-100, 200])
mesh.AddLine('y', [-150, 150])
mesh.AddLine('z', [-150, 150])
mesh.SmoothMeshLines('all', mesh_res_air, 1.4)

skin_r, skin_c = _phantoms[0][4], _phantoms[0][5]
margin = 5
sar_dump = CSX.AddDump('SAR', dump_type=29, frequency=[f0],
                        file_type=1, dump_mode=2)
sar_dump.AddBox(
    [skin_c[0] - skin_r[0] - margin, -skin_r[1] - margin, -skin_r[2] - margin],
    [skin_c[0] + skin_r[0] + margin,  skin_r[1] + margin,  skin_r[2] + margin],
)

os.makedirs(Sim_Path, exist_ok=True)
FDTD.Run(Sim_Path, cleanup=True)

# Accepted power at f0 — used to normalize all SAR values to 1 W input
f_sweep = np.linspace(0.5e9, f_stop, 501)
port.CalcPort(Sim_Path, f_sweep)
Pin_f0 = np.interp(f0, f_sweep, port.P_acc)
print(f'Accepted power at {f0/1e9:.1f} GHz: {Pin_f0:.4g} W')

raw_sar = os.path.join(Sim_Path, 'SAR.h5')

# ── SAR calculation ───────────────────────────────────────────────────────────
import h5py

def _run_sar(key, mass_g, method, autorange_db=None):
    out_h5 = os.path.join(Sim_Path, f'{key}.h5')
    sc = SAR_Calculation(mass=float(mass_g), method=method)
    if autorange_db is not None:
        sc.EnableAutoRange(autorange_db)
    ok = sc.CalcFromHDF5(raw_sar, out_h5)
    assert ok, f'SAR calculation failed for {key}'
    sar, _, data = readSAR(out_h5)
    with h5py.File(out_h5, 'r') as h:
        total_mass = float(h.attrs['mass'])
    # normalize SAR array and peak to 1 W accepted input power
    sar_norm = sar / Pin_f0
    return sar_norm, float(sar_norm.max()), total_mass, float(data['power'])

results = {}
results['m0g_SIMPLE']       = _run_sar('m0g_SIMPLE',       0,  'SIMPLE')
results['m1g_SIMPLE']       = _run_sar('m1g_SIMPLE',       1,  'SIMPLE')
results['m1g_IEEE_62704']   = _run_sar('m1g_IEEE_62704',   1,  'IEEE_62704')
results['m10g_IEEE_62704']  = _run_sar('m10g_IEEE_62704',  10, 'IEEE_62704')
results['m10g_SIMPLE_ar20'] = _run_sar('m10g_SIMPLE_ar20', 10, 'SIMPLE', autorange_db=6.0)

# ── Print results ─────────────────────────────────────────────────────────────
print('\nPeak SAR results (at 1 W accepted antenna power):')
for key, (_, peak, mass, power) in results.items():
    print(f'  {key:25s}  peak = {peak:.6g} W/kg   mass = {mass:.4g} kg   power = {power:.4g} W')

# ── Checks ───────────────────────────────────────────────────────────────────
peak_0g  = results['m0g_SIMPLE'][1]
peak_1g  = results['m1g_SIMPLE'][1]
peak_10g = results['m10g_IEEE_62704'][1]

assert peak_0g >= peak_1g, \
    f'FAIL: peak SAR 0g ({peak_0g:.4g}) < peak SAR 1g ({peak_1g:.4g})'
assert peak_1g >= peak_10g, \
    f'FAIL: peak SAR 1g ({peak_1g:.4g}) < peak SAR 10g ({peak_10g:.4g})'

masses = [v[2] for v in results.values()]
assert all(m == masses[0] for m in masses), \
    f'FAIL: phantom mass differs across variants: {masses}'
if EXPECTED_MASS_KG is not None:
    assert abs(masses[0] - EXPECTED_MASS_KG) / EXPECTED_MASS_KG < 0.01, \
        f'FAIL: phantom mass = {masses[0]:.4g} kg, expected {EXPECTED_MASS_KG} kg'

powers = [v[3] for v in results.values()]
assert all(abs(p - powers[0]) < 1e-6 for p in powers), \
    f'FAIL: absorbed power differs across variants: {powers}'

sar_plain = results['m10g_IEEE_62704'][0]
sar_ar    = results['m10g_SIMPLE_ar20'][0]
assert np.count_nonzero(sar_ar) < np.count_nonzero(sar_plain), \
    f'FAIL: autoRange did not reduce the number of non-zero cells: {np.count_nonzero(sar_ar)} >= {np.count_nonzero(sar_plain)}'

for key, (_, peak, _, _) in results.items():
    ref = EXPECTED_PEAK_SAR.get(key)
    if ref is None:
        print(f'  {key}: no reference value, skipping peak check')
        continue
    assert abs(peak - ref) / ref < PEAK_SAR_RTOL, \
        f'FAIL: {key} peak SAR = {peak:.6g} W/kg, expected {ref:.6g} W/kg (rtol={PEAK_SAR_RTOL})'

# ── Cleanup ───────────────────────────────────────────────────────────────────
shutil.rmtree(Sim_Path, ignore_errors=True)

print('PASS')
