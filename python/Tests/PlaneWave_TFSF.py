# -*- coding: utf-8 -*-
"""
 Plane-wave excitation test — total-field / scattered-field (TF/SF) box

 An empty TF/SF plane-wave box (excitation type 10) in free space.  Inside
 the box the incident plane wave must appear with the requested amplitude,
 polarisation and propagation direction; outside the box (scattered-field
 region) the field must vanish, since there is no scatterer.

 The test is repeated for normal incidence and for an oblique direction
 of incidence (all three wave-vector components non-zero).

 Pass criteria (per direction)
   inside:  peak |E| within 5 % of the excitation amplitude (1 V/m)
            polarisation: |E . e_dir| / |E| > 0.99 at the peak
            |E| / |H| within 2 % of Z0
            delay between two probes = (d . k_dir)/c within 2 %  (at 2 GHz)
   outside: peak |E| < -60 dB of the incident field (leakage)

 Tested with
  - python 3.13
  - openEMS v0.0.37+

 (c) 2026 Sean Mollet <sean@malmoset.com>

"""

import os, tempfile
import numpy as np

from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import C0, Z0
from openEMS.utilities import DFT_time2freq

unit     = 1e-3    # drawing unit: mm
res      = 2       # mesh resolution (mm), lambda/50 at f_stop
sim_box  = 60      # half size of the simulation domain (mm), excluding PML
pw_box   = 36      # half size of the TF/SF box (mm)

f_start = 0.5e9
f_stop  = 3.5e9
f_eval  = 2e9      # frequency for the delay check

# probe positions (mm)
p_in1 = np.array([-12, -12, -12])
p_in2 = np.array([ 12,  12,  12])
p_out = {'x-': [-48,   0,   0], 'x+': [ 48,  0,  0],
         'y-': [  0, -48,   0], 'y+': [  0, 48,  0],
         'z-': [  0,   0, -48], 'z+': [  0,  0, 48]}


def load_probe(Sim_Path, name):
    """ returns (t, field[3, Nt]) of a time-domain field probe """
    data = np.loadtxt(os.path.join(Sim_Path, name), comments='%')
    return data[:, 0], data[:, 1:4].T


def run_direction(label, k_dir, E_dir):
    k_dir = np.array(k_dir, dtype=float) / np.linalg.norm(k_dir)
    E_dir = np.array(E_dir, dtype=float) / np.linalg.norm(E_dir)
    assert abs(np.dot(k_dir, E_dir)) < 1e-12, 'E_dir must be perpendicular to k_dir'

    Sim_Path = os.path.join(tempfile.gettempdir(), 'PlaneWave_TFSF_' + label)

    FDTD = openEMS(EndCriteria=1e-5)
    FDTD.SetGaussExcite(0.5*(f_start + f_stop), 0.5*(f_stop - f_start))
    FDTD.SetBoundaryCond(['PML_8'] * 6)

    CSX = ContinuousStructure()
    FDTD.SetCSX(CSX)
    mesh = CSX.GetGrid()
    mesh.SetDeltaUnit(unit)
    for ax in 'xyz':
        mesh.AddLine(ax, np.arange(-sim_box - 8*res, sim_box + 8*res + res/2, res))

    pw = CSX.AddExcitation('plane_wave', exc_type=10, exc_val=E_dir)
    pw.SetPropagationDir(k_dir)
    pw.SetFrequency(f_eval)
    pw.AddBox([-pw_box]*3, [pw_box]*3)

    for name, pos in (('in1', p_in1), ('in2', p_in2)):
        CSX.AddProbe('et_' + name, p_type=2).AddPoint(pos)
        CSX.AddProbe('ht_' + name, p_type=3).AddPoint(pos)
    for name, pos in p_out.items():
        CSX.AddProbe('et_out_' + name, p_type=2).AddPoint(pos)

    FDTD.Run(Sim_Path, cleanup=True)

    ### inside the TF/SF box
    t, E1 = load_probe(Sim_Path, 'et_in1')
    _, H1 = load_probe(Sim_Path, 'ht_in1')
    _, E2 = load_probe(Sim_Path, 'et_in2')

    E_mag = np.linalg.norm(E1, axis=0)
    H_mag = np.linalg.norm(H1, axis=0)
    i_pk  = np.argmax(E_mag)
    E_pk  = E_mag[i_pk]
    pol   = abs(np.dot(E1[:, i_pk], E_dir)) / E_pk
    Z     = E_pk / np.max(H_mag)

    # delay between the two probes from the phase of the projected field
    e1 = DFT_time2freq(t, np.dot(E_dir, E1), [f_eval])[0]
    e2 = DFT_time2freq(t, np.dot(E_dir, E2), [f_eval])[0]
    delay     = -np.angle(e2 / e1) / (2*np.pi*f_eval)
    delay_ana = np.dot(p_in2 - p_in1, k_dir) * unit / C0

    print(f'  peak |E| = {E_pk:.4f} V/m,  polarisation = {pol:.4f},  |E|/|H| = {Z:.1f} Ohm')
    print(f'  delay in2-in1 = {delay*1e12:.1f} ps  (analytic {delay_ana*1e12:.1f} ps)')

    assert abs(E_pk - 1) < 0.05, \
        f'FAIL [{label}]: peak |E| = {E_pk:.4f} V/m inside the TF/SF box, expected 1 V/m (+/- 5 %)'
    assert pol > 0.99, \
        f'FAIL [{label}]: field polarisation {pol:.4f} along E_dir, expected > 0.99'
    assert abs(Z/Z0 - 1) < 0.02, \
        f'FAIL [{label}]: |E|/|H| = {Z:.1f} Ohm, expected Z0 = {Z0:.1f} Ohm (+/- 2 %)'
    assert abs(delay/delay_ana - 1) < 0.02, \
        f'FAIL [{label}]: delay {delay*1e12:.1f} ps, expected {delay_ana*1e12:.1f} ps (+/- 2 %)'

    ### outside the TF/SF box: no scatterer, no field
    leak = {name: np.max(np.linalg.norm(load_probe(Sim_Path, 'et_out_' + name)[1], axis=0)) / E_pk
            for name in p_out}
    worst = max(leak, key=leak.get)
    print(f'  max leakage outside the box: {20*np.log10(leak[worst]):.1f} dB ({worst})')
    assert leak[worst] < 1e-3, \
        f'FAIL [{label}]: field outside the TF/SF box at {worst}: {20*np.log10(leak[worst]):.1f} dB, expected < -60 dB'

    return t, E1, E2


results = {}
for label, k_dir, E_dir in (('normal',  [1, 0, 0], [0, 0, 1]),
                            ('oblique', [1, 2, 3], [2, -1, 0])):
    print('Testing {} incidence, k = {}'.format(label, k_dir))
    results[label] = run_direction(label, k_dir, E_dir)
    print('PASS [{}]'.format(label))

print('PASS')

if 0:  # set to 1 for debugging plots
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(len(results), 1, num='TF/SF probes', tight_layout=True, sharex=True)
    for ax, (label, (t, E1, E2)) in zip(axes, results.items()):
        ax.plot(t*1e9, np.linalg.norm(E1, axis=0), 'k-',  linewidth=2, label='|E| in1')
        ax.plot(t*1e9, np.linalg.norm(E2, axis=0), 'r--', linewidth=2, label='|E| in2')
        ax.set_title(label)
        ax.grid()
        ax.legend()
    axes[-1].set_xlabel('Time (ns)')
    plt.show()
