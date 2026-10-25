# -*- coding: utf-8 -*-
"""
 Lossy metal test — surface impedance of LossyMetal bodies

 1. Plane-wave transmission through a thin copper layer (1D channel of PEC and
    PMC walls, PML at both z-ends, as in Conducting_Sheet.py). The layer is a
    LossyMetal body thinner than the mesh, which is modelled as a sheet:

        T(f) = 2 / (2 + Z0 Y(f)),   Y = 2 sigma / gamma * tanh(gamma t / 2)

    once as a 2D primitive with the Thickness of the property and once as a
    3D box of thickness t collapsed onto one mesh line.

 2. A rectangular TE101 cavity with LossyMetal walls (3D boxes, the inner
    faces on mesh lines). The quality factor from the ring-down is compared
    with the closed form (Pozar, Microwave Engineering, eq. 6.46)

        Qc = (k a d)^3 b eta / (2 pi^2 Rs) / (2 a^3 b + 2 b d^3 + a^3 d + a d^3)

    smooth and with the gradient roughness model, for which Rs is the real
    part of the rough surface impedance from a 1D solution of the model. The
    conductivity is lowered to keep the ring-down short.

 Pass criteria
   1: |dB(T_fdtd) - dB(T_analytic)| < 0.4 dB, |arg(T_fdtd / T_analytic)| < 3 deg
      (1 .. 10 GHz)
   2: Q within 5 % of the closed form, smooth and rough

 Tested with
  - python 3.12
  - openEMS v0.0.37+

"""

import os, tempfile
import numpy as np
from math import erf

from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import Z0, MUE0, C0
from openEMS.utilities import DFT_time2freq

unit = 1e-3   # drawing unit: mm


### 1. Plane wave through a thin copper layer
dz, width, length = 0.25, 1.0, 150
z_src, z_sheet, z_probe = 20, 60, 110
freq = np.linspace(1e9, 10e9, 91)


def run_layer(name, kind=None, sigma=0, t=0):
    FDTD = openEMS(NrTS=40000, EndCriteria=1e-5)
    FDTD.SetGaussExcite(5.5e9, 4.5e9)
    FDTD.SetBoundaryCond(['PEC', 'PEC', 'PMC', 'PMC', 'PML_8', 'PML_8'])
    CSX = ContinuousStructure()
    FDTD.SetCSX(CSX)
    mesh = CSX.GetGrid()
    mesh.SetDeltaUnit(unit)
    mesh.AddLine('x', [0, width/2, width])
    mesh.AddLine('y', [0, width/2, width])
    mesh.AddLine('z', np.arange(0, length + dz, dz))
    CSX.AddExcitation('plane', exc_type=0, exc_val=[1, 0, 0]).AddBox([0, 0, z_src], [width, width, z_src])
    CSX.AddProbe('et', p_type=2).AddPoint([width/2, width/2, z_probe])
    if kind == '2D':
        lm = CSX.AddLossyMetal('layer', conductivity=sigma, thickness=t)
        lm.AddBox([0, 0, z_sheet], [width, width, z_sheet], priority=10)
    elif kind == '3D':
        lm = CSX.AddLossyMetal('layer', conductivity=sigma)
        lm.AddBox([0, 0, z_sheet - t/2/unit], [width, width, z_sheet + t/2/unit], priority=10)
    Sim_Path = os.path.join(tempfile.gettempdir(), 'LossyMetal_' + name)
    FDTD.Run(Sim_Path, cleanup=True, exact_endcriteria=True)
    data = np.loadtxt(os.path.join(Sim_Path, 'et'), comments='%')
    return DFT_time2freq(data[:, 0], data[:, 1], freq)


print('Running reference (no layer)...')
E_ref = run_layer('Ref')
for sigma, t in ((5.8e7, 1e-6), (5.8e7, 18e-6)):
    gamma = np.sqrt(1j*2*np.pi*freq*MUE0*sigma)
    T_ana = 2/(2 + Z0*2*sigma/gamma*np.tanh(gamma*t/2))
    for kind in ('2D', '3D'):
        name = 'layer_{}_{:g}um'.format(kind, t*1e6)
        T = run_layer(name, kind, sigma, t)/E_ref
        dB_err = np.max(np.abs(20*np.log10(np.abs(T/T_ana))))
        ph_err = np.max(np.abs(np.angle(T/T_ana, deg=True)))
        print(f'  {name}: max |dB error| = {dB_err:.3f} dB, max |phase error| = {ph_err:.2f} deg')
        assert dB_err < 0.4, f'FAIL [{name}]: |T| deviates {dB_err:.3f} dB from analytic, expected < 0.4 dB'
        assert ph_err < 3, f'FAIL [{name}]: arg(T) deviates {ph_err:.2f} deg from analytic, expected < 3 deg'
        print(f'PASS [{name}]')


### 2. TE101 cavity with lossy (and rough) walls
a, b, d, wall, res = 22.86, 10.16, 20.0, 1.0, 0.5   # mm
sigma_cav = 5.8e6         # lowered conductivity: Q ~ 2500, short ring-down
rq_cav = 6.3e-6           # Rq / skin depth ~ 3 at 10 GHz


def zs_gradient(f, sigma, rq):
    """ surface impedance of the gradient roughness model (mean plane reference),
        integrating Y' = j w mu0 Y^2 - sigma(x) out of the metal """
    w = 2*np.pi*f
    g = np.sqrt(1j*w*MUE0*sigma)
    if rq <= 0:
        return g/sigma
    delta = np.sqrt(2/(w*MUE0*sigma))
    x0, x = -6*rq, 6*rq
    n = int(np.ceil((x - x0)/(min(rq, delta)/40)))
    h = (x - x0)/n
    Y = sigma/g
    s = lambda x: sigma*0.5*(1 + erf(x/rq/np.sqrt(2)))
    fy = lambda x, Y: 1j*w*MUE0*Y*Y - s(x)
    for _ in range(n):
        k1 = fy(x, Y); k2 = fy(x - h/2, Y - h/2*k1); k3 = fy(x - h/2, Y - h/2*k2); k4 = fy(x - h, Y - h*k3)
        Y -= h/6*(k1 + 2*k2 + 2*k3 + k4)
        x -= h
    return 1/Y - 1j*w*MUE0*abs(x0)


def run_cavity(name, rq):
    FDTD = openEMS(NrTS=60000, EndCriteria=0)
    FDTD.SetGaussExcite(10e9, 3e9)
    FDTD.SetBoundaryCond(['PEC']*6)
    CSX = ContinuousStructure()
    FDTD.SetCSX(CSX)
    mesh = CSX.GetGrid()
    mesh.SetDeltaUnit(unit)
    for ax, L in (('x', a), ('y', b), ('z', d)):
        mesh.AddLine(ax, np.r_[-wall, -wall/2, np.linspace(0, L, int(np.ceil(L/res)) + 1), L + wall/2, L + wall])
    lm = CSX.AddLossyMetal('wall', conductivity=sigma_cav, roughness_model='gradient' if rq > 0 else 'none', rq=rq)
    W = wall
    for p1, p2 in (([-W, -W, -W], [a+W, b+W, 0]), ([-W, -W, d], [a+W, b+W, d+W]),
                   ([-W, -W, 0], [0, b+W, d]), ([a, -W, 0], [a+W, b+W, d]),
                   ([0, -W, 0], [a, 0, d]), ([0, b, 0], [a, b+W, d])):
        lm.AddBox(p1, p2, priority=10)
    x_src, z_src = mesh.GetLine('x', 2 + 15), mesh.GetLine('z', 2 + 13)
    CSX.AddExcitation('src', exc_type=0, exc_val=[0, 1, 0]).AddBox([x_src, 0, z_src], [x_src, b, z_src])
    CSX.AddProbe('et', p_type=2).AddPoint([0.6*a, b/2, 0.6*d])
    Sim_Path = os.path.join(tempfile.gettempdir(), 'LossyMetal_' + name)
    FDTD.Run(Sim_Path, cleanup=True)

    data = np.loadtxt(os.path.join(Sim_Path, 'et'), comments='%')
    t, ey = data[:, 0], data[:, 2]
    N, dt = len(t), t[1] - t[0]
    f = np.fft.rfftfreq(N, dt)
    S = np.fft.rfft(ey)
    f_res = f[np.argmax(np.abs(S)*(np.abs(f - 10e9) < 2e9))]
    S[np.abs(f - f_res) > 1e9] = 0          # band pass around the TE101 resonance
    env = np.abs(np.fft.ifft(np.r_[S, np.zeros(N - len(S))]*2))   # analytic signal envelope
    sel = (t > 0.15*t[-1]) & (t < 0.85*t[-1])
    tau = -1/np.polyfit(t[sel], np.log(env[sel]), 1)[0]
    Q = np.pi*f_res*tau

    k = 2*np.pi*f_res/C0
    A, B, D = a*unit, b*unit, d*unit
    Rs = np.real(zs_gradient(f_res, sigma_cav, rq))
    Qc = (k*A*D)**3*B*Z0/(2*np.pi**2*Rs)/(2*A**3*B + 2*B*D**3 + A**3*D + A*D**3)
    return f_res, Q, Qc


for name, rq in (('cavity_smooth', 0), ('cavity_rough', rq_cav)):
    f_res, Q, Qc = run_cavity(name, rq)
    err = Q/Qc - 1
    print(f'  {name}: f = {f_res/1e9:.4f} GHz, Q = {Q:.0f}, closed form {Qc:.0f}, error {err*100:+.1f} %')
    assert abs(err) < 0.05, f'FAIL [{name}]: Q = {Q:.0f}, expected {Qc:.0f} +/- 5 %'
    print(f'PASS [{name}]')

print('PASS')
