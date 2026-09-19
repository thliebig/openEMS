# -*- coding: utf-8 -*-
"""
 Free-space benchmark — speed of the FDTD updates on a large uniform mesh

 A soft z-dipole in the center of an N^3 uniform mesh (1 mm cells, default
 N = 300, i.e. 27 million cells) radiates a Gaussian pulse (5 +- 4 GHz) for a
 fixed number of timesteps (default 800) without an end criterion. The run is
 made twice: with PML_8 on all sides and with PEC walls. One field probe, no
 dumps, so the time is spent in the field updates (and, with PML, the UPML
 update). This is the "Free space, 300^3 cells" benchmark of benchmarks.md.

 Usage: FreeSpace_Benchmark.py [engine] [N] [timesteps]
   engine: gpu (default), gpu-reference, multithreaded, sse-compressed, basic, ...

 Pass criteria (per boundary)
   the run did all timesteps and reported its speed (MCells/s)
   engine=gpu: the GPU engine was created
 and
   the probe signals of both runs agree to 1e-4 of the peak value until the
     first wall reflection of the PEC run can reach the probe, and differ by
     more than 1e-2 afterwards

 Tested with
  - python 3.13
  - openEMS v0.0.37+

 (c) 2026 Sean Mollet <sean@malmoset.com>

"""

import os, re, sys, tempfile, ctypes, ctypes.util
import numpy as np

from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import C0

engine = sys.argv[1] if len(sys.argv) > 1 else 'gpu'
N      = int(sys.argv[2]) if len(sys.argv) > 2 else 300
NrTS   = int(sys.argv[3]) if len(sys.argv) > 3 else 800

unit    = 1e-3                   # drawing unit: mm, one cell per mm
f0, fc  = 5e9, 4e9               # Gaussian excitation
src     = [N/2, N/2, N/2]        # dipole center
probe   = [N/2 + 10, N/2, N/2]   # 10 mm from the dipole
RTOL    = 1e-4


def run(Sim_Path, bc):
    """ Run the benchmark, return the openEMS console output and the probe signal (t, Ez) """
    FDTD = openEMS(NrTS=NrTS, EndCriteria=0)
    FDTD.SetGaussExcite(f0, fc)
    FDTD.SetBoundaryCond([bc]*6)

    CSX = ContinuousStructure()
    FDTD.SetCSX(CSX)
    mesh = CSX.GetGrid()
    mesh.SetDeltaUnit(unit)
    for ax in 'xyz':
        mesh.AddLine(ax, np.arange(N + 1)*1.0)

    exc = CSX.AddExcitation('dipole', exc_type=0, exc_val=[0, 0, 1])
    exc.AddBox([src[0], src[1], src[2] - 1], [src[0], src[1], src[2] + 1])
    CSX.AddProbe('et', p_type=2).AddPoint(probe)

    # capture the console output of the openEMS library
    log_file = Sim_Path + '.log'
    sys.stdout.flush()
    saved = os.dup(1)
    with open(log_file, 'w') as log:
        os.dup2(log.fileno(), 1)
        try:
            FDTD.Run(Sim_Path, cleanup=True, engine=engine)
        finally:
            try:   # flush the C stdio buffer of the openEMS library
                ctypes.CDLL(ctypes.util.find_library('c')).fflush(None)
            except (OSError, AttributeError, TypeError):
                pass
            os.dup2(saved, 1)
            os.close(saved)
    with open(log_file) as log:
        output = log.read()

    # field probe file columns: t/s, Ex, Ey, Ez
    data = np.loadtxt(os.path.join(Sim_Path, 'et'), comments='%')
    return output, data[:, 0], data[:, 3]


print(f'Free-space benchmark: {N}^3 = {N**3/1e6:.1f} million cells, {NrTS} timesteps, engine={engine}')
results = {}
for bc in ('PML_8', 'PEC'):
    output, t, Ez = run(os.path.join(tempfile.gettempdir(), f'FreeSpace_Benchmark_{bc}'), bc)

    created = re.search(r'Create FDTD engine \((.*)\)\s*$', output, re.M)
    speed   = re.search(r'Speed: *([0-9.]+) MCells/s', output)
    steps   = re.search(r'Time for (\d+) iterations', output)
    assert speed and steps, f'FAIL [{bc}]: the run reported no speed'
    assert int(steps.group(1)) >= NrTS, \
        f'FAIL [{bc}]: the run ended after {steps.group(1)} of {NrTS} timesteps'
    if engine == 'gpu':
        assert created and created.group(1).startswith('GPU'), f'FAIL [{bc}]: engine=gpu did not create the GPU engine'
    print(f'  {bc:6s} {float(speed.group(1)):10.1f} MCells/s  ({created.group(1) if created else engine})')
    results[bc] = (t, Ez)

### the boundaries may only change the probe signal once the first reflection can arrive
# shortest path source -> wall -> probe: to the mirror image of the source in
# the inner surface of the PML (8 cells from the border), the nearer wall for both runs
images = [src[:i] + [2*w - src[i]] + src[i+1:] for i in range(3) for w in (8, N - 8)]
path   = min(np.linalg.norm(np.subtract(probe, image)) for image in images)
t_refl = 0.95 * path*unit / C0

(t_pml, E_pml), (t_pec, E_pec) = results['PML_8'], results['PEC']
assert np.array_equal(t_pml, t_pec), 'FAIL: the probe times of both runs differ'
m = t_pml < t_refl
assert t_pml[-1] > t_refl, f'FAIL: the run ends before the first reflection ({t_pml[-1]*1e9:.2f} ns)'
peak = np.max(np.abs(E_pml[m]))
dev  = np.max(np.abs(E_pml[m] - E_pec[m])) / peak
print(f'  PML vs PEC before the first reflection ({t_refl*1e9:.2f} ns): max. deviation {dev:.1e} of the peak value')
assert dev < RTOL, f'FAIL: PML and PEC runs deviate {dev:.1e} before the first reflection, expected < {RTOL:.0e}'
dev_after = np.max(np.abs(E_pml[~m] - E_pec[~m])) / peak
print(f'  PML vs PEC after the first reflection: max. deviation {dev_after:.1e} of the peak value')
assert dev_after > 100*RTOL, f'FAIL: PML and PEC runs agree after the first reflection, the boundaries had no effect'

print('PASS')
