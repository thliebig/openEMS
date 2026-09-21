#!/usr/bin/env python3
"""Check that the installed openEMS modules and libraries contain a GPU backend.

usage: check_gpu_engine.py hip|metal

Runs a tiny simulation with engine='gpu' and checks the console output of the
engine: it must have created the backend, or reported that there is no device
for it. openEMS only prints that for a backend it was built with: a build without
it silently uses the reference backend on the CPU, which fails this check. So the
check works on machines without a GPU, e.g. CI runners.
"""

import ctypes
import ctypes.util
import os
import sys
import tempfile

from CSXCAD import ContinuousStructure
from openEMS import openEMS

EXPECTED = {'hip': ('backend: HIP', 'no HIP device found'),
            'metal': ('backend: Metal', 'no Metal device found')}


def run_captured(sim_path):
    """run a tiny simulation with the GPU engine, return the console output of openEMS"""
    FDTD = openEMS(NrTS=50, EndCriteria=0)
    FDTD.SetGaussExcite(5e9, 4e9)
    FDTD.SetBoundaryCond(['PEC'] * 6)
    CSX = ContinuousStructure()
    FDTD.SetCSX(CSX)
    mesh = CSX.GetGrid()
    mesh.SetDeltaUnit(1e-3)
    for ax in 'xyz':
        mesh.AddLine(ax, [0, 1, 2, 3, 4, 5, 6])
    CSX.AddExcitation('d', exc_type=0, exc_val=[0, 0, 1]).AddBox([3, 3, 2], [3, 3, 4])
    CSX.AddProbe('et', p_type=2).AddPoint([2, 2, 3])

    log_file = sim_path + '.log'
    sys.stdout.flush()
    sys.stderr.flush()
    saved = [os.dup(1), os.dup(2)]
    with open(log_file, 'w') as log:
        os.dup2(log.fileno(), 1)
        os.dup2(log.fileno(), 2)
        try:
            FDTD.Run(sim_path, cleanup=True, engine='gpu')
        finally:
            try:   # flush the C stdio buffers of the openEMS library
                ctypes.CDLL(ctypes.util.find_library('c') or 'ucrtbase').fflush(None)
            except (OSError, AttributeError, TypeError):
                pass
            os.dup2(saved[0], 1)
            os.dup2(saved[1], 2)
            for fd in saved:
                os.close(fd)
    with open(log_file) as log:
        return log.read()


def main():
    backend = sys.argv[1].lower()
    created, no_device = EXPECTED[backend]
    output = run_captured(os.path.join(tempfile.gettempdir(), 'check_gpu_engine'))
    lines = [l for l in output.splitlines() if 'GPU' in l or 'backend' in l]
    print('\n'.join(lines))
    if created in output:
        print(f'OK: the {backend} backend runs the simulation')
    elif no_device in output:
        print(f'OK: the {backend} backend is built in (no device on this machine)')
    else:
        sys.exit(f'FAIL: no {backend} backend in this build (engine=gpu used another backend)')


if __name__ == '__main__':
    main()
