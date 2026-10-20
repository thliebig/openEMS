# SPDX-License-Identifier: GPL-3.0-or-later
"""Synthetic UPML allocation/decay fixture; optional, outside the test runner."""
import argparse
import json
import os
from pathlib import Path

import numpy as np
from CSXCAD import ContinuousStructure
from openEMS import openEMS

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('output', type=Path, help='new directory for this run')
parser.add_argument('--mode', choices=('fused', 'regions', 'separate', 'control'),
                    default='fused')
args = parser.parse_args()
args.output.mkdir(parents=True, exist_ok=False)

for key in ('OPENEMS_HIP_FUSED_STEP', 'OPENEMS_HIP_UPML_FUSION'):
    os.environ.pop(key, None)
if args.mode == 'regions':
    os.environ['OPENEMS_HIP_FUSED_STEP'] = '0'
elif args.mode == 'separate':
    os.environ['OPENEMS_HIP_UPML_FUSION'] = '0'

fdtd = openEMS(NrTS=5000, EndCriteria=1e-6)
fdtd.SetGaussExcite(5.5e9, 4.5e9)
fdtd.SetBoundaryCond(['PEC' if args.mode == 'control' else 'PML_8'] * 6)
csx = ContinuousStructure()
fdtd.SetCSX(csx)
mesh = csx.GetGrid()
mesh.SetDeltaUnit(1e-3)
for axis in 'xyz':
    mesh.AddLine(axis, np.arange(-40, 41))
# Dampen the residual static field as well as the propagating pulse.
csx.AddMaterial('loss', epsilon=1, kappa=0.02).AddBox([-40]*3, [40]*3)
csx.AddExcitation('dipole', exc_type=0, exc_val=[0, 0, 1]).AddBox([0, 0, -1], [0, 0, 1])
csx.AddProbe('et', p_type=2).AddPoint([4, 3, 2])
csx.AddProbe('ht', p_type=3).AddPoint([-3, 5, 0])
fdtd.Write2XML(str(args.output / 'model.xml'))
fdtd.Run(str(args.output), engine='gpu', exact_endcriteria=True, verbose=2,
         numThreads=2)

probes = {}
for name in ('et', 'ht'):
    values = np.loadtxt(args.output / name, comments='%')[:, 1:]
    peak = float(np.max(np.abs(values)))
    assert np.isfinite(values).all() and peak > 0, f'Invalid {name} probe'
    probes[name] = dict(samples=len(values), peak=peak)
record = dict(mode=args.mode, probes=probes)
(args.output / 'summary.json').write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps(record))
