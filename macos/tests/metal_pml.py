#!/usr/bin/env python3
"""Compare GPU UPML with CPU UPML on the same Metal fields, and with SSE."""

import argparse
import os
from pathlib import Path
import shutil
import re
import tempfile

from metal_fields import compare, h5_arrays, make_model, np, run


def check_fields(reference, result, rtol, atol, l2_limit, pointwise=True):
    for field in ('Et.h5', 'Ht.h5'):
        for directory in (reference, result):
            arrays = h5_arrays(directory / field)
            samples = [a for name, a in arrays.items() if name.startswith('FieldData/')]
            if not samples or not any(np.any(a) for a in samples):
                raise AssertionError('Missing or zero field data: ' + str(directory / field))
            if not all(np.isfinite(a).all() for a in arrays.values()):
                raise AssertionError('Non-finite data: ' + str(directory / field))
        stats = compare(reference / field, result / field, rtol, atol)
        print('    {}: max abs {:.6g}, relative L2 {:.6g}, pointwise failures {}'.format(
            field, stats['max_abs'], stats['relative_l2'], stats['failures']))
        if (pointwise and stats['failures']) or stats['relative_l2'] > l2_limit:
            raise AssertionError('Field tolerance exceeded: ' + field)


def run_case(args, label, cells, steps, boundaries, nonuniform=False):
    root = Path(tempfile.mkdtemp(prefix='openems-metal-pml-'))
    try:
        model = root / 'model.xml'
        # A short pulse exercises propagation into PML and post-pulse decay.
        make_model(model, cells, steps, nonuniform, boundaries, frequency=5e9)
        sse_time, _ = run(args.openems, model, 'sse', root / 'sse')
        cpu_time, cpu_log = run(args.openems, model, 'metal', root / 'cpu',
                                args.fp64_reference, compress=True, pml=False)
        gpu_time, gpu_log = run(args.openems, model, 'metal', root / 'gpu',
                                args.fp64_reference, compress=True, pml=True)
        _, dense_log = run(args.openems, model, 'metal', root / 'dense',
                           args.fp64_reference, compress=False, pml=True)
        old_layout = os.environ.get('OPENEMS_METAL_PML_LAYOUT')
        try:
            os.environ['OPENEMS_METAL_PML_LAYOUT'] = 'scalar'
            _, scalar_log = run(args.openems, model, 'metal', root / 'scalar',
                                args.fp64_reference, compress=True, pml=True)
        finally:
            if old_layout is None:
                os.environ.pop('OPENEMS_METAL_PML_LAYOUT', None)
            else:
                os.environ['OPENEMS_METAL_PML_LAYOUT'] = old_layout
        # Isolate each new optimization against the identical arithmetic path.
        variants = [('duplicate-storage', 'OPENEMS_METAL_PML_REUSE'),
                    ('dense-pml', 'OPENEMS_METAL_PML_COMPRESS'),
                    ('late-ec-free', 'OPENEMS_METAL_EARLY_EC_FREE'),
                    ('serial-coefficients', 'OPENEMS_METAL_SERIAL_COEFFICIENTS')]
        variant_logs = {}
        for label_variant, setting in variants:
            previous = os.environ.get(setting)
            try:
                os.environ[setting] = '1' if setting == 'OPENEMS_METAL_SERIAL_COEFFICIENTS' else '0'
                _, variant_logs[label_variant] = run(args.openems, model, 'metal', root / label_variant,
                                                     args.fp64_reference, compress=True, pml=True)
            finally:
                if previous is None: os.environ.pop(setting, None)
                else: os.environ[setting] = previous
        expected = sum(b.startswith('PML') for b in boundaries)
        marker = 'Metal: GPU UPML conditioning: {} regions'.format(expected)
        for log in (gpu_log, dense_log):
            if expected and (marker not in log or 'Metal: UPML layout: indexed' not in log):
                raise AssertionError('Indexed GPU UPML not exercised: ' + log)
            if not expected and 'Metal: GPU UPML conditioning:' in log:
                raise AssertionError('Unexpected UPML regions')
        if expected and 'Metal: UPML layout: scalar' not in scalar_log:
            raise AssertionError('Scalar GPU UPML fallback not exercised')
        if 'Metal: CPU UPML conditioning selected' not in cpu_log:
            raise AssertionError('CPU UPML fallback not exercised')
        print('{}: {}, {} steps; SSE/CPU-PML/GPU-PML wall {:.3f}/{:.3f}/{:.3f}s'.format(
            label, cells, steps, sse_time, cpu_time, gpu_time))
        print('  CPU/GPU UPML:')
        check_fields(root / 'cpu', root / 'gpu', args.rtol, args.atol, args.l2)
        print('  SSE/GPU UPML:')
        # SSE/Metal already differ near cancellation zeros in long runs.
        # Gate that comparison on global relative L2; CPU/GPU PML above
        # also enforces pointwise tolerances, isolating this implementation.
        check_fields(root / 'sse', root / 'gpu', args.rtol, args.atol, args.l2,
                     pointwise=False)
        for field in ('Et.h5', 'Ht.h5'):
            reference = h5_arrays(root / 'gpu' / field)
            for label_variant, _ in variants:
                result = h5_arrays(root / label_variant / field)
                if result.keys() != reference.keys():
                    raise AssertionError('Optimization datasets differ: ' + label_variant)
                for key, a in reference.items():
                    b = result[key]
                    if a.shape != b.shape or a.dtype != b.dtype or a.tobytes() != b.tobytes():
                        raise AssertionError('Optimization bits differ: ' + label_variant + '/' + field + '/' + key)
            dense = h5_arrays(root / 'dense' / field)
            packed = h5_arrays(root / 'gpu' / field)
            scalar = h5_arrays(root / 'scalar' / field)
            if scalar.keys() != packed.keys():
                raise AssertionError('Scalar/indexed datasets differ')
            for name, a in scalar.items():
                b = packed[name]
                if a.shape != b.shape or a.dtype != b.dtype or a.tobytes() != b.tobytes():
                    raise AssertionError('Scalar/indexed bits differ: ' + field + '/' + name)
            if dense.keys() != packed.keys():
                raise AssertionError('Dense/compressed datasets differ')
            for name, a in dense.items():
                b = packed[name]
                if a.shape != b.shape or a.dtype != b.dtype or a.tobytes() != b.tobytes():
                    raise AssertionError('Dense/compressed bits differ: ' + field + '/' + name)
        if args.fp64_reference:
            prefix = 'Metal FP64 update reference:'
            diagnostic = [s for s in gpu_log.splitlines() if s.startswith(prefix)]
            if not diagnostic or diagnostic != [s for s in dense_log.splitlines() if s.startswith(prefix)]:
                raise AssertionError('Missing or inconsistent FP64 diagnostics')
            if re.search(r'\b(?:nan|inf)\b', diagnostic[-1], re.IGNORECASE):
                raise AssertionError('Non-finite FP64 diagnostic')
            print('  ' + diagnostic[-1])
        if expected:
            match = re.search(r'Metal: UPML duplicate bytes avoided: (\d+)', gpu_log)
            if not match or int(match[1]) <= 0:
                raise AssertionError('UPML storage reuse not exercised')
            if 'Metal: UPML duplicate bytes avoided: 0' not in variant_logs['duplicate-storage']:
                raise AssertionError('Duplicate-storage control not exercised')
            if not any(s in gpu_log for s in ('Metal: lossless UPML coefficients:', 'Metal: UPML coefficient dictionary fallback')):
                raise AssertionError('UPML compression/fallback not exercised')
        print('  Dense/compressed, scalar/indexed and reuse/copy: bit-identical')
    finally:
        if args.keep:
            print('  kept ' + str(root))
        else:
            shutil.rmtree(root)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--openems', default=os.environ.get('OPENEMS_BIN', 'openEMS'))
    parser.add_argument('--rtol', type=float, default=5e-4)
    parser.add_argument('--atol', type=float, default=2e-6)
    parser.add_argument('--l2', type=float, default=1e-4)
    parser.add_argument('--fp64-reference', action='store_true')
    parser.add_argument('--keep', action='store_true')
    args = parser.parse_args()
    # Metal optimizations are defaults; per-case controls explicitly turn them off.
    if os.environ.get('OPENEMS_METAL_PML_LAYOUT', 'indexed') != 'indexed':
        parser.error('Unset OPENEMS_METAL_PML_LAYOUT to test the default indexed path')
    cases = [
        ('no-pml', (17, 19, 21), 100, ['PEC'] * 6, False),
        ('six-faces-corners', (24, 23, 22), 400, ['PML_4'] * 6, False),
        ('overlapping-slabs', (8, 15, 14), 400,
         ['PML_4', 'PML_4', 'PEC', 'PEC', 'PEC', 'PEC'], False),
        ('nonuniform-asymmetric', (29, 27, 25), 400,
         ['PML_3', 'PML_5', 'PML_4', 'PML_3', 'PML_5', 'PML_4'], True),
        ('mixed-mur-pec-pml', (25, 24, 23), 400,
         ['PML_4', 'MUR', 'PEC', 'PML_3', 'PML_4', 'PEC'], False),
        ('long-run', (32, 31, 30), 1200, ['PML_8'] * 6, False),
        ('thin-board', (65, 66, 32), 600, ['PML_8'] * 6, True),
    ]
    # Each physical face alone: includes all four z-line-count remainders.
    for face in range(6):
        boundaries = ['PEC'] * 6
        boundaries[face] = 'PML_3'
        cases.append(('face-' + str(face), (21, 22, 20 + face % 4), 250,
                      boundaries, False))
    for case in cases:
        run_case(args, *case)


if __name__ == '__main__':
    main()
