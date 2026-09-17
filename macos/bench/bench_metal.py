#!/usr/bin/env python3
"""Benchmark the openEMS Metal engine against the fastest SSE + multithread engine.

The Metal engine is compared against the multithreaded engine (compressed SSE +
boost threads), sweeping thread counts to find the fastest CPU configuration.
Runs are interleaved across configurations and repetitions, and medians are
reported, so thermal/scheduler drift affects every configuration equally.

This measures performance only. Correctness (SSE vs Metal field equality) is
covered by the other metal_*.py scripts under macos/tests/.

See ../doc/metal-benchmark.rst for the exact methodology and representative
results.

Usage:
    python bench_metal.py --openems /path/to/openEMS
    python bench_metal.py --openems /path/to/openEMS --cells 256 256 256 --steps 1000 --reps 3
"""

import argparse
import json
import os
import platform
import re
import statistics
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path

HERE = Path(__file__).resolve().parent          # macos/bench
REPO = HERE.parents[1]                           # repository root
MACOS_TESTS = REPO / 'macos' / 'tests'
sys.path.insert(0, str(MACOS_TESTS))
from metal_fields import make_model  # noqa: E402


# /usr/bin/time reports peak resident set size per run: -l on macOS (bytes),
# -v on Linux (kbytes). If it is missing, memory is reported as "-".
_TIME_BIN = Path('/usr/bin/time')
_TIME_ARGS = ['-l'] if platform.system() == 'Darwin' else ['-v']


def parse_peak_rss(text):
    if not _TIME_BIN.exists():
        return None
    if platform.system() == 'Darwin':
        m = re.search(r'^\s*(\d+)\s+maximum resident set size', text, re.M)
        if m:
            return int(m.group(1)) / (1024.0 * 1024.0)
    else:
        m = re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)', text)
        if m:
            return int(m.group(1)) / 1024.0
    return None


def find_binary(explicit):
    if explicit:
        return Path(explicit).resolve()
    env = os.environ.get('OPENEMS_BIN')
    if env:
        return Path(env).resolve()
    candidates = [
        REPO / 'build' / 'openEMS',
        REPO / 'install' / 'bin' / 'openEMS',
        HERE / 'openEMS',
    ]
    for c in candidates:
        if c.exists():
            return c.resolve()
    raise SystemExit('openEMS binary not found; pass --openems /path/to/openEMS '
                     'or set OPENEMS_BIN')


def strip_dumps(model):
    tree = ET.parse(model)
    props = tree.getroot().find('.//Properties')
    for child in list(props):
        if child.tag == 'DumpBox':
            props.remove(child)
    tree.write(model)


def make_model_file(path, cells, steps, boundaries):
    make_model(path, cells, steps, boundaries=boundaries, frequency=5e9)
    strip_dumps(path)


def parse_output(stdout):
    """Extract the reported cell count, iteration count and stepping time."""
    info = {'step': None, 'cells': None, 'iters': None, 'speed': None, 'notes': []}
    for line in stdout.splitlines():
        if line.startswith('Time for '):
            left, right = line.split(':', 1)          # ... : 1.80285 sec
            info['step'] = float(right.split()[0])
            tok = left.split()                        # Time for 600 iterations with 3.01e+06 cells
            info['iters'] = int(tok[2])
            info['cells'] = float(tok[5])
        elif line.startswith('Speed:'):
            info['speed'] = float(line.split()[1])
        if 'Metal: UPML layout:' in line:
            info['notes'].append(line.split(':', 1)[1].strip())
        if 'Metal: in-place diamond E/H pipeline:' in line:
            info['notes'].append('diamond' if line.rstrip().endswith('enabled') else 'legacy')
        if line.startswith('Metal: in-place diamond update:'):
            info['notes'].append(line.split(':', 2)[2].strip())
        if 'coefficient dictionary limit' in line:
            info['notes'].append('dense-fallback')
    info['note'] = '+'.join(dict.fromkeys(info['notes']))
    return info


def run_once(binary, model, engine, threads, outdir, env_overrides=None):
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [str(binary), str(model), '--engine=' + engine]
    if threads is not None:
        cmd.append('--numThreads=%d' % threads)
    if _TIME_BIN.exists():
        cmd = [str(_TIME_BIN)] + _TIME_ARGS + cmd
    start = time.perf_counter()
    env = os.environ.copy()
    if env_overrides:
        env.update(env_overrides)
    proc = subprocess.run(cmd, cwd=outdir, env=env, text=True,
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    wall = time.perf_counter() - start
    if proc.returncode:
        sys.stderr.write(proc.stdout[-4000:])
        raise SystemExit('%s --numThreads=%s failed (rc=%d)'
                         % (engine, threads, proc.returncode))
    info = parse_output(proc.stdout)
    if info['step'] is None:
        raise SystemExit('could not parse stepping time from output:\n' + proc.stdout[-4000:])
    info['peak_rss_mb'] = parse_peak_rss(proc.stdout)
    info['wall'] = wall
    info['setup'] = wall - info['step']
    return info


def median_config(infos, threads):
    step = statistics.median(i['step'] for i in infos)
    wall = statistics.median(i['wall'] for i in infos)
    setup = statistics.median(i['setup'] for i in infos)
    cells = infos[-1]['cells']
    iters = infos[-1]['iters']
    rss = [i['peak_rss_mb'] for i in infos if i.get('peak_rss_mb') is not None]
    return {
        'engine': 'metal' if threads is None else 'multithreaded',
        'threads': threads,
        'step_s': round(step, 4),
        'wall_s': round(wall, 4),
        'setup_s': round(setup, 4),
        'mcells_per_s': round(cells * iters / step / 1e6, 1),
        'peak_rss_mb': round(statistics.median(rss), 1) if rss else None,
        'note': infos[-1]['note'],
        'runs': [{'wall_s': round(i['wall'], 4), 'step_s': round(i['step'], 4),
                  'peak_rss_mb': (round(i['peak_rss_mb'], 1)
                                  if i.get('peak_rss_mb') is not None else None)}
                 for i in infos],
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--openems', help='path to the openEMS binary')
    ap.add_argument('--cells', type=int, nargs=3, default=[192, 160, 96],
                    help='grid cell counts in x y z (default: 192 160 96)')
    ap.add_argument('--steps', type=int, default=600, help='number of timesteps')
    ap.add_argument('--reps', type=int, default=5,
                    help='repetitions per configuration (default: 5)')
    ap.add_argument('--boundaries', default='PEC',
                    help="boundary list or single value, e.g. 'PEC' or "
                         "'PML_8,PML_8,PEC,PEC,PML_8,PML_8' (default: PEC)")
    ap.add_argument('--mt-threads', default='8,10,14',
                    help='comma-separated thread counts for the multithreaded engine '
                         '(default: 8,10,14)')
    ap.add_argument('--engines', choices=['both', 'metal', 'mt'], default='both',
                    help="which engines to run (default: both); 'mt' is useful for a "
                         'CPU-only baseline or a harness smoke test')
    ap.add_argument('--include-auto', action='store_true',
                    help='also run multithreaded without --numThreads '
                         '(note: the auto path runs single-threaded, see README)')
    ap.add_argument('--metal-legacy', action='store_true',
                    help='explicitly benchmark the legacy two-kernel Metal path')
    ap.add_argument('--compare-metal-legacy', action='store_true',
                    help='benchmark both in-place diamond and legacy Metal paths')
    ap.add_argument('--outdir', default=str(HERE / 'out'),
                    help='output directory (default: ./out)')
    ap.add_argument('--json', help='write the summary as JSON to this path')
    ap.add_argument('--tag', default='run', help='label used in output subdirectories')
    args = ap.parse_args()

    binary = find_binary(args.openems)
    if not binary.exists():
        raise SystemExit('openEMS binary not found: %s' % binary)

    boundaries = args.boundaries.split(',')
    if len(boundaries) == 1:
        boundaries *= 6
    if len(boundaries) != 6:
        raise SystemExit('--boundaries needs 1 or 6 values')

    out = Path(args.outdir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    model = out / 'model.xml'
    make_model_file(model, tuple(args.cells), args.steps, boundaries)

    threads = [int(t) for t in args.mt_threads.split(',') if t.strip()]
    configs = []
    if args.engines in ('both', 'metal'):
        if not args.metal_legacy:
            configs.append(('metal', None, {}))
        if args.metal_legacy or args.compare_metal_legacy:
            configs.append(('metal-legacy', None,
                            {'OPENEMS_METAL_FUSED_PIPELINE': '0'}))
    if args.engines in ('both', 'mt'):
        configs += [('mt-%d' % t, t, {}) for t in threads]
        if args.include_auto:
            configs.append(('mt-auto', 'auto', {}))
    if not configs:
        raise SystemExit('no engines selected')

    print('openEMS : %s' % binary)
    print('model   : %s cells, %s boundaries, %d steps'
          % ('x'.join(map(str, args.cells)), ','.join(boundaries), args.steps))
    print('configs : %s' % ', '.join(name for name, _, _ in configs))
    print('reps    : %d (interleaved, medians reported)' % args.reps)

    samples = {name: [] for name, _, _ in configs}
    for rep in range(args.reps):
        order = configs if rep % 2 == 0 else list(reversed(configs))
        for name, thr, extra_env in order:
            samples[name].append(
                run_once(binary, model, 'metal' if thr is None else 'multithreaded',
                         None if thr == 'auto' else thr, out / ('%s-%s-r%d' % (args.tag, name, rep)),
                         extra_env))

    results = []
    for name, thr, _ in configs:
        r = median_config(samples[name], None if thr is None else thr)
        r['name'] = name
        results.append(r)

    metal = next((r for r in results if r['name'] in ('metal', 'metal-legacy')), None)
    mt = [r for r in results if r['name'].startswith('mt-') and r['threads']]
    best = min(mt, key=lambda r: r['step_s']) if mt else None

    hdr = '%-9s %9s %9s %9s %9s %9s   %s' % ('config', 'wall[s]', 'step[s]', 'setup[s]',
                                             'MCells/s', 'RSS[MB]', 'note')
    print()
    print(hdr)
    print('-' * (len(hdr) + 6))
    for r in results:
        rss = '%9s' % '-' if r['peak_rss_mb'] is None else '%9.1f' % r['peak_rss_mb']
        print('%-9s %9.3f %9.3f %9.3f %9.1f %s   %s'
              % (r['name'], r['wall_s'], r['step_s'], r['setup_s'],
                 r['mcells_per_s'], rss, r['note']))
    speedup = None
    if best:
        print()
        print('fastest non-metal: %s (step %.3f s, wall %.3f s)'
              % (best['name'], best['step_s'], best['wall_s']))
        if metal:
            speedup = {
                'step': round(best['step_s'] / metal['step_s'], 3),
                'wall': round(best['wall_s'] / metal['wall_s'], 3),
                'baseline': best['name'],
            }
            print('Metal vs %s: stepping %.2fx, wall %.2fx'
                  % (best['name'], speedup['step'], speedup['wall']))

    if args.json:
        report = {
            'binary': str(binary),
            'cells': list(args.cells),
            'boundaries': boundaries,
            'steps': args.steps,
            'reps': args.reps,
            'results': results,
            'fastest_non_metal': best,
            'metal_speedup': speedup,
        }
        Path(args.json).write_text(json.dumps(report, indent=2) + '\n')
        print('wrote %s' % args.json)


if __name__ == '__main__':
    main()
