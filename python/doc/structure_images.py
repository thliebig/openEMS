# -*- coding: utf-8 -*-
"""
Render the structure of every Python tutorial, without simulating.

Each script runs in its own process with SHOW_STRUCTURE=1 and
CSXCAD_VIEWER_OUTPUT=<out>, so its Show() call writes <name>_view.png,
<name>_views.png and <name>_describe.json into <out> instead of opening a
window. openEMS.Run is replaced by an exit, so no simulation starts.

This generates the structure images for the documentation and checks the
structure viewer on real structures; --examples adds the (undocumented)
examples to the check.

usage: python structure_images.py [-o OUT_DIR] [--timeout SEC] [--examples] [script.py ...]
"""

import argparse
import glob
import json
import os
import subprocess
import sys

PY_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def run_one(script):
    """ Run `script` (child process) with openEMS.Run stopping it. """
    import runpy
    import openEMS
    from openEMS import openEMS as _openEMS

    class _NoRun(_openEMS):
        def Run(self, *args, **kw):
            sys.exit(0)
    openEMS.openEMS = _NoRun
    sys.argv = [script]
    sys.path.insert(0, os.path.dirname(script))   # as "python script.py" does
    runpy.run_path(script, run_name='__main__')


def main():
    p = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    p.add_argument('scripts', nargs='*', help='scripts to run (default: all tutorials)')
    p.add_argument('--examples', action='store_true', help='also run the examples')
    p.add_argument('-o', '--out', default='structure_images', help='output directory')
    p.add_argument('--timeout', type=float, default=600, help='seconds per script')
    p.add_argument('--run-one', help=argparse.SUPPRESS)
    args = p.parse_args()
    if args.run_one:
        return run_one(args.run_one)

    dirs = ['Tutorials'] + (['Examples'] if args.examples else [])
    scripts = args.scripts or sorted(f for d in dirs for f in glob.glob(os.path.join(PY_DIR, d, '*.py')))
    out = os.path.abspath(args.out)
    os.makedirs(out, exist_ok=True)
    env = dict(os.environ, SHOW_STRUCTURE='1', CSXCAD_VIEWER_OUTPUT=out, MPLBACKEND='Agg')
    failed = []
    for script in scripts:
        name = os.path.splitext(os.path.basename(script))[0]
        for old in glob.glob(os.path.join(out, name + '_*')):
            os.remove(old)
        try:
            r = subprocess.run([sys.executable, os.path.abspath(__file__), '--run-one', os.path.abspath(script)],
                               env=env, cwd=os.path.dirname(os.path.abspath(script)), capture_output=True,
                               text=True, timeout=args.timeout)
            code, log = r.returncode, r.stdout + r.stderr
        except subprocess.TimeoutExpired:
            code, log = 'timeout', ''
        problems = check(out, name)
        if code != 0:
            problems.insert(0, 'exit {}: {}'.format(code, log.strip().splitlines()[-1] if log.strip() else ''))
        print('{:<30} {}'.format(name, 'ok' if not problems else '; '.join(problems)))
        if problems:
            failed.append(name)
    print('{} of {} scripts ok'.format(len(scripts) - len(failed), len(scripts)))
    return 1 if failed else 0


def check(out, name):
    """ Problems with the outputs of one script (empty list if fine). """
    files = sorted(glob.glob(os.path.join(out, name + '*_describe.json')))
    if not files:
        return ['no structure shown']
    problems = []
    for fn in files:
        for suffix in ('_view.png', '_views.png'):
            if not os.path.isfile(fn.replace('_describe.json', suffix)):
                problems.append('missing ' + os.path.basename(fn.replace('_describe.json', suffix)))
        with open(fn) as f:
            d = json.load(f)
        for r in d['unrendered']:
            if r['reason'] != 'unsupported':
                problems.append('unrendered {prop}: {type} #{primitive_id} ({reason})'.format(**r))
        for prop in d['properties']:
            if prop['shown'] and prop['bbox'] is None:
                problems.append('visible property without geometry: ' + prop['name'])
        if not d['grid']:
            problems.append('no grid')
    return problems


if __name__ == '__main__':
    sys.exit(main())
