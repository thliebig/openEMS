#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright (C) 2026 Thorsten Liebig (Thorsten.Liebig@gmx.de)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
"""openEMS speed benchmark, see README.md.

    bench.py run [--openems PATH [--openems PATH2]] [options] [-o result.json]
    bench.py compare result.json[:LABEL] [other.json[:LABEL]]

`run` writes the cases of cases.py to XML and times the openEMS binary on
them; given two binaries it interleaves their runs. `compare` compares the
speed and the probe results of two binaries, from one or two result files.
"""

import os
import re
import sys
import json
import time
import zlib
import base64
import shutil
import argparse
import platform
import datetime
import statistics
import subprocess

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cases  # noqa: E402

FORMAT_VERSION = 1
ENGINES = ('multithreaded', 'gpu', 'basic', 'sse', 'sse-compressed')
SINGLE_THREAD_ENGINES = ('basic', 'sse', 'sse-compressed')
MISMATCH = 1e-4     # probe deviation relative to its peak that is reported as a mismatch
NOISE = 0.03        # a speed ratio within this (or the spread of the runs) is reported as noise


# --- machine ------------------------------------------------------------------

def _read(fn):
    try:
        with open(fn) as f:
            return f.read().strip()
    except OSError:
        return None


def _cpu_list(text):
    """ the CPUs of a Linux cpu list, e.g. '0-3,8-11' """
    cpus = set()
    for part in text.split(','):
        if '-' in part:
            a, b = part.split('-')
            cpus.update(range(int(a), int(b) + 1))
        elif part:
            cpus.add(int(part))
    return cpus


def _sysctl(name):
    try:
        return subprocess.run(['sysctl', '-n', name], capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def _powershell(cmd):
    try:
        return subprocess.run(['powershell', '-NoProfile', '-Command', cmd],
                              capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def physical_cores():
    """ The physical cores this process may run on, the performance cores only on a
        hybrid CPU (Intel P/E cores on Linux, Apple). None if unknown. """
    system = platform.system()
    if system == 'Linux':
        cpus = os.sched_getaffinity(0)
        p_cores = _read('/sys/devices/cpu_core/cpus')   # only present on a hybrid Intel CPU
        if p_cores:
            cpus &= _cpu_list(p_cores)
        cores = set()
        for c in cpus:
            topo = '/sys/devices/system/cpu/cpu{}/topology/'.format(c)
            pkg, core = _read(topo + 'physical_package_id'), _read(topo + 'core_id')
            if pkg is None or core is None:
                return None
            cores.add((pkg, core))
        return len(cores) or None
    if system == 'Darwin':
        n = _sysctl('hw.perflevel0.physicalcpu') or _sysctl('hw.physicalcpu')
        return int(n) if n else None
    if system == 'Windows':
        try:
            import psutil
            return psutil.cpu_count(logical=False)
        except ImportError:
            n = _powershell('(Get-CimInstance Win32_Processor | Measure-Object NumberOfCores -Sum).Sum')
            return int(n) if n and n.isdigit() else None
    return None


def machine_info():
    system = platform.system()
    cpu = ram = None
    if system == 'Linux':
        m = re.search(r'^model name\s*:\s*(.*)$', _read('/proc/cpuinfo') or '', re.M)
        cpu = m.group(1) if m else None
        m = re.search(r'^MemTotal:\s*(\d+) kB', _read('/proc/meminfo') or '', re.M)
        ram = int(m.group(1)) * 1024 if m else None
    elif system == 'Darwin':
        cpu = _sysctl('machdep.cpu.brand_string')
        ram = int(_sysctl('hw.memsize') or 0) or None
    elif system == 'Windows':
        cpu = _powershell('(Get-CimInstance Win32_Processor).Name')
        n = _powershell('(Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory')
        ram = int(n) if n and n.isdigit() else None
    return {
        'cpu': cpu or platform.processor() or platform.machine(),
        'physical_cores': physical_cores(),
        'logical_cpus': os.cpu_count(),
        'ram_gib': round(ram / 2**30, 1) if ram else None,
        'os': platform.platform(),
        'python': platform.python_version(),
        'env': {k: v for k, v in os.environ.items() if k.startswith('OPENEMS_')},
    }


# --- openEMS ------------------------------------------------------------------

def binary_info(path):
    out = subprocess.run([path, '--help'], capture_output=True, text=True).stdout
    version = re.search(r'version (\S+)', out)
    csxcad = re.search(r'CSXCAD -- Version: (\S+)', out)
    return {'path': path,
            'version': version.group(1) if version else None,
            'csxcad': csxcad.group(1) if csxcad else None}


def gpu_devices(path, xml, rundir, gpu_device):
    """ The GPU devices the binary lists, empty if its GPU engine cannot run """
    cmd = [path, xml, '--engine=gpu', '--no-simulation', '--verbose=1', '--numThreads=1']
    if gpu_device is not None:
        cmd.append('--gpu-device={}'.format(gpu_device))
    out = subprocess.run(cmd, cwd=rundir, capture_output=True, text=True)
    text = out.stdout + out.stderr
    if 'Create FDTD operator (GPU)' not in text:
        return []
    m = re.search(r'^GPU devices:\n((?:  .*\n)*)', text, re.M)
    return [l.strip() for l in m.group(1).splitlines()] if m else ['unknown device']


def _encode(a):
    return {'shape': list(a.shape), 'data': base64.b64encode(zlib.compress(a.astype('<f8').tobytes())).decode()}


def _decode(d):
    return np.frombuffer(zlib.decompress(base64.b64decode(d['data'])), '<f8').reshape(d['shape'])


def read_probes(rundir):
    """ the probe files of a run: the files without an extension """
    probes = {}
    for f in sorted(os.listdir(rundir)):
        if '.' not in f and os.path.isfile(os.path.join(rundir, f)):
            probes[f] = _encode(np.loadtxt(os.path.join(rundir, f), comments='%', ndmin=2))
    return probes


def _reap(p, block):
    """ (exit code, peak RSS in MiB) of the finished process `p`, None if it still runs """
    if hasattr(os, 'wait4'):
        # the resource usage of this very child, RUSAGE_CHILDREN would be the maximum of all
        pid, status, ru = os.wait4(p.pid, 0 if block else os.WNOHANG)
        if pid == 0:
            return None
        p.returncode = os.WEXITSTATUS(status) if os.WIFEXITED(status) else -os.WTERMSIG(status)
        return p.returncode, ru.ru_maxrss / (2**20 if sys.platform == 'darwin' else 2**10)
    if block:
        p.wait()
    return None if p.poll() is None else (p.returncode, None)


def _remove(fn):
    try:
        os.remove(fn)
    except FileNotFoundError:
        pass


def run_once(path, xml, rundir, engine, threads, gpu_device, duration, log):
    """ Run openEMS and stop it `duration` seconds after its time-stepping started, by the
        file ABORT in its working directory: the same time for every engine and machine """
    cmd = [path, xml, '--engine={}'.format(engine), '--numThreads={}'.format(threads)]
    if engine == 'gpu' and gpu_device is not None:
        cmd.append('--gpu-device={}'.format(gpu_device))
    abort = os.path.join(rundir, 'ABORT')
    _remove(abort)   # a left-over one would stop the run at once
    done = None
    with open(log, 'w') as logf:
        t0 = time.perf_counter()
        p = subprocess.Popen(cmd, cwd=rundir, stdout=logf, stderr=subprocess.STDOUT)
        try:
            # the setup, until the time-stepping starts
            started = None
            while done is None and started is None:
                time.sleep(0.02)
                done = _reap(p, False)
                with open(log) as f:
                    if 'Running FDTD engine' in f.read():
                        started = time.perf_counter()
            while done is None and time.perf_counter() - started < duration:
                time.sleep(0.02)
                done = _reap(p, False)
            if done is None:
                with open(abort, 'w'):
                    pass
                done = _reap(p, True)
        except BaseException:
            p.kill()
            raise
        finally:
            _remove(abort)
        wall = time.perf_counter() - t0
    code, rss = done
    with open(log) as f:
        text = f.read()
    t = re.search(r'^Time for (\d+) iterations with (\S+) cells : (\S+) sec', text, re.M)
    s = re.search(r'^Speed: (\S+) MCells/s', text, re.M)
    if code != 0 or not t or not s:
        tail = '\n'.join(text.splitlines()[-15:])
        return {'error': 'exit code {}, see {}\n{}'.format(code, log, tail)}
    if 'Found file "ABORT"' not in text:
        return {'error': 'stopped on its own after {} timesteps, see {}'.format(t.group(1), log)}
    run_s = float(t.group(3))
    r = {'mcells_s': float(s.group(1)), 'timesteps': int(t.group(1)), 'cells': float(t.group(2)), 'run_s': run_s,
         'wall_s': wall, 'setup_s': wall - run_s, 'rss_mib': rss}
    if engine == 'gpu':
        backend = re.search(r'^Create FDTD engine \(GPU, backend: (.*)\)$', text, re.M)
        r['gpu_backend'] = backend.group(1) if backend else None
        r['gpu_on_device'] = '(host fallback)' not in text
    return r


def _stats(values):
    return {'median': statistics.median(values), 'min': min(values), 'max': max(values)}


# --- run ----------------------------------------------------------------------

def _list(text, choices, what):
    items = [x.strip() for x in text.split(',') if x.strip()]
    for x in items:
        if x not in choices:
            sys.exit('bench.py: unknown {} "{}", choose from {}'.format(what, x, ', '.join(choices)))
    return items


def cmd_run(args):
    names = _list(args.cases, list(cases.CASES), 'case') if args.cases else list(cases.CASES)
    sizes = _list(args.size, cases.SIZES, 'size') if args.size else list(cases.DEFAULT_SIZES)
    engines = _list(args.engine, ENGINES, 'engine') if args.engine else None

    paths = args.openems or [shutil.which('openEMS') or 'openEMS']
    paths = [os.path.abspath(shutil.which(p) or p) for p in paths]   # the runs change the directory
    labels = args.label or [chr(ord('A') + n) for n in range(len(paths))]
    if len(labels) != len(paths) or len(set(labels)) != len(labels):
        sys.exit('bench.py: give one distinct --label per --openems')

    machine = machine_info()
    if args.threads:
        threads = [int(t) for t in args.threads.split(',')]
    else:
        threads = [machine['physical_cores'] or os.cpu_count()]

    workdir = os.path.abspath(args.workdir or datetime.datetime.now().strftime('bench-%Y%m%d-%H%M%S'))
    os.makedirs(os.path.join(workdir, 'xml'), exist_ok=True)
    output = os.path.abspath(args.output or os.path.join(workdir, 'result.json'))

    # the XS ideal case: GPU check and warm-up
    probe_xml = os.path.join(workdir, 'xml', 'ideal_XS.xml')
    cases.write_case('ideal', 'XS', probe_xml)

    binaries = {}
    for label, path in zip(labels, paths):
        if not os.path.isfile(path):
            sys.exit('bench.py: no openEMS binary at {}'.format(path))
        info = binary_info(path)
        gpudir = os.path.join(workdir, label, 'gpu_check')
        os.makedirs(gpudir, exist_ok=True)
        info['gpu_devices'] = gpu_devices(path, probe_xml, gpudir, args.gpu_device)
        binaries[label] = info
        print('{}: {} (openEMS {}, CSXCAD {}){}'.format(
            label, path, info['version'], info['csxcad'],
            ', GPU: ' + '; '.join(info['gpu_devices']) if info['gpu_devices'] else ''))
    print('machine: {}, {} physical cores / {} logical CPUs, {} GiB RAM, {}'.format(
        machine['cpu'], machine['physical_cores'], machine['logical_cpus'], machine['ram_gib'], machine['os']))
    if machine['env']:
        print('environment: ' + ' '.join('{}={}'.format(k, v) for k, v in machine['env'].items()))

    if engines is None:
        engines = ['multithreaded'] + (['gpu'] if any(b['gpu_devices'] for b in binaries.values()) else [])

    result = {'format': FORMAT_VERSION, 'cases_version': cases.CASES_VERSION,
              'date': datetime.datetime.now().isoformat(timespec='seconds'),
              'command': ' '.join(sys.argv), 'machine': machine, 'binaries': binaries,
              'repeat': args.repeat, 'duration': args.duration, 'points': []}

    # one discarded run per binary: loads the libraries and wakes the CPU up
    for label, info in binaries.items():
        rundir = os.path.join(workdir, label, 'warmup')
        os.makedirs(rundir, exist_ok=True)
        run_once(info['path'], probe_xml, rundir, 'multithreaded', threads[0], None, min(1.0, args.duration),
                 os.path.join(rundir, 'openEMS.log'))

    # all points first, for the progress count
    todo = []
    for name in names:
        for size in sizes:
            for engine in engines:
                if name == 'cylindrical' and engine in cases.CARTESIAN_ONLY_ENGINES:
                    continue    # same as multithreaded on a cylindrical mesh
                for nt in ([1] if engine in SINGLE_THREAD_ENGINES else threads):
                    todo.append((name, size, engine, nt))
    if any(e in cases.CARTESIAN_ONLY_ENGINES for e in engines) and 'cylindrical' in names:
        print('note: cylindrical runs multithreaded only, the other CPU engines are the same there')
    if 'gpu' in engines:
        for label, info in binaries.items():
            if not info['gpu_devices']:
                print('note: {} has no usable GPU device, no gpu runs for it'.format(label))

    xml_done = set()
    nruns = sum(args.repeat * sum(1 for l in binaries if e != 'gpu' or binaries[l]['gpu_devices'])
                for _, _, e, _ in todo)
    count = 0
    try:
        for name, size, engine, nt in todo:
            xml = os.path.join(workdir, 'xml', '{}_{}.xml'.format(name, size))
            if xml not in xml_done:
                cases.write_case(name, size, xml)
                xml_done.add(xml)
            run_labels = [l for l in binaries if engine != 'gpu' or binaries[l]['gpu_devices']]
            points = {}
            for label in run_labels:
                rundir = os.path.join(workdir, label, '{}_{}_{}_t{}'.format(name, size, engine, nt))
                os.makedirs(rundir, exist_ok=True)
                points[label] = {'case': name, 'size': size, 'engine': engine, 'threads': nt,
                                 'label': label, 'runs': [], 'error': None, 'dir': rundir}
            # interleaved: A B A B ...
            for rep in range(args.repeat):
                for label in run_labels:
                    pt = points[label]
                    if pt['error']:
                        continue
                    count += 1
                    r = run_once(binaries[label]['path'], xml, pt['dir'], engine, nt, args.gpu_device, args.duration,
                                 os.path.join(pt['dir'], 'openEMS_{}.log'.format(rep)))
                    tag = '[{:3d}/{}] {:14s} {:2s} {:14s} t{:<3d} {}'.format(
                        count, nruns, name, size, engine, nt, label)
                    if 'error' in r:
                        pt['error'] = r['error']
                        print(tag + ' FAILED: ' + r['error'], flush=True)
                        continue
                    print(tag + ' {:9.1f} MC/s  setup {:6.2f} s'.format(r['mcells_s'], r['setup_s']), flush=True)
                    if rep == 0:
                        pt['probes'] = read_probes(pt['dir'])
                    pt['runs'].append(r)
            for label in run_labels:
                pt = points[label]
                if pt['runs']:
                    pt['mcells_s'] = _stats([r['mcells_s'] for r in pt['runs']])
                    pt['setup_s'] = _stats([r['setup_s'] for r in pt['runs']])
                    rss = [r['rss_mib'] for r in pt['runs'] if r['rss_mib'] is not None]
                    pt['rss_mib'] = max(rss) if rss else None
                    if engine == 'gpu':
                        pt['gpu_backend'] = pt['runs'][0]['gpu_backend']
                        pt['gpu_on_device'] = all(r['gpu_on_device'] for r in pt['runs'])
                result['points'].append(pt)
    except KeyboardInterrupt:
        print('\ninterrupted, writing the results so far')
        result['interrupted'] = True

    with open(output, 'w') as f:
        json.dump(result, f, indent=1)

    print()
    for label in binaries:
        print_table(result, label)
    if len(binaries) > 1:
        a, b = list(binaries)[:2]
        print_compare((result, a), (result, b))
    print('results: ' + output)


# --- report -------------------------------------------------------------------

def _spread_text(*points):
    """ the spread of the speed of `points` (the largest), '-' if one of them is a single run """
    if any(len(pt['runs']) < 2 for pt in points):
        return '-'
    return '{:.1f}%'.format(100 * max(_spread(pt['mcells_s']) for pt in points))


def _spread(s):
    return (s['max'] - s['min']) / s['median'] if s['median'] else 0.0


def print_table(result, label):
    info = result['binaries'][label]
    print('{}: openEMS {} ({})'.format(label, info['version'], info['path']))
    print('{:14s} {:4s} {:14s} {:>4s} {:>7s} {:>9s} {:>9s} {:>6s} {:>8s} {:>8s}'.format(
        'case', 'size', 'engine', 'thr', 'MCells', 'timesteps', 'MC/s', 'spread', 'setup s', 'RSS MiB'))
    for pt in result['points']:
        if pt['label'] != label:
            continue
        head = '{:14s} {:4s} {:14s} {:4d}'.format(pt['case'], pt['size'], pt['engine'], pt['threads'])
        if not pt['runs']:
            print(head + '  FAILED')
            continue
        note = '  (extensions on the host)' if pt.get('gpu_on_device') is False else ''
        print(head + ' {:7.2f} {:9d} {:9.1f} {:>6s} {:8.2f} {:>8s}{}'.format(
            pt['runs'][0]['cells'] / 1e6, int(statistics.median(r['timesteps'] for r in pt['runs'])),
            pt['mcells_s']['median'], _spread_text(pt),
            pt['setup_s']['median'], '{:.0f}'.format(pt['rss_mib']) if pt['rss_mib'] else '-', note))
    print()


def probe_deviation(pa, pb):
    """ max. deviation of the probes of b from a, relative to each probe's peak, over the
        samples both have: the runs stop after a time, not after a number of timesteps """
    from openEMS.compare import deviation
    dev = 0.0
    for name in set(pa) | set(pb):
        if name not in pa or name not in pb:
            return np.inf
        a, b = _decode(pa[name]), _decode(pb[name])
        n = min(len(a), len(b))
        a, b = a[:n], b[:n]
        if a.shape != b.shape or not np.allclose(a[:, 0], b[:, 0], rtol=1e-9, atol=0):
            return np.inf   # the time axis has to agree
        dev = max(dev, deviation(a[:, 1:], b[:, 1:]))
    return dev


def _machine_key(r):
    m = r['machine']
    return (m['cpu'], m['physical_cores'], m['logical_cpus'], m['ram_gib'], m['os'])


def print_compare(side_a, side_b):
    """ compare (result, label) b with a; returns the number of probe mismatches """
    (ra, la), (rb, lb) = side_a, side_b
    print('compare {} (A) with {} (B)'.format(la, lb))
    for tag, r, l in (('A', ra, la), ('B', rb, lb)):
        m, b = r['machine'], r['binaries'][l]
        print('  {}: openEMS {}, {}, {} cores, {}, {}'.format(tag, b['version'], m['cpu'],
                                                            m['physical_cores'], m['os'], r['date']))
    if _machine_key(ra) != _machine_key(rb) or ra['machine']['env'] != rb['machine']['env']:
        print('  different machines (or environment): the ratio is the hardware, not only the code')
    elif ra is not rb:
        print('  separate runs: interleaved runs (two --openems in one run) are more reliable')
    if min(ra['repeat'], rb['repeat']) < 2:
        print('  single runs: no spread, a change within {:.0f}% is noise; --repeat 3 for a reliable comparison'.format(100 * NOISE))
    probes = ra['cases_version'] == rb['cases_version']
    if not probes:
        print('  different case versions ({} / {}): no probe comparison'.format(ra['cases_version'], rb['cases_version']))

    def points(r, l):
        d = {}
        for pt in r['points']:
            if pt['label'] == l and pt['runs']:
                d.setdefault((pt['case'], pt['size'], pt['engine']), {})[pt['threads']] = pt
        return d
    pa, pb = points(ra, la), points(rb, lb)

    print('{:14s} {:4s} {:14s} {:>7s} {:>9s} {:>9s} {:>6s} {:>7s} {:>8s} {:>9s}  {}'.format(
        'case', 'size', 'engine', 'thr A/B', 'A MC/s', 'B MC/s', 'B/A', 'spread', 'change', 'setup B/A', 'probes'))
    mismatches = 0
    for key in [k for k in pa if k in pb]:
        ta, tb = pa[key], pb[key]
        if len(ta) == 1 and len(tb) == 1:
            pairs = [(list(ta.values())[0], list(tb.values())[0])]
        else:
            pairs = [(ta[t], tb[t]) for t in sorted(set(ta) & set(tb))]
        for a, b in pairs:
            ratio = b['mcells_s']['median'] / a['mcells_s']['median']
            spread = max(_spread(a['mcells_s']), _spread(b['mcells_s']))   # 0 for single runs
            change = 'noise' if abs(ratio - 1) <= max(spread, NOISE) else '{:+.0f}%'.format(100 * (ratio - 1))
            setup = b['setup_s']['median'] / a['setup_s']['median'] if a['setup_s']['median'] > 0 else np.nan
            if not probes or 'probes' not in a or 'probes' not in b:
                pr = '-'
            else:
                dev = probe_deviation(a['probes'], b['probes'])
                pr = 'identical' if dev == 0 else '{:.1e}'.format(dev)
                if dev > MISMATCH:
                    pr += ' MISMATCH'
                    mismatches += 1
            print('{:14s} {:4s} {:14s} {:>7s} {:9.1f} {:9.1f} {:6.3f} {:>7s} {:>8s} {:9.2f}  {}'.format(
                key[0], key[1], key[2], '{}/{}'.format(a['threads'], b['threads']),
                a['mcells_s']['median'], b['mcells_s']['median'], ratio, _spread_text(a, b), change, setup, pr))
    only = sorted(set(pa) ^ set(pb))
    if only:
        print('  only in one of them: ' + ', '.join(' '.join(k) for k in only))
    print()
    return mismatches


def _load(spec, which):
    """ result file and label from 'file.json[:LABEL]'; `which` picks a label if none is given """
    fn, label = spec, None
    if not os.path.isfile(spec) and ':' in spec:
        fn, label = spec.rsplit(':', 1)
    with open(fn) as f:
        r = json.load(f)
    labels = list(r['binaries'])
    if label is None:
        if which >= len(labels):
            sys.exit('bench.py: {} holds {} only, give a second result file'.format(fn, ', '.join(labels)))
        label = labels[which]
    elif label not in labels:
        sys.exit('bench.py: no binary {} in {}, it holds {}'.format(label, fn, ', '.join(labels)))
    return r, label


def cmd_compare(args):
    if args.other is None:
        a, b = _load(args.base, 0), _load(args.base, 1)
    else:
        a, b = _load(args.base, 0), _load(args.other, 0)
    return 1 if print_compare(a, b) else 0


def main(argv=None):
    parser = argparse.ArgumentParser(prog='bench.py', description='openEMS speed benchmark')
    sub = parser.add_subparsers(dest='command', required=True)

    p = sub.add_parser('run', help='run the benchmark')
    p.add_argument('--openems', action='append', metavar='PATH',
                   help='openEMS binary, twice for an interleaved A/B comparison (default: openEMS in PATH)')
    p.add_argument('--label', action='append', help='name of each binary (default: A, B)')
    p.add_argument('--cases', help='comma separated, default: all of ' + ', '.join(cases.CASES))
    p.add_argument('--size', help='comma separated of {}, default: {}'.format(
        ', '.join(cases.SIZES), ','.join(cases.DEFAULT_SIZES)))
    p.add_argument('--engine', help='comma separated of {}, default: multithreaded, and gpu with a device'.format(
        ', '.join(ENGINES)))
    p.add_argument('--threads', help='comma separated thread counts, default: the physical cores')
    p.add_argument('--gpu-device', type=int, help='GPU device of the gpu engine (default: the best)')
    p.add_argument('--repeat', type=int,
                   help='timed runs per point (default: 1, and 3 with two binaries, a comparison needs the spread)')
    p.add_argument('--duration', type=float, default=3.0,
                   help='seconds of time-stepping per run, the setup comes on top (default: 3)')
    p.add_argument('--workdir', help='for the XML, the logs and the probes (default: ./bench-<date>-<time>)')
    p.add_argument('-o', '--output', help='result file (default: <workdir>/result.json)')

    p = sub.add_parser('compare', help='compare two binaries')
    p.add_argument('base', metavar='BASE.json[:LABEL]',
                   help='reference result; alone: compares its first two binaries')
    p.add_argument('other', metavar='OTHER.json[:LABEL]', nargs='?', help='result compared with the reference')

    args = parser.parse_args(argv)
    if args.command == 'run':
        if args.repeat is None:
            args.repeat = 3 if args.openems and len(args.openems) > 1 else 1
        if args.repeat < 1:
            parser.error('--repeat must be at least 1')
        return cmd_run(args)
    return cmd_compare(args)


if __name__ == '__main__':
    sys.exit(main())
