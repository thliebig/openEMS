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
    bench.py verify [--openems PATH [--openems PATH2]] [options]

`run` writes the cases of cases.py to XML and times the openEMS binary on
them; given two binaries it interleaves their runs. `compare` compares the
speed and the probe results of two binaries, from one or two result files.
`verify` runs the cases for a fixed number of timesteps on every binary and
engine and compares all their probes and field dumps.
"""

import os
import re
import sys
import json
import time
import zlib
import atexit
import base64
import shutil
import argparse
import platform
import datetime
import tempfile
import statistics
import subprocess

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cases  # noqa: E402

FORMAT_VERSION = 1
ENGINES = ('multithreaded', 'gpu', 'basic', 'cpu')
MISMATCH = 1e-4     # probe or dump deviation relative to its peak that is reported as a mismatch
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


# --- scratch space ------------------------------------------------------------

def _fs_type(path):
    """ the file system type of `path`, None if unknown (only known on Linux) """
    path = os.path.realpath(path)
    mount, fstype = '', None
    for line in (_read('/proc/mounts') or '').splitlines():
        f = line.split()
        if len(f) > 2 and (path == f[1] or path.startswith(f[1].rstrip('/') + '/')) and len(f[1]) > len(mount):
            mount, fstype = f[1], f[2]
    return fstype


def _mem_available():
    m = re.search(r'^MemAvailable:\s*(\d+) kB', _read('/proc/meminfo') or '', re.M)
    return int(m.group(1)) * 1024 if m else None


def _files(path, dumps):
    """ the files of a run: its field dumps (with an extension), or its probes (without) """
    return [e for e in os.scandir(path) if e.is_file() and ('.' in e.name) == dumps]


def _bytes(path):
    return sum(e.stat().st_size for e in os.scandir(path) if e.is_file())


def default_scratch():
    """ /dev/shm, a RAM disk, where there is one """
    return '/dev/shm' if os.path.isdir('/dev/shm') and os.access('/dev/shm', os.W_OK) else None


class Scratch:
    """ Where openEMS runs, with its field dumps: a directory per run, removed as soon as
        the run is done, under `base`. A run is stopped early when its files grow beyond
        `limit_gib` (default: a quarter of the free space, on a tmpfs of the free RAM) or
        when the file system, on a tmpfs the RAM, gets down to RESERVE. """
    RESERVE = 2**30

    def __init__(self, base, limit_gib=None):
        self.ram = _fs_type(base) == 'tmpfs'
        self.dir = tempfile.mkdtemp(prefix='openems-bench-', dir=base)
        atexit.register(shutil.rmtree, self.dir, ignore_errors=True)
        self.limit = limit_gib * 2**30 if limit_gib else self._available() // 4

    def _available(self):
        free = shutil.disk_usage(self.dir).free
        mem = _mem_available() if self.ram else None
        return min(free, mem) if mem is not None else free

    def describe(self):
        return '{}{}, at most {:.3g} GiB per run'.format(self.dir, ' (RAM)' if self.ram else '', self.limit / 2**30)

    def rundir(self, name):
        path = os.path.join(self.dir, name)
        os.makedirs(path)
        return path

    def guard(self, rundir):
        """ for run_once: the reason to stop the run in `rundir`, None to go on """
        last = [0.0]

        def check():
            now = time.perf_counter()
            if now - last[0] < 0.1:
                return None
            last[0] = now
            if _bytes(rundir) > self.limit:
                return 'its files reached {:.3g} GiB'.format(self.limit / 2**30)
            if self._available() < self.RESERVE:
                return 'low on ' + ('RAM' if self.ram else 'disk space')
            return None
        return check

    def keep_probes(self, rundir, dest):
        """ copy the probes of a run to `dest` """
        for e in _files(rundir, False):
            shutil.copyfile(e.path, os.path.join(dest, e.name))

    def remove(self, rundir):
        shutil.rmtree(rundir, ignore_errors=True)


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
    # the operator line ends in the engine type on every mesh; it is printed only once a
    # backend with a usable device exists (--no-simulation creates no engine)
    if not re.search(r'^Create .*FDTD operator.* \(GPU\)$', text, re.M):
        return []
    m = re.search(r'^GPU devices:\n((?:  .*\n)*)', text, re.M)
    return [l.strip() for l in m.group(1).splitlines()] if m else ['unknown device']


def _encode(a):
    return {'shape': list(a.shape), 'data': base64.b64encode(zlib.compress(a.astype('<f8').tobytes())).decode()}


def _decode(d):
    return np.frombuffer(zlib.decompress(base64.b64decode(d['data'])), '<f8').reshape(d['shape'])


def read_probes(rundir):
    """ the probe files of a run: the files without an extension, except et and ht,
        which hold the excitation signal """
    probes = {}
    for f in sorted(os.listdir(rundir)):
        if '.' not in f and f not in ('et', 'ht') and os.path.isfile(os.path.join(rundir, f)):
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


# what the CPU engine chose, from its log (-v)
ENGINE_INFO = (
    ('cache', r'^Engine_CPU: cache (\S+ MiB, .*)$'),
    ('kernels', r'^Engine_CPU: kernels for (\S+)$'),
    ('layout', r'^Engine_CPU: fast axis (\S+), layout (\S+)$'),
    ('blocking', r'^Engine_CPU: temporal blocking, up to (\d+) timesteps per block, tiles of (\d+) lines'),
    ('flat', r'^Engine_CPU: no temporal blocking: (.*)$'),
    ('trials', r'^Multithreaded engine: tried (.*); running (.*)$'),
)


def engine_info(text):
    """ the settings the CPU engine chose for a run (empty for the other engines) """
    info = {}
    for key, rx in ENGINE_INFO:
        m = re.search(rx, text, re.M)
        if m:
            info[key] = m.group(1) if len(m.groups()) == 1 else list(m.groups())
    return info


def engine_note(info):
    """ a short line of engine_info() for the tables """
    if not info:
        return ''
    parts = [info.get('kernels', '?')]
    if 'layout' in info:
        parts.append('fast ' + info['layout'][0])
    if 'blocking' in info:
        parts.append('k{} W{}'.format(*info['blocking']))
    elif 'flat' in info:
        parts.append('flat (' + info['flat'] + ')')
    if 'trials' in info:
        parts.append('trials kept ' + info['trials'][1])
    if 'cache' in info:
        parts.append('cache ' + info['cache'].split(',')[0])
    return ', '.join(parts)


def _option(arg):
    """ the name of an openEMS option; -v, -vv and -vvv are --verbose """
    return '--verbose' if re.match(r'-v+$', arg) else arg.split('=')[0]


def run_once(path, xml, rundir, engine, threads, gpu_device, duration, log, env=None, guard=None, args=None):
    """ Run openEMS and stop it `duration` seconds after its time-stepping started, by the
        file ABORT in its working directory: the same time for every engine and machine.
        With `duration` None it runs to its end. `env`: variables added to the environment,
        `args`: openEMS options (a variant); an option given here replaces the same one of
        bench.py, openEMS takes an option only once. `guard` is called during the
        time-stepping and stops the run early when it returns a reason, see Scratch. """
    own = ['--engine={}'.format(engine), '--numThreads={}'.format(threads), '-v']
    if engine == 'gpu' and gpu_device is not None:
        own.append('--gpu-device={}'.format(gpu_device))
    args = list(args or [])
    given = {_option(a) for a in args}
    cmd = [path, xml] + [a for a in own if _option(a) not in given] + args
    abort = os.path.join(rundir, 'ABORT')
    _remove(abort)   # a left-over one would stop the run at once
    done = stopped = None
    with open(log, 'w') as logf:
        t0 = time.perf_counter()
        run_env = dict(os.environ, **(env or {}))
        # the CPU engine tries its tiles during its first seconds (OPENEMS_CPU_TRIAL), a start-up cost
        # like the setup: off unless asked for, so a run of a few seconds times the steady state
        if '--engine=cpu' in cmd and 'OPENEMS_CPU_TRIAL' not in run_env:
            run_env['OPENEMS_CPU_TRIAL'] = '0'
        p = subprocess.Popen(cmd, cwd=rundir, stdout=logf, stderr=subprocess.STDOUT, env=run_env)
        try:
            # the setup, until the time-stepping starts
            started = None
            while done is None and started is None:
                time.sleep(0.02)
                done = _reap(p, False)
                with open(log) as f:
                    if 'Running FDTD engine' in f.read():
                        started = time.perf_counter()
            while done is None and (duration is None or time.perf_counter() - started < duration):
                time.sleep(0.02)
                done = _reap(p, False)
                stopped = guard() if guard and done is None else None
                if stopped:
                    break
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
    if stopped and duration is None:
        return {'error': 'stopped after {} timesteps, {}'.format(t.group(1), stopped)}
    if duration is not None and 'Found file "ABORT"' not in text:
        return {'error': 'stopped on its own after {} timesteps, see {}'.format(t.group(1), log)}
    run_s = float(t.group(3))
    dump_bytes = sum(e.stat().st_size for e in _files(rundir, True))
    r = {'mcells_s': float(s.group(1)), 'timesteps': int(t.group(1)), 'cells': float(t.group(2)), 'run_s': run_s,
         'wall_s': wall, 'setup_s': wall - run_s, 'rss_mib': rss, 'engine_info': engine_info(text),
         'dump_mb_s': dump_bytes / 1e6 / run_s if run_s > 0 else None, 'stopped': stopped}
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


def _variants(specs):
    """ [(name, {VAR: value}, [openEMS options])] from --variant NAME:VAR=VALUE,--OPTION,...;
        [(None, {}, [])] without """
    if not specs:
        return [(None, {}, [])]
    out = []
    for spec in specs:
        name, _, rest = spec.partition(':')
        env, opts = {}, []
        for item in [x for x in rest.split(',') if x]:
            if item.startswith('-'):
                opts.append(item)
                continue
            var, eq, value = item.partition('=')
            if not eq or not var:
                sys.exit('bench.py: --variant {}: give VAR=VALUE pairs or openEMS options after the name'.format(spec))
            env[var] = value
        if not name or name in [n for n, _, _ in out]:
            sys.exit('bench.py: --variant {}: give each variant a distinct name'.format(spec))
        out.append((name, env, opts))
    return out


def _variant_text(info):
    """ the environment and the openEMS options of a variant, for the headers """
    return ' '.join(['{}={}'.format(k, v) for k, v in (info.get('env') or {}).items()] + (info.get('args') or []))


def _dumps(text):
    if text == 'all':
        return list(cases.DUMPS)
    return _list(text, list(cases.DUMPS), 'dump')


def _point_name(name, size, dumps):
    return '{}_{}'.format(name, size) + ('' if dumps == 'none' else '_' + dumps)


def setup(args):
    """ The contestants (every binary with every variant), the machine, the thread counts,
        the workdir and the scratch space of a run or verify """
    paths = args.openems or [shutil.which('openEMS') or 'openEMS']
    paths = [os.path.abspath(shutil.which(p) or p) for p in paths]   # the runs change the directory
    bin_labels = args.label or [chr(ord('A') + n) for n in range(len(paths))]
    if len(bin_labels) != len(paths) or len(set(bin_labels)) != len(bin_labels):
        sys.exit('bench.py: give one distinct --label per --openems')
    # the contestants: every binary with every variant, named by the variant with one binary
    variants = _variants(args.variant)
    contestants = []   # (label, binary label, path, variant, env, openEMS options)
    for bl, path in zip(bin_labels, paths):
        for vname, env, opts in variants:
            if vname is None:
                label = bl
            else:
                label = vname if len(paths) == 1 else '{}-{}'.format(bl, vname)
            contestants.append((label, bl, path, vname, env, opts))

    machine = machine_info()
    if args.threads:
        threads = [int(t) for t in args.threads.split(',')]
    else:
        threads = [machine['physical_cores'] or os.cpu_count()]

    workdir = os.path.abspath(args.workdir or datetime.datetime.now().strftime(
        '{}-%Y%m%d-%H%M%S'.format('bench' if args.command == 'run' else 'verify')))
    os.makedirs(os.path.join(workdir, 'xml'), exist_ok=True)
    scratch = Scratch(args.scratch or default_scratch() or workdir, args.scratch_limit)

    # the XS ideal case: GPU check and warm-up
    probe_xml = os.path.join(workdir, 'xml', 'ideal_XS.xml')
    cases.write_case('ideal', 'XS', probe_xml)

    binaries = {}
    checked = {}   # per binary
    for label, bl, path, vname, env, opts in contestants:
        if not os.path.isfile(path):
            sys.exit('bench.py: no openEMS binary at {}'.format(path))
        if path not in checked:
            info = binary_info(path)
            gpudir = os.path.join(workdir, bl, 'gpu_check')
            os.makedirs(gpudir, exist_ok=True)
            info['gpu_devices'] = gpu_devices(path, probe_xml, gpudir, args.gpu_device)
            checked[path] = info
            print('{}: {} (openEMS {}, CSXCAD {}){}'.format(
                bl, path, info['version'], info['csxcad'],
                ', GPU: ' + '; '.join(info['gpu_devices']) if info['gpu_devices'] else ''))
        binaries[label] = dict(checked[path], variant=vname, env=env, args=opts)
        if vname is not None:
            print('  variant {}: {}'.format(label, _variant_text(binaries[label]) or 'as is'))
    print('machine: {}, {} physical cores / {} logical CPUs, {} GiB RAM, {}'.format(
        machine['cpu'], machine['physical_cores'], machine['logical_cpus'], machine['ram_gib'], machine['os']))
    if machine['env']:
        print('environment: ' + ' '.join('{}={}'.format(k, v) for k, v in machine['env'].items()))
    print('scratch: ' + scratch.describe())
    return binaries, machine, threads, workdir, scratch, probe_xml


def cmd_run(args):
    names = _list(args.cases, list(cases.CASES), 'case') if args.cases else list(cases.CASES)
    sizes = _list(args.size, cases.SIZES, 'size') if args.size else list(cases.DEFAULT_SIZES)
    dumps = _dumps(args.dumps) if args.dumps else ['none']
    engines = _list(args.engine, ENGINES, 'engine') if args.engine else None

    binaries, machine, threads, workdir, scratch, probe_xml = setup(args)
    output = os.path.abspath(args.output or os.path.join(workdir, 'result.json'))
    if engines is None:
        engines = ['multithreaded'] + (['gpu'] if any(b['gpu_devices'] for b in binaries.values()) else [])

    result = {'format': FORMAT_VERSION, 'cases_version': cases.CASES_VERSION,
              'date': datetime.datetime.now().isoformat(timespec='seconds'),
              'command': ' '.join(sys.argv), 'machine': machine, 'binaries': binaries,
              'repeat': args.repeat, 'duration': args.duration,
              'scratch': {'ram': scratch.ram, 'limit_gib': scratch.limit / 2**30}, 'points': []}

    # one discarded run per binary: loads the libraries and wakes the CPU up
    for label, info in binaries.items():
        logdir = os.path.join(workdir, label, 'warmup')
        os.makedirs(logdir, exist_ok=True)
        rundir = scratch.rundir('warmup_' + label)
        run_once(info['path'], probe_xml, rundir, 'multithreaded', threads[0], None, min(1.0, args.duration),
                 os.path.join(logdir, 'openEMS.log'), info['env'], args=info['args'])
        scratch.remove(rundir)

    # all points first, for the progress count
    todo = []
    for name in names:
        for size in sizes:
            for dump in dumps:
                for engine in engines:
                    for nt in threads:
                        todo.append((name, size, dump, engine, nt))
    if 'gpu' in engines:
        for label, info in binaries.items():
            if not info['gpu_devices']:
                print('note: {} has no usable GPU device, no gpu runs for it'.format(label))

    xml_done = set()
    nruns = sum(args.repeat * sum(1 for l in binaries if e != 'gpu' or binaries[l]['gpu_devices'])
                for _, _, _, e, _ in todo)
    count = 0
    try:
        for name, size, dump, engine, nt in todo:
            xml = os.path.join(workdir, 'xml', _point_name(name, size, dump) + '.xml')
            if xml not in xml_done:
                cases.write_case(name, size, xml, dump)
                xml_done.add(xml)
            run_labels = [l for l in binaries if engine != 'gpu' or binaries[l]['gpu_devices']]
            points = {}
            for label in run_labels:
                pointdir = os.path.join(workdir, label, '{}_{}_t{}'.format(_point_name(name, size, dump), engine, nt))
                os.makedirs(pointdir, exist_ok=True)
                points[label] = {'case': name, 'size': size, 'dumps': dump, 'engine': engine, 'threads': nt,
                                 'label': label, 'runs': [], 'error': None, 'dir': pointdir}
            # interleaved: A B A B ...
            for rep in range(args.repeat):
                for label in run_labels:
                    pt = points[label]
                    if pt['error']:
                        continue
                    count += 1
                    rundir = scratch.rundir('{}_{}'.format(os.path.basename(pt['dir']), label))
                    try:
                        info = binaries[label]
                        r = run_once(info['path'], xml, rundir, engine, nt, args.gpu_device, args.duration,
                                     os.path.join(pt['dir'], 'openEMS_{}.log'.format(rep)), info['env'],
                                     scratch.guard(rundir), info['args'])
                        if 'error' not in r:
                            if rep == 0:
                                pt['probes'] = read_probes(rundir)
                            scratch.keep_probes(rundir, pt['dir'])
                    finally:
                        scratch.remove(rundir)
                    tag = '[{:3d}/{}] {:14s} {:2s} {}{:14s} t{:<3d} {}'.format(
                        count, nruns, name, size, '' if dumps == ['none'] else '{:7s} '.format(dump), engine, nt, label)
                    if 'error' in r:
                        pt['error'] = r['error']
                        print(tag + ' FAILED: ' + r['error'], flush=True)
                        continue
                    print(tag + ' {:9.1f} MC/s  setup {:6.2f} s{}{}'.format(
                        r['mcells_s'], r['setup_s'],
                        '  dumps {:7.1f} MB/s'.format(r['dump_mb_s']) if dump != 'none' else '',
                        '  stopped after {:.1f} s: {}'.format(r['run_s'], r['stopped']) if r['stopped'] else ''),
                        flush=True)
                    pt['runs'].append(r)
            for label in run_labels:
                pt = points[label]
                if pt['runs']:
                    pt['mcells_s'] = _stats([r['mcells_s'] for r in pt['runs']])
                    pt['setup_s'] = _stats([r['setup_s'] for r in pt['runs']])
                    pt['dump_mb_s'] = statistics.median(r['dump_mb_s'] or 0 for r in pt['runs'])
                    rss = [r['rss_mib'] for r in pt['runs'] if r['rss_mib'] is not None]
                    pt['rss_mib'] = max(rss) if rss else None
                    pt['engine_info'] = pt['runs'][0].get('engine_info', {})
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
    # each against the first: the reference binary, or the first variant
    labels = list(binaries)
    for b in labels[1:]:
        print_compare((result, labels[0]), (result, b))
    print('results: ' + output)


# --- report -------------------------------------------------------------------

def _spread_text(*points):
    """ the spread of the speed of `points` (the largest), '-' if one of them is a single run """
    if any(len(pt['runs']) < 2 for pt in points):
        return '-'
    return '{:.1f}%'.format(100 * max(_spread(pt['mcells_s']) for pt in points))


def _spread(s):
    return (s['max'] - s['min']) / s['median'] if s['median'] else 0.0


def _with_dumps(*results):
    """ whether the tables have the dumps columns: only for a run with --dumps """
    return any(pt.get('dumps', 'none') != 'none' for r in results for pt in r['points'])


def print_table(result, label):
    info = result['binaries'][label]
    variant = _variant_text(info)
    print('{}: openEMS {} ({}){}'.format(label, info['version'], info['path'], ', ' + variant if variant else ''))
    dumps = _with_dumps(result)
    print('{:14s} {:4s} {}{:14s} {:>4s} {:>7s} {:>9s} {:>9s} {:>6s} {:>8s} {:>8s}{}'.format(
        'case', 'size', 'dumps   ' if dumps else '', 'engine', 'thr', 'MCells', 'timesteps', 'MC/s', 'spread',
        'setup s', 'RSS MiB', ' dump MB/s' if dumps else ''))
    for pt in result['points']:
        if pt['label'] != label:
            continue
        dump = pt.get('dumps', 'none')
        head = '{:14s} {:4s} {}{:14s} {:4d}'.format(pt['case'], pt['size'], '{:7s} '.format(dump) if dumps else '',
                                                   pt['engine'], pt['threads'])
        if not pt['runs']:
            print(head + '  FAILED')
            continue
        notes = []
        if pt.get('gpu_on_device') is False:
            notes.append('(extensions on the host)')
        if engine_note(pt.get('engine_info')):
            notes.append(engine_note(pt.get('engine_info')))
        stopped = [r['stopped'] for r in pt['runs'] if r.get('stopped')]
        if stopped:
            notes.append('(stopped early: {})'.format(stopped[0]))
        rate = pt.get('dump_mb_s')
        print(head + ' {:7.2f} {:9d} {:9.1f} {:>6s} {:8.2f} {:>8s}{}{}'.format(
            pt['runs'][0]['cells'] / 1e6, int(statistics.median(r['timesteps'] for r in pt['runs'])),
            pt['mcells_s']['median'], _spread_text(pt),
            pt['setup_s']['median'], '{:.0f}'.format(pt['rss_mib']) if pt['rss_mib'] else '-',
            ' {:>9s}'.format('{:.1f}'.format(rate) if rate or dump != 'none' else '-') if dumps else '',
            ''.join('  ' + n for n in notes)))
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
        variant = _variant_text(b)
        print('  {}: openEMS {}, {}, {} cores, {}, {}{}'.format(tag, b['version'], m['cpu'],
                                                              m['physical_cores'], m['os'], r['date'],
                                                              ', ' + variant if variant else ''))
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
                d.setdefault((pt['case'], pt['size'], pt.get('dumps', 'none'), pt['engine']), {})[pt['threads']] = pt
        return d
    pa, pb = points(ra, la), points(rb, lb)

    dumps = _with_dumps(ra, rb)
    print('{:14s} {:4s} {}{:14s} {:>7s} {:>9s} {:>9s} {:>6s} {:>7s} {:>8s} {:>9s}  {}'.format(
        'case', 'size', 'dumps   ' if dumps else '', 'engine', 'thr A/B', 'A MC/s', 'B MC/s', 'B/A', 'spread', 'change',
        'setup B/A', 'probes'))
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
            print('{:14s} {:4s} {}{:14s} {:>7s} {:9.1f} {:9.1f} {:6.3f} {:>7s} {:>8s} {:9.2f}  {}'.format(
                key[0], key[1], '{:7s} '.format(key[2]) if dumps else '', key[3], '{}/{}'.format(a['threads'], b['threads']),
                a['mcells_s']['median'], b['mcells_s']['median'], ratio, _spread_text(a, b), change, setup, pr))
    only = sorted(set(pa) ^ set(pb))
    if only:
        print('  only in one of them: ' + ', '.join(' '.join(k[:2] + k[2 if dumps else 3:]) for k in only))
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


# --- verify -------------------------------------------------------------------

def cmd_verify(args):
    """ Run every case for a fixed number of timesteps on every binary and engine and
        compare all probes and field dumps with the run of the first binary's first engine """
    from openEMS.compare import compare_sim_results
    names = _list(args.cases, list(cases.CASES), 'case') if args.cases else list(cases.CASES)
    sizes = _list(args.size, cases.SIZES, 'size') if args.size else ['XS']
    dumps = _dumps(args.dumps) if args.dumps else list(cases.DUMPS)
    engines = _list(args.engine, ENGINES, 'engine') if args.engine else None

    binaries, machine, threads, workdir, scratch, _ = setup(args)
    if engines is None:
        engines = ['multithreaded'] + (['gpu'] if any(b['gpu_devices'] for b in binaries.values()) else [])
    runs = [(l, e) for l in binaries for e in engines if e != 'gpu' or binaries[l]['gpu_devices']]
    width = max(6, max(len(l) for l in binaries))
    print('reference: {} {}, {} timesteps, deviation relative to the peak, mismatch above {:g}'.format(
        runs[0][0], runs[0][1], args.timesteps, args.rtol))
    print()
    print('{:14s} {:4s} {:7s} {:{w}s} {:14s} {:>9s} {:>7s} {:>9s}  {}'.format(
        'case', 'size', 'dumps', 'binary', 'engine', 'timesteps', 'outputs', 'deviation', 'worst', w=width))

    rows = []
    for name in names:
        for size in sizes:
            for dump in dumps:
                xml = os.path.join(workdir, 'xml', _point_name(name, size, dump) + '.xml')
                cases.write_case(name, size, xml, dump, args.timesteps)
                ref = None
                try:
                    for n, (label, engine) in enumerate(runs):
                        info = binaries[label]
                        logdir = os.path.join(workdir, label, '{}_{}'.format(_point_name(name, size, dump), engine))
                        os.makedirs(logdir, exist_ok=True)
                        rundir = scratch.rundir('{}_{}'.format(os.path.basename(logdir), label))
                        row = {'case': name, 'size': size, 'dumps': dump, 'label': label, 'engine': engine}
                        try:
                            r = run_once(info['path'], xml, rundir, engine, threads[0], args.gpu_device, None,
                                         os.path.join(logdir, 'openEMS.log'), info['env'], scratch.guard(rundir),
                                         info['args'])
                            if 'error' in r:
                                row['error'] = r['error']
                            else:
                                row['timesteps'] = r['timesteps']
                                scratch.keep_probes(rundir, logdir)
                                if n == 0:
                                    ref, rundir = rundir, None   # kept until the case is done
                                elif ref is not None:
                                    devs = compare_sim_results(ref, rundir)
                                    row['outputs'] = len(devs)
                                    row['worst'], dev = max(devs, key=lambda d: d[1])
                                    row['deviation'] = float(dev)
                        finally:
                            if rundir:
                                scratch.remove(rundir)
                        rows.append(row)
                        head = '{:14s} {:4s} {:7s} {:{w}s} {:14s}'.format(name, size, dump, label, engine, w=width)
                        if 'error' in row:
                            lines = row['error'].splitlines()
                            print(head + ' FAILED: ' + lines[0] + (': ' + lines[-1].strip() if len(lines) > 1 else ''),
                                  flush=True)
                        elif n == 0:
                            print(head + ' {:9d} {:>7s} {:>9s}'.format(row['timesteps'], '-', 'reference'), flush=True)
                        elif 'deviation' not in row:
                            print(head + ' {:9d}  no reference'.format(row['timesteps']), flush=True)
                        else:
                            dev = row['deviation']
                            print(head + ' {:9d} {:7d} {:>9s}  {}{}'.format(
                                row['timesteps'], row['outputs'], 'identical' if dev == 0 else '{:.1e}'.format(dev),
                                row['worst'] if dev > 0 else '', '  MISMATCH' if dev > args.rtol else ''),
                                flush=True)
                finally:
                    if ref:
                        scratch.remove(ref)

    failed = sum(1 for r in rows if 'error' in r)
    mismatches = sum(1 for r in rows if r.get('deviation', 0) > args.rtol)
    output = os.path.join(workdir, 'verify.json')
    with open(output, 'w') as f:
        json.dump({'format': FORMAT_VERSION, 'cases_version': cases.CASES_VERSION,
                   'date': datetime.datetime.now().isoformat(timespec='seconds'), 'command': ' '.join(sys.argv),
                   'machine': machine, 'binaries': binaries, 'timesteps': args.timesteps, 'rtol': args.rtol,
                   'rows': rows}, f, indent=1)
    print()
    print('{} runs, {} failed, {} mismatches; results: {}'.format(len(rows), failed, mismatches, output))
    return 1 if failed or mismatches else 0


def main(argv=None):
    parser = argparse.ArgumentParser(prog='bench.py', description='openEMS speed benchmark')
    sub = parser.add_subparsers(dest='command', required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument('--openems', action='append', metavar='PATH',
                        help='openEMS binary, twice for an interleaved A/B comparison (default: openEMS in PATH)')
    common.add_argument('--label', action='append', help='name of each binary (default: A, B)')
    common.add_argument('--variant', action='append', metavar='NAME:VAR=VALUE|--OPTION[,...]',
                        help='run every point once more with these environment variables and openEMS options, '
                             'e.g. flat:OPENEMS_CPU_TBLOCK=0 or nodumps:--disable-dumps; an option replaces the '
                             'same one of bench.py (--engine, --numThreads, --gpu-device, -v); several, '
                             'interleaved like binaries; "NAME:" runs as is')
    common.add_argument('--cases', help='comma separated, default: all of ' + ', '.join(cases.CASES))
    common.add_argument('--engine', help='comma separated of {}, default: multithreaded, and gpu with a device'.format(
        ', '.join(ENGINES)))
    common.add_argument('--threads', help='comma separated thread counts, default: the physical cores')
    common.add_argument('--gpu-device', type=int, help='GPU device of the gpu engine (default: the best)')
    common.add_argument('--workdir', help='for the XML, the logs and the probes '
                                          '(default: ./bench-<date>-<time>, ./verify-<date>-<time>)')
    common.add_argument('--scratch', metavar='DIR',
                        help='where openEMS runs and writes its dumps, a directory per run that is removed right '
                             'after it (default: /dev/shm, a RAM disk, if there is one, else the workdir)')
    common.add_argument('--scratch-limit', type=float, metavar='GIB',
                        help='stop a run when its files reach this size (default: a quarter of the free space, '
                             'on a RAM disk of the free RAM)')

    p = sub.add_parser('run', parents=[common], help='run the benchmark')
    p.add_argument('--size', help='comma separated of {}, default: {}'.format(
        ', '.join(cases.SIZES), ','.join(cases.DEFAULT_SIZES)))
    p.add_argument('--dumps', help='comma separated of {}, or all, default: none'.format(', '.join(cases.DUMPS)))
    p.add_argument('--repeat', type=int,
                   help='timed runs per point (default: 1, and 3 with two binaries or variants, a comparison needs the spread)')
    p.add_argument('--duration', type=float, default=3.0,
                   help='seconds of time-stepping per run, the setup comes on top (default: 3)')
    p.add_argument('-o', '--output', help='result file (default: <workdir>/result.json)')

    p = sub.add_parser('compare', help='compare two binaries')
    p.add_argument('base', metavar='BASE.json[:LABEL]',
                   help='reference result; alone: compares its first two binaries')
    p.add_argument('other', metavar='OTHER.json[:LABEL]', nargs='?', help='result compared with the reference')

    p = sub.add_parser('verify', parents=[common],
                       help='compare the probes and dumps of all binaries and engines, after a fixed number of '
                            'timesteps')
    p.add_argument('--size', help='comma separated of {}, default: XS'.format(', '.join(cases.SIZES)))
    p.add_argument('--dumps', help='comma separated of {}, or all, default: all'.format(', '.join(cases.DUMPS)))
    p.add_argument('--timesteps', type=int, default=2000, help='timesteps of every run (default: %(default)s)')
    p.add_argument('--rtol', type=float, default=MISMATCH,
                   help='max. deviation relative to the peak (default: %(default)g)')

    args = parser.parse_args(argv)
    if args.command == 'verify':
        return cmd_verify(args)
    if args.command == 'run':
        if args.repeat is None:
            n = len(args.openems or [None]) * len(args.variant or [None])
            args.repeat = 3 if n > 1 else 1
        if args.repeat < 1:
            parser.error('--repeat must be at least 1')
        return cmd_run(args)
    return cmd_compare(args)


if __name__ == '__main__':
    sys.exit(main())
