"""Estimate the run time, and the price, of a job on a given GPU.

openEMS --dry-run works out the scope of a simulation without running it: the
number of cells, the dumps and the maximum number of timesteps. Together with
the measured speed of a machine (assets/gpu_performance.json, from the
benchmarks of the GPU engine) that gives a run time, and with the price per
hour of a rented machine, a price.

    from estimate import Scope, estimate, MACHINES
    scope = Scope.from_file("sim/dry_run.json")
    print(estimate(scope, "RTX 5080", timesteps=20000))

or from the command line:

    python3 estimate.py sim/dry_run.json --timesteps 20000

What it cannot know is how many timesteps the run really takes. openEMS stops
on its end criteria, i.e. when the energy in the domain has decayed, and how
long that takes depends on the structure (a resonant one rings for longer).
Thus every number here is "for N timesteps", and the default N is the maximum
of the model, which is an upper bound, usually a generous one.
"""

import json
import os

_HERE = os.path.dirname(os.path.abspath(__file__))

# The write path of the dumps in the time domain: HDF5 writes them serially, and
# the rented machines sustained between 350 and 900 MB/s. Only a job that dumps
# the time domain fields is affected; a frequency domain dump writes once, at the
# end, and a run with probes alone writes nothing worth counting.
DISK_MB_S = 500.0
DISK_MB_S_RANGE = (350.0, 900.0)

# openEMS builds the operator before it starts, and post-processes afterwards, which
# takes a few seconds and depends on the host, not on the GPU. The median of the
# rented machines was 2.1 s per million cells (0.7 on the fastest host, 7.4 on the
# slowest), measured as the difference between the whole run and the timestepping.
SETUP_S_PER_MCELL = 2.1

# A rented machine bills from the moment it starts, and it has to pull the image and
# start the container before it can run anything. The instances of Vast.ai of the
# benchmarks needed between two and five minutes. This is added to the price, not to
# the run time.
RENT_OVERHEAD_S = 300.0


def _load_machines():
    with open(os.path.join(_HERE, 'gpu_performance.json')) as fh:
        return json.load(fh)['machines']


MACHINES = _load_machines()


class Scope:
    """what openEMS --dry-run found out about a simulation"""

    def __init__(self, data):
        self.cells = int(data['cells'])
        self.lines = tuple(data.get('lines', ()))
        self.max_timesteps = int(data['max_timesteps'])
        self.timestep_s = float(data.get('timestep_s', 0.0))
        self.end_criteria = float(data.get('end_criteria', 0.0))
        self.td_bytes_per_ts = float(data.get('time_domain_dump_bytes_per_timestep', 0))
        self.fd_bytes = float(data.get('frequency_domain_dump_bytes', 0))
        self.engine = data.get('engine', '?')
        self.processings = data.get('processings', [])

    @classmethod
    def from_file(cls, path):
        with open(path) as fh:
            return cls(json.load(fh))


def resolve(name):
    """Give the full name of a machine, from a part of it ("5080", "Metal").

    The names hold a comma ("Apple M5 Max, GPU (Metal)"), so a list of them on
    the command line cannot be split on commas alone.
    """
    if name in MACHINES:
        return name
    hits = [m for m in MACHINES if name.lower() in m.lower()]
    if len(hits) == 1:
        return hits[0]
    if not hits:
        raise KeyError('no machine matches %r, known: %s'
                       % (name, ', '.join(sorted(MACHINES))))
    raise KeyError('%r matches several machines: %s' % (name, ', '.join(sorted(hits))))


def rate_mcells_s(machine, cells):
    """Give the speed of a machine on a mesh of this size, in MCells/s.

    A GPU does not reach its full speed on a small mesh, because it cannot fill
    all of its cores. rate(N) = rate_inf * N / (N + half_cells) goes through the
    two measured points of the benchmarks.
    """
    m = MACHINES[resolve(machine)]
    half = m['half_cells']
    if not half:
        return m['rate_inf_mcells_s']
    return m['rate_inf_mcells_s'] * cells / (cells + half)


def estimate(scope, machine, timesteps=None, excitations=1, price_per_hour=None,
             disk_mb_s=DISK_MB_S, rent_overhead_s=RENT_OVERHEAD_S):
    """Estimate one job: a dict of seconds, and of dollars if a price is given.

    :param timesteps: timesteps of one excitation, default the maximum of the model
    :param excitations: number of runs, one per excited port
    """
    machine = resolve(machine)
    steps = int(timesteps or scope.max_timesteps)
    rate = rate_mcells_s(machine, scope.cells)
    gpu_s = scope.cells * steps / (rate * 1e6)

    # the dumps go through one thread, in parallel with the timestepping: a run is
    # as slow as the slower of the two
    write_bytes = scope.td_bytes_per_ts * steps + scope.fd_bytes
    write_s = write_bytes / (disk_mb_s * 1e6)
    run_s = max(gpu_s, write_s)

    setup_s = SETUP_S_PER_MCELL * scope.cells / 1e6
    total_s = (run_s + setup_s) * excitations

    out = {
        'machine': machine,
        'timesteps': steps,
        'timesteps_are_upper_bound': timesteps is None,
        'cells': scope.cells,
        'rate_mcells_s': round(rate, 1),
        'gpu_s': round(gpu_s, 1),
        'write_s': round(write_s, 1),
        'write_bound': write_s > gpu_s,
        'setup_s': round(setup_s, 1),
        'excitations': excitations,
        'total_s': round(total_s, 1),
        'write_gb': round(write_bytes / 1e9, 2),
    }
    if price_per_hour:
        out['price_per_hour'] = price_per_hour
        out['rent_overhead_s'] = rent_overhead_s
        out['cost_usd'] = round(price_per_hour * (total_s + rent_overhead_s) / 3600.0, 4)
    return out


def format_table(rows):
    """format the estimates of several machines as a table"""
    w = max([len(r['machine']) for r in rows] + [7])
    head = '%-*s %9s %8s %8s %8s %9s' % (w, 'machine', 'MCells/s', 'GPU', 'total',
                                         '$/h', 'cost')
    out = [head, '-' * len(head)]
    for r in rows:
        price = '%.3f' % r['price_per_hour'] if 'price_per_hour' in r else '-'
        cost = '$%.3f' % r['cost_usd'] if 'cost_usd' in r else '-'
        mark = ' (writes)' if r['write_bound'] else ''
        out.append('%-*s %9.0f %7.1fs %7.1fs %8s %9s%s'
                   % (w, r['machine'], r['rate_mcells_s'], r['gpu_s'], r['total_s'],
                      price, cost, mark))
    return '\n'.join(out)


def vast_price(machine, cli=None):
    """Ask the Vast.ai CLI for the cheapest offer of this GPU, $/h, or None.

    The names of the GPUs of Vast are the ones of the table, without the memory
    size ("RTX 5060 Ti 16 GB" is "RTX 5060 Ti"), and Super is "S".
    """
    import shutil
    import subprocess

    cli = cli or shutil.which('vastai') or os.path.expanduser(
        '~/.local/share/vastai/bin/vastai')
    if not os.path.exists(cli):
        return None
    name = machine.split(',')[0]
    for a, b in ((' 16 GB', ''), (' 10 GB', ''), (' 40 GB', ''),
                 (' Ti Super', 'S Ti'), (' Super', 'S')):
        name = name.replace(a, b)
    try:
        raw = subprocess.run(
            [cli, 'search', 'offers',
             'num_gpus=1 gpu_name="%s" rentable=true direct_port_count>=1' % name,
             '-o', 'dph_total', '--limit', '20', '--raw'],
            capture_output=True, text=True, timeout=60).stdout
        offers = json.loads(raw)
    except Exception:
        return None
    ok = [o for o in offers
          if (o.get('geolocation') or '').split(',')[-1].strip() != 'CN']
    return min(o['dph_total'] for o in ok) if ok else None


def main(argv=None):
    import argparse

    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('dry_run_json', help='the dry_run.json of openEMS --dry-run')
    p.add_argument('--timesteps', type=int,
                   help='timesteps per excitation (default: the maximum of the '
                        'model, an upper bound)')
    p.add_argument('--excitations', type=int, default=1,
                   help='number of runs, one per excited port (default 1)')
    p.add_argument('--machines',
                   help='comma separated; a part of each name is enough, e.g. '
                        '"5080,Metal". Default: every machine of the table')
    p.add_argument('--price', type=float, help='$/h of the machine, for one machine')
    p.add_argument('--vast', action='store_true',
                   help='look up the cheapest price of each GPU on Vast.ai')
    p.add_argument('--rent-overhead-s', type=float, default=RENT_OVERHEAD_S,
                   help='seconds a rented machine bills before it can run, for the '
                        'price (default %d)' % RENT_OVERHEAD_S)
    p.add_argument('--disk-mb-s', type=float, default=DISK_MB_S,
                   help='write speed for the dumps of the time domain (default %d)'
                        % DISK_MB_S)
    args = p.parse_args(argv)

    scope = Scope.from_file(args.dry_run_json)
    names = ([m.strip() for m in args.machines.split(',') if m.strip()]
             if args.machines else sorted(MACHINES))
    rows = []
    for name in names:
        try:
            name = resolve(name)
        except KeyError as err:
            raise SystemExit(str(err))
        price = args.price if len(names) == 1 else None
        if args.vast:
            price = vast_price(name) or price
        rows.append(estimate(scope, name, args.timesteps, args.excitations, price,
                             args.disk_mb_s, args.rent_overhead_s))
    rows.sort(key=lambda r: r['total_s'])

    print('%d cells, %s timesteps x %d excitation(s)'
          % (scope.cells, rows[0]['timesteps'],
             rows[0]['excitations']))
    if rows[0]['timesteps_are_upper_bound']:
        print('the timesteps are the maximum of the model: openEMS stops on its '
              'end criteria (%g),\nusually well before that, so these are upper '
              'bounds' % scope.end_criteria)
    if scope.td_bytes_per_ts:
        print('dumps of the time domain: %.1f GB, at %d MB/s'
              % (rows[0]['write_gb'], args.disk_mb_s))
    print()
    print(format_table(rows))


if __name__ == '__main__':
    main()
