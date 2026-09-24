# Run time estimates

How long a simulation takes on a given GPU, and what it costs to rent one.

AI disclosure: measured and written up with Claude Opus 5 (Claude Code).

## The dry run

`openEMS --dry-run` builds the model as a real run does, reports what the run
would cost and stops before the timestepping:

    openEMS simulation.xml --dry-run

or from Python, `FDTD.Run(sim_path, dry_run=True)`. It writes `dry_run.json`
into the simulation directory:

    {
      "cells": 2416581,
      "lines": [123, 111, 177],
      "max_timesteps": 1000000000,
      "timestep_s": 9.823618e-13,
      "end_criteria": 0.0001,
      "time_domain_dump_bytes_per_timestep": 329976,
      "frequency_domain_dump_bytes": 0,
      "processings": [ ... one entry per probe and dump ... ]
    }

The setup itself runs, so these are the numbers of the operator, not a guess:
the mesh is the real mesh and the dumps are the real dumps.

## The estimate

    python3 estimate.py sim/dry_run.json --timesteps 20000
    python3 estimate.py sim/dry_run.json --timesteps 20000 --machines "5080,Metal"
    python3 estimate.py sim/dry_run.json --timesteps 20000 --vast

gives

    machine      MCells/s      GPU    total      $/h      cost
    ----------------------------------------------------------
    RTX 5090        13168     2.7s     7.8s    0.446    $0.038
    RTX 5070 Ti      7369     4.9s    10.0s    0.201    $0.017
    RTX 4070         4414     8.2s    13.2s    0.135    $0.012

`--vast` asks the Vast.ai CLI for the cheapest offer of each GPU, and the cost
counts the five minutes a rented machine bills before it can run anything
(`--rent-overhead-s`).

## What it knows, and what it cannot know

`gpu_performance.json` comes from the measurements in [benchmarks.md](../benchmarks.md),
through `gen_gpu_table.py`. Every machine there ran the same two tests, which
give two points on its curve: the horn (2.42 million cells) and free space
(27 million cells). A GPU is slower per cell on a small mesh, because it cannot
fill all of its cores, so the rate follows

    rate(N) = rate_inf * N / (N + half_cells)

and the two points fix the two parameters. `half_cells` lands near one million
for most of the GPUs, and near two million for an H200.

**The number of timesteps is the one thing no dry run can tell you.** openEMS
stops when the energy in the domain has decayed to the end criteria, and how
long that takes depends on the structure: a resonant one rings for much longer
than a matched line. `max_timesteps` is only the ceiling of the model, and it
is usually far above the run. Thus pass `--timesteps` when you know the figure
(a previous run of the same board prints it), and read the default as an upper
bound.

Two further limits:

- **Writing the dumps of the time domain can be slower than the simulation.**
  The estimator assumes 500 MB/s (`--disk-mb-s`), and the rented machines ran
  between 239 and 1171 MB/s. A far field box that records frequencies instead
  of the time domain fields avoids the question, and openEMS writes 12 MB
  instead of 4.7 GB.
- **The host matters for everything outside the timestepping.** Building the
  operator and post-processing took 0.7 s per million cells on the fastest
  machine and 7.4 s on the slowest; the estimator uses the median, 2.1 s.

`test_estimate.py` checks the model against every row of benchmarks.md: the
median error is -1 % for a run that records frequencies (31 of 33 machines
within 30 %) and -13 % for one that dumps the time domain (28 of 33 within
40 %, the spread of the write speeds).
