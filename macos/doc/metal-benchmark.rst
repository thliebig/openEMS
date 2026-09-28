Metal engine benchmark
======================

``macos/bench/bench_metal.py`` measures how much faster ``--engine=metal`` is
than the fastest CPU engine on the same model, plus the peak process memory of
each run. It measures performance only; correctness is covered by
``macos/tests/metal_*.py``.

What is compared
----------------

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Config
     - Engine
   * - ``metal``
     - ``--engine=metal``: in-place diamond wavefront, 4 timesteps per block
   * - ``metal-legacy``
     - ``OPENEMS_METAL_FUSED_PIPELINE=0``: separate E/H dispatches per
       half-step (diagnostic path, kept for A/B)
   * - ``mt-N``
     - ``--engine=multithreaded --numThreads=N`` (compressed SSE); the thread
       sweep picks the fastest CPU configuration

The model is a uniform grid with one dielectric box and a Gaussian excitation,
built by ``make_model()`` from ``macos/tests/metal_fields.py``. Dumps are
stripped so I/O does not distort timing. Operator setup runs on the CPU for
every engine.

Running
-------

.. code-block:: sh

   # quick run: 3.0M cells, 600 steps
   python macos/bench/bench_metal.py --openems build/openEMS --reps 3

   # large run: 17.0M cells, 1000 steps, with the legacy path for comparison
   python macos/bench/bench_metal.py --openems build/openEMS \
       --cells 256 256 256 --steps 1000 --reps 3 --mt-threads 8,10,14 \
       --compare-metal-legacy --json out/pec17m.json

   # legacy UPML conditioners
   python macos/bench/bench_metal.py --metal-legacy --boundaries PML_8

``--openems`` may be omitted when ``OPENEMS_BIN`` is set or ``build/openEMS``
exists. ``--engines mt`` runs only the CPU baseline (no GPU needed); see
``--help`` for all options.

Configurations are interleaved across repetitions and the median of each is
reported, so thermal throttling and background load hit every engine equally.

* ``wall`` — whole process: startup, operator build, stepping.
* ``step`` — the solver's ``Time for N iterations`` line; ``setup = wall - step``.
* ``MCells/s`` — cells × timesteps / ``step``.
* ``RSS`` — peak resident set size of the whole process (``/usr/bin/time``).
* speedup — fastest ``mt-N`` time / Metal time, for ``step`` and for ``wall``.

``--json`` also writes every per-repetition sample.

Results
-------

Medians of three interleaved repetitions, PEC boundaries, Release build.
*Speedup* is the diamond Metal engine against that row (step / wall).

Apple M4 Pro (14 cores, macOS 26.3)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The sweep over 8, 10 and 14 threads picked **mt-10** at 3.0M cells and **mt-8**
at 17.0M cells.

.. list-table::
   :header-rows: 1
   :widths: 22 16 10 10 11 10 21

   * - Workload
     - Engine
     - wall [s]
     - step [s]
     - MCells/s
     - RSS [MB]
     - Speedup
   * - 3.0M cells, 600 steps
     - **metal**
     - 1.358
     - 0.314
     - 5752
     - 463
     - —
   * -
     - metal-legacy
     - 1.862
     - 0.816
     - 2216
     - 460
     - 2.60x / 1.37x
   * -
     - mt-10
     - 1.866
     - 0.737
     - 2455
     - 320
     - **2.34x / 1.37x**
   * - 17.0M cells, 1000 steps
     - **metal**
     - 9.585
     - 4.200
     - 4042
     - 2408
     - —
   * -
     - metal-legacy
     - 10.949
     - 5.580
     - 3042
     - 2398
     - 1.33x / 1.14x
   * -
     - mt-8
     - 15.983
     - 7.930
     - 2141
     - 1656
     - **1.89x / 1.67x**

Apple M5 Max (18 cores, macOS 26.6)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The sweep over 6, 10, 14 and 18 threads picked **mt-14** at 3.0M cells and
**mt-10** at 17.0M cells. Measured before shaders were built with
``-fno-fast-math``; on the M4 Pro that change moved diamond stepping by under 3%.

.. list-table::
   :header-rows: 1
   :widths: 22 16 10 10 11 10 21

   * - Workload
     - Engine
     - wall [s]
     - step [s]
     - MCells/s
     - RSS [MB]
     - Speedup
   * - 3.0M cells, 600 steps
     - **metal**
     - 1.063
     - 0.123
     - 14734
     - 462
     - —
   * -
     - metal-legacy
     - 1.452
     - 0.511
     - 3536
     - 460
     - 4.15x / 1.37x
   * -
     - mt-14
     - 1.775
     - 0.774
     - 2336
     - 320
     - **6.29x / 1.67x**
   * - 17.0M cells, 1000 steps
     - **metal**
     - 6.342
     - 1.472
     - 11532
     - 2406
     - —
   * -
     - metal-legacy
     - 7.594
     - 2.734
     - 6208
     - 2398
     - 1.86x / 1.20x
   * -
     - mt-10
     - 13.067
     - 6.510
     - 2608
     - 1656
     - **4.42x / 2.06x**

Analysis
--------

* **Against the CPU**, diamond Metal steps 1.9-2.3x faster on the M4 Pro and
  4.4-6.3x faster on the M5 Max. Wall-clock gains are smaller (1.4-2.1x)
  because operator setup stays on the CPU and is the same for every engine.
* **Against the legacy path**, the diamond kernel steps 1.3-4.2x faster. Four
  timesteps reuse each tile's cache working set, and four phase dispatches
  replace eight whole-grid dispatches per block. The gain is largest on small
  grids, where the legacy path is dispatch-bound.
* **Memory**: both Metal paths keep one E/H field pair; the diamond schedule and
  tile-local source records add under 10 MB. Metal RSS is about 45% above the
  CPU engine at both sizes. Coefficient dictionaries shrink the GPU working
  set, not the host process.
* **Tiling**: a shortest span of two cells leaves enough threadgroups on large
  grids; wider tiles reduced occupancy and lost the large-grid gain.

Caveats
-------

* **Always pass** ``--numThreads``. Without it ``--engine=multithreaded`` (and
  therefore ``--engine=fastest``) runs a single worker
  (``Engine_Multithread::Init``); ``--include-auto`` reproduces this.
* **This is Metal's best case.** A regular grid compresses extremely well; real
  PCB geometry with dense coefficients can lose much of the speedup (see
  ``metal-engine.rst``).
* **Scope.** ADE models abort on the default path; UPML runs in the diamond, and
  ``--metal-legacy --boundaries PML_8`` measures the old conditioners. These
  are finite-run timings, not accuracy or stability checks.
