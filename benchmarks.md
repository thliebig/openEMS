# GPU engine benchmarks

Wall-clock benchmarks of the GPU engine (`--engine=gpu`) against the
multithreaded CPU engine, on a realistic antenna simulation.

AI disclosure: measured and written up with Claude Opus 5 (Claude Code).

## Horn antenna with coaxial pin feed

The simulation is the horn antenna tutorial by Thorsten Liebig (`Horn_Antenna.py`,
openEMS v0.37 tutorials):

- a pyramidal horn fed by a coaxial pin (lumped port) in a closed waveguide,
- 10-20 GHz Gaussian excitation,
- PML_8 on all sides,
- a NF2FF box with time-domain dumps,
- end criterion -40 dB.

The mesh has 123 x 111 x 177 = **2.42 million cells**. Every run stopped after
**14900 timesteps** and gave the same results: directivity 16.9 dBi, aperture
efficiency 64.6 %.

The benchmark ran the script unchanged, except that the AppCSXCAD geometry
viewer is not launched. That includes the post-processing: port calculation,
far field at 15 GHz, and plots written with the non-interactive matplotlib
backend.

### Results

Commit 9bdf4b1: the fused CUDA step and the background HDF5 dumps (see the
notes below).

| Machine | Host CPU | Engine | Total run | Timestepping | Speed (MCells/s) | Peak host memory | GPU memory |
|---|---|---|---|---|---|---|---|
| Apple M5 Max (CPU) | Apple M5 Max | multithreaded | 123.2 s | 118.7 s | 303 | 716 MiB | - |
| Apple M5 Max (GPU) | Apple M5 Max | Metal | **16.1 s** | **10.9 s** | **3302** | 1218 MiB (1) | 343 MiB (1) |
| RTX 2080 Ti | Xeon E5-2673 v4 | CUDA | 41.1 s | 28.4 s | 1270 | 805 MiB | 367 MiB |
| RTX 3090 Ti | Threadripper PRO 5955WX | CUDA | 21.6 s | 14.9 s | 2414 | 802 MiB | 478 MiB |
| RTX 4090 | EPYC 7K62 | CUDA | 33.1 s | 20.9 s | 1721 | 775 MiB | 603 MiB |
| RTX 5070 Ti | Ryzen 7 5700X | CUDA | 27.5 s | 20.6 s | 1746 | 810 MiB | 438 MiB |
| RTX 5090 | EPYC 7742 | CUDA | 22.3 s | 11.8 s | 3057 | 814 MiB | 714 MiB |
| A100 SXM4 40 GB | EPYC 7K62 | CUDA | 28.6 s | 17.5 s | 2058 | 776 MiB | 633 MiB |
| H200 | Xeon Platinum 8488C | CUDA | 24.5 s | 16.9 s | 2127 | 900 MiB | 735 MiB |

(1) Unified memory: the Metal buffers (343 MiB) are part of the host memory
figure. The peak physical footprint was 1259 MiB.

Column definitions:

- **Total run**: the wall-clock time of the whole script, including setup and
  post-processing.
- **Timestepping**: openEMS's "Time for N iterations" figure, which includes
  the field processing during the run (probes, dumps).
- **Speed**: cells times timesteps per second of timestepping, in millions
  (openEMS's "Speed" figure): 2.42 million cells x 14900 timesteps divided by
  the timestepping time.
- **Peak host memory**: the peak resident set size of the process
  (`/usr/bin/time -l` on macOS, `/usr/bin/time -v` on Linux).
- **GPU memory**:
  - CUDA: the peak `nvidia-smi` memory in use during the run minus the idle
    baseline. It includes the HIP context, which accounts for most of it and
    varies by GPU and driver. The simulation's own buffers are ~100 MiB.
  - Metal: the peak of the "graphics" categories of `footprint`.

### Earlier runs

Timestepping time by commit:

- e64e076: the CUDA and Metal optimizations of `Optimizations-HIP.md` and
  `Optimizations-Metal.md`,
- 119949d: the field dumps of the NF2FF box computed on several threads,
- 9bdf4b1: the table above.

| Machine | e64e076 | 119949d | 9bdf4b1 |
|---|---|---|---|
| Apple M5 Max (CPU) | 142.7 s | 122.4 s | 118.7 s |
| Apple M5 Max (Metal) | 51.1 s | 28.9 s | 10.9 s |
| RTX 2080 Ti | 138.0 s | 72.9 s | 28.4 s |
| RTX 3090 Ti | 87.1 s | 39.4 s | 14.9 s |
| RTX 4090 (1) | 122.4 s | 63.7 s | 20.9 s |
| RTX 5070 Ti (1) | 143.6 s | 81.4 s | 20.6 s |

(1) Different hosts: the earlier runs were on an EPYC 7542 (RTX 4090) and a
Ryzen 9 7945HX (RTX 5070 Ti).

## Free space, 300^3 cells

The GPU kernels alone, without the host work of the horn example:

- 300 x 300 x 300 = 27 million cells, uniform mesh,
- a soft dipole source in the center, one field probe,
- 800 timesteps, no end criterion, no dumps,
- PML_8 on all sides, or PEC walls instead.

The script is `python/Tests/FreeSpace_Benchmark.py` (`python
FreeSpace_Benchmark.py [engine] [N] [timesteps]`, default `gpu 300 800`). It
also checks that the PML and PEC runs agree at the probe until the first wall
reflection can reach it.

| Machine | Engine | PML_8 (MCells/s) | PEC (MCells/s) | PML_8 before the fused step (53f1954) |
|---|---|---|---|---|
| Apple M5 Max (CPU) | multithreaded | 485 | 1106 | - |
| Apple M5 Max (GPU) | Metal | 4193 | 4712 | - |
| RTX 2080 Ti | CUDA | 5844 | 8410 | 5574 |
| RTX 3090 Ti | CUDA | 9931 | 14130 | 8134 |
| RTX 4090 | CUDA | 12744 | 18421 | 8385 |
| RTX 5070 Ti | CUDA | 10321 | 14174 | 7014 |
| RTX 5090 | CUDA | 20072 | 27641 | 12933 |
| A100 SXM4 40 GB | CUDA | 11536 | 18929 | ~7480 (1) |
| H200 | CUDA | **25197** | **41030** | 19813 |

(1) Measured at 250^3.

The fused step (commits d2867d7 to 9bdf4b1) computes the voltage and the
current update of a timestep in one pass over the fields. Each block works on
a tile of 31 x 7 (z, y) lines and marches along 4 x lines. The UPML regions
are updated in the same kernel. The results are bit-identical to the separate
kernels.

### Machines

| Machine | GPU | Host CPU | RAM | OS, driver |
|---|---|---|---|---|
| MacBook Pro | Apple M5 Max (Metal) | Apple M5 Max | unified | macOS 26 |
| Vast.ai container | RTX 2080 Ti (22 GB) | Xeon E5-2673 v4, 80 threads, 2.3 GHz | 251 GiB | Ubuntu 24.04, driver 580.126, CUDA 12.8 |
| Vast.ai container | RTX 3090 Ti (24 GB) | Threadripper PRO 5955WX, 32 threads | 125 GiB | Ubuntu 24.04, driver 595.71, CUDA 12.8 |
| Vast.ai container | RTX 4090 (24 GB) | EPYC 7K62, 192 threads | 503 GiB | Ubuntu 24.04, driver 580.126, CUDA 12.8 |
| Vast.ai container | RTX 5070 Ti (16 GB) | Ryzen 7 5700X, 16 threads | 62 GiB | Ubuntu 24.04, driver 595.91, CUDA 12.8 |
| Vast.ai container | RTX 5090 (32 GB) | EPYC 7742, 256 threads | 503 GiB | Ubuntu 24.04, driver 580.142, CUDA 12.8 |
| Vast.ai container | A100 SXM4 (40 GB) | EPYC 7K62, 192 threads | 503 GiB | Ubuntu 24.04, driver 595.84, CUDA 12.8 |
| Vast.ai container | H200 (141 GB) | Xeon Platinum 8488C, 192 threads | 1999 GiB | Ubuntu 24.04, driver 615.71, CUDA 12.8 |

The thread counts are the threads visible in the container. openEMS uses the
CPUs the container's quota allows (commit 2d51fe3), which can be fewer.

Build and run details:

- Everything was built from this branch in release mode, with the CUDA
  architecture of each GPU (75, 86, 89, 120, 80, 90).
- The Python bindings were built against each build.
- The engine was selected with `engine='gpu'` or `'multithreaded'` in
  `openEMS.Run()`.

## Notes

- **The horn example is limited by the host, not by the GPU.**
  - The CUDA runs reach 8-24 % of each GPU's 300^3 speed. At its PML_8
    speed, the RTX 5090 would need ~1.8 s for the 14900 timesteps of 2.42
    million cells, against 11.8 s measured.
  - The rest is the field processing between batches (NF2FF dumps every 25
    timesteps, probes) and the launch and synchronization overhead of a small
    grid.
  - The M5 Max GPU is the fastest machine here although its kernels are the
    slowest: with unified memory, the dumps are computed from a field
    snapshot on a background thread while the GPU continues (commits 7c335b7
    to 4641181).

- **NF2FF time-domain dumps.**
  - `CreateNF2FFBox()` without a frequency dumps E and H on the six box
    surfaces every Nyquist interval (25 timesteps here).
  - Before commit 119949d, the interpolation of the dumped nodes
    (`ProcessFields::CalcField()`) ran on one thread while the GPU idled. That
    was most of the GPU runs' time.
  - The dumps are now split over up to one thread per thousand nodes. They are
    bit-identical to the single-threaded dumps (all 29850 datasets of the
    horn example).
  - Since 958205f and 7c335b7, the dump file stays open and all HDF5 writes of
    the time-domain dumps run on a background thread, on every engine.

- **Next step for CUDA: field snapshots on the device.** On CUDA, the dumps
  still copy the fields to the host and compute the dumped values while the
  GPU waits. A copy into a second device buffer, downloaded on a separate
  stream and processed on the background thread, would overlap that work with
  the next batch, as on Metal.

- **Outside the timestepping** (total run minus timestepping: 5.2 s on the
  M5 Max GPU, 6.7 s on the RTX 3090 Ti, up to 12.7 s on the RTX 2080 Ti) is
  the setup and the post-processing. The setup runs on the host and depends
  on its single-thread speed. Since 79a89e0, the far-field calculation lists
  the time-domain datasets once instead of looking each one up by index.

- **Run-to-run variation.** The end criterion is checked at wall-clock
  intervals, so the stop timestep can differ between engines and machines.
  Here all runs stopped at 14900 timesteps. Earlier runs on the M5 Max CPU with
  an older build stopped at 15480 and 15654.
