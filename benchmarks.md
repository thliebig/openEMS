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

Commit b8ee1b7: the fused CUDA step, and the CUDA field dumps evaluated on the
device (see the notes below). The CPU row is from commit 9bdf4b1; the changes
since do not affect the CPU engine.

| Machine | Host CPU | Engine | Total run | Timestepping | Speed (MCells/s) | Peak host memory | GPU memory |
|---|---|---|---|---|---|---|---|
| Apple M5 Max (CPU) | Apple M5 Max | multithreaded | 123.2 s | 118.7 s | 303 | 716 MiB | - |
| Apple M5 Max (GPU) | Apple M5 Max | Metal | 15.4 s | 10.9 s | 3315 | 1218 MiB (1) | 343 MiB (1) |
| RTX 2080 Ti | Xeon E5-2673 v4 | CUDA | 21.9 s | 8.5 s | 4244 | 815 MiB | 457 MiB |
| RTX 3090 Ti | Threadripper PRO 5955WX | CUDA | 11.9 s | 5.2 s | 6930 | 808 MiB | 568 MiB |
| RTX 4090 | EPYC 7K62 | CUDA | 20.4 s | 8.5 s | 4234 | 804 MiB | 693 MiB |
| RTX 5070 Ti | Ryzen 7 5700X | CUDA | **11.8 s** | **4.9 s** | **7360** | 816 MiB | 530 MiB |
| RTX 5090 | EPYC 7742 | CUDA | 20.0 s | 8.9 s | 4049 | 821 MiB | 806 MiB |
| A100 SXM4 40 GB | EPYC 7K62 | CUDA | 19.0 s | 8.1 s | 4444 | 802 MiB | 723 MiB |
| H200 | Xeon Platinum 8488C | CUDA | 22.0 s | 13.6 s | 2657 | 903 MiB | 827 MiB |

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
- 9bdf4b1: the background HDF5 dumps, the Metal field snapshots and the
  fused CUDA step,
- b8ee1b7: the table above.

| Machine | e64e076 | 119949d | 9bdf4b1 | b8ee1b7 |
|---|---|---|---|---|
| Apple M5 Max (CPU) | 142.7 s | 122.4 s | 118.7 s | - |
| Apple M5 Max (Metal) | 51.1 s | 28.9 s | 10.9 s | 10.9 s |
| RTX 2080 Ti | 138.0 s | 72.9 s | 28.4 s | 8.5 s |
| RTX 3090 Ti | 87.1 s | 39.4 s | 14.9 s | 5.2 s |
| RTX 4090 (1) | 122.4 s | 63.7 s | 20.9 s | 8.5 s |
| RTX 5070 Ti (1) | 143.6 s | 81.4 s | 20.6 s | 4.9 s |
| RTX 5090 | - | - | 11.8 s | 8.9 s |
| A100 SXM4 40 GB | - | - | 17.5 s | 8.1 s |
| H200 | - | - | 16.9 s | 13.6 s |

(1) Different hosts: the first two runs were on an EPYC 7542 (RTX 4090) and a
Ryzen 9 7945HX (RTX 5070 Ti).

### NF2FF box with frequencies

The example only evaluates the far field at 15 GHz. With
`CreateNF2FFBox(frequency=[10e9, 15e9, 20e9])`, the box records those
frequencies during the run instead of dumping the time domain fields. On CUDA
the sums are kept on the device (commit 3ba4b92) and downloaded once at the
end. Since commit d2e0162 the recording is sampled with the `OverSampling`
factor, like the time domain dumps; at the Nyquist rate, as before, the
pattern at the 20 GHz band edge was up to 18 dB off.

RTX 5090 (EPYC 7742), commit d2e0162:

| NF2FF box | Total run | Timestepping | Speed (MCells/s) | Dump files |
|---|---|---|---|---|
| time domain dumps (the example) | 25.6 s | 9.1 s | 3956 | 4.7 GB |
| frequencies 10, 15, 20 GHz | **9.5 s** | **2.8 s** | **12783** | 12 MB |

The far-field patterns of both agree within 0.001 dB at 10, 15 and 20 GHz,
and the directivity to 4 digits (14.014, 16.874, 16.908 dBi). The total run
also drops because the far-field calculation no longer reads 4.7 GB of dumps.

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

- **NF2FF time-domain dumps.**
  - `CreateNF2FFBox()` without a frequency dumps E and H on the six box
    surfaces every 6 timesteps here (a quarter of the Nyquist interval of 25
    timesteps): 2484 times, 29808 datasets and 4.7 GB of HDF5 files in the
    run.
  - Before commit 119949d, the interpolation of the dumped nodes
    (`ProcessFields::CalcField()`) ran on one thread while the GPU idled. That
    was most of the GPU runs' time. It is now split over up to one thread per
    thousand nodes, bit-identical to the single-threaded dumps.
  - Since 958205f and 7c335b7, the dump files stay open and all HDF5 writes of
    the time-domain dumps run on a background thread, on every engine.
  - Metal (unified memory): at a dump timestep the fields are copied on the
    device into a snapshot, and the background thread computes the dumps from
    it while the GPU continues (commits 7c335b7 to 4641181).
  - CUDA (commit 115b334): the dumps register their node interpolation before
    the run, and a snapshot is one kernel that evaluates all dumped values on
    the device into a packed buffer (~2 MB here). It is downloaded on a second
    stream while the next timesteps run. Copying the fields instead (58 MB per
    snapshot) was limited by PCIe: 11.0 s on the RTX 5070 Ti (PCIe 4.0 x8),
    against 6.6 s. The dumps are bit-identical.

- **The fused CUDA step** needs all extensions of the simulation to support
  it. The lumped port of this example creates a lumped RLC extension (without
  elements), which blocked it until commit b8ee1b7: 6.6 s -> 4.9 s on the RTX
  5070 Ti.

- **What limits the CUDA runs now.**
  - The RTX 5070 Ti and 3090 Ti runs are within ~1.5 s of their GPU time: at
    their 300^3 speed the kernels take ~3.5 s.
  - The others are limited by writing the dumps on their host: the HDF5
    writes are serial (the library is not thread-safe) and write 4.7 GB into
    the container's overlay file system. On the RTX 5090 host (EPYC 7742) the
    dump thread was busy 9.0 s of the 9.1 s run, 6.2 s of it in the HDF5
    writes, while the GPU needed 2.7 s. On the RTX 5070 Ti host (Ryzen 7
    5700X) the same writes took 2.4 s.
  - A NF2FF box with frequencies writes no time-domain dumps, see "NF2FF box
    with frequencies" above.

- **Outside the timestepping** (total run minus timestepping: 4.5 s on the
  M5 Max GPU, 6.7 s on the RTX 3090 Ti, up to 13.4 s on the RTX 2080 Ti) is
  the setup and the post-processing. The setup runs on the host and depends
  on its single-thread speed. Since 79a89e0, the far-field calculation lists
  the time-domain datasets once instead of looking each one up by index.

- **Run-to-run variation.** The end criterion is checked at wall-clock
  intervals, so the stop timestep can differ between engines and machines.
  Here all runs stopped at 14900 timesteps. Earlier runs on the M5 Max CPU with
  an older build stopped at 15480 and 15654.
