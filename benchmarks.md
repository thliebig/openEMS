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

The runs were made twice:

- **before**: commit e64e076, the CUDA and Metal optimizations of
  `Optimizations-HIP.md` and `Optimizations-Metal.md`,
- **after**: commit 119949d, the field dumps of the NF2FF box computed on
  several threads (see the notes below).

Memory figures are from the "after" runs. The "before" runs were within 60 MiB.

| Machine | Engine | Total run (before -> after) | Timestepping (before -> after) | Speed after | Peak host memory | GPU memory |
|---|---|---|---|---|---|---|
| Apple M5 Max (CPU) | multithreaded | 169.3 -> 149.3 s | 142.7 -> 122.4 s | 294 MCells/s | 553 MiB | - |
| Apple M5 Max (GPU) | Metal | 77.4 -> **55.6 s** | 51.1 -> **28.9 s** | **1248 MCells/s** | 852 MiB (1) | 232 MiB (1) |
| RTX 2080 Ti | CUDA | 157.1 -> 91.4 s | 138.0 -> 72.9 s | 494 MCells/s | 636 MiB | 351 MiB |
| RTX 3090 Ti | CUDA | 95.7 -> **48.2 s** | 87.1 -> **39.4 s** | 914 MCells/s | 623 MiB | 462 MiB |
| RTX 4090 | CUDA | 136.5 -> 78.5 s | 122.4 -> 63.7 s | 565 MCells/s | 624 MiB | 587 MiB |
| RTX 5070 Ti | CUDA | 175.9 -> 111.0 s | 143.6 -> 81.4 s | 442 MCells/s | 632 MiB | 422 MiB |

(1) Unified memory: the Metal buffers (232 MiB) are part of the host memory
figure. The peak physical footprint was about 880 MiB.

Column definitions:

- **Total run**: the wall-clock time of the whole script, including setup and
  post-processing.
- **Timestepping**: openEMS's "Time for N iterations" figure, which includes
  the field processing during the run (probes, dumps).
- **Peak host memory**: the peak resident set size of the process
  (`/usr/bin/time -l` on macOS, `/usr/bin/time -v` on Linux).
- **GPU memory**:
  - CUDA: the peak `nvidia-smi` memory in use during the run minus the idle
    baseline. It includes the HIP context, which accounts for most of it and
    varies by GPU and driver. The simulation's own buffers are ~100 MiB.
  - Metal: the peak of the "graphics" categories of `footprint`.

### Machines

| Machine | GPU | Host CPU | RAM | OS, driver |
|---|---|---|---|---|
| MacBook Pro | Apple M5 Max (Metal) | Apple M5 Max | unified | macOS 26 |
| Vast.ai container | RTX 2080 Ti (22 GB) | Xeon E5-2673 v4, 80 threads, 2.3 GHz | 251 GB | Ubuntu 24.04, driver 580.126, CUDA 12.8 |
| Vast.ai container | RTX 3090 Ti (24 GB) | Threadripper PRO 5955WX, 32 threads | 125 GB | Ubuntu 22.04, driver 595.71, CUDA 12.8 |
| Vast.ai container | RTX 4090 (24 GB) | EPYC 7542, 128 threads, 2.9 GHz | 251 GB | Ubuntu 22.04, driver 580.159, CUDA 12.8 |
| Vast.ai container | RTX 5070 Ti (16 GB) | Ryzen 9 7945HX, 32 threads | 15 GB | Ubuntu 22.04, driver 580.178, CUDA 12.8 |

Build and run details:

- Everything was built from this branch in release mode, with the CUDA
  architecture of each GPU (75, 86, 89, 120).
- The Python bindings were built against each build.
- The engine was selected with `engine='gpu'` or `'multithreaded'` in
  `openEMS.Run()`.

## Notes

- **This simulation is limited by the host, not by the GPU.** The ranking of
  the CUDA machines follows the single-thread speed of their host CPUs: the
  RTX 4090 host (EPYC 7542) is slower than the RTX 3090 Ti host (Threadripper
  PRO 5955WX). For comparison, the Metal engine runs a 200^3 free-space
  benchmark without dumps at 3810 MCells/s, and the RTX 2080 Ti at 5210
  MCells/s.

- **NF2FF time-domain dumps.**
  - `CreateNF2FFBox()` without a frequency dumps E and H on the six box
    surfaces every Nyquist interval (25 timesteps here).
  - Before commit 119949d, the interpolation of the dumped nodes
    (`ProcessFields::CalcField()`) ran on one thread while the GPU idled. That
    was most of the GPU runs' time.
  - The dumps are now split over up to one thread per thousand nodes. They are
    bit-identical to the single-threaded dumps (all 29850 datasets of the
    horn example).

- **Remaining host bottleneck.** On the M5 Max GPU run, about two thirds of the
  timestepping time is still host processing and one third is waiting for the
  GPU. Two ideas are left:
  - Overlap the processing with the GPU: snapshot the fields on the device at
    the end of a batch and process the snapshot while the next batch runs. The
    time would then approach the larger of processing and GPU instead of their
    sum. This needs a change of the main loop (`openEMS::RunFDTD()`).
  - Keep the HDF5 dump files open. `HDF5_File_Writer` opens and closes the file
    for every dataset, which is about 20 % of the processing time.

- **Outside the timestepping** (total run minus timestepping: 8.8 s on the
  RTX 3090 Ti up to 29.6 s on the RTX 5070 Ti) is the setup and the
  post-processing. In a profile of the M5 Max run, the far-field calculation
  (`CalcNF2FF`, reading the time-domain dumps) spent most of its time looking
  up the HDF5 datasets by index (`HDF5_File_Reader::GetDataSetNameByIndex`).

- **Run-to-run variation.** The end criterion is checked at wall-clock
  intervals, so the stop timestep can differ between engines and machines.
  Here all runs stopped at 14900 timesteps. Earlier runs on the M5 Max CPU with
  an older build stopped at 15480 and 15654.
