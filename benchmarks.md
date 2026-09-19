# GPU engine benchmarks

The GPU engine (`engine='gpu'`: Metal on macOS, CUDA elsewhere) against the
multithreaded CPU engine.

AI disclosure: measured and written up with Claude Opus 5 (Claude Code).

## Tests

**Horn antenna**: the horn antenna tutorial by Thorsten Liebig
(`Horn_Antenna.py`, openEMS v0.37 tutorials), run unchanged including its
post-processing (port, far field at 15 GHz, plots), without the geometry
viewer.

- a pyramidal horn fed by a coaxial pin (lumped port), 10-20 GHz excitation,
- 123 x 111 x 177 = 2.42 million cells, PML_8 on all sides,
- a NF2FF box with time-domain dumps (4.7 GB of HDF5 in the run),
- every run stopped after 14900 timesteps with the same result: directivity
  16.9 dBi, aperture efficiency 64.6 %.

**Free space**: `python/Tests/FreeSpace_Benchmark.py`, the field updates alone.

- 300 x 300 x 300 = 27 million cells, uniform mesh,
- a soft dipole in the center, one probe, 800 timesteps, no dumps,
- PML_8 on all sides, or PEC walls instead.

## Results

| Machine | Host CPU | Horn: total run | Horn: timestepping | Horn: MCells/s | Horn: host memory | Horn: GPU memory | Free space PML_8: MCells/s | Free space PEC: MCells/s |
|---|---|---|---|---|---|---|---|---|
| Apple M5 Max, CPU (multithreaded) | Apple M5 Max | 123.2 s | 118.7 s | 303 | 716 MiB | - | 485 | 1106 |
| Apple M5 Max, GPU (Metal) | Apple M5 Max | 15.4 s | 10.9 s | 3315 | 1218 MiB (1) | 343 MiB (1) | 4193 | 4712 |
| RTX 2080 Ti | Xeon E5-2673 v4 | 21.9 s | 8.5 s | 4244 | 815 MiB | 457 MiB | 5844 | 8410 |
| RTX 3090 Ti | Threadripper PRO 5955WX | 11.9 s | 5.2 s | 6930 | 808 MiB | 568 MiB | 9931 | 14130 |
| RTX 4090 | EPYC 7K62 | 20.4 s | 8.5 s | 4234 | 804 MiB | 693 MiB | 12744 | 18421 |
| RTX 5070 Ti | Ryzen 7 5700X | **11.8 s** | **4.9 s** | **7360** | 816 MiB | 530 MiB | 10321 | 14174 |
| RTX 5090 | EPYC 7742 | 20.0 s | 8.9 s | 4049 | 821 MiB | 806 MiB | 20072 | 27641 |
| A100 SXM4 40 GB | EPYC 7K62 | 19.0 s | 8.1 s | 4444 | 802 MiB | 723 MiB | 11536 | 18929 |
| H200 | Xeon Platinum 8488C | 22.0 s | 13.6 s | 2657 | 903 MiB | 827 MiB | **25197** | **41030** |

(1) Unified memory: the Metal buffers are part of the host memory figure.

- **Total run**: wall-clock time of the whole script, including setup and
  post-processing.
- **Timestepping**: openEMS's "Time for N iterations", including the field
  processing during the run (probes, dumps).
- **MCells/s**: cells x timesteps per second of timestepping.
- **Host memory**: peak resident set size (`/usr/bin/time`).
- **GPU memory**: CUDA: peak `nvidia-smi` memory in use minus the idle
  baseline, mostly the HIP context (the simulation's buffers are ~100 MiB).
  Metal: peak of the "graphics" categories of `footprint`.

The CUDA machines are Vast.ai containers (Ubuntu 24.04, CUDA 12.8), built in
release mode for the architecture of each GPU.

## Notes

- **The horn is limited by its NF2FF dumps on most hosts.** The dumps are
  evaluated on the GPU and written by a background thread, but HDF5 writes
  serially: 4.7 GB into the container file system. On the RTX 5090 host that
  thread was busy for the whole run (6.2 s of HDF5 writes) while the GPU
  needed 2.7 s. The fastest single-thread hosts (RTX 5070 Ti, 3090 Ti) come
  closest to their GPU time.

- **NF2FF box with frequencies.** The tutorial only evaluates the far field at
  15 GHz. With `CreateNF2FFBox(frequency=[10e9, 15e9, 20e9])` the box records
  the frequencies on the GPU during the run instead of dumping time-domain
  fields: on the RTX 5090, 9.5 s total and 2.8 s timestepping (12783
  MCells/s), 12 MB instead of 4.7 GB, with the same far field (within 0.001
  dB at all three frequencies).
