# GPU engine benchmarks

The GPU engine (`engine='gpu'`: Metal on macOS, CUDA elsewhere) against the
multithreaded CPU engine.

AI disclosure: measured and written up with Claude Opus 5 (Claude Code).

## Tests

**Horn antenna, time-domain NF2FF**: the horn antenna tutorial by Thorsten
Liebig (`Horn_Antenna.py`, openEMS v0.37 tutorials), run unchanged including
its post-processing (port, far field at 15 GHz, plots), without the geometry
viewer.

- a pyramidal horn fed by a coaxial pin (lumped port), 10-20 GHz excitation,
- 123 x 111 x 177 = 2.42 million cells, PML_8 on all sides,
- a NF2FF box with time-domain dumps (4.7 GB of HDF5 in the run),
- every run stopped after 14900 timesteps with the same result: directivity
  16.9 dBi, aperture efficiency 64.6 %.

**Horn antenna, frequency-domain NF2FF**: the same script with
`CreateNF2FFBox(frequency=[f0])`: the box records the 15 GHz fields during the
run instead of dumping the time-domain fields. Same stop timestep and result.

**Free space**: `python/Tests/FreeSpace_Benchmark.py`, the field updates alone.

- 300 x 300 x 300 = 27 million cells, uniform mesh,
- a soft dipole in the center, one probe, 800 timesteps, no dumps,
- PML_8 on all sides, or PEC walls instead.

## Results

### Horn antenna, time-domain NF2FF

| Machine | Host CPU | Total run | Timestepping | MCells/s | Host memory | GPU memory |
|---|---|---|---|---|---|---|
| Apple M5 Max, CPU (multithreaded) | Apple M5 Max | 123.2 s | 118.7 s | 303 | 716 MiB | - |
| Apple M5 Max, GPU (Metal) | Apple M5 Max | 15.4 s | 10.9 s | 3307 | 1216 MiB (1) | 343 MiB (1) |
| GTX 1080 Ti | EPYC 7551 | 32.4 s | 19.0 s | 1900 | 838 MiB | 436 MiB |
| GTX 1660 Ti | Ryzen 9 3900X | 27.1 s | 19.5 s | 1846 | 817 MiB | 373 MiB |
| RTX 2060 | i9-9900K | 24.0 s | 15.9 s | 2265 | 817 MiB | 385 MiB |
| RTX 2060 Super | Ryzen 5 3600 | 20.3 s | 11.3 s | 3187 | 816 MiB | 391 MiB |
| RTX 2070 | Xeon E3-1225 V2 | 28.9 s | 12.8 s | 2815 | 811 MiB | 397 MiB |
| RTX 2070 Super | Xeon E5-2673 v3 | 26.3 s | 14.4 s | 2506 | 816 MiB | 403 MiB |
| RTX 2080 Ti | Xeon E5-2673 v4 | 22.1 s | 8.4 s | 4277 | 813 MiB | 457 MiB |
| RTX 3050 | Xeon E5-2660 v3 | 32.3 s | 20.6 s | 1744 | 810 MiB | 388 MiB |
| RTX 3060 | Ryzen 9 5900X | 21.2 s | 13.0 s | 2766 | 810 MiB | 410 MiB |
| RTX 3060 Ti | Xeon E5-2690 v3 | 26.4 s | 12.5 s | 2888 | 811 MiB | 438 MiB |
| RTX 3070 | Ryzen 5 5600X | 16.6 s | 9.6 s | 3739 | 811 MiB | 460 MiB |
| RTX 3070 Ti | i7-8700 | 15.8 s | 8.1 s | 4430 | 810 MiB | 466 MiB |
| RTX 3080 10 GB | Ryzen 9 5950X | 15.2 s | 6.7 s | 5388 | 889 MiB | 522 MiB |
| RTX 3080 Ti | EPYC 7352 | 19.1 s | 8.5 s | 4239 | 810 MiB | 557 MiB |
| RTX 3090 | Ryzen 5 5600G | 12.5 s | 5.5 s | 6571 | 807 MiB | 558 MiB |
| RTX 3090 Ti | Threadripper PRO 5955WX | 11.9 s | 5.2 s | 6939 | 808 MiB | 568 MiB |
| RTX 4060 | Xeon E5-2673 v4 | 27.3 s | 13.5 s | 2674 | 810 MiB | 400 MiB |
| RTX 4060 Ti 16 GB | i7-12700 | 17.8 s | 11.8 s | 3042 | 809 MiB | 426 MiB |
| RTX 4070 | Ryzen 9 5900X | 17.1 s | 8.2 s | 4417 | 808 MiB | 460 MiB |
| RTX 4070 Super | Xeon E5-2673 v4 | 22.9 s | 10.7 s | 3366 | 810 MiB | 488 MiB |
| RTX 4070 Ti Super | EPYC 7702P | 16.8 s | 6.8 s | 5298 | 802 MiB | 516 MiB |
| RTX 4080 | Threadripper PRO 3975WX | 13.9 s | 5.7 s | 6332 | 807 MiB | 544 MiB |
| RTX 4080 Super | Xeon E5-2673 v4 | 25.8 s | 10.1 s | 3557 | 806 MiB | 555 MiB |
| RTX 4090 | EPYC 7K62 | 20.9 s | 8.8 s | 4116 | 795 MiB | 693 MiB |
| RTX 5060 | Xeon E5-2683 v4 | 28.5 s | 12.4 s | 2892 | 818 MiB | 420 MiB |
| RTX 5060 Ti 16 GB | Ryzen 9 5950X | 15.1 s | 8.9 s | 4066 | 820 MiB | 436 MiB |
| RTX 5070 | Threadripper PRO 3955WX | 14.4 s | 6.1 s | 5906 | 816 MiB | 470 MiB |
| RTX 5070 Ti | Ryzen 7 5700X | 11.7 s | 4.9 s | 7361 | 818 MiB | 530 MiB |
| RTX 5080 | Ryzen 9 9900X | **8.0 s** | **4.2 s** | **8478** | 818 MiB | 569 MiB |
| RTX 5090 | EPYC 7742 | 20.3 s | 9.3 s | 3894 | 821 MiB | 806 MiB |
| A100 SXM4 40 GB | EPYC 7K62 | 18.9 s | 8.0 s | 4530 | 801 MiB | 723 MiB |
| H200 | Xeon Platinum 8488C | 20.3 s | 13.0 s | 2775 | 904 MiB | 827 MiB |

### Horn antenna, frequency-domain NF2FF

| Machine | Host CPU | Total run | Timestepping | MCells/s | Host memory | GPU memory |
|---|---|---|---|---|---|---|
| Apple M5 Max, CPU (multithreaded) | Apple M5 Max | 121.8 s | 120.1 s | 300 | 496 MiB | - |
| Apple M5 Max, GPU (Metal) | Apple M5 Max | 12.4 s | 10.7 s | 3377 | 943 MiB (1) | 343 MiB (1) |
| GTX 1080 Ti | EPYC 7551 | 37.1 s | 19.2 s | 1879 | 659 MiB | 442 MiB |
| GTX 1660 Ti | Ryzen 9 3900X | 24.1 s | 19.6 s | 1838 | 637 MiB | 379 MiB |
| RTX 2060 | i9-9900K | 20.8 s | 16.2 s | 2219 | 637 MiB | 391 MiB |
| RTX 2060 Super | Ryzen 5 3600 | 16.9 s | 11.4 s | 3162 | 637 MiB | 397 MiB |
| RTX 2070 | Xeon E3-1225 V2 | 21.4 s | 12.0 s | 2998 | 636 MiB | 403 MiB |
| RTX 2070 Super | Xeon E5-2673 v3 | 18.1 s | 11.6 s | 3113 | 636 MiB | 409 MiB |
| RTX 2080 Ti | Xeon E5-2673 v4 | 15.6 s | 8.3 s | 4342 | 636 MiB | 463 MiB |
| RTX 3050 | Xeon E5-2660 v3 | 27.3 s | 20.7 s | 1737 | 630 MiB | 394 MiB |
| RTX 3060 | Ryzen 9 5900X | 17.9 s | 13.1 s | 2751 | 630 MiB | 416 MiB |
| RTX 3060 Ti | Xeon E5-2690 v3 | 17.6 s | 10.7 s | 3353 | 631 MiB | 444 MiB |
| RTX 3070 | Ryzen 5 5600X | 13.8 s | 9.6 s | 3743 | 631 MiB | 466 MiB |
| RTX 3070 Ti | i7-8700 | 12.6 s | 8.1 s | 4422 | 631 MiB | 472 MiB |
| RTX 3080 10 GB | Ryzen 9 5950X | 11.8 s | 6.7 s | 5416 | 711 MiB | 528 MiB |
| RTX 3080 Ti | EPYC 7352 | 12.7 s | 6.3 s | 5675 | 631 MiB | 563 MiB |
| RTX 3090 | Ryzen 5 5600G | 10.4 s | 5.5 s | 6498 | 629 MiB | 564 MiB |
| RTX 3090 Ti | Threadripper PRO 5955WX | 9.0 s | 5.2 s | 6911 | 630 MiB | 574 MiB |
| RTX 4060 | Xeon E5-2673 v4 | 20.5 s | 13.7 s | 2629 | 630 MiB | 406 MiB |
| RTX 4060 Ti 16 GB | i7-12700 | 15.5 s | 12.0 s | 3006 | 629 MiB | 432 MiB |
| RTX 4070 | Ryzen 9 5900X | 13.4 s | 8.2 s | 4414 | 629 MiB | 466 MiB |
| RTX 4070 Super | Xeon E5-2673 v4 | 14.1 s | 7.9 s | 4570 | 630 MiB | 494 MiB |
| RTX 4070 Ti Super | EPYC 7702P | 12.0 s | 6.1 s | 5880 | 627 MiB | 522 MiB |
| RTX 4080 | Threadripper PRO 3975WX | 10.4 s | 5.6 s | 6388 | 629 MiB | 550 MiB |
| RTX 4080 Super | Xeon E5-2673 v4 | 13.7 s | 5.5 s | 6563 | 627 MiB | 561 MiB |
| RTX 4090 | EPYC 7K62 | 11.0 s | 4.0 s | 8925 | 623 MiB | 699 MiB |
| RTX 5060 | Xeon E5-2683 v4 | 18.2 s | 9.7 s | 3728 | 639 MiB | 426 MiB |
| RTX 5060 Ti 16 GB | Ryzen 9 5950X | 12.8 s | 8.9 s | 4048 | 640 MiB | 442 MiB |
| RTX 5070 | Threadripper PRO 3955WX | 11.1 s | 6.1 s | 5898 | 637 MiB | 476 MiB |
| RTX 5070 Ti | Ryzen 7 5700X | 8.9 s | 4.9 s | 7369 | 638 MiB | 536 MiB |
| RTX 5080 | Ryzen 9 9900X | **6.7 s** | 4.3 s | 8412 | 638 MiB | 575 MiB |
| RTX 5090 | EPYC 7742 | 8.9 s | 2.7 s | 13168 | 640 MiB | 812 MiB |
| A100 SXM4 40 GB | EPYC 7K62 | 11.3 s | 4.8 s | 7448 | 631 MiB | 729 MiB |
| H200 | Xeon Platinum 8488C | 7.1 s | **2.3 s** | **15362** | 729 MiB | 833 MiB |

### Free space

| Machine | Host CPU | PML_8: MCells/s | PEC: MCells/s |
|---|---|---|---|
| Apple M5 Max, CPU (multithreaded) | Apple M5 Max | 485 | 1106 |
| Apple M5 Max, GPU (Metal) | Apple M5 Max | 4192 | 4687 |
| GTX 1080 Ti | EPYC 7551 | 2619 | 4029 |
| GTX 1660 Ti | Ryzen 9 3900X | 2268 | 3399 |
| RTX 2060 | i9-9900K | 2674 | 4300 |
| RTX 2060 Super | Ryzen 5 3600 | 3913 | 6084 |
| RTX 2070 | Xeon E3-1225 V2 | 3772 | 5997 |
| RTX 2070 Super | Xeon E5-2673 v3 | 3948 | 6151 |
| RTX 2080 Ti | Xeon E5-2673 v4 | 5803 | 8429 |
| RTX 3050 | Xeon E5-2660 v3 | 2205 | 3156 |
| RTX 3060 | Ryzen 9 5900X | 3649 | 5245 |
| RTX 3060 Ti | Xeon E5-2690 v3 | 4397 | 6476 |
| RTX 3070 | Ryzen 5 5600X | 4941 | 6505 |
| RTX 3070 Ti | i7-8700 | 5858 | 8740 |
| RTX 3080 10 GB | Ryzen 9 5950X | 7458 | 10404 |
| RTX 3080 Ti | EPYC 7352 | 7486 | 9905 |
| RTX 3090 | Ryzen 5 5600G | 9439 | 13229 |
| RTX 3090 Ti | Threadripper PRO 5955WX | 9939 | 14109 |
| RTX 4060 | Xeon E5-2673 v4 | 3536 | 4796 |
| RTX 4060 Ti 16 GB | i7-12700 | 4067 | 5021 |
| RTX 4070 | Ryzen 9 5900X | 5894 | 9059 |
| RTX 4070 Super | Xeon E5-2673 v4 | 6119 | 8195 |
| RTX 4070 Ti Super | EPYC 7702P | 7970 | 12029 |
| RTX 4080 | Threadripper PRO 3975WX | 8898 | 12066 |
| RTX 4080 Super | Xeon E5-2673 v4 | 9142 | 12744 |
| RTX 4090 | EPYC 7K62 | 12695 | 18416 |
| RTX 5060 | Xeon E5-2683 v4 | 5050 | 7017 |
| RTX 5060 Ti 16 GB | Ryzen 9 5950X | 5516 | 7046 |
| RTX 5070 | Threadripper PRO 3955WX | 8160 | 10822 |
| RTX 5070 Ti | Ryzen 7 5700X | 10346 | 14179 |
| RTX 5080 | Ryzen 9 9900X | 11706 | 15124 |
| RTX 5090 | EPYC 7742 | 19774 | 27644 |
| A100 SXM4 40 GB | EPYC 7K62 | 11473 | 18938 |
| H200 | Xeon Platinum 8488C | **25178** | **41076** |

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
release mode for the architecture of each GPU, on drivers 550 to 610. The GTX
1080 Ti (Pascal) and the GTX 1660 Ti (Turing) ran a build for all architectures
from Pascal on, the one of the packages.

## Notes

- **Time-domain NF2FF is limited by its dumps on most hosts.** The dumps are
  evaluated on the GPU and written by a background thread, but HDF5 writes
  serially: 4.7 GB into the container file system. On the RTX 5090 host that
  thread was busy for the whole run (6.2 s of HDF5 writes) while the GPU
  needed 2.7 s. The fastest single-thread hosts (RTX 5070 Ti, 3090 Ti) come
  closest to their GPU time.

- **Frequency-domain NF2FF** writes 12 MB instead of 4.7 GB and gives the same
  far field. On CUDA the frequencies are summed on the GPU. On Metal a
  background thread sums them from field snapshots while the GPU continues.
