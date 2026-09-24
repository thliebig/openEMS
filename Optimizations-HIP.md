# HIP GPU engine: optimization notes

Performance work on the HIP backend of the GPU engine (`--engine=gpu`,
`FDTD/hip/`): what was done, what it gained, and what is left. Many of the
changes are ports of the Metal work (see `Optimizations-Metal.md`).

AI disclosure: measured, implemented and written up with Claude Opus 5 (Claude Code).

## Benchmark

Release build, RTX 2080 Ti (616 GB/s), Xeon E5-2673 v4 host with 80 threads.
The mesh is free space, 200^3 = 8 million cells, with an 8-cell PML on all
sides, a soft dipole source, one field probe and 4000 timesteps. The graded
variant uses mesh spacings that vary with a period of 97 lines in all three
directions (1.47 million distinct coefficient sets). Profiles were taken with
`nvprof`.

## Results

| Step (commit) | Uniform mesh | PEC walls instead of PML | Graded mesh |
|---|---|---|---|
| Baseline | 1211 MCells/s | 1558 | - |
| End-criterion energy on the GPU | 1292 | - | - |
| Fields read from the device on demand | 2456 | 4544 | - |
| Compressed update coefficients | 2990 | 6766 | - |
| UPML fused with the main updates, compressed UPML coefficients | 4418 | - | 3480 |
| UPML regions along z in the main kernels | **5210** | - | **3988** |

For comparison, the multithreaded CPU engine runs the benchmark at 186 MCells/s
on the same machine (80 threads).

Every step is bit-exact: `python/Tests/GPU_Engine.py` shows the same
deviations before and after. The coefficient steps were also run with each
layout forced (full arrays, 16 bit and 32 bit indices). All physics tests in
`python/Tests/` pass on the HIP engine.

### What was done

1. **Fields read from the device on demand** (`Engine_GPU`, shared with all
   backends without shared memory). This was the largest item.
   - Before, both fields were copied to the host after every batch: 47 % of
     the GPU time, since the benchmark's batches are only ~6 timesteps long.
   - Now the host copy is marked out of date at the end of a batch.
     `GetVolt`/`GetCurr` read the z-line of a value from the device
     (`GPU_Backend::DownloadRange()`), which is what probes need.
   - After 256 lines in a batch, the whole field is copied once, which is what
     dumps need. A batch that needed the whole field copies it at its first
     read in the next batch.
2. **End-criterion energy on the GPU.** This is the Metal kernel: sums in float
   along x per (y,z) line, and the host sums the lines in double. Without it,
   the energy check alone would force a full copy.
3. **Compressed update coefficients.** A 16 or 32 bit set index per node and the
   distinct sets. The set search is shared with Metal (`gpu_coeff_sets.cpp`).
4. **Fused and compressed UPML**, as on Metal. The fusion conditions
   (`GPU_UPMLFusionBox()`) are shared by both backends.
5. **UPML regions along z in the main kernels.**
   - Those regions are only 9 nodes deep along z, the contiguous direction.
     Their rows touch 36 bytes of the fields, which wastes most of every memory
     transaction: 124 us per launch against 49 us for the other regions.
   - They cover exactly the x/y range of the main updates, so the main kernels
     now run along full z lines and do the fused UPML update for those nodes.
   - This could be ported to Metal too, where the same layout effect should
     exist.

## Where the time goes now

- The main kernels take ~72 % of the GPU time and the x/y UPML regions ~24 %,
  both at ~450-510 GB/s. That is close to what the card delivers in practice.
- Gaps between the kernels and host work are negligible for this mesh size.
- Further gains need fewer bytes per timestep, e.g. temporal blocking
  (several timesteps per pass over a tile). That is a large change.

## Remaining ideas

- **Kernel launch overhead** (CUDA Graphs, fewer launches). A timestep is ~11
  launches; for 8 million cells that is ~2 % of the time. It matters for small
  meshes of ~100k cells, where the GPU work per timestep is only tens of
  microseconds.
- **Multi-grid levels on separate streams.** All grids of a cylindrical
  multi-grid share one stream. The sub-grids are small and do not fill the GPU.
- **FMA contraction** (`--fmad=false` keeps the results exact). This is probably
  irrelevant, since the kernels are memory bound.
- **Setup** (not a priority: real jobs are dominated by the computation).
  - At 200^3 the host setup takes ~17 s on this Xeon, against ~3.7 s on an M5
    Max: `Calc_EC`, `CalcTimestep`, the UPML build and the per-node
    coefficients.
  - `Calc_EC` scales poorly beyond ~10 threads on this machine: 1.7 s with 10
    threads, 4.9 s with 80.

## Known limits (not performance, but related)

- Float precision only.
- 32-bit indexing in the kernels: meshes up to ~1.4 billion cells.
- Device memory: ~24 bytes per node for the fields plus the compressed
  coefficients, or 72 bytes per node with the full arrays. The host keeps a
  pinned mirror of the fields (24 bytes per node).
- `CMAKE_HIP_ARCHITECTURES` is set to `native`, so a build only runs on GPUs of
  the build machine's architecture. Packaged builds need an explicit list (for
  example `75;80;86;89;90`) plus PTX for newer GPUs.
- The host fallback, used only for an extension without a CUDA implementation,
  moves both fields across PCIe twice per half-step.
