# Metal GPU engine: optimization notes

Possible performance improvements for the Metal backend of the GPU engine
(`--engine=gpu`, `FDTD/metal/`), with the measurements they are based on.

AI disclosure: measured and written up with Claude Opus 5 (Claude Code).

## Baseline measurement

Release build, Apple M5 Max, free-space mesh of 200^3 = 8 million cells with
8-cell PML on all sides, soft dipole source, one field probe, 4000 timesteps.

- Engine speed: **1.8 GCells/s** (24 s wall-clock).
- Main thread (macOS `sample`, 8 s during the run):

| Main thread activity | Share |
|---|---|
| Waiting for the GPU (`waitUntilCompleted`) | 96.2 % |
| End-criterion energy on the host (`Engine_Interface_FDTD::CalcFastEnergy`) | 3.6 % |
| Encoding kernels | ~0.2 % |

For comparison, on the small meshes of the simulation tests the multithreaded
CPU engine ran at 60-380 MCells/s and the Metal engine at 330-2800 MCells/s.

## 1. Compressed update coefficients (largest expected gain)

The main update kernels read four coefficient arrays (`vv`, `vi`, `ii`, `iv`),
3 floats each per cell. With the fields, one timestep moves roughly
120 bytes per cell. At ~550 GB/s memory bandwidth that bounds the engine at
about 4-5 GCells/s, so the current kernels reach ~40 % of the bound.

Most cells share one of a few coefficient sets (same material, same mesh
spacing). `Operator_SSE_Compressed` already exploits this on the CPU: it stores
the distinct sets once and a small index per cell. Doing the same on the GPU
replaces 48 bytes of coefficients per cell by a 2-4 byte index, which should cut
the memory traffic by about 40 %.

Notes:
- Build the table from the final coefficients in `GPU_Backend_Metal::Init()`
  (or take it from an `Operator_SSE_Compressed`-derived GPU operator).
- Fall back to the full arrays if the number of distinct sets is too large.
- Coefficients differ per direction; the table entry holds all 12 values.

## 2. Main update kernels: memory access

- Try other threadgroup shapes than the current (execution width x 4 x 1),
  e.g. longer runs along z (the contiguous direction) and more rows per group.
- Each thread reads the neighbouring currents/voltages of the previous line
  in y and x; check whether staging a tile in threadgroup memory helps, or
  whether the GPU caches already catch these reads.
- Measure with Xcode's Metal System Trace / GPU counters (bandwidth,
  occupancy) before and after each change.

## 3. Fewer, larger dispatches

A timestep currently encodes about ten dispatches, several of them tiny (the six
UPML regions, the excitation, the Mur planes). Each dispatch has a fixed cost
and a serial encoder waits for the previous one to finish.

- Merge the UPML regions of one hook into a single dispatch (one kernel over a
  list of regions).
- Fold the excitation into the main update kernel, or into the UPML dispatch.
- Consider concurrent dispatch (`MTLDispatchTypeConcurrent`) with explicit
  memory barriers for independent extension kernels.

## 4. End-criterion energy on the GPU (small, certain gain)

Every 100 timesteps (at least, see `RunFDTD`) the host sums the field energy
over the whole grid, single-threaded, while the GPU is idle: ~3.5 % of the run
time. The share should stay about the same for larger meshes, since both the sum
and the GPU work between checks scale with the number of cells.

- Compute the energy with a parallel reduction kernel at the end of the batch.
- Better: do not stall for it. Encode the reduction into the batch and read the
  result of the previous batch, so the check is one batch late but never blocks.
- The engine interface would need a way to ask the engine for the energy
  (e.g. a virtual `Engine::CalcFastEnergy()` that `Engine_GPU` overrides).

## 5. Field processing without full synchronization

Each batch ends with a full device synchronization so the host can run the
probes and dumps on the shared memory. Batches are at most a Nyquist period long
(28 timesteps in the baseline), so the GPU idles for the host processing every
batch.

- Probes only read a few values: record them on the device every timestep (like
  the steady-state extension does) and copy them out asynchronously.
- For dumps the synchronization is needed, but their host work could overlap
  with the next batch if the fields of the dump step are first copied to a
  staging buffer on the device.

## Known limits (not performance, but related)

- Float precision only (Metal has no double).
- 32-bit indexing in the kernels: meshes up to ~1.4 billion cells.
- The Metal backend keeps all fields and coefficients in GPU-visible memory;
  with compression (1.) the coefficient memory shrinks accordingly.
