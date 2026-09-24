# Metal GPU engine: optimization notes

Performance work on the Metal backend of the GPU engine (`--engine=gpu`,
`FDTD/metal/`): what was done, what it gained, and what is left.

AI disclosure: measured, implemented and written up with Claude Opus 5 (Claude Code).

## Benchmark

Release build, Apple M5 Max, free-space mesh of 200^3 = 8 million cells with
8-cell PML on all sides, soft dipole source, one field probe, 4000 timesteps.
The graded variant uses mesh spacings that vary with a period of 97 lines in all
three directions (1.47 million distinct coefficient sets).

## Results

| Step (commit) | Uniform mesh | Graded mesh |
|---|---|---|
| Baseline | 1957 MCells/s | - |
| Compressed update coefficients | 2496 | - |
| End-criterion energy on the GPU | 2631 | - |
| Concurrent UPML region dispatches | 2690 | - |
| UPML fused with the main updates | 3269 | - |
| 32 bit set indices, fusion without compression | 3255 | 2577 (full arrays: 2131) |
| Compressed UPML coefficients | **3810** | **2910** |

Setup (operator build and upload) went from 6.3 s to 3.7 s, the same as the
multithreaded CPU engine. The whole benchmark run went from 22.7 s to 12.0 s.

Every step is bit-exact: `python/Tests/GPU_Engine.py` shows the same
deviations from the basic engine before and after, and was run with each
coefficient layout forced (full arrays, 16 bit and 32 bit indices).

### What was done

1. **Compressed update coefficients.** Each node stores a set index (16 bit, or
   32 bit for more than 65536 sets). The distinct sets of 12 coefficients are
   stored once. The full arrays remain the fallback when index and sets would be
   more than half their size. The set search (`Metal_FindSets()`) is shared
   with the UPML.
2. **End-criterion energy on the GPU.** `GPU_Backend::CalcFastEnergy()` sums the
   energy along x per (y,z) line in float, and the host sums the lines in
   double. The result agrees with the host sum to ~3e-7. Before, the host
   summed single-threaded while the GPU idled (~3.5 % of the run).
3. **Concurrent dispatch.** The encoder is concurrent, with an explicit barrier
   before every dispatch except between dispatches of the same group. The UPML
   regions of a hook are disjoint and form a group.
4. **UPML fused with the main updates.** This was the largest single item.
   - With PEC walls instead of PML, the benchmark ran at 5187 MCells/s, against
     2657 with the PML. The PML took half the time for a quarter of the cells.
   - The UPML hooks run directly before and after the main update. Only the
     steady-state extension, which has no pre/post hooks, has a higher priority.
   - So each region now does pre-update, main update and post-update in one
     kernel, with the same operations. The main kernels only cover the box
     inside the regions.
   - The conditions are checked at the first UPML hook: the extension order,
     that all UPML extensions run on the device, and that the regions tile the
     outside of a box. The base grid of a cylindrical multi-grid fails the last
     condition and keeps the separate kernels.
5. **Compressed UPML coefficients**, the same scheme with sets of 18 coefficients
   per cell, used by the fused kernels.
6. **Parallel operator build** (`Operator_GPU`, shared with all GPU backends).
   The material and PEC evaluation is split over threads like
   `Operator_Multithread`.

### Tried without gain

- **Threadgroup shapes** for the 3D dispatches: 12 shapes from 32x1 to 256x1 and
  32x16, some 3D. All were within 1.5 % of the current 32x4x1, which was the
  fastest.
- **Packing short rows** into threadgroups for the thin z-slab UPML regions
  (9 nodes along z): no change.

## Where the time goes now

- The main thread spends over 99 % of the run waiting for the GPU. Kernel
  encoding, the batch synchronization and the host processing no longer
  register.
- The main updates alone (PEC walls) move ~76 bytes per node and timestep. At
  5187 MCells/s that is ~390 GB/s of memory traffic.
- The fused PML regions still read and write the UPML flux (24 bytes per
  field, cell and half-step). In the z-slabs (9 nodes along z) each row touches
  only 36 bytes of the fields, so most of every cache line is wasted. That is a
  consequence of the field layout.

## Remaining ideas

- **March along x inside a thread** in the main kernels. Keeping the previous
  plane in registers saves the neighbour reads through the cache. Untested; the
  gain depends on how well the caches already serve these reads.
- **Probes on the device.** Record probe values every timestep, as the
  steady-state extension does, and copy them out without a full
  synchronization. This only matters for small meshes with many short batches.
  For large meshes the batch synchronization no longer shows.
- **Fold the excitation into the fused/main kernels.** It is one tiny dispatch
  per timestep, so the gain is small.
- **Setup:** the remaining 3.7 s at 200^3 is shared with the CPU engines
  (geometry, timestep, extensions). For short simulations it is now the larger
  part of the run.

## Known limits (not performance, but related)

- Float precision only (Metal has no double).
- 32-bit indexing in the kernels: meshes up to ~1.4 billion cells.
- Metal cannot flush denormals to zero without its fast-math mode, which would
  also break the exact agreement elsewhere. The CPU SSE/multithreaded engines
  flush them, so on the cylindrical multi-grid wedge test the Metal results
  deviate from those engines by ~4e-6 of the peak value.
