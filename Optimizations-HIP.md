# HIP GPU engine: optimization notes

Possible performance improvements for the HIP backend of the GPU engine
(`--engine=gpu`, `FDTD/hip/`). Many ideas are shared with the Metal backend
(see `Optimizations-Metal.md`). This file ranks them for a discrete GPU and adds
the ones specific to CUDA.

AI disclosure: written up with Claude Opus 5 (Claude Code).

## Status of the measurements

There is no dedicated benchmark of the HIP backend yet. The figures below are
estimates from the code and the hardware specifications, unless they are marked
as measured. Before acting on any item, run the 200^3 free-space benchmark used
for the Metal notes (8-cell PML, soft dipole, one probe, 4000 timesteps).
Profile it with Nsight Systems for the timeline and Nsight Compute (`ncu`) for
the kernels.

Measured: the simulation tests in `python/Tests/` on an RTX 2080 Ti against the
multithreaded engine on the same machine (80 CPU cores). These are wall-clock
times per test script, including Python and post-processing:

| | CUDA | CPU (multithreaded) |
|---|---|---|
| All 18 physics tests | 96 s | 315 s |
| Largest gains (LumpedRLC, MSL_With_Local_Absorbers) | 11.3 s, 3.9 s | 63.3 s, 26.4 s |
| Smallest tests (Conducting_Sheet, SteadyState) | 1.4 s, 1.4 s | 0.9 s, 0.6 s |

The smallest tests are slower on the GPU. Item 3 explains why.

Reference numbers used below: the RTX 2080 Ti has ~616 GB/s of memory
bandwidth. We assume PCIe 3.0 x16 at ~12 GB/s for pinned transfers.

## 1. Do not copy the fields to the host after every batch (CUDA specific, likely the largest gain)

`Engine_GPU::FinishBatch()` downloads both full fields to the host mirror at
the end of every batch, so the host can run probes, dumps and the end-criterion
energy. On Metal this is only a synchronization, because the memory is shared.
On a discrete GPU it is a PCIe transfer of 24 bytes per cell.

Estimate for 200^3 cells:
- Each batch downloads 192 MB, which takes ~16 ms at 12 GB/s.
- Batches are at most a Nyquist period long (28 timesteps in the Metal
  baseline).
- At an assumed 3 GCells/s, a batch computes in ~75 ms.
- The download therefore adds roughly 20 %, and the GPU idles during it.

- Download only what the processing of this batch needs. Most batches only feed
  probes, which read a few values. Dumps need full fields, but only at their
  own interval.
- Record the probe values on the device every timestep, as the steady-state
  extension does, and copy them out asynchronously (see 5.).
- Compute the end-criterion energy on the device (see 4.). Otherwise the energy
  check alone forces a full download every time it runs.
- The host mirror then becomes stale between downloads. Anything that reads it
  must request the download first, for example through `Engine_Interface_FDTD`
  or a flag set by the processing classes.

## 2. Compressed update coefficients

This is the same idea as item 1 of the Metal notes. The main kernels read four
coefficient arrays of 3 floats per cell. With the fields, one timestep moves
about 120 bytes per cell, 48 of them coefficients. At 616 GB/s that bounds the
engine at ~5 GCells/s.

Storing each distinct coefficient set once, plus a 2-4 byte index per cell (as
`Operator_SSE_Compressed` does on the CPU), should cut the memory traffic by
about 40 %. It would also cut the device memory from ~72 to ~28 bytes per cell
(see the known limits below). Keep the full arrays as a fallback when there are
too many distinct sets.

## 3. Kernel launch overhead: fewer launches, CUDA Graphs

Every extension hook is its own kernel launch. With a PML on all six sides, the
UPML alone adds 6 regions x 4 hooks = 24 launches per timestep. The two main
updates, the excitation and the other extensions come on top. At a few
microseconds of host time per launch this costs ~100 µs per timestep:

- On large meshes it is hidden: a 200^3 timestep takes ~2.7 ms of GPU time.
- On small meshes it dominates: a 100k-cell timestep takes ~33 µs of GPU work,
  so the GPU mostly waits for the host. Together with the HIP context creation
  (see 7.), this is the likely reason the smallest tests run slower than on the
  CPU.

- Merge the UPML regions of one hook into a single launch (one kernel over a
  list of regions). Do the same for the TF/SF planes where they do not share
  edges.
- Capture a batch of timesteps in a CUDA Graph and launch the graph instead.
  This removes nearly all launch overhead. The difficulty is that several
  extensions pass per-timestep values by value from the host:
  - the timestep number (excitation, TF/SF, steady-state),
  - the rotated slot indices (lumped RLC),
  - the start timestep test (Mur ABC `IsActive()`).

  These values would have to come from a device-side timestep counter,
  incremented by a tiny kernel at the end of each timestep. The steady-state
  period check runs on the host and would split the graph at those timesteps.

## 4. End-criterion energy on the GPU

This is the same idea as item 4 of the Metal notes, but more important here,
because the host computation needs the full fields on the host (see 1.).

- Compute the energy with a parallel reduction on the device, in double
  precision for the accumulation, to stay close to the host result.
- Better: run the reduction asynchronously and read the result of the previous
  batch. The check is then one batch late but never stalls the GPU.
- The engine interface needs a way to ask the engine for the energy, for
  example a virtual `Engine::CalcFastEnergy()` that `Engine_GPU` overrides.

## 5. Overlap transfers and host processing with computation

All transfers currently run on the single work stream and end with a stream
synchronization, so the GPU idles while the host processes the fields.

- Copy the data a batch needs into a device staging buffer, then download the
  staging buffer on a second stream while the next batch computes. Record an
  event after the staging copy and wait on it before the host processing.
- The host mirror is already pinned (`PinHostMemory`), which asynchronous
  copies require.

## 6. Main update kernels

- **March along x inside a thread.** Currently each block handles one x plane
  (`blockIdx.z`), and the neighbour at x-1 (voltages) or x+1 (currents) is
  reloaded through the L2 cache. Letting each thread loop over x keeps the
  previous plane in registers. This is the usual 2.5D blocking of GPU FDTD
  codes and saves one field read per component.
- **Block shape.** `Impl::Block()` uses 32x8 threads along (z, y). For meshes
  with a z line count just above a multiple of 32 (for example nz=40), a large
  share of the threads is idle. Try picking the shape from nz, or mapping a
  flattened (y, z) index onto the threads.
- **FMA contraction.** The kernels are built with `--fmad=false`, so the
  results match the CPU engines exactly. The kernels are memory bound, so this
  probably costs little. Measure it once with FMA allowed. If the gain matters,
  offer a build option that trades exact agreement for speed.
- Check register use and occupancy with `ncu` (or `-Xptxas -v`) before and
  after each change.

## 7. Startup cost

- HIP context creation takes a significant fraction of a second, and
  `GPU_Backend_HIP::New()` triggers it only when the engine is created, after
  the operator setup. Starting the context earlier on a background thread (for
  example with `hipFree(0)` when the GPU engine is selected) would overlap it
  with the operator setup.
- `GPU_Backend_HIP::Init()` gathers the coefficients through four virtual
  calls per cell and component, single-threaded. For large meshes, gather them
  in parallel, or read the operator arrays directly when the operator type is
  known.
- `Alloc()` waits for the stream after each buffer initialization. That is
  cheap for the few dozen buffers of a simulation, but it could be batched into
  one synchronization at the end of the setup.

## 8. Multi-grid levels on separate streams

All grids of a cylindrical multi-grid simulation share one stream, so the base
grid and the sub-grids update one after the other. Sub-grids are small and do
not fill the GPU. Giving each level its own stream, with events at the coupling
points (`SyncVoltages`, `SyncCurrents`, `InterpolateToBase`), would let the
updates of different levels overlap.

## Known limits (not performance, but related)

- Float precision only. The kernels could be templated on the float type, but
  doubles are slow on consumer GPUs.
- 32-bit indexing in the kernels: meshes up to ~1.4 billion cells.
- Device memory: ~72 bytes per cell for the fields and the coefficients, plus
  the extension buffers. That is ~300 million cells on a 22 GB card, ~780
  million with compressed coefficients (2.). The host keeps a pinned mirror of
  the fields (24 bytes per cell).
- `CMAKE_HIP_ARCHITECTURES` is set to `native`, so a build only runs on GPUs of
  the build machine's architecture. Packaged builds need an explicit list (for
  example `75;80;86;89;90`) plus PTX for newer GPUs.
- The host fallback, used only for an extension without a CUDA implementation,
  moves both fields across PCIe twice per half-step. That is orders of
  magnitude slower than the device path, and much worse than on Metal.
