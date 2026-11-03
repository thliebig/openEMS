# On-demand accumulated E/H exports

## Why this is useful

Frequency-domain recording surfaces already accumulate complex E/H DFT sums
throughout FDTD, but their normal HDF5 field data is written at postprocessing.
For a long run, users can inspect ports and energy decay while the radiation
pattern remains unavailable. Stopping the run to obtain those fields also
ends the simulation whose convergence they wanted to inspect.

An on-demand export makes a coherent prefix available to the existing NF2FF
tool in a second process, then lets FDTD continue accumulating. This enables
radiation-pattern progress views and diagnostics without time-domain field
volume dumps or a new excitation. It does not accelerate NF2FF or change the
solver's convergence criterion. No NF2FF backend patch is required.

## Usage and protocol

Start a new run with `--dft-snapshots` (Python: `Run(..., dft_snapshots=True)`;
Matlab: pass the command-line option to `RunOpenEMS`). The default does no
request-file polling. Existing binaries and running processes cannot acquire
this feature without a new run.

From a second process, using the same simulation directory:

```python
from openEMS.snapshots import request_dft_snapshot
snapshot = request_dft_snapshot('/path/to/simulation', 'preview-001')
# Wait for snapshot / 'snapshot.xml', or dft-snapshots/preview-001.error.
# Then run the normal NF2FF transform on this directory's E/H files.
```

Matlab/Octave: `snapshot = RequestDFTSnapshot(sim_path, 'preview-001')`.
The Matlab helper assumes one producer; the Python helper also rejects
simultaneous pending-file publication atomically using a hard link.

The request file is `DFT_SNAPSHOT`, containing one ID: 1–64 ASCII letters,
digits, hyphens or underscores. Publish it atomically, never write it in place.
There may be only one outstanding request. IDs must not be reused. The solver
polls at most once per wall-clock second, after all processings of a timestep.
It renames the request to `DFT_SNAPSHOT.active` while exporting, drains all
asynchronous processing, reads device DFT buffers without resetting them and
writes existing HDF5 E/H dumps to `dft-snapshots/.<id>.tmp/`. Once every write
has succeeded and closed, it atomically renames that directory to `<id>/`.
Only the published directory is readable input; never read the staging files.

`snapshot.xml` records the ID, build version, common timestep, simulation time,
last monitored energy estimates, end criterion and exported files. Each HDF5
file records its sample count, interval, last sampled timestep/time and E/H
time staggering. E and H retain their own Yee half-step conventions; they are
not relabelled as samples at the same physical time. The ordinary final field
files and accumulators are not changed. Disabled dumps and non-HDF5/non-E/H
processings are excluded. Missing/failed samples fail the snapshot; this does
not turn an incomplete recording surface into a valid NF2FF box.

A published snapshot is an **unconverged prefix**, regardless of whether its
pattern looks stable. Check its age, prefix length and field-energy decay.
It is not a restart checkpoint and does not contain FDTD/material/PML states.
Absolute gain/efficiency calculations also need incident port power extracted
from the same prefix; do not combine these fields with a later/full-run input
power. Normalized patterns do not require that port normalization.
The existing NF2FF tool may run independently after publication. Use a
separate result filename for each snapshot and angular grid; do not reuse
a result cached from another prefix. No NF2FF implementation is changed.

The synchronous copy and disk write temporarily pause stepping. They allocate
no second set of full-domain fields, but NF2FF itself has separate CPU/GPU
memory requirements. Sparse requests are recommended. File-system failure or
a duplicate filename/ID leaves diagnostic staging data and an `<id>.error`
file; partial results are not published. Review stale `.active` and staging
files after a crash before resubmitting. Requests are trusted local files,
not a remote-control or authentication service. No socket or port is opened.

With a fixed step cap or `--exact-endcriteria`, snapshots do not change the
numerical stopping schedule. The default end criterion can be evaluated on a
wall-clock schedule in long runs, so extra I/O can change the final overshoot
in timesteps. This is a stopping-schedule effect, not accumulator corruption.

Run `tests/manual/dft_snapshots.py --binary /path/to/openEMS --engine basic
--output /new/test/directory` for stopped-prefix and uninterrupted controls.
Use `--engine gpu` with the selected backend to check device DFT readback,
`--legacy` for legacy HDF5, and increase `--steps` for very fast devices so the
run lasts beyond the one-second request polling interval.

Pass `--nf2ff` to the manual test to consume the second snapshot with the
existing CPU NF2FF tool while FDTD is still running. This tests the complete
export/consumer workflow without an optional GPU postprocessing backend.
The same-engine field comparisons allow `rtol=1e-6` plus `1e-7` times the
reference array peak as an absolute floor (about float32 rounding precision);
probe comparisons use `rtol=1e-6`, `atol=0`. Maximum errors are recorded,
including when the measured difference is exactly zero.
