/*
*	Copyright (C) 2026 Sean Mollet (sean@malmoset.com)
*
*	This program is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
*	(at your option) any later version.
*
*	This program is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*
*	You should have received a copy of the GNU General Public License
*	along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef GPU_BACKEND_H
#define GPU_BACKEND_H

#include <complex>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "openems_global.h"
#include "tools/constants.h"
#include "tools/global.h"
#include "tools/arraylib/array_nijk.h"
#include "engine_field_gather.h"
#include "multigrid_interpolation.h"

class Operator;
struct GPU_CellModel;
struct Engine_Box;

//! The hooks of a half step around the main update, in this order: pre (the extensions in reverse order), post, apply (see Engine)
enum GPU_Hook {GPU_PRE, GPU_POST, GPU_APPLY};

//! Device implementation of an engine extension, see GPU_CreateExtension() (gpu_ext.h)
/*!
  Engine_GPU calls it in the same order as the engine extension it replaces,
  while the fields stay on the device. It works on the device data of the
  backend that created it, so Engine_GPU deletes it before that backend.
  */
class GPU_Extension
{
public:
	virtual ~GPU_Extension() {}

	//! Hook \a h of the voltage (\a currents false) or current half step of timestep \a ts, as Engine_Extension::DoPreVoltageUpdates() to Apply2Current()
	/*!
	  On the nodes with x in [box.start[0], box.stop[0]): a tile of the temporal
	  blocking (the boxes span y and z) or the whole domain.
	  */
	virtual void Hook(bool currents, GPU_Hook h, const Engine_Box& box, unsigned int ts) = 0;

	//! Whether timestep \a ts needs the main updates of the whole domain, one after the other: no temporal blocking and no fused step then (e.g. the steady-state check reads the energy of all fields between them)
	virtual bool NeedsFlatStep(unsigned int ts) const {(void)ts; return false;}
	//! For the fused step (GPU_Backend::PrepareFused()): the nodes (x*ny + y)*nz + z whose voltages the hooks between the main voltage and current update write, and whose currents they read or write
	/*!
	  These hooks are the post and apply hooks of the voltages and the pre hook of the currents.
	  Returns false if the extension does not know them: no fused step then.
	  */
	virtual bool MidStepNodes(std::vector<uint64_t>& volt, std::vector<uint64_t>& curr) const {(void)volt; (void)curr; return false;}
};

//! The per-point operations of the device extensions (GPU_Ext_*, kernels/ext_ops.h) a backend runs with GPU_Backend::RunPoints(), as CM_POINT_OPS lists them
enum GPU_PointOp
{
	GPU_OP_EXCITATION,
	GPU_OP_SHEET,
	GPU_OP_LUMPED_RLC,
	GPU_OP_TFSF,
	GPU_OP_STEADY_STATE,
	GPU_OP_CYLINDER,
	GPU_OP_MULTIGRID,
	GPU_OP_DISPERSIVE,
};

//! Abstract interface to the device used by Engine_GPU
/*!
  A backend owns the device memory: the voltage and current fields and the cell
  model of the operator. It performs the main FDTD updates on the device, and runs
  the device extensions of all backends (GPU_Ext_*, gpu_ext.h) through
  DeviceAlloc(), DeviceField() and RunPoints(). What it must provide is pure
  virtual; the rest (snapshots, device DFT sums, the fused step, ...) only
  speeds things up, and Engine_GPU takes another path without it.

  Engine_GPU keeps a host mirror of the fields in the basic Engine layout (ArrayNIJK)
  for the field processing, and brings it up to date with Download(). If the device shares its
  memory with the host (unified memory), the host mirror is the device memory
  itself (see GetShared()) and synchronizing means waiting for the device.

  Init() is called once and before everything else; the rest may follow in any
  order, except where a method says otherwise. The device work of all calls is
  queued in the order they are made and all of them run on the engine thread,
  only WaitSnapshot() may be called from another thread.

  The backend outlives what works on its device data (sub-grid backend, device
  extensions, multi-grid coupling), and the caller frees those first.
  */
class GPU_Backend
{
public:
	//! Create a backend by name: "reference" or "auto" (the best device backend available), as chosen by --engine=gpu-reference / --engine=gpu. It's the responsibility of the caller to free it.
	/*!
	  Returns NULL if "auto" finds no device backend: there is no silent fallback to the
	  (slow) reference backend.
	  */
	static GPU_Backend* New(const std::string& name);

	virtual ~GPU_Backend() {}

	virtual std::string GetName() const = 0;

	//! The bytes of the device's last-level cache, 0 if the backend cannot tell; \a source says where they come from (temporal blocking, see Engine_GPU::PrepareBlocking())
	virtual size_t CacheBytes(std::string& source) const {source.clear(); return 0;}

	//! Allocate the device fields (all zero) and upload \a model, the cell model of the grid of \a op
	/*!
	  The model fuses the UPML and (where it can) the dispersive materials into
	  the main updates (kernels/cell_update.h); Engine_GPU drops their engine
	  extensions. The results are within single precision of the CPU engine on
	  the same model (Engine_CPU with the flat sweep, its generic kernels and
	  the fast axis z): the coefficients in float per node (CM_FACTOR).

	  Returns false if the backend cannot run this grid. After a successful Init(),
	  a device failure (a failed launch, transfer or required allocation, a lost
	  device) throws std::runtime_error; it ends the run with an error message, see
	  main() and openEMS::RunFDTD(). Memory that only speeds things up (the fused
	  step, snapshots) is optional: without it the backend takes the slower path.
	  */
	virtual bool Init(const Operator* op, const GPU_CellModel& model) = 0;

	//! Device memory of \a bytes, initialized with \a host or zero, freed with the backend
	virtual void* DeviceAlloc(size_t bytes, const void* host=NULL) = 0;
	//! Copy \a bytes of the device memory \a dev to \a host, after the work so far (waits for it)
	virtual void DeviceDownload(void* host, const void* dev, size_t bytes) = 0;
	//! The device voltages (false) or currents (true), ArrayNIJK layout
	virtual float* DeviceField(bool currents) = 0;
	//! The main voltage or current update on the cell model of the nodes with x in [box.start[0], box.stop[0]): a tile of the temporal blocking (the boxes span y and z), or the whole domain
	virtual void UpdateVoltages(const Engine_Box& box) = 0;
	virtual void UpdateCurrents(const Engine_Box& box) = 0;
	//! Prepare the fused step: the main voltage and current update in one pass over the fields, see FusedUpdate()
	/*!
	  \a late are the nodes (x*ny + y)*nz + z, sorted, whose currents must wait for
	  the hooks between the two main updates: the hooks read them, or their curl
	  reads a voltage the hooks write. Returns false if the backend has no fused
	  step, or not the memory for it (it is only a speed-up); then the main
	  updates stay separate.
	  */
	virtual bool PrepareFused(const std::vector<uint64_t>& late) {UNUSED(late); return false;}
	//! Free what PrepareFused() allocated: the main updates stay separate from here on
	virtual void ReleaseFused() {}
	//! The new voltages and the currents except the late ones; the hooks between the main updates see the new voltages and the old currents
	virtual void FusedUpdate() {throw std::logic_error("GPU_Backend: no fused step");}
	//! The late currents, after the hooks between the main updates: the currents are new from here on
	virtual void SeamUpdate() {throw std::logic_error("GPU_Backend: no fused step");}
	//! Run the point operation \a op (of type \a kind, see kernels/ext_ops.h) for its points [begin, end), in parallel or, if \a serial, in order
	virtual void RunPoints(GPU_PointOp kind, const void* op, uint64_t begin, uint64_t end, bool serial) = 0;

	//! Copy the device voltages (\a currents false) or currents to the host mirror \a field
	virtual void Download(bool currents, ArrayLib::ArrayNIJK<FDTD_FLOAT>& field) = 0;

	//! Host pointer to the device voltages (\a currents false) or currents in the ArrayNIJK layout, if the device shares memory with the host, else NULL. Valid after Init().
	virtual FDTD_FLOAT* GetShared(bool currents) const {UNUSED(currents); return NULL;}
	//! Wait until the device finished all work, required before the host accesses shared memory
	virtual void Synchronize() {}

	//! Copy the current fields into the snapshot and return host pointers to it (ArrayNIJK layout), valid until the next snapshot. Returns false if the backend cannot.
	/*!
	  The copy may still run on the device: WaitSnapshot() before reading it.
	  */
	virtual bool SnapshotFields(const FDTD_FLOAT* &volt, const FDTD_FLOAT* &curr) {UNUSED(volt); UNUSED(curr); return false;}
	//! Wait until the copy of the snapshot is done, may be called from another thread
	virtual void WaitSnapshot() {}
	//! Whether the snapshots can copy only the raw values the dumps read, see SetSnapshotGather()
	virtual bool CanSnapshotGather() const {return false;}
	//! Snapshots copy the raw values of these regions instead of the fields
	/*!
	  SnapshotFields() then returns the voltages of \a volt_regions in \a volt and the currents
	  of \a curr_regions in \a curr, region after region (see Engine_GatherRegion); the host
	  interpolates them. Called before the first snapshot. Returns false if the backend
	  cannot (no snapshots then).
	  */
	virtual bool SetSnapshotGather(const std::vector<Engine_GatherRegion>& volt_regions, const std::vector<Engine_GatherRegion>& curr_regions) {UNUSED(volt_regions); UNUSED(curr_regions); return false;}

	//! Sums over timesteps of the weighted raw values of a dump on the device (frequency domain dumps), see Engine_Interface_Base::CreateFieldDFT()
	/*!
	  \param h_field the raw values \a region of the currents, else of the voltages
	  \param count number of sums (weights per AccumulateFieldDFT())
	  The host interpolates the sums, see Field_Gather_FDTD::EvaluateSums().
	  Returns an id, or -1 if the backend cannot.
	  */
	virtual int AddFieldDFT(bool h_field, const Engine_GatherRegion& region, unsigned int count) {UNUSED(h_field); UNUSED(region); UNUSED(count); return -1;}
	//! Add \a weights[n] times the current raw values to sum n of \a id (in order with the updates)
	virtual void AccumulateFieldDFT(int id, const std::vector<std::complex<float>>& weights) {UNUSED(id); UNUSED(weights);}
	//! Sum \a k of \a id, one value per raw value of its region
	virtual bool ReadFieldDFT(int id, unsigned int k, std::vector<std::complex<float>>& sums) {UNUSED(id); UNUSED(k); UNUSED(sums); return false;}

	//! Sums of the squared voltages and currents of the first numNodes[n] nodes in each direction, see Engine_Interface_FDTD::CalcFastEnergy(). Returns false if the backend cannot compute them.
	virtual bool CalcFastEnergy(const unsigned int numNodes[3], double& E_energy, double& H_energy) {UNUSED(numNodes); UNUSED(E_energy); UNUSED(H_energy); return false;}

	//! Create a backend for a sub-grid (e.g. a cylindrical multi-grid level) whose work is ordered with the work of this backend. It's the responsibility of the caller to free it.
	virtual GPU_Backend* NewSubGridBackend() = 0;
};

#endif // GPU_BACKEND_H
