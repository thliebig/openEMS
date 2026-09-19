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

#include <string>
#include <vector>

#include "tools/constants.h"
#include "tools/global.h"
#include "tools/arraylib/array_nijk.h"

class Operator;
class Engine;
class Engine_Extension;

//! Device implementation of an engine extension, see GPU_Backend::CreateExtension()
/*!
  Engine_GPU calls it through the same hooks and in the same order as the engine
  extension it replaces, while the fields stay on the device. It works on the
  device data of the backend that created it.
  */
class GPU_Extension
{
public:
	virtual ~GPU_Extension() {}

	virtual void DoPreVoltageUpdates() {}
	virtual void DoPostVoltageUpdates() {}
	virtual void Apply2Voltages() {}

	virtual void DoPreCurrentUpdates() {}
	virtual void DoPostCurrentUpdates() {}
	virtual void Apply2Current() {}

	//! Called after each IterateTS() of the engine, with the host mirror up to date
	virtual void Synchronize() {}
};

//! Coupling of a cylindrical multi-grid level and its sub-grid, see Engine_GPU_CylinderMultiGrid
/*!
  Interpolation data of Operator_CylinderMultiGrid, independent of the operator
  classes. The sub-grid has the same r and z lines as the base grid (up to its
  size) and every other alpha line. For each alpha line a of the base grid and
  direction n (0: r and z components, 1: alpha component), a field of the base
  grid is interpolated from the sub-grid alpha lines pos_2p[n][a] and pos_2pp[n][a]
  with the weights w_2p[n][a] and w_2pp[n][a] (_v: voltages, _i: currents).
  */
struct GPU_MultiGridInterpolation
{
	unsigned int split_pos;   //!< first base grid r line not covered by the sub-grid
	std::vector<unsigned int> pos_v_2p[2], pos_v_2pp[2], pos_i_2p[2], pos_i_2pp[2];
	std::vector<float> w_v_2p[2], w_v_2pp[2], w_i_2p[2], w_i_2pp[2];
};

//! Device implementation of the coupling of a multi-grid level and its sub-grid, see GPU_Backend::CreateMultiGridLink()
class GPU_MultiGridLink
{
public:
	virtual ~GPU_MultiGridLink() {}

	//! Base grid voltages at r line split_pos-1 to the sub-grid, after the voltage updates of both grids
	virtual void SyncVoltages() = 0;
	//! Sub-grid currents to the base grid at r line split_pos-2, after the current updates of both grids
	virtual void SyncCurrents() = 0;
	//! Interpolate the sub-grid fields to the base grid inside the sub-grid (r lines below split_pos-1 / split_pos-2)
	virtual void InterpolateToBase() = 0;
};

//! Nodes of a UPML region on the device: start node and number of nodes per direction
struct GPU_UPMLRegion
{
	unsigned int start[3];
	unsigned int size[3];
};

//! Whether the UPML of \a eng can run fused with the main updates, and the node box [start, stop) the main updates then cover
/*!
  Fused kernels do the UPML pre-update, the main update and the UPML post-update of
  the nodes of a region in one pass. This is the same as the separate steps if:
  - the UPML hooks run directly before and after the main update: only extensions
    without device pre/post hooks (steady-state) come before the UPML extensions,
  - all UPML extensions of the engine run on the device (\a regions has one entry per UPML extension),
  - the regions cover exactly the nodes outside a box, which the main updates then cover.
  A voltage/current update only reads its own node of the updated field, so the
  regions and the box can be updated in any order.
  */
bool GPU_UPMLFusionBox(Engine* eng, const std::vector<GPU_UPMLRegion>& regions, const unsigned int numLines[3],
                       unsigned int start[3], unsigned int stop[3]);

//! Abstract interface to the device used by Engine_GPU
/*!
  A backend owns the device memory: the voltage and current fields and the update
  coefficients of the operator. It performs the main FDTD updates on the device.

  Engine_GPU keeps a host mirror of the fields in the basic Engine layout (ArrayNIJK),
  used by the field processing and by extensions without a device implementation,
  and synchronizes it with the Download/Upload methods. If the device shares its
  memory with the host (unified memory), the host mirror is the device memory
  itself (see GetSharedVoltages()) and synchronizing means waiting for the device.
  */
class GPU_Backend
{
public:
	//! Create a backend by name: "reference" or "auto" (the best device backend available). It's the responsibility of the caller to free it.
	static GPU_Backend* New(const std::string& name);

	virtual ~GPU_Backend() {}

	virtual std::string GetName() const = 0;

	//! Allocate the device fields (all zero) and upload the update coefficients of \a op
	virtual bool Init(const Operator* op) = 0;

	//! Main voltage update over the whole domain, equivalent to Engine::UpdateVoltages
	virtual void UpdateVoltages() = 0;
	//! Main current update over the whole domain, equivalent to Engine::UpdateCurrents
	virtual void UpdateCurrents() = 0;

	//! Copy the device voltages to the host mirror
	virtual void DownloadVoltages(ArrayLib::ArrayNIJK<FDTD_FLOAT>& volt) = 0;
	//! Copy the device currents to the host mirror
	virtual void DownloadCurrents(ArrayLib::ArrayNIJK<FDTD_FLOAT>& curr) = 0;
	//! Copy the host mirror voltages to the device
	virtual void UploadVoltages(const ArrayLib::ArrayNIJK<FDTD_FLOAT>& volt) = 0;
	//! Copy the host mirror currents to the device
	virtual void UploadCurrents(const ArrayLib::ArrayNIJK<FDTD_FLOAT>& curr) = 0;
	//! Copy \a count values of the device voltages or \a currents from the flat ArrayNIJK index \a offset to \a dst, returns false if the backend cannot
	virtual bool DownloadRange(bool currents, size_t offset, size_t count, FDTD_FLOAT* dst) {UNUSED(currents); UNUSED(offset); UNUSED(count); UNUSED(dst); return false;}

	//! Host pointer to the device voltages in the ArrayNIJK layout, if the device shares memory with the host, else NULL. Valid after Init().
	virtual FDTD_FLOAT* GetSharedVoltages() const {return NULL;}
	//! Host pointer to the device currents in the ArrayNIJK layout, if the device shares memory with the host, else NULL. Valid after Init().
	virtual FDTD_FLOAT* GetSharedCurrents() const {return NULL;}
	//! Wait until the device finished all work, required before the host accesses shared memory
	virtual void Synchronize() {}

	//! Copy the current fields into snapshot \a slot (0 or 1) and return host pointers to them (ArrayNIJK layout), valid until the next snapshot into the same slot. Returns false if the backend cannot.
	virtual bool SnapshotFields(unsigned int slot, const FDTD_FLOAT* &volt, const FDTD_FLOAT* &curr) {UNUSED(slot); UNUSED(volt); UNUSED(curr); return false;}

	//! Sums of the squared voltages and currents of the first numNodes[n] nodes in each direction, see Engine_Interface_FDTD::CalcFastEnergy(). Returns false if the backend cannot compute them.
	virtual bool CalcFastEnergy(const unsigned int numNodes[3], double& E_energy, double& H_energy) {UNUSED(numNodes); UNUSED(E_energy); UNUSED(H_energy); return false;}

	//! Create a backend for a sub-grid (e.g. a cylindrical multi-grid level) whose work is ordered with the work of this backend. It's the responsibility of the caller to free it.
	virtual GPU_Backend* NewSubGridBackend() = 0;
	//! Create the device coupling of this (base grid) backend and \a sub_grid (created by NewSubGridBackend()), or NULL if this backend has none
	virtual GPU_MultiGridLink* CreateMultiGridLink(GPU_Backend* sub_grid, const GPU_MultiGridInterpolation& interpol) {UNUSED(sub_grid); UNUSED(interpol); return NULL;}

	//! Create the device implementation of the engine extension \a eng_ext of engine \a eng, or NULL if this backend has none
	virtual GPU_Extension* CreateExtension(Engine_Extension* eng_ext, Engine* eng) {UNUSED(eng_ext); UNUSED(eng); return NULL;}
};

#endif // GPU_BACKEND_H
