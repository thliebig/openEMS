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

#ifndef ENGINE_GPU_H
#define ENGINE_GPU_H

#include "engine.h"
#include "gpu_backend.h"

#include <string>
#include <unordered_map>
#include <vector>

class GPU_Backend;
class GPU_Extension;

//! FDTD engine running the main updates on a GPU backend
/*!
  The fields live on the device (see GPU_Backend). The base class arrays
  (volt_ptr, curr_ptr) serve as host mirror in the basic engine layout, read by
  the field processing through GetVolt/GetCurr and by the engine extensions.
  If the device shares its memory with the host, the host mirror is a view of
  the device memory and no copies are needed.

  An extension runs on the device if its operator extension IsGPUSave() and the
  backend provides a device implementation (GPU_Backend::CreateExtension()).
  Every concrete operator extension of openEMS reports IsGPUSave(), so a fallback
  means the backend has no implementation of one of them (the reference backend
  has none at all).
  If all extensions do, the fields stay on the device for a whole IterateTS()
  call. Without shared memory, the host mirror is then not updated at the end of
  a batch: GetVolt/GetCurr read the z-line of a value from the device (e.g. for
  a probe), or the whole field if a batch reads many lines (e.g. for a dump).

  Otherwise all extensions run on the host mirror, and the engine uploads both fields before and downloads the updated field
  after each main update (host fallback, correct but slow). Note that the pre-
  and post-update hooks may write the fields too, despite the Engine_Extension
  documentation: the UPML in all four, the TF/SF and the cylinder extension in
  the post-update hooks.

  The engine type is Engine::GPU, which the extension dispatcher and
  Engine_Interface_FDTD handle like the basic engine, because that is the layout
  of the host mirror. The engine owns the backend and the device extensions and
  deletes the extensions first.
  */
class Engine_GPU : public Engine
{
public:
	//! Create the engine for \a op on the GPU backend \a backend (see GPU_Backend::New), or for a sub-grid on a backend sharing the work stream of \a parent
	static Engine_GPU* New(const Operator* op, const std::string& backend, GPU_Backend* parent=NULL);
	virtual ~Engine_GPU();

	//! Create the backend, the device extensions and decide device or host fallback
	virtual void Init();
	//! Free the device extensions before the backend they belong to
	virtual void Reset();

	//! A batch of timesteps, the host mirror is brought up to date at its end
	virtual bool IterateTS(unsigned int iterTS);

	//! One voltage half-step: extensions and main update (with sub-grids, see Engine_GPU_CylinderMultiGrid)
	virtual void VoltageHalfStep();
	//! One current half-step: extensions and main update
	virtual void CurrentHalfStep();
	//! Count the completed timestep
	virtual void NextTimestep();
	//! After a batch of timesteps: update the host mirror for the field processing
	virtual void FinishBatch();

	//! Run all extensions on the host mirror, even if they have a device implementation
	virtual void ForceHostFallback();

	//! The device backend holding the fields
	GPU_Backend* GetBackend() const {return m_Backend;}

	//! True if the extensions run on the host mirror, which is then kept in sync during IterateTS()
	bool FieldsOnHost() const {return m_FieldsOnHost;}

	//! Field energy sums on the device (see GPU_Backend::CalcFastEnergy()), false if the host has to compute them
	bool CalcFastEnergy(const unsigned int numNodes[3], double& E_energy, double& H_energy) const;

	//! Host mirror of the fields (basic engine layout)
	ArrayLib::ArrayNIJK<FDTD_FLOAT>& HostVoltages() const {return *volt_ptr;}
	ArrayLib::ArrayNIJK<FDTD_FLOAT>& HostCurrents() const {return *curr_ptr;}

	//! Make the device results visible in the host mirror
	void VoltagesToHost();
	void CurrentsToHost();
	//! Update the host mirror if it is out of date (see GetVolt())
	void UpdateHostMirror();

	//! Copy of the current fields in snapshot \a slot (0 or 1), see GPU_Backend::SnapshotFields(); false if not available
	bool SnapshotFields(unsigned int slot, const FDTD_FLOAT* &volt, const FDTD_FLOAT* &curr);
	//! Wait until snapshot \a slot can be read, may be called from another thread
	void WaitSnapshot(unsigned int slot) const;
	//! Whether the snapshots evaluate the dumped nodes on the device (see AddSnapshotGather()), instead of copying the fields
	bool CanSnapshotGather() const;
	//! Evaluate the dump entries \a entries of the E (\a h_field false) or H field at every snapshot
	/*!
	  The snapshot then holds their values from \a offset on, in the voltages (E) or currents (H)
	  returned by SnapshotFields(). Only before the first snapshot, returns false otherwise.
	  */
	bool AddSnapshotGather(bool h_field, const std::vector<GPU_GatherEntry>& entries, size_t& offset);

	//! Sums of weighted dump entries on the device, see GPU_Backend::AddFieldDFT(); -1 if not available
	int AddFieldDFT(bool h_field, const std::vector<GPU_GatherEntry>& entries, unsigned int count);
	void AccumulateFieldDFT(int id, const std::vector<std::complex<float>>& weights) {m_Backend->AccumulateFieldDFT(id, weights);}
	bool ReadFieldDFT(int id, std::vector<std::complex<float>>& sums) {return m_Backend->ReadFieldDFT(id, sums);}

	//! Field values, read from the device if the host mirror is out of date (see class description)
	virtual FDTD_FLOAT GetVolt(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const;
	virtual FDTD_FLOAT GetVolt(unsigned int n, const unsigned int pos[3]) const {return GetVolt(n, pos[0], pos[1], pos[2]);}
	virtual FDTD_FLOAT GetCurr(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const;
	virtual FDTD_FLOAT GetCurr(unsigned int n, const unsigned int pos[3]) const {return GetCurr(n, pos[0], pos[1], pos[2]);}
	//! Make the host mirror changes visible to the device
	void VoltagesToDevice();
	void CurrentsToDevice();

protected:
	Engine_GPU(const Operator* op, const std::string& backend, GPU_Backend* parent);

	std::string m_BackendName;
	GPU_Backend* m_ParentBackend;
	GPU_Backend* m_Backend;
	bool m_FieldsOnHost;
	bool m_SharedMemory; //!< the host mirror is the device memory
	bool m_SnapshotsUsed; //!< SnapshotFields() succeeded once
	std::vector<GPU_GatherEntry> m_SnapshotGather[2];   //!< dump entries of the snapshots: E, H (see AddSnapshotGather())
	bool m_SnapshotGatherSet;                           //!< m_SnapshotGather was passed to the backend
	bool m_FieldDFTUsed;                                //!< AddFieldDFT() succeeded once

	//! device implementations of the engine extensions, same order as m_Eng_exts (fast path only)
	std::vector<GPU_Extension*> m_GPU_exts;

	//! Out-of-date host mirror of one field (see GetVolt()): the z-lines read from the device in this batch
	struct StaleField
	{
		bool stale;
		bool full_this_batch;   //!< this batch read the whole field
		bool full_last_batch;   //!< the last batch read the whole field, so will this one
		std::unordered_map<unsigned int, std::vector<FDTD_FLOAT>> lines;
	};
	mutable StaleField m_StaleVolt, m_StaleCurr;
	//! Read a value of a stale field from the device, keeping its z-line for the rest of the batch
	FDTD_FLOAT ReadStale(bool currents, unsigned int n, unsigned int x, unsigned int y, unsigned int z) const;
	//! The device updated the field, the host mirror is out of date from here on
	void MarkStale(StaleField& field);

	void ClearGPUExtensions();
};

#endif // ENGINE_GPU_H
