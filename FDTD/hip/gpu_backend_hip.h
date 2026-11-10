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

#ifndef GPU_BACKEND_HIP_H
#define GPU_BACKEND_HIP_H

#include "FDTD/gpu_backend.h"
#include "FDTD/gpu_plugin.h"

//! GPU backend using HIP, for AMD GPUs and for NVIDIA GPUs through CUDA
/*!
  The fields and coefficients live in device memory. The host mirror of
  Engine_GPU is a separate copy, synchronized with explicit transfers (the host
  arrays are page-locked on first use for full transfer speed). All kernels of
  the grids of a simulation run in order on one HIP stream (see
  NewSubGridBackend()), which is only synchronized when the host needs the fields.
  Where the extensions allow it, a timestep is a single fused kernel with the UPML
  regions folded into it (see Impl::DecideFusedStep()).
  The field dumps are evaluated on the device at a snapshot and downloaded on a
  second stream, overlapping the work (see SnapshotFields()). Frequency domain
  dumps are summed on the device (see AddFieldDFT()).

  The grid can be split along x into parts on several devices (or several parts on
  one), see New(const int*, int) and Init(). Each part has its own state (Impl) and
  streams, and the parts exchange the field lines at their cuts (see ExchangeHalos()).
  With more than one part there is no fused step and no multi-grid, and only the
  extensions with a device implementation per part run on the device (see
  CreateExtension()). The snapshots and sums of the dumps are taken per part, of the
  x lines it owns, and joined on the host (see SetSnapshotGatherParts()).

  This header is plain C++, the HIP state is kept in Impl (see hip_internal.h).
  */
class GPU_Backend_HIP : public GPU_Backend
{
public:
	//! Create the backend on HIP device \a device, or return NULL if it cannot run there
	static GPU_Backend_HIP* New(int device);
	//! Create the backend with the grid split into \a count parts on the HIP \a devices (a device may repeat), or return NULL if it cannot run there
	static GPU_Backend_HIP* New(const int* devices, int count);
	//! Number of HIP devices, 0 without a runtime or driver
	static int DeviceCount();
	//! Describe HIP device \a n
	static bool DeviceInfo(int n, GPU_DeviceInfo* info);
	//! This backend and its devices, built into libopenEMS or as a plugin (see FDTD/gpu_plugin.h)
	static const GPU_Plugin* Plugin();
	virtual ~GPU_Backend_HIP();

	virtual std::string GetName() const;

	virtual bool Init(const Operator* op);

	virtual void UpdateVoltages();
	virtual void UpdateCurrents();

	virtual void DownloadVoltages(ArrayLib::ArrayNIJK<FDTD_FLOAT>& volt);
	virtual void DownloadCurrents(ArrayLib::ArrayNIJK<FDTD_FLOAT>& curr);
	virtual void UploadVoltages(const ArrayLib::ArrayNIJK<FDTD_FLOAT>& volt);
	virtual void UploadCurrents(const ArrayLib::ArrayNIJK<FDTD_FLOAT>& curr);

	virtual bool DownloadRange(bool currents, size_t offset, size_t count, FDTD_FLOAT* dst);
	virtual void Synchronize();
	virtual bool SnapshotFields(const FDTD_FLOAT* &volt, const FDTD_FLOAT* &curr);
	virtual void WaitSnapshot();
	virtual bool CanSnapshotGather() const {return true;}
	virtual bool SetSnapshotGather(const std::vector<Engine_GatherRegion>& volt_regions, const std::vector<Engine_GatherRegion>& curr_regions);
	virtual int AddFieldDFT(bool h_field, const Engine_GatherRegion& region, unsigned int count);
	virtual void AccumulateFieldDFT(int id, const std::vector<std::complex<float>>& weights);
	virtual bool ReadFieldDFT(int id, unsigned int k, std::vector<std::complex<float>>& sums);
	virtual bool CalcFastEnergy(const unsigned int numNodes[3], double& E_energy, double& H_energy);

	virtual GPU_Backend* NewSubGridBackend();
	virtual GPU_MultiGridLink* CreateMultiGridLink(GPU_Backend* sub_grid, const MultiGridInterpolation& interpol);

	virtual GPU_Extension* CreateExtension(Engine_Extension* eng_ext, Engine* eng);

	struct Impl;
	Impl* GetImpl() const {return d;}

protected:
	GPU_Backend_HIP(const std::vector<Impl*>& parts);

	//! Allocate and upload the part \a p, owning the global x lines [x0, x1)
	bool InitPart(Impl* p, const Operator* op, const unsigned int numLines[3], unsigned int x0, unsigned int x1);
	//! Copy the field lines at the cuts between the parts that the next update reads
	void ExchangeHalos(bool currents);
	//! Make the sending parts of ExchangeHalos() wait for the copies before they write the lines again
	void ReleaseHalos(bool currents);
	//! Copy a field of all parts from or to the host mirror \a host of \a size values (several parts only)
	void CopyParts(bool currents, bool upload, float* host, size_t size);
	//! SetSnapshotGather() with several parts: each part gathers the raw values on the x lines it owns
	bool SetSnapshotGatherParts(const std::vector<Engine_GatherRegion>& volt_regions, const std::vector<Engine_GatherRegion>& curr_regions);
	//! SnapshotFields() with several parts
	bool SnapshotFieldsParts(const FDTD_FLOAT* &volt, const FDTD_FLOAT* &curr);

	Impl* d;                     //!< the first part
	std::vector<Impl*> m_Parts;  //!< the parts along x, at least one
	bool m_SnapshotParts;        //!< several parts: the snapshots are set up (see SetSnapshotGatherParts())
	size_t m_SnapshotNV;         //!< several parts: the raw voltages of a snapshot (of all parts)
};

#endif // GPU_BACKEND_HIP_H
