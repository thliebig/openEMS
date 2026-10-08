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

  This header is plain C++, the HIP state is kept in Impl (see hip_internal.h).
  */
class GPU_Backend_HIP : public GPU_Backend
{
public:
	//! Create the backend on the current HIP device, or return NULL if there is none
	static GPU_Backend_HIP* New();
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
	virtual bool SnapshotFields(unsigned int slot, const FDTD_FLOAT* &volt, const FDTD_FLOAT* &curr);
	virtual void WaitSnapshot(unsigned int slot);
	virtual bool CanSnapshotGather() const {return true;}
	virtual bool SetSnapshotGather(const std::vector<GPU_GatherEntry>& volt_entries, const std::vector<GPU_GatherEntry>& curr_entries);
	virtual int AddFieldDFT(const std::vector<GPU_GatherEntry>& entries, bool currents, unsigned int count);
	virtual void AccumulateFieldDFT(int id, const std::vector<std::complex<float>>& weights);
	virtual bool ReadFieldDFT(int id, std::vector<std::complex<float>>& sums);
	virtual bool CalcFastEnergy(const unsigned int numNodes[3], double& E_energy, double& H_energy);

	virtual GPU_Backend* NewSubGridBackend();
	virtual GPU_MultiGridLink* CreateMultiGridLink(GPU_Backend* sub_grid, const GPU_MultiGridInterpolation& interpol);

	virtual GPU_Extension* CreateExtension(Engine_Extension* eng_ext, Engine* eng);

	struct Impl;
	Impl* GetImpl() const {return d;}

protected:
	GPU_Backend_HIP(Impl* impl);

	Impl* d;
};

#endif // GPU_BACKEND_HIP_H
