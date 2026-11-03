/*
*	Copyright (C) 2026 Sean Mollet (sean@malmoset.com)
*	Copyright (C) 2026 Alexandru Iercosan-Lucaci (alex.lucaci@alacrity.ro)
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

#ifndef GPU_BACKEND_SYCL_H
#define GPU_BACKEND_SYCL_H

#include "FDTD/gpu_backend.h"
#include "FDTD/gpu_plugin.h"

//! GPU backend using SYCL (Intel oneAPI DPC++), for Intel GPUs through Level Zero
/*!
  The fields and coefficients live in device memory (USM device allocations). The
  host mirror of Engine_GPU is a separate copy, synchronized with explicit transfers
  (SYCL cannot page-lock memory it did not allocate, so it stays pageable). All
  kernels of the grids of a simulation run in order on one in-order SYCL queue (see
  NewSubGridBackend()), which is only synchronized when the host needs the fields.
  A second in-order queue in the same context is meant for the snapshot downloads
  (see SYCL_Context::copy_q).

  This header is plain C++, the SYCL state is kept in Impl (see sycl_internal.h).
  */
class GPU_Backend_SYCL : public GPU_Backend
{
public:
	//! Create the backend on SYCL device \a device (a Level Zero GPU), or return NULL if it cannot run there
	static GPU_Backend_SYCL* New(int device);
	//! Number of Level Zero GPU devices, 0 without a runtime or driver
	static int DeviceCount();
	//! Describe SYCL device \a n
	static bool DeviceInfo(int n, GPU_DeviceInfo* info);
	//! This backend and its devices, built into libopenEMS or as a plugin (see FDTD/gpu_plugin.h)
	static const GPU_Plugin* Plugin();
	virtual ~GPU_Backend_SYCL();

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
	virtual GPU_MultiGridLink* CreateMultiGridLink(GPU_Backend* sub_grid, const GPU_MultiGridInterpolation& interpol);

	virtual GPU_Extension* CreateExtension(Engine_Extension* eng_ext, Engine* eng);

	struct Impl;
	Impl* GetImpl() const {return d;}

protected:
	GPU_Backend_SYCL(Impl* impl);

	Impl* d;
};

#endif // GPU_BACKEND_SYCL_H
