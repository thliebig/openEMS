/*
*	Copyright (C) 2026 Thorsten Liebig (Thorsten.Liebig@gmx.de)
*	Copyright (C) 2026 Alexandru Iercosan-Lucaci (alex.lucaci@alacrity.ro)
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

#ifndef GPU_BACKEND_SYCL_H
#define GPU_BACKEND_SYCL_H

#include "FDTD/gpu_backend.h"
#include "FDTD/gpu_plugin.h"

//! GPU backend using SYCL, for Intel GPUs (Level Zero); a port of the HIP backend (FDTD/hip)
/*!
  The fields and the cell model live in device memory. The host mirror of
  Engine_GPU is a separate copy, synchronized with explicit transfers. All
  kernels of the grids of a simulation run in order on one in-order queue (see
  NewSubGridBackend()), which is only waited for when the host needs the fields.
  Where the extensions allow it, the voltages and currents of a timestep can be
  one fused kernel (see PrepareFused(), sycl_cell_model.cpp).
  The field dumps are evaluated on the device at a snapshot and downloaded on a
  second queue, overlapping the work (see SnapshotFields()). Frequency domain
  dumps are summed on the device (see AddFieldDFT()).

  The devices are the Level Zero GPUs; with ONEAPI_DEVICE_SELECTOR set (SYCL's
  own variable) every device it selects instead, e.g. opencl:cpu to run the
  kernels on the CPU for tests.

  This header is plain C++, the SYCL state is kept in Impl (see sycl_internal.h).
  */
class GPU_Backend_SYCL : public GPU_Backend
{
public:
	//! Create the backend on SYCL device \a device (see DeviceCount()), or return NULL if it cannot run there
	static GPU_Backend_SYCL* New(int device);
	//! Number of devices, 0 without a runtime or driver
	static int DeviceCount();
	//! Describe device \a n
	static bool DeviceInfo(int n, GPU_DeviceInfo* info);
	//! This backend and its devices, as a plugin (see FDTD/gpu_plugin.h)
	static const GPU_Plugin* Plugin();
	virtual ~GPU_Backend_SYCL();

	virtual std::string GetName() const;
	virtual size_t CacheBytes(std::string& source) const;

	virtual bool Init(const Operator* op, const GPU_CellModel& model);
	virtual void* DeviceAlloc(size_t bytes, const void* host=NULL);
	virtual void DeviceDownload(void* host, const void* dev, size_t bytes);
	virtual float* DeviceField(bool currents);
	virtual void RunPoints(GPU_PointOp kind, const void* op, uint64_t begin, uint64_t end, bool serial, const CM_PointFilter& filter);

	virtual void UpdateVoltages(const Engine_Box& box);
	virtual void UpdateCurrents(const Engine_Box& box);
	virtual bool PrepareFused(const std::vector<uint64_t>& late);
	virtual void ReleaseFused();
	virtual void FusedUpdate();
	virtual void SeamUpdate();

	virtual void Download(bool currents, ArrayLib::ArrayNIJK<FDTD_FLOAT>& field);

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

	struct Impl;
	Impl* GetImpl() const {return d;}

protected:
	GPU_Backend_SYCL(Impl* impl);

	//! The grid dimensions and the device fields (all zero), the first part of Init()
	bool InitFields(const Operator* op);

	Impl* d;
};

#endif // GPU_BACKEND_SYCL_H
