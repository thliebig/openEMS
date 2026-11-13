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

#ifndef GPU_BACKEND_METAL_H
#define GPU_BACKEND_METAL_H

#include <vector>

#include "FDTD/gpu_backend.h"
#include "FDTD/gpu_plugin.h"

//! GPU backend using Apple Metal
/*!
  The fields and the cell model live in shared-storage Metal buffers (unified
  memory), which Engine_GPU uses directly as its host mirror: there is no
  download, Download() only waits for the device. All kernels of a batch of
  timesteps are encoded into one compute encoder, which is only committed when
  the host needs the fields.

  The main updates and the device extensions are the code of every backend on
  the cell model (kernels/cell_update.h, kernels/ext_ops.h, gpu_ext.h), compiled
  at run time with the entry points of metal/metal_kernels.metal. Its structs
  hold device pointers: the gpuAddress of the buffers (Metal 3, macOS 13), which
  every encoder declares resident.

  The field dumps are copied on the device at a snapshot (see SnapshotFields())
  and evaluated on the dump thread while the device keeps stepping; the frequency
  domain dumps are summed there too, where HIP sums them on the device.

  This header is plain C++, the Metal state is kept in Impl (see metal_internal.h).
  */
class GPU_Backend_Metal : public GPU_Backend
{
public:
	//! Create the backend on the default Metal device, or return NULL if there is none
	static GPU_Backend_Metal* New();
	//! This backend and its device (the default one), see FDTD/gpu_plugin.h
	static const GPU_Plugin* Plugin();
	virtual ~GPU_Backend_Metal();

	virtual std::string GetName() const;

	virtual bool Init(const Operator* op, const GPU_CellModel& model);

	virtual void UpdateVoltages(const Engine_Box& box);
	virtual void UpdateCurrents(const Engine_Box& box);

	virtual void* DeviceAlloc(size_t bytes, const void* host=NULL);
	virtual void DeviceDownload(void* host, const void* dev, size_t bytes);
	virtual float* DeviceField(bool currents);
	virtual void RunPoints(GPU_PointOp kind, const void* op, uint64_t begin, uint64_t end, bool serial);

	virtual void Download(bool currents, ArrayLib::ArrayNIJK<FDTD_FLOAT>& field);

	virtual FDTD_FLOAT* GetShared(bool currents) const;
	virtual void Synchronize();
	virtual bool CalcFastEnergy(const unsigned int numNodes[3], double& E_energy, double& H_energy);
	virtual bool SnapshotFields(const FDTD_FLOAT* &volt, const FDTD_FLOAT* &curr);
	virtual void WaitSnapshot();

	virtual GPU_Backend* NewSubGridBackend();

	struct Impl;
	Impl* GetImpl() const {return d;}

protected:
	GPU_Backend_Metal(Impl* impl);

	//! The main voltage or current update of the x lines [x0, x1)
	void Update(bool voltages, unsigned int x0, unsigned int x1);

	Impl* d;
};

#endif // GPU_BACKEND_METAL_H
