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

#ifndef GPU_BACKEND_REFERENCE_H
#define GPU_BACKEND_REFERENCE_H

#include <memory>
#include <vector>

#include "gpu_backend.h"
#include "gpu_cell_model.h"

//! Reference GPU backend running on the CPU
/*!
  Keeps its own "device" copies of the fields and of the cell model, separate
  from the host mirror of Engine_GPU, so that all synchronization paths are
  exercised. It runs the kernels of the device backends (kernels/cell_update.h)
  on the host threads, with the same data, as a device would; within single
  precision of the CPU engine on the same model (its flat sweep, generic
  kernels, fast axis z), whose coefficients are in double per run (CM_FACTOR).

  It needs no GPU and no GPU SDK and is built unconditionally, so Engine_GPU can
  be built and tested anywhere (--engine=gpu-reference). It is never chosen
  by itself: --engine=gpu without a device is an error. It is the oracle of
  the GPU engine test: python/Tests/Engine_Compare.py compares it with the CPU
  engine, and HIP with it bit for bit (Metal within single precision), so a
  difference between the two is in a device backend, not in the engine or
  the model.

  It runs the device extensions and the multi-grid coupling of every backend
  (kernels/ext_ops.h, through DeviceAlloc() and RunPoints()) and sums the energy;
  the field dumps read the host mirror of Engine_GPU (no snapshots, no device
  DFT sums), which is exercised so.
  */
class GPU_Backend_Reference : public GPU_Backend
{
public:
	GPU_Backend_Reference();
	virtual ~GPU_Backend_Reference() {}

	virtual std::string GetName() const {return "reference (CPU)";}

	virtual bool Init(const Operator* op, const GPU_CellModel& model);

	virtual void* DeviceAlloc(size_t bytes, const void* host=NULL);
	virtual void DeviceDownload(void* host, const void* dev, size_t bytes);
	virtual bool CalcFastEnergy(const unsigned int numNodes[3], double& E_energy, double& H_energy);
	virtual float* DeviceField(bool currents) {return currents ? m_curr.data() : m_volt.data();}
	virtual void RunPoints(GPU_PointOp kind, const void* op, uint64_t begin, uint64_t end, bool serial);

	virtual void UpdateVoltages(const Engine_Box& box);
	virtual void UpdateCurrents(const Engine_Box& box);

	virtual void Download(bool currents, ArrayLib::ArrayNIJK<FDTD_FLOAT>& field);

	//! On the CPU there is no work stream to share, a sub-grid just gets its own backend
	virtual GPU_Backend* NewSubGridBackend() {return new GPU_Backend_Reference();}

protected:
	void Update(bool voltages, unsigned int x0, unsigned int x1);

	unsigned int numLines[3];

	// "device" memory, never handed to the engine (no GetShared())
	ArrayLib::ArrayNIJK<FDTD_FLOAT> m_volt;
	ArrayLib::ArrayNIJK<FDTD_FLOAT> m_curr;
	std::unique_ptr<GPU_CellModel> m_Model;
	std::vector<float> m_State[2];   //!< the state of the UPML and dispersive nodes, see CM_Field
	CM_Field m_Field[2];
	std::vector< std::unique_ptr<char[]> > m_Memory;   //!< of DeviceAlloc()
};

#endif // GPU_BACKEND_REFERENCE_H
