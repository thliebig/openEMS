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

// Objective-C++ only: the Metal state of the backend

#ifndef METAL_INTERNAL_H
#define METAL_INTERNAL_H

#import <Metal/Metal.h>

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "gpu_backend_metal.h"
#include "FDTD/kernels/cell_update.h"

//! Grid size (number of mesh lines)
struct Metal_GridDim
{
	uint32_t nx, ny, nz;
};

//! Device, work stream and kernels, shared by the backends of all grids of a simulation (see NewSubGridBackend())
struct Metal_Context
{
	id<MTLDevice> device;
	id<MTLCommandQueue> queue;
	id<MTLCommandBuffer> cmd;
	id<MTLComputeCommandEncoder> enc;

	std::map<std::string, id<MTLComputePipelineState>> pipelines;
	id<MTLLibrary> library;   //!< all kernels, see Pipeline()

	//! the buffers the kernels reach through pointers in their arguments (the cell model, the device extensions):
	//! every encoder declares them resident (useResources), see Encoder()
	std::vector< id<MTLBuffer> > resources;
	//! the gpuAddress of every buffer of DeviceAlloc() and its buffer, for DeviceDownload()
	std::map<uint64_t, id<MTLBuffer>> addresses;

	//! The pipeline of kernel \a function, the kernels compiled at the first call (metal_cm_source.h, see FDTD/CMakeLists.txt)
	id<MTLComputePipelineState> Pipeline(const char* function);

	//! The encoder of the current batch, created on demand
	id<MTLComputeCommandEncoder> Encoder();

	//! Commit the current batch and wait for it to finish
	void Flush();

	//! Commit the current batch without waiting, returns its command buffer (nil if empty)
	id<MTLCommandBuffer> Commit();
};


//! State of one grid
/*!
  ARC owns the Metal objects below: they are released with the Impl, which the backend deletes.
  */
struct GPU_Backend_Metal::Impl
{
	std::shared_ptr<Metal_Context> ctx;   //!< shared with the sub-grid backends: one queue, so all grids run in order
	id<MTLDevice> device;

	Metal_GridDim dim;
	size_t numCells;    //!< nx*ny*nz, the field buffers hold 3*numCells values

	id<MTLBuffer> volt, curr;
	//! The main updates on the cell model (see Init()): per field the node index (16 or 32 bit) and the kernel arguments
	bool cm_index16 = true;
	id<MTLBuffer> cm_index[2];
	CM_Field cm_field[2];

	id<MTLBuffer> energy;           //!< per-line energy sums, see GPU_Backend_Metal::CalcFastEnergy()
	id<MTLBuffer> snap_volt, snap_curr;         //!< the field snapshot, see GPU_Backend_Metal::SnapshotFields()
	id<MTLCommandBuffer> snap_cmd;              //!< the command buffer copying it
	bool snap_failed;                           //!< not enough memory for the snapshots: the dumps read the fields


	id<MTLComputePipelineState> Pipeline(const char* function) {return ctx->Pipeline(function);}
	id<MTLComputeCommandEncoder> Encoder() {return ctx->Encoder();}
	void Flush() {ctx->Flush();}

	//! Buffer of \a bytes in shared storage (host and device see the same memory), initialized with \a data or zero
	id<MTLBuffer> NewBuffer(size_t bytes, const void* data=NULL);
	//! NewBuffer() that the kernels reach through a pointer: resident in every encoder, its gpuAddress returned
	uint64_t NewResource(size_t bytes, const void* data=NULL);

	//! Dispatch \a pso with one thread per (i,j,k), i fastest; the pipeline and its arguments must be set on Encoder()
	/*!
	  A dispatch sees the results of all previous dispatches (a barrier before it).
	  */
	void Dispatch(id<MTLComputePipelineState> pso, size_t ni, size_t nj=1, size_t nk=1);

};

#endif // METAL_INTERNAL_H
