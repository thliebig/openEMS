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

// Objective-C++ only: shared by the Metal backend and the Metal extensions

#ifndef METAL_INTERNAL_H
#define METAL_INTERNAL_H

#import <Metal/Metal.h>

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "gpu_backend_metal.h"

//! Grid size as passed to the kernels (number of mesh lines), see METAL_COMMON_SOURCE
struct Metal_GridDim
{
	uint32_t nx, ny, nz;
};

//! Common kernel code, prepended to every kernel source
/*!
  Fields and coefficients use the ArrayNIJK layout of the host:
  index(n,x,y,z) = n*nx*ny*nz + x*ny*nz + y*nz + z
  */
#define METAL_COMMON_SOURCE R"MSL(
#include <metal_stdlib>
using namespace metal;

struct GridDim { uint nx; uint ny; uint nz; };

inline uint nijk(constant GridDim& N, uint n, uint x, uint y, uint z)
{
	return ((n*N.nx + x)*N.ny + y)*N.nz + z;
}
)MSL"

//! Device, work stream and kernels, shared by the backends of all grids of a simulation (see NewSubGridBackend())
struct Metal_Context
{
	id<MTLDevice> device;
	id<MTLCommandQueue> queue;
	id<MTLCommandBuffer> cmd;
	id<MTLComputeCommandEncoder> enc;

	std::map<std::string, id<MTLComputePipelineState>> pipelines;

	//! Compile (once) and return the pipeline of kernel \a function in \a source (without METAL_COMMON_SOURCE)
	id<MTLComputePipelineState> Pipeline(const char* source, const char* function);

	//! The encoder of the current batch, created on demand
	id<MTLComputeCommandEncoder> Encoder();

	//! Commit the current batch and wait for it to finish
	void Flush();
};

//! State of one grid
struct GPU_Backend_Metal::Impl
{
	std::shared_ptr<Metal_Context> ctx;
	id<MTLDevice> device;

	Metal_GridDim dim;
	size_t numCells;    //!< nx*ny*nz, the field buffers hold 3*numCells values

	id<MTLBuffer> volt, curr;
	id<MTLBuffer> vv, vi, ii, iv;   //!< full coefficients, if not compressed
	id<MTLBuffer> index, coeff;     //!< compressed coefficients: set index per node, sets (see update_voltages_c)
	id<MTLBuffer> energy;           //!< per-line energy sums, see GPU_Backend_Metal::CalcFastEnergy()

	id<MTLComputePipelineState> Pipeline(const char* source, const char* function) {return ctx->Pipeline(source, function);}
	id<MTLComputeCommandEncoder> Encoder() {return ctx->Encoder();}
	void Flush() {ctx->Flush();}

	//! Shared buffer of \a bytes, initialized with \a data or zero
	id<MTLBuffer> NewBuffer(size_t bytes, const void* data=NULL);

	//! Dispatch \a pso with one thread per (i,j,k), i fastest; the pipeline and its arguments must be set on Encoder()
	void Dispatch(id<MTLComputePipelineState> pso, size_t ni, size_t nj=1, size_t nk=1);

	//! Bind the grid dimension to \a index
	void SetGridDim(unsigned int index);
};

//! Factory of a Metal extension: the device implementation of \a eng_ext, or NULL if \a eng_ext is not of its type
typedef GPU_Extension* (*Metal_ExtensionFactory)(GPU_Backend_Metal::Impl* d, Engine_Extension* eng_ext, Engine* eng);

GPU_Extension* Metal_CreateExt_Excitation(GPU_Backend_Metal::Impl* d, Engine_Extension* eng_ext, Engine* eng);
GPU_Extension* Metal_CreateExt_UPML(GPU_Backend_Metal::Impl* d, Engine_Extension* eng_ext, Engine* eng);
GPU_Extension* Metal_CreateExt_Mur_ABC(GPU_Backend_Metal::Impl* d, Engine_Extension* eng_ext, Engine* eng);
GPU_Extension* Metal_CreateExt_LorentzMaterial(GPU_Backend_Metal::Impl* d, Engine_Extension* eng_ext, Engine* eng);
GPU_Extension* Metal_CreateExt_LumpedRLC(GPU_Backend_Metal::Impl* d, Engine_Extension* eng_ext, Engine* eng);
GPU_Extension* Metal_CreateExt_TFSF(GPU_Backend_Metal::Impl* d, Engine_Extension* eng_ext, Engine* eng);
GPU_Extension* Metal_CreateExt_Absorbing_BC(GPU_Backend_Metal::Impl* d, Engine_Extension* eng_ext, Engine* eng);
GPU_Extension* Metal_CreateExt_SteadyState(GPU_Backend_Metal::Impl* d, Engine_Extension* eng_ext, Engine* eng);
GPU_Extension* Metal_CreateExt_Cylinder(GPU_Backend_Metal::Impl* d, Engine_Extension* eng_ext, Engine* eng);

#endif // METAL_INTERNAL_H
