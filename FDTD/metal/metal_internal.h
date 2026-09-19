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

#include <functional>
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

// Main update coefficients of node i, see GPU_Backend_Metal::Impl::coeff_mode:
// a[n*s] and b[n*s] are vv/vi (voltages, set_offset 0) or ii/iv (currents, set_offset 6)
// of direction n. Mode 0: ca/cb are the full arrays, else ca is the table of sets.
struct MainCoeff { const device float* a; const device float* b; uint s; };

inline MainCoeff main_coeff(const device void* index, const device float* ca, const device float* cb,
                            uint mode, uint sn, uint i, uint set_offset)
{
	MainCoeff c;
	if (mode==0)
	{
		c.a = ca + i; c.b = cb + i; c.s = sn;
		return c;
	}
	const uint set = (mode==1) ? uint(((const device ushort*)index)[i]) : ((const device uint*)index)[i];
	c.a = ca + 12*set + set_offset; c.b = c.a + 3; c.s = 1;
	return c;
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

	//! group of the last dispatch, see GPU_Backend_Metal::Impl::Dispatch()
	const void* group;

	//! Compile (once) and return the pipeline of kernel \a function in \a source (without METAL_COMMON_SOURCE)
	id<MTLComputePipelineState> Pipeline(const char* source, const char* function);

	//! The encoder of the current batch, created on demand
	id<MTLComputeCommandEncoder> Encoder();

	//! Commit the current batch and wait for it to finish
	void Flush();
};

//! Coefficients stored as distinct sets and a set index per item, see Metal_FindSets()
#define METAL_MAX_SET_WIDTH 18
struct Metal_CoeffSets
{
	uint32_t mode;                //!< 1: 16 bit index, 2: 32 bit index
	size_t count;                 //!< number of sets
	std::vector<uint32_t> index;  //!< set of each item
	std::vector<float> table;     //!< the sets, \a width values each
};

//! Find the distinct sets of \a width values of \a count items, \a get(i, values) returns the values of item i
/*!
  Returns false if a set index plus the sets would be larger than half the full
  arrays (\a width floats per item). 16 bit indices for up to 65536 sets, else 32 bit.
  */
bool Metal_FindSets(size_t count, unsigned int width, const std::function<void(size_t, float*)>& get, Metal_CoeffSets& sets);

//! State of one grid
struct GPU_Backend_Metal::Impl
{
	std::shared_ptr<Metal_Context> ctx;
	id<MTLDevice> device;

	Metal_GridDim dim;
	size_t numCells;    //!< nx*ny*nz, the field buffers hold 3*numCells values

	id<MTLBuffer> volt, curr;
	//! Main update coefficients: 0: full arrays vv, vi, ii, iv; 1/2: a 16/32 bit set index per node and the sets (see main_coeff())
	uint32_t coeff_mode;
	id<MTLBuffer> vv, vi, ii, iv;
	id<MTLBuffer> index, coeff;
	id<MTLBuffer> energy;           //!< per-line energy sums, see GPU_Backend_Metal::CalcFastEnergy()

	//! The main updates cover the nodes in [main_start, main_stop), the fused UPML kernels the others (see metal_ext_upml.mm)
	Metal_GridDim main_start, main_stop;
	//! UPML extensions of this grid and whether they run fused with the main updates (-1: not decided yet), see metal_ext_upml.mm
	std::vector<GPU_Extension*> upml;
	int upml_fused;

	id<MTLComputePipelineState> Pipeline(const char* source, const char* function) {return ctx->Pipeline(source, function);}
	id<MTLComputeCommandEncoder> Encoder() {return ctx->Encoder();}
	void Flush() {ctx->Flush();}

	//! Shared buffer of \a bytes, initialized with \a data or zero
	id<MTLBuffer> NewBuffer(size_t bytes, const void* data=NULL);

	//! Dispatch \a pso with one thread per (i,j,k), i fastest; the pipeline and its arguments must be set on Encoder()
	/*!
	  A dispatch sees the results of all previous dispatches, unless it and the
	  previous dispatch belong to the same \a group: dispatches of one group may run
	  concurrently, they must not access the same memory (e.g. the disjoint UPML regions).
	  */
	void Dispatch(id<MTLComputePipelineState> pso, size_t ni, size_t nj=1, size_t nk=1, const void* group=NULL);

	//! Buffer of the set indices of \a sets (16 or 32 bit)
	id<MTLBuffer> NewIndexBuffer(const Metal_CoeffSets& sets);

	//! Bind the grid dimension to \a index
	void SetGridDim(unsigned int index);

	//! Bind the main update coefficients for main_coeff(): the full arrays \a full_a and \a full_b (vv/vi or ii/iv) or the compressed ones
	void SetCoefficients(unsigned int index_idx, unsigned int a_idx, unsigned int b_idx, unsigned int mode_idx,
	                     id<MTLBuffer> full_a, id<MTLBuffer> full_b);
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
