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
#include "FDTD/gpu_coeff_sets.h"

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

// ---- UPML (see metal_ext_upml.mm), shared with the main updates (UPML regions along z)
struct UPMLParam { uint sx, sy, sz; uint lx, ly, lz; };

inline uint upml_local(constant UPMLParam& P, uint n, uint3 gid)
{
	return ((n*P.lx + gid.z)*P.ly + gid.y)*P.lz + gid.x;
}

// UPML coefficients of the region cell l0 (see Metal_Ext_UPML::m_CoeffMode): old[n*s],
// fo[n*s] and fn[n*s] are vv/vvfo/vvfn (set_offset 0) or ii/iifo/iifn (set_offset 9) of
// direction n. Mode 0: the full arrays, else a 16/32 bit set index and the sets of 18.
struct UPMLCoeff { const device float* old; const device float* fo; const device float* fn; uint s; };

inline UPMLCoeff upml_coeff(uint mode, const device void* index, const device float* sets,
                            const device float* c_old, const device float* c_fo, const device float* c_fn,
                            uint l0, uint cells, uint set_offset)
{
	UPMLCoeff c;
	if (mode==0)
	{
		c.old = c_old + l0; c.fo = c_fo + l0; c.fn = c_fn + l0; c.s = cells;
		return c;
	}
	const uint set = (mode==1) ? uint(((const device ushort*)index)[l0]) : ((const device uint*)index)[l0];
	c.old = sets + 18*set + set_offset; c.fo = c.old + 3; c.fn = c.old + 6; c.s = 1;
	return c;
}

// Fused: upml_pre, the main update (update_voltages) and upml_post of the region node
// gid = local (z, y, x), operation by operation the same. Valid because the UPML hooks run
// directly before and after the main update (see GPU_UPMLFusionBox()), and the voltage
// update only reads the voltage of its own node.
inline void upml_fused_volt_node(device float* volt, const device float* curr, device float* flux,
                                 const device float* c_old, const device float* c_fo, const device float* c_fn,
                                 const device void* index, const device float* ca, const device float* cb, uint mode,
                                 const device void* u_index, const device float* u_sets, uint u_mode,
                                 constant GridDim& N, constant UPMLParam& P, uint3 gid)
{
	const uint x = P.sx + gid.z, y = P.sy + gid.y, z = P.sz + gid.x;
	const uint sn = N.nx*N.ny*N.nz;
	const uint i  = nijk(N, 0, x, y, z);
	const uint xm = (x>0) ? N.ny*N.nz : 0;
	const uint ym = (y>0) ? N.nz : 0;
	const uint zm = (z>0) ? 1 : 0;
	const MainCoeff C = main_coeff(index, ca, cb, mode, sn, i, 0);
	const uint cells = P.lx*P.ly*P.lz;
	const UPMLCoeff U = upml_coeff(u_mode, u_index, u_sets, c_old, c_fo, c_fn, (gid.z*P.ly + gid.y)*P.lz + gid.x, cells, 0);
	float curl[3];
	curl[0] = (curr[2*sn+i] - curr[2*sn+i-ym] - curr[sn+i] + curr[sn+i-zm]);
	curl[1] = (curr[i] - curr[i-zm] - curr[2*sn+i] + curr[2*sn+i-xm]);
	curl[2] = (curr[sn+i] - curr[sn+i-xm] - curr[i] + curr[i-ym]);
	for (uint n=0; n<3; ++n)
	{
		const uint l = upml_local(P, n, gid);
		const uint g = n*sn + i;
		const float f_help = U.old[n*U.s]*volt[g] - U.fo[n*U.s]*flux[l];
		float v;
		v  = flux[l] * C.a[n*C.s];
		v += C.b[n*C.s] * curl[n];
		flux[l] = v;
		volt[g] = f_help + U.fn[n*U.s]*v;
	}
}

// same for the currents (update_currents), which are not updated on the last mesh lines
inline void upml_fused_curr_node(device float* curr, const device float* volt, device float* flux,
                                 const device float* c_old, const device float* c_fo, const device float* c_fn,
                                 const device void* index, const device float* ca, const device float* cb, uint mode,
                                 const device void* u_index, const device float* u_sets, uint u_mode,
                                 constant GridDim& N, constant UPMLParam& P, uint3 gid)
{
	const uint x = P.sx + gid.z, y = P.sy + gid.y, z = P.sz + gid.x;
	const uint sn = N.nx*N.ny*N.nz;
	const uint i  = nijk(N, 0, x, y, z);
	const bool update = (x+1<N.nx) && (y+1<N.ny) && (z+1<N.nz);
	const uint cells = P.lx*P.ly*P.lz;
	const UPMLCoeff U = upml_coeff(u_mode, u_index, u_sets, c_old, c_fo, c_fn, (gid.z*P.ly + gid.y)*P.lz + gid.x, cells, 9);
	float curl[3] = {0, 0, 0};
	MainCoeff C = {ca, cb, 0};
	if (update)
	{
		const uint xp = N.ny*N.nz;
		const uint yp = N.nz;
		C = main_coeff(index, ca, cb, mode, sn, i, 6);
		curl[0] = (volt[2*sn+i] - volt[2*sn+i+yp] - volt[sn+i] + volt[sn+i+1]);
		curl[1] = (volt[i] - volt[i+1] - volt[2*sn+i] + volt[2*sn+i+xp]);
		curl[2] = (volt[sn+i] - volt[sn+i+xp] - volt[i] + volt[i+yp]);
	}
	for (uint n=0; n<3; ++n)
	{
		const uint l = upml_local(P, n, gid);
		const uint g = n*sn + i;
		const float f_help = U.old[n*U.s]*curr[g] - U.fo[n*U.s]*flux[l];
		float c = flux[l];
		if (update)
		{
			c  = c * C.a[n*C.s];
			c += C.b[n*C.s] * curl[n];
		}
		flux[l] = c;
		curr[g] = f_help + U.fn[n*U.s]*c;
	}
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

	//! Commit the current batch without waiting, returns its command buffer (nil if empty)
	id<MTLCommandBuffer> Commit();
};


//! A UPML region along z, updated by the main kernels (see update_voltages), [0]: voltages, [1]: currents
struct Metal_ZSlab
{
	bool active;
	struct {uint32_t sx, sy, sz, lx, ly, lz;} region;
	id<MTLBuffer> flux[2], c_old[2], c_fo[2], c_fn[2];
	uint32_t u_mode;
	id<MTLBuffer> u_index, u_sets;
};

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
	id<MTLBuffer> snap_volt[2], snap_curr[2];   //!< field snapshots, see GPU_Backend_Metal::SnapshotFields()
	id<MTLCommandBuffer> snap_cmd[2];           //!< the command buffers copying them

	//! The main updates cover the nodes in [main_start, main_stop), the fused UPML kernels the others (see metal_ext_upml.mm)
	Metal_GridDim main_start, main_stop;
	//! UPML extensions of this grid and whether they run fused with the main updates (-1: not decided yet), see metal_ext_upml.mm
	std::vector<GPU_Extension*> upml;
	int upml_fused;
	//! UPML regions along z updated by the main kernels (lo: below main_start.nz, hi: from main_stop.nz on)
	Metal_ZSlab zslab[2];

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
	id<MTLBuffer> NewIndexBuffer(const GPU_CoeffSets& sets);

	//! Bind the grid dimension to \a index
	void SetGridDim(unsigned int index);

	//! No UPML regions updated by the main kernels
	void ResetZSlabs();
	//! Nodes [b, e) of the main kernels: main_start/main_stop and the UPML regions along z they update
	void MainRange(Metal_GridDim& b, Metal_GridDim& e) const;
	//! Bind the UPML regions along z for the main update of the voltages (\a field 0) or currents (1)
	void SetZSlabs(int field);

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
