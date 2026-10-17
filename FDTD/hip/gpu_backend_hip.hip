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

#include "hip/hip_runtime.h"
#include <algorithm>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <set>
#include <iostream>
#include <limits>
#include <stdexcept>

#include "hip_internal.h"
#include "FDTD/operator.h"
#include "FDTD/extensions/engine_ext_upml.h"
#include "FDTD/extensions/engine_ext_excitation.h"
#include "FDTD/extensions/engine_ext_lumpedRLC.h"
#include "FDTD/extensions/engine_ext_steadystate.h"

// Main FDTD updates, operation by operation the same as Engine::UpdateVoltages/UpdateCurrents,
// over the nodes in [B, E). One thread per mesh node: z along x-threads (contiguous), y along
// y-threads, x along the grid z-dimension. The coefficients are read with main_coeff().
// Nodes below Z.z_lo or from Z.z_hi on belong to UPML regions along z (see HIP_ZSlabs).
__global__ void update_voltages(float* __restrict__ volt, const float* __restrict__ curr,
                                const void* index, const float* ca, const float* cb, unsigned int mode,
                                HIP_GridDim N, HIP_GridDim B, HIP_GridDim E, HIP_ZSlabs Z)
{
	// with UPML regions along z: slab lo, main nodes and slab hi in whole warps each (no divergence)
	const unsigned int tz = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int z = B.nz + tz;
	if (Z.t_int)
	{
		if (tz<Z.t_lo)
		{
			if (B.nz+tz>=Z.z_lo)
				return;
		}
		else if (tz<Z.t_lo+Z.t_int)
		{
			z = Z.z_lo + (tz-Z.t_lo);
			if (z>=Z.z_hi)
				return;
		}
		else
			z = Z.z_hi + (tz-Z.t_lo-Z.t_int);
	}
	const unsigned int y = B.ny + blockIdx.y*blockDim.y + threadIdx.y;
	const unsigned int x = B.nx + blockIdx.z;
	if (z>=E.nz || y>=E.ny || x>=E.nx)
		return;
	if (z<Z.z_lo)
	{
		upml_fused_volt_node(volt, volt, curr, Z.volt_lo, N, Z.R_lo, x-Z.R_lo.sx, y-Z.R_lo.sy, z-Z.R_lo.sz);
		return;
	}
	if (z>=Z.z_hi)
	{
		upml_fused_volt_node(volt, volt, curr, Z.volt_hi, N, Z.R_hi, x-Z.R_hi.sx, y-Z.R_hi.sy, z-Z.R_hi.sz);
		return;
	}
	const size_t sn = (size_t)N.nx*N.ny*N.nz;
	const size_t i  = nijk(N, 0, x, y, z);
	const size_t xm = (x>0) ? (size_t)N.ny*N.nz : 0;   // shift to the previous line, none on the first
	const size_t ym = (y>0) ? N.nz : 0;
	const size_t zm = (z>0) ? 1 : 0;
	const HIP_MainCoeff C = main_coeff(index, ca, cb, mode, sn, i, 0);
	const float* vv = C.a;
	const float* vi = C.b;
	const size_t s = C.s;
	float v;

	//for x
	v  = volt[i] * vv[0];
	v += vi[0] * (curr[2*sn+i] - curr[2*sn+i-ym] - curr[sn+i] + curr[sn+i-zm]);
	volt[i] = v;

	//for y
	v  = volt[sn+i] * vv[s];
	v += vi[s] * (curr[i] - curr[i-zm] - curr[2*sn+i] + curr[2*sn+i-xm]);
	volt[sn+i] = v;

	//for z
	v  = volt[2*sn+i] * vv[2*s];
	v += vi[2*s] * (curr[sn+i] - curr[sn+i-xm] - curr[i] + curr[i-ym]);
	volt[2*sn+i] = v;
}

// the currents on the last mesh line of each direction are not updated
__global__ void update_currents(float* __restrict__ curr, const float* __restrict__ volt,
                                const void* index, const float* ca, const float* cb, unsigned int mode,
                                HIP_GridDim N, HIP_GridDim B, HIP_GridDim E, HIP_ZSlabs Z)
{
	// with UPML regions along z: slab lo, main nodes and slab hi in whole warps each (no divergence)
	const unsigned int tz = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int z = B.nz + tz;
	if (Z.t_int)
	{
		if (tz<Z.t_lo)
		{
			if (B.nz+tz>=Z.z_lo)
				return;
		}
		else if (tz<Z.t_lo+Z.t_int)
		{
			z = Z.z_lo + (tz-Z.t_lo);
			if (z>=Z.z_hi)
				return;
		}
		else
			z = Z.z_hi + (tz-Z.t_lo-Z.t_int);
	}
	const unsigned int y = B.ny + blockIdx.y*blockDim.y + threadIdx.y;
	const unsigned int x = B.nx + blockIdx.z;
	if (z>=E.nz || y>=E.ny || x>=E.nx)
		return;
	if (z<Z.z_lo)
	{
		upml_fused_curr_node(curr, curr, volt, Z.curr_lo, N, Z.R_lo, x-Z.R_lo.sx, y-Z.R_lo.sy, z-Z.R_lo.sz);
		return;
	}
	if (z>=Z.z_hi)
	{
		upml_fused_curr_node(curr, curr, volt, Z.curr_hi, N, Z.R_hi, x-Z.R_hi.sx, y-Z.R_hi.sy, z-Z.R_hi.sz);
		return;
	}
	if (z+1>=N.nz || y+1>=N.ny || x+1>=N.nx)
		return;
	const size_t sn = (size_t)N.nx*N.ny*N.nz;
	const size_t i  = nijk(N, 0, x, y, z);
	const size_t xp = (size_t)N.ny*N.nz;
	const size_t yp = N.nz;
	const HIP_MainCoeff C = main_coeff(index, ca, cb, mode, sn, i, 6);
	const float* ii = C.a;
	const float* iv = C.b;
	const size_t s = C.s;
	float c;

	//for x
	c  = curr[i] * ii[0];
	c += iv[0] * (volt[2*sn+i] - volt[2*sn+i+yp] - volt[sn+i] + volt[sn+i+1]);
	curr[i] = c;

	//for y
	c  = curr[sn+i] * ii[s];
	c += iv[s] * (volt[i] - volt[i+1] - volt[2*sn+i] + volt[2*sn+i+xp]);
	curr[sn+i] = c;

	//for z
	c  = curr[2*sn+i] * ii[2*s];
	c += iv[2*s] * (volt[sn+i] - volt[sn+i+xp] - volt[i] + volt[i+yp]);
	curr[2*sn+i] = c;
}

// ---- fused step

// the voltages of node (x,y,z), the same operations as update_voltages
template<typename I>
__device__ __forceinline__ void fused_volt(const float* __restrict__ volt, const float* __restrict__ curr,
                                           const void* index, const float* va, const float* vb, unsigned int mode,
                                           const HIP_GridDim& N, unsigned int x, unsigned int y, unsigned int z, float v[3])
{
	const I sn = (I)N.nx*N.ny*N.nz;
	const I i  = nijk<I>(N, 0, x, y, z);
	const I xm = (x>0) ? (I)N.ny*N.nz : 0;
	const I ym = (y>0) ? N.nz : 0;
	const I zm = (z>0) ? 1 : 0;
	const HIP_MainCoeffT<I> C = main_coeff<I>(index, va, vb, mode, sn, i, 0);
	float t;
	t  = volt[i] * C.a[0];
	t += C.b[0] * (curr[2*sn+i] - curr[2*sn+i-ym] - curr[sn+i] + curr[sn+i-zm]);
	v[0] = t;
	t  = volt[sn+i] * C.a[C.s];
	t += C.b[C.s] * (curr[i] - curr[i-zm] - curr[2*sn+i] + curr[2*sn+i-xm]);
	v[1] = t;
	t  = volt[2*sn+i] * C.a[2*C.s];
	t += C.b[2*C.s] * (curr[sn+i] - curr[sn+i-xm] - curr[i] + curr[i-ym]);
	v[2] = t;
}

// The tile is as wide as a wavefront, so that a wavefront is one row of it: 32 threads on
// NVIDIA and on the RDNA GPUs of AMD, 64 on the CDNA ones. The kernel is compiled for both
// widths and the backend picks by the warp size of the device (see UpdateVoltages), so one
// binary serves either.
// OPENEMS_HIP_FUSED_TZ and _TY override the shape, to measure it on a new GPU.
#define FUSED_TZ 32   // threads along z, the last one computes the border voltages only
#define FUSED_TY 8    // threads along y, the last row computes the border voltages only
#define FUSED_XC 4    // x lines per block (default): short blocks are faster, despite the voltages of one more x line

// Fused step: the voltages and currents of the main nodes [B, E) in one pass, from the current
// fields (volt_in, curr_in) into the next buffers (volt_out, curr_out); every value is computed
// with the same operations as update_voltages/update_currents.
// One kernel for the whole timestep: the current update reads the new voltages from shared
// memory, instead of a second pass reading them back from global memory.
// A block marches along FUSED_XC x lines of a (FUSED_TZ-1)x(FUSED_TY-1) tile in (z,y). At each
// x it computes the new voltages of the tile and its border (+z, +y) into shared memory and of
// the next x line (registers), then the new currents of the tile. The border voltages of other
// tiles are computed again (they only read the current fields), those outside the main nodes
// (UPML regions) come from volt_out, where the UPML regions stored them before.
// With the UPML regions in the kernel (F.count>0, [B, E) is the whole grid), the region nodes
// are computed like the main nodes, from the current flux (the new voltage flux goes to a
// second buffer, written by the node's own block; the current flux in place).
template<typename I, unsigned int TZ, unsigned int TY>
__global__ void __launch_bounds__(TZ*TY) update_fused(
	const float* __restrict__ volt_in, const float* __restrict__ curr_in, float* volt_out, float* __restrict__ curr_out,
	const void* index, const float* va, const float* vb, const float* ia, const float* ib, unsigned int mode,
	HIP_GridDim N, HIP_GridDim B, HIP_GridDim E, HIP_FusedRegions F, unsigned int xc)
{
	__shared__ float sV[3][TY][TZ];
	const unsigned int tz = threadIdx.x, ty = threadIdx.y;
	const unsigned int z = B.nz + blockIdx.x*(TZ-1) + tz;
	const unsigned int y = B.ny + blockIdx.y*(TY-1) + ty;
	const unsigned int xs = B.nx + blockIdx.z*xc;
	const unsigned int xe = min(xs+xc, E.nx);
	const bool main_yz = (z<E.nz) && (y<E.ny);                          // a main node line
	const bool owner = main_yz && (tz<TZ-1) && (ty<TY-1);  // written by this block
	const bool in_grid = (z<N.nz) && (y<N.ny);
	const I sn = (I)N.nx*N.ny*N.nz;

	// the new voltages at x: computed on main nodes, else from the UPML regions (or unused)
	auto volt_at = [&](unsigned int x, float v[3])
	{
		if (main_yz && (x<E.nx))
		{
			const int r = F.count ? fused_region_of(F, x, y, z) : -1;
			if (r<0)
				fused_volt<I>(volt_in, curr_in, index, va, vb, mode, N, x, y, z, v);
			else
				upml_volt_value<I>(volt_in, curr_in, F.volt[r], N, F.R[r], x, y, z, owner && (x<xe), v);
			if (owner && (x<xe))
			{
				const I i = nijk<I>(N, 0, x, y, z);
				volt_out[i] = v[0]; volt_out[sn+i] = v[1]; volt_out[2*sn+i] = v[2];
			}
		}
		else if (in_grid && (x<N.nx))
		{
			const I i = nijk<I>(N, 0, x, y, z);
			v[0] = volt_out[i]; v[1] = volt_out[sn+i]; v[2] = volt_out[2*sn+i];
		}
		else
			v[0] = v[1] = v[2] = 0;
	};

	float v[3];
	volt_at(xs, v);
	sV[0][ty][tz] = v[0]; sV[1][ty][tz] = v[1]; sV[2][ty][tz] = v[2];
	__syncthreads();
	for (unsigned int x=xs; x<xe; ++x)
	{
		float v1[3];   // the new voltages at x+1
		volt_at(x+1, v1);
		// the currents at x, not updated on the last mesh lines (see update_currents)
		if (owner)
		{
			const bool update = (x+1<N.nx) && (y+1<N.ny) && (z+1<N.nz);
			float curl[3] = {0, 0, 0};
			if (update)
			{
				curl[0] = (sV[2][ty][tz] - sV[2][ty+1][tz] - sV[1][ty][tz] + sV[1][ty][tz+1]);
				curl[1] = (sV[0][ty][tz] - sV[0][ty][tz+1] - sV[2][ty][tz] + v1[2]);
				curl[2] = (sV[1][ty][tz] - v1[1] - sV[0][ty][tz] + sV[0][ty+1][tz]);
			}
			const int r = F.count ? fused_region_of(F, x, y, z) : -1;
			if (r>=0)
				upml_curr_value<I>(curr_in, curr_out, F.curr[r], N, F.R[r], x, y, z, update, curl);
			else if (update)
			{
				const I i = nijk<I>(N, 0, x, y, z);
				const HIP_MainCoeffT<I> C = main_coeff<I>(index, ia, ib, mode, sn, i, 6);
				float c;
				//for x
				c  = curr_in[i] * C.a[0];
				c += C.b[0] * curl[0];
				curr_out[i] = c;
				//for y
				c  = curr_in[sn+i] * C.a[C.s];
				c += C.b[C.s] * curl[1];
				curr_out[sn+i] = c;
				//for z
				c  = curr_in[2*sn+i] * C.a[2*C.s];
				c += C.b[2*C.s] * curl[2];
				curr_out[2*sn+i] = c;
			}
		}
		__syncthreads();
		sV[0][ty][tz] = v1[0]; sV[1][ty][tz] = v1[1]; sV[2][ty][tz] = v1[2];
		__syncthreads();
	}
}

// The arguments of update_fused, so that a launch is written once for every tile shape
struct HIP_FusedArgs
{
	const float* volt_in;
	const float* curr_in;
	float* volt_out;
	float* curr_out;
	const void* index;
	const float* va;
	const float* vb;
	const float* ia;
	const float* ib;
	unsigned int mode;
	HIP_GridDim N, B, E;
	HIP_FusedRegions F;
	unsigned int xc;
};

template<typename I, unsigned int TZ, unsigned int TY>
static void launch_fused(hipStream_t stream, const HIP_FusedArgs& a)
{
	const dim3 block(TZ, TY, 1);
	const dim3 grid((a.E.nz-a.B.nz+TZ-2)/(TZ-1), (a.E.ny-a.B.ny+TY-2)/(TY-1), (a.E.nx-a.B.nx+a.xc-1)/a.xc);
	update_fused<I, TZ, TY><<<grid, block, 0, stream>>>(a.volt_in, a.curr_in, a.volt_out, a.curr_out, a.index,
	                                                 a.va, a.vb, a.ia, a.ib, a.mode, a.N, a.B, a.E, a.F, a.xc);
}

template<typename I, unsigned int TZ>
static void launch_fused_ty(hipStream_t stream, const HIP_FusedArgs& a, unsigned int ty)
{
	switch (ty)
	{
		case 4:  launch_fused<I, TZ, 4>(stream, a); break;
		case 16: launch_fused<I, TZ, 16>(stream, a); break;
		default: launch_fused<I, TZ, 8>(stream, a); break;
	}
}

// the index type of the fused step, see GPU_Backend_HIP::Impl::index32
template<unsigned int TZ>
static void launch_fused_tz(hipStream_t stream, const HIP_FusedArgs& a, unsigned int ty, bool index32)
{
	if (index32)
		launch_fused_ty<unsigned int, TZ>(stream, a, ty);
	else
		launch_fused_ty<size_t, TZ>(stream, a, ty);
}

// the raw values of a dump region (see HIP_GatherRegion) for a snapshot: only these go over
// PCIe, not the whole fields; the host interpolates them
__global__ void gather_region(const float* __restrict__ f, HIP_GatherRegion R, HIP_GridDim N, float* __restrict__ out)
{
	const size_t e = (size_t)blockIdx.x*blockDim.x + threadIdx.x;
	if (e<R.total)
		out[e] = f[region_field_index(R, N, e)];
}

// frequency domain dumps: sums[k*total + e] += raw value(e) * w[k] for the weights of one timestep,
// a float times a complex float, then a complex sum; the host interpolates the sums
#define DFT_WEIGHTS 32   // weights per launch, passed by value in the kernel parameters
struct HIP_DFTWeights { float2 w[DFT_WEIGHTS]; };

__global__ void field_dft(const float* __restrict__ f, HIP_GatherRegion R, HIP_GridDim N,
                          HIP_DFTWeights W, unsigned int nw, float2* __restrict__ sums)
{
	const size_t e = (size_t)blockIdx.x*blockDim.x + threadIdx.x;
	if (e>=R.total)
		return;
	const float v = f[region_field_index(R, N, e)];
	for (unsigned int k=0; k<nw; ++k)
	{
		float2 s = sums[k*R.total + e];
		s.x += v*W.w[k].x;
		s.y += v*W.w[k].y;
		sums[k*R.total + e] = s;
	}
}

// the currents of the given main nodes again (update_currents), after the voltage extensions
// changed voltages they read (see GPU_Backend_HIP::Impl::fixup)
__global__ void update_currents_nodes(const float* __restrict__ curr_in, float* __restrict__ curr_out, const float* __restrict__ volt,
                                      const size_t* nodes, unsigned int count,
                                      const void* index, const float* ia, const float* ib, unsigned int mode, HIP_GridDim N)
{
	const unsigned int k = blockIdx.x*blockDim.x + threadIdx.x;
	if (k>=count)
		return;
	const size_t sn = (size_t)N.nx*N.ny*N.nz;
	const size_t i  = nodes[k];
	const size_t xp = (size_t)N.ny*N.nz;
	const size_t yp = N.nz;
	const HIP_MainCoeff C = main_coeff(index, ia, ib, mode, sn, i, 6);
	float c;
	c  = curr_in[i] * C.a[0];
	c += C.b[0] * (volt[2*sn+i] - volt[2*sn+i+yp] - volt[sn+i] + volt[sn+i+1]);
	curr_out[i] = c;
	c  = curr_in[sn+i] * C.a[C.s];
	c += C.b[C.s] * (volt[i] - volt[i+1] - volt[2*sn+i] + volt[2*sn+i+xp]);
	curr_out[sn+i] = c;
	c  = curr_in[2*sn+i] * C.a[2*C.s];
	c += C.b[2*C.s] * (volt[sn+i] - volt[sn+i+xp] - volt[i] + volt[i+yp]);
	curr_out[2*sn+i] = c;
}

// Field energy: the squared voltages and currents of the nodes below L, summed along x
// in float, one thread per (y,z) line. The host sums the lines in double.
__global__ void energy_lines(const float* __restrict__ volt, const float* __restrict__ curr, float2* partial,
                             HIP_GridDim N, HIP_GridDim L)
{
	const unsigned int z = blockIdx.x*blockDim.x + threadIdx.x;
	const unsigned int y = blockIdx.y*blockDim.y + threadIdx.y;
	if (z>=L.nz || y>=L.ny)
		return;
	const size_t sn = (size_t)N.nx*N.ny*N.nz;
	float E = 0, H = 0;
	for (unsigned int x=0; x<L.nx; ++x)
	{
		const size_t i = nijk(N, 0, x, y, z);
		E += volt[i]*volt[i];
		E += volt[sn+i]*volt[sn+i];
		E += volt[2*sn+i]*volt[2*sn+i];
		H += curr[i]*curr[i];
		H += curr[sn+i]*curr[sn+i];
		H += curr[2*sn+i]*curr[2*sn+i];
	}
	partial[y*L.nz + z] = make_float2(E, H);
}

/***************************** helpers *****************************/

namespace
{
struct KernelTimeTable
{
	hipEvent_t start = NULL, stop = NULL;
	std::map<std::string, std::pair<double, size_t>> times;   // total ms, launches
	~KernelTimeTable()
	{
		if (times.empty())
			return;
		double total = 0;
		for (const auto& t : times)
			total += t.second.first;
		std::vector<std::pair<double, std::string>> sorted;
		for (const auto& t : times)
			sorted.push_back(std::make_pair(t.second.first, t.first));
		std::sort(sorted.rbegin(), sorted.rend());
		fprintf(stderr, "HIP kernel times (%.1f ms total):\n", total);
		for (const auto& t : sorted)
		{
			const size_t n = times[t.second].second;
			fprintf(stderr, "  %6.1f %%  %10.1f ms  %8zu x %8.1f us  %s\n", 100*t.first/total, t.first, n, 1000*t.first/n, t.second.c_str());
		}
	}
};
KernelTimeTable& KernelTimes() {static KernelTimeTable table; return table;}
}

bool HIP_KernelTimes::Enabled()
{
	static const bool enabled = getenv("OPENEMS_HIP_KERNEL_TIMES") && (atoi(getenv("OPENEMS_HIP_KERNEL_TIMES"))!=0);
	return enabled;
}

void HIP_KernelTimes::Begin(hipStream_t stream)
{
	KernelTimeTable& t = KernelTimes();
	if (!t.start)
	{
		hipEventCreate(&t.start);
		hipEventCreate(&t.stop);
	}
	hipEventRecord(t.start, stream);
}

void HIP_KernelTimes::End(hipStream_t stream, const char* name)
{
	KernelTimeTable& t = KernelTimes();
	hipEventRecord(t.stop, stream);
	hipEventSynchronize(t.stop);
	float ms = 0;
	hipEventElapsedTime(&ms, t.start, t.stop);
	std::pair<double, size_t>& entry = t.times[name];
	entry.first += ms;
	++entry.second;
}

void HIP_Check(hipError_t err, const char* what)
{
	if (err!=hipSuccess)
		throw std::runtime_error(std::string("GPU_Backend_HIP: ") + what + " failed: " + hipGetErrorString(err));
}

HIP_Context::HIP_Context()
{
	device = 0;
	stream = 0;
	copy_stream = 0;
	warp_size = 32;
}

HIP_Context::~HIP_Context()
{
	if (stream)
	{
		hipStreamSynchronize(stream);
		hipStreamDestroy(stream);
	}
	if (copy_stream)
	{
		hipStreamSynchronize(copy_stream);
		hipStreamDestroy(copy_stream);
	}
}

GPU_Backend_HIP::Impl::Impl()
{
	dim.nx = dim.ny = dim.nz = 0;
	numCells = 0;
	volt = curr = vv = vi = ii = iv = NULL;
	coeff_mode = 0;
	index = NULL;
	coeff = NULL;
	energy = NULL;
	energy_count = 0;
	fused_step = -1;
	fused_blockers = 0;
	fregions.count = 0;
	volt_next = curr_next = NULL;
	fixup = NULL;
	index32 = false;
	fixup_count = 0;
	snap_nv = snap_n = 0;
	snap_dev = NULL;
	snap_host = NULL;
	snap_taken = false;
	HIP_Check(hipEventCreateWithFlags(&snap_evaluated, hipEventDisableTiming), "hipEventCreate");
	HIP_Check(hipEventCreateWithFlags(&snap_done, hipEventDisableTiming), "hipEventCreate");
}

GPU_Backend_HIP::Impl::~Impl()
{
	if (ctx)
	{
		hipStreamSynchronize(ctx->stream);
		hipStreamSynchronize(ctx->copy_stream);
	}
	FreeSnapshots();
	hipEventDestroy(snap_evaluated);
	hipEventDestroy(snap_done);
	for (size_t n=0; n<m_Allocations.size(); ++n)
		hipFree(m_Allocations.at(n));
	for (std::set<void*>::iterator it=m_Pinned.begin(); it!=m_Pinned.end(); ++it)
		hipHostUnregister(*it);
}

bool GPU_Backend_HIP::Impl::DecideFusedStep()
{
	if (fused_step>=0)
		return fused_step>0;
	fused_step = 0;
	const char* env = getenv("OPENEMS_HIP_FUSED_STEP");
	if ((fused_blockers>0) || (env && (atoi(env)==0)))
		return false;

	// The currents of the nodes that read a voltage changed between the half-steps are computed
	// again. With the UPML regions in the kernel, region nodes cannot (the current flux is
	// updated in place): then the region kernels update the regions.
	const bool folded = fregions.count>0;
	std::set<size_t> nodes;
	for (size_t k=0; k<volt_modified.size(); ++k)
	{
		const size_t node = volt_modified[k] % numCells;
		const unsigned int p[3] = {(unsigned int)(node/((size_t)dim.ny*dim.nz)), (unsigned int)((node/dim.nz)%dim.ny), (unsigned int)(node%dim.nz)};
		const unsigned int start[3] = {main_start.nx, main_start.ny, main_start.nz};
		const unsigned int stop[3] = {main_stop.nx, main_stop.ny, main_stop.nz};
		const unsigned int n_lines[3] = {dim.nx, dim.ny, dim.nz};
		for (int m=-1; m<3; ++m)   // the node itself and its lower neighbours
		{
			unsigned int c[3] = {p[0], p[1], p[2]};
			if (m>=0)
			{
				if (c[m]==0)
					continue;
				--c[m];
			}
			bool ok = true;
			for (int n=0; n<3; ++n)
				ok &= (folded || ((c[n]>=start[n]) && (c[n]<stop[n]))) && (c[n]+1<n_lines[n]);
			if (!ok)
				continue;
			if (folded && (fused_region_of(fregions, c[0], c[1], c[2])>=0))
				fregions.count = 0;   // a region node: no regions in the kernel
			nodes.insert(((size_t)c[0]*dim.ny + c[1])*dim.nz + c[2]);
		}
	}
	std::vector<size_t> list;
	for (std::set<size_t>::const_iterator it=nodes.begin(); it!=nodes.end(); ++it)
	{
		const unsigned int c[3] = {(unsigned int)(*it/((size_t)dim.ny*dim.nz)), (unsigned int)((*it/dim.nz)%dim.ny), (unsigned int)(*it%dim.nz)};
		if ((c[0]>=main_start.nx) && (c[0]<main_stop.nx) && (c[1]>=main_start.ny) && (c[1]<main_stop.ny) && (c[2]>=main_start.nz) && (c[2]<main_stop.nz))
			list.push_back(*it);   // main nodes only, the region kernels update the others after the extensions
		else if (fregions.count)
			list.push_back(*it);
	}
	// The buffers are allocated last, after the field dumps: without them the run is only slower.
	volt_next = TryAlloc<float>(3*numCells);
	curr_next = TryAlloc<float>(3*numCells);
	if ((volt_next==NULL) || (curr_next==NULL))
	{
		Free(volt_next);
		Free(curr_next);
		volt_next = curr_next = NULL;
		std::cout << "GPU_Backend_HIP: not enough device memory for the fused step, using separate updates" << std::endl;
		return false;
	}
	fixup_count = list.size();
	fixup = Alloc<size_t>(fixup_count, list.data());

	// both buffers start with the current fields
	HIP_Check(hipMemcpyAsync(volt_next, volt, 3*numCells*sizeof(float), hipMemcpyDeviceToDevice, Stream()), "fused step buffers");
	HIP_Check(hipMemcpyAsync(curr_next, curr, 3*numCells*sizeof(float), hipMemcpyDeviceToDevice, Stream()), "fused step buffers");
	Flush();
	fused_step = 1;
	if (g_settings.GetVerboseLevel()>1)
	{
		std::cout << "GPU_Backend_HIP: fused voltage and current updates" << (index32 ? ", 32 bit indices" : "");
		if (fregions.count)
			std::cout << " (" << fregions.count << " UPML regions in the kernel)";
		std::cout << std::endl;
	}
	return true;
}

void GPU_Backend_HIP::Impl::ResetZSlabs()
{
	zslab_lo = zslab_hi = false;
	zslabs.z_lo = main_start.nz;
	zslabs.z_hi = main_stop.nz;
}

void GPU_Backend_HIP::Impl::MainRange(HIP_GridDim& B, HIP_GridDim& E, unsigned int& threads_z)
{
	B = main_start;
	E = main_stop;
	if (zslab_lo)
		B.nz = zslabs.R_lo.sz;
	if (zslab_hi)
		E.nz = zslabs.R_hi.sz + zslabs.R_hi.lz;
	threads_z = E.nz - B.nz;
	zslabs.t_lo = zslabs.t_int = 0;
	// Slab lo, main nodes and slab hi in whole warps each: a warp shared by UPML and main
	// nodes runs both code paths (on a RTX 5070 Ti 3 % slower for 300^3 cells)
	if (zslab_lo || zslab_hi)
	{
		const unsigned int w = ctx->warp_size;
		auto warps = [w](unsigned int n) {return (n+w-1)/w*w;};
		zslabs.t_lo = warps(zslabs.z_lo - B.nz);
		zslabs.t_int = warps(zslabs.z_hi - zslabs.z_lo);
		threads_z = zslabs.t_lo + zslabs.t_int + warps(E.nz - zslabs.z_hi);
	}
}

void GPU_Backend_HIP::Impl::FreeSnapshots()
{
	for (size_t r=0; r<snap_lines.size(); ++r)
		Free(snap_lines[r]);
	snap_lines.clear();
	snap_regions[0].clear();
	snap_regions[1].clear();
	hipFree(snap_dev);
	hipHostFree(snap_host);
	snap_dev = snap_host = NULL;
	snap_nv = snap_n = 0;
	snap_taken = false;
}

void GPU_Backend_HIP::Impl::Flush()
{
	HIP_Check(hipStreamSynchronize(Stream()), "stream synchronization");
}

void GPU_Backend_HIP::Impl::CheckLaunch(const char* kernel)
{
	HIP_Check(hipGetLastError(), kernel);
}

void GPU_Backend_HIP::Impl::PinHostMemory(void* ptr, size_t bytes)
{
	if (m_Pinned.count(ptr))
		return;
	// not fatal: unpinned memory only makes the transfers slower
	if (hipHostRegister(ptr, bytes, hipHostRegisterDefault)==hipSuccess)
		m_Pinned.insert(ptr);
	else
		hipGetLastError();
}

dim3 GPU_Backend_HIP::Impl::Block(size_t ni, size_t nj)
{
	if (nj<=1)
		return dim3(256, 1, 1);
	if (ni<=16)
		return dim3(16, 16, 1);
	return dim3(32, 8, 1);
}

dim3 GPU_Backend_HIP::Impl::Grid(dim3 block, size_t ni, size_t nj, size_t nk)
{
	return dim3((ni+block.x-1)/block.x, (nj+block.y-1)/block.y, nk);
}

/***************************** GPU_Backend_HIP *****************************/

//! Smallest kernel there is, to see whether this build has code for this device (see New())
__global__ void probe_device()
{
}

GPU_Backend_HIP* GPU_Backend_HIP::New()
{
	int count = 0;
	if ((hipGetDeviceCount(&count)!=hipSuccess) || (count==0))
	{
		hipGetLastError();
		return NULL;
	}
	Impl* impl = new Impl();
	impl->ctx = std::make_shared<HIP_Context>();
	HIP_Check(hipGetDevice(&impl->ctx->device), "hipGetDevice");
	hipDeviceProp_t prop;
	HIP_Check(hipGetDeviceProperties(&prop, impl->ctx->device), "hipGetDeviceProperties");
	impl->ctx->name = prop.name;
	impl->ctx->warp_size = prop.warpSize > 0 ? prop.warpSize : 32;
	HIP_Check(hipStreamCreateWithFlags(&impl->ctx->stream, hipStreamNonBlocking), "hipStreamCreate");
	HIP_Check(hipStreamCreateWithFlags(&impl->ctx->copy_stream, hipStreamNonBlocking), "hipStreamCreate");

	// A GPU whose architecture this build has no code for is found and described like any
	// other, and launching on it takes the process down with it rather than returning an
	// error. Ask about a kernel instead of running one: this reports the missing code and
	// leaves the CPU to fall back to.
	hipFuncAttributes attr;
	hipError_t err = hipFuncGetAttributes(&attr, reinterpret_cast<const void*>(probe_device));
	if (err!=hipSuccess)
	{
		// the architecture is the name to add to the build, but only the AMD platform
		// fills it in: the NVIDIA one leaves gcnArchName of the properties untouched
		std::string arch;
#ifdef __HIP_PLATFORM_AMD__
		arch = std::string(" (") + prop.gcnArchName + ")";
#endif
		std::cerr << "GPU_Backend_HIP: no kernels for " << impl->ctx->name << arch
		          << " in this build: " << hipGetErrorString(err) << std::endl;
		hipGetLastError();
		delete impl;
		return NULL;
	}
	return new GPU_Backend_HIP(impl);
}

GPU_Backend* GPU_Backend_HIP::NewSubGridBackend()
{
	// no fused step with sub-grids (cylindrical multi-grid)
	d->fused_blockers++;
	Impl* impl = new Impl();
	impl->fused_blockers++;
	impl->ctx = d->ctx;
	return new GPU_Backend_HIP(impl);
}

GPU_Backend_HIP::GPU_Backend_HIP(Impl* impl)
{
	d = impl;
}

GPU_Backend_HIP::~GPU_Backend_HIP()
{
	delete d;
}

std::string GPU_Backend_HIP::GetName() const
{
	return "HIP (" + d->ctx->name + ")";
}

bool GPU_Backend_HIP::Init(const Operator* op)
{
	unsigned int numLines[3];
	for (int n=0; n<3; ++n)
		numLines[n] = op->GetNumberOfLines(n, true);
	d->dim.nx = numLines[0];
	d->dim.ny = numLines[1];
	d->dim.nz = numLines[2];
	d->numCells = (size_t)numLines[0]*numLines[1]*numLines[2];
	d->main_start.nx = d->main_start.ny = d->main_start.nz = 0;
	d->main_stop = d->dim;
	d->upml_fused = -1;
	d->ResetZSlabs();

	// the main updates launch one block per x line along z, which the device limits (64k on AMD and NVIDIA)
	int max_grid_z = 0;
	HIP_Check(hipDeviceGetAttribute(&max_grid_z, hipDeviceAttributeMaxGridDimZ, d->ctx->device), "hipDeviceGetAttribute");
	if ((max_grid_z>0) && (d->dim.nx>(unsigned int)max_grid_z))
	{
		std::cerr << "GPU_Backend_HIP::Init: Error: " << d->dim.nx << " mesh lines in x, the device allows at most "
		          << max_grid_z << "; swap the axes, or use a CPU engine" << std::endl;
		return false;
	}

	// The fused step indexes with 32 bit on NVIDIA while the fields fit (faster there), with
	// 64 bit on AMD; OPENEMS_HIP_INDEX=32 or 64 overrides.
	const bool fits32 = 3*d->numCells <= UINT_MAX;
#ifdef __HIP_PLATFORM_AMD__
	d->index32 = false;
#else
	d->index32 = fits32;
#endif
	if (const char* env = getenv("OPENEMS_HIP_INDEX"))
		d->index32 = (atoi(env)==32);
	if (d->index32 && !fits32)
	{
		std::cerr << "GPU_Backend_HIP: the mesh is too large for 32 bit indices, the fused step uses 64 bit" << std::endl;
		d->index32 = false;
	}

	const size_t count = 3*d->numCells;
	d->volt = d->Alloc<float>(count);
	d->curr = d->Alloc<float>(count);

	// upload the final operator coefficients, including all changes by operator extensions.
	// Most nodes share one of a few coefficient sets (same material and mesh spacing): store
	// every distinct set of the 12 coefficients of a node once and a set index per node
	GPU_CoeffSets sets;
	if (GPU_FindSets(d->numCells, 12, [&](size_t i, float* values) {GPU_NodeCoefficients(op, numLines, i, values);}, sets))
	{
		d->coeff_mode = sets.mode;
		d->index = d->AllocIndex(sets);
		d->coeff = d->Alloc<float>(sets.table.size(), sets.table.data());
		if (g_settings.GetVerboseLevel()>1)
			std::cout << "GPU_Backend_HIP: " << sets.count << " distinct coefficient sets, compressed update coefficients ("
			          << (sets.mode==1 ? 16 : 32) << " bit index)" << std::endl;
		return true;
	}

	// the full arrays, one host copy at a time
	d->coeff_mode = 0;
	std::vector<float> values;
	auto upload = [&](int c) {GPU_CoefficientArray(op, numLines, c, values); return d->Alloc<float>(count, values.data());};
	d->vv = upload(0);
	d->vi = upload(1);
	d->ii = upload(2);
	d->iv = upload(3);
	return true;
}

void GPU_Backend_HIP::UpdateVoltages()
{
	const bool c = d->coeff_mode>0;
	if (d->DecideFusedStep())
	{
		// the main nodes, or the whole grid with the UPML regions in the kernel (see update_fused)
		HIP_GridDim B = d->main_start;
		HIP_GridDim E = d->main_stop;
		if (d->fregions.count)
		{
			B.nx = B.ny = B.nz = 0;
			E = d->dim;
		}
		static const unsigned int xc = getenv("OPENEMS_HIP_FUSED_XC") ? std::max(1, atoi(getenv("OPENEMS_HIP_FUSED_XC"))) : FUSED_XC;
		const bool timed = HIP_KernelTimes::Enabled();
		if (timed)
			HIP_KernelTimes::Begin(d->Stream());
		// the tile is a wavefront wide (see FUSED_TZ); OPENEMS_HIP_FUSED_TY sets the other side
		static const unsigned int ty = getenv("OPENEMS_HIP_FUSED_TY") ? std::max(2, atoi(getenv("OPENEMS_HIP_FUSED_TY"))) : FUSED_TY;
		const HIP_FusedArgs a = { (const float*)d->volt, (const float*)d->curr, d->volt_next, d->curr_next, (const void*)d->index,
		                          (const float*)(c ? d->coeff : d->vv), (const float*)(c ? d->coeff : d->vi),
		                          (const float*)(c ? d->coeff : d->ii), (const float*)(c ? d->coeff : d->iv),
		                          d->coeff_mode, d->dim, B, E, d->fregions, xc };
#ifdef __HIP_PLATFORM_AMD__
		// only AMD has wavefronts of 64; compiling them for NVIDIA would be dead code
		static const unsigned int tz = getenv("OPENEMS_HIP_FUSED_TZ") ? atoi(getenv("OPENEMS_HIP_FUSED_TZ")) : 0;
		if (tz ? (tz>=64) : (d->ctx->warp_size>=64))
			launch_fused_tz<64>(d->Stream(), a, ty, d->index32);
		else
#endif
			launch_fused_tz<32>(d->Stream(), a, ty, d->index32);
		d->CheckLaunch("update_fused");
		if (timed)
			HIP_KernelTimes::End(d->Stream(), "update_fused");
		// only the voltages: UpdateCurrents swaps the currents, once the fix-up has read the old ones
		std::swap(d->volt, d->volt_next);
		return;
	}
	HIP_GridDim B, E;
	unsigned int threads_z;
	d->MainRange(B, E, threads_z);
	HIP_Launch(d, "update_voltages", update_voltages, threads_z, E.ny-B.ny, E.nx-B.nx,
	            d->volt, (const float*)d->curr, (const void*)d->index, (const float*)(c ? d->coeff : d->vv), (const float*)(c ? d->coeff : d->vi), d->coeff_mode, d->dim, B, E, d->zslabs);
}

void GPU_Backend_HIP::UpdateCurrents()
{
	const bool c = d->coeff_mode>0;
	if (d->fused_step>0)
	{
		// done by update_fused, but for the currents next to voltages changed by extensions since
		if (d->fixup_count)
			HIP_Launch(d, "update_currents_nodes", update_currents_nodes, d->fixup_count, 1, 1, (const float*)d->curr, d->curr_next, (const float*)d->volt,
			            (const size_t*)d->fixup, d->fixup_count, (const void*)d->index,
			            (const float*)(c ? d->coeff : d->ii), (const float*)(c ? d->coeff : d->iv), d->coeff_mode, d->dim);
		std::swap(d->curr, d->curr_next);
		return;
	}
	HIP_GridDim B, E;
	unsigned int threads_z;
	d->MainRange(B, E, threads_z);
	HIP_Launch(d, "update_currents", update_currents, threads_z, E.ny-B.ny, E.nx-B.nx,
	            d->curr, (const float*)d->volt, (const void*)d->index, (const float*)(c ? d->coeff : d->ii), (const float*)(c ? d->coeff : d->iv), d->coeff_mode, d->dim, B, E, d->zslabs);
}

static void CheckSize(size_t host_count, size_t device_count)
{
	if (host_count!=device_count)
		throw std::runtime_error("GPU_Backend_HIP: host and device field size mismatch");
}

void GPU_Backend_HIP::DownloadVoltages(ArrayLib::ArrayNIJK<FDTD_FLOAT>& volt)
{
	CheckSize(volt.size(), 3*d->numCells);
	d->PinHostMemory(volt.data(), volt.size()*sizeof(FDTD_FLOAT));
	HIP_Check(hipMemcpyAsync(volt.data(), d->volt, volt.size()*sizeof(float), hipMemcpyDeviceToHost, d->Stream()), "download");
	d->Flush();
}

void GPU_Backend_HIP::DownloadCurrents(ArrayLib::ArrayNIJK<FDTD_FLOAT>& curr)
{
	CheckSize(curr.size(), 3*d->numCells);
	d->PinHostMemory(curr.data(), curr.size()*sizeof(FDTD_FLOAT));
	HIP_Check(hipMemcpyAsync(curr.data(), d->curr, curr.size()*sizeof(float), hipMemcpyDeviceToHost, d->Stream()), "download");
	d->Flush();
}

void GPU_Backend_HIP::UploadVoltages(const ArrayLib::ArrayNIJK<FDTD_FLOAT>& volt)
{
	CheckSize(volt.size(), 3*d->numCells);
	d->PinHostMemory(volt.data(), volt.size()*sizeof(FDTD_FLOAT));
	HIP_Check(hipMemcpyAsync(d->volt, volt.data(), volt.size()*sizeof(float), hipMemcpyHostToDevice, d->Stream()), "upload");
	d->Flush();
}

void GPU_Backend_HIP::UploadCurrents(const ArrayLib::ArrayNIJK<FDTD_FLOAT>& curr)
{
	CheckSize(curr.size(), 3*d->numCells);
	d->PinHostMemory(curr.data(), curr.size()*sizeof(FDTD_FLOAT));
	HIP_Check(hipMemcpyAsync(d->curr, curr.data(), curr.size()*sizeof(float), hipMemcpyHostToDevice, d->Stream()), "upload");
	d->Flush();
}

bool GPU_Backend_HIP::DownloadRange(bool currents, size_t offset, size_t count, FDTD_FLOAT* dst)
{
	if (offset+count > 3*d->numCells)
		return false;
	const float* src = (currents ? d->curr : d->volt) + offset;
	HIP_Check(hipMemcpyAsync(dst, src, count*sizeof(float), hipMemcpyDeviceToHost, d->Stream()), "range download");
	d->Flush();
	return true;
}

void GPU_Backend_HIP::Synchronize()
{
	d->Flush();
}

void* GPU_Backend_HIP::Impl::UploadRegion(const Engine_GatherRegion& region, HIP_GatherRegion& dev)
{
	std::vector<unsigned int> lines;
	for (int n=0; n<3; ++n)
		for (int k=0; k<3; ++k)
			lines.insert(lines.end(), region.lines[n][k].begin(), region.lines[n][k].end());
	unsigned int* buf = TryAlloc<unsigned int>(lines.size());
	if (!buf)
		return NULL;
	HIP_Check(hipMemcpy(buf, lines.data(), lines.size()*sizeof(unsigned int), hipMemcpyHostToDevice), "hipMemcpy");
	size_t pos = 0, start = 0;
	for (int n=0; n<3; ++n)
	{
		for (int k=0; k<3; ++k)
		{
			dev.lines[n][k] = buf + pos;
			dev.size[n][k] = region.lines[n][k].size();
			pos += dev.size[n][k];
		}
		dev.start[n] = start;
		start += region.Size(n);
	}
	dev.total = start;
	return buf;
}

bool GPU_Backend_HIP::SetSnapshotGather(const std::vector<Engine_GatherRegion>& volt_regions, const std::vector<Engine_GatherRegion>& curr_regions)
{
	HIP_Check(hipStreamSynchronize(d->ctx->copy_stream), "stream synchronization");
	d->FreeSnapshots();
	size_t nv = 0, n = 0;
	for (size_t r=0; r<volt_regions.size(); ++r)
		nv += volt_regions[r].Size();
	n = nv;
	for (size_t r=0; r<curr_regions.size(); ++r)
		n += curr_regions[r].Size();
	if ((n==0) || (n>UINT_MAX))
		return false;
	// not fatal: without snapshots the dumps read the fields from the host mirror
	bool ok = (hipMalloc(&d->snap_dev, n*sizeof(float))==hipSuccess)
	          && (hipHostMalloc(&d->snap_host, n*sizeof(float))==hipSuccess);
	const std::vector<Engine_GatherRegion>* regions[2] = {&volt_regions, &curr_regions};
	for (int f=0; f<2 && ok; ++f)
		for (size_t r=0; r<regions[f]->size() && ok; ++r)
		{
			HIP_GatherRegion dev;
			void* buf = d->UploadRegion(regions[f]->at(r), dev);
			ok = (buf!=NULL);
			if (ok)
			{
				d->snap_lines.push_back(buf);
				d->snap_regions[f].push_back(dev);
			}
		}
	if (!ok)
	{
		hipGetLastError();
		d->FreeSnapshots();
		std::cerr << "GPU_Backend_HIP: not enough memory for field snapshots, the field dumps wait for the device" << std::endl;
		return false;
	}
	d->snap_nv = nv;
	d->snap_n = n;
	return true;
}

bool GPU_Backend_HIP::SnapshotFields(const FDTD_FLOAT* &volt, const FDTD_FLOAT* &curr)
{
	if (!d->snap_n)
		return false;
	// the raw values on the work stream, once the last snapshot was downloaded from snap_dev
	if (d->snap_taken)
		HIP_Check(hipStreamWaitEvent(d->Stream(), d->snap_done, 0), "hipStreamWaitEvent");
	size_t offset = 0;
	for (int f=0; f<2; ++f)
		for (size_t r=0; r<d->snap_regions[f].size(); ++r)
		{
			const HIP_GatherRegion& R = d->snap_regions[f][r];
			HIP_Launch(d, "gather_region", gather_region, R.total, 1, 1, (const float*)(f ? d->curr : d->volt), R, d->dim, d->snap_dev + offset);
			offset += R.total;
		}
	HIP_Check(hipEventRecord(d->snap_evaluated, d->Stream()), "hipEventRecord");
	// the download on the copy stream, while the work stream continues with the next timesteps
	hipStream_t copy = d->ctx->copy_stream;
	HIP_Check(hipStreamWaitEvent(copy, d->snap_evaluated, 0), "hipStreamWaitEvent");
	HIP_Check(hipMemcpyAsync(d->snap_host, d->snap_dev, d->snap_n*sizeof(float), hipMemcpyDeviceToHost, copy), "snapshot download");
	HIP_Check(hipEventRecord(d->snap_done, copy), "hipEventRecord");
	d->snap_taken = true;
	volt = d->snap_host;
	curr = d->snap_host + d->snap_nv;
	return true;
}

int GPU_Backend_HIP::AddFieldDFT(bool h_field, const Engine_GatherRegion& region, unsigned int count)
{
	const size_t n = region.Size();
	if ((n==0) || (count==0))
		return -1;
	Impl::FieldDFT dft;
	dft.count = count;
	dft.currents = h_field;
	void* lines = d->UploadRegion(region, dft.region);
	dft.sums = lines ? d->TryAlloc<float2>(n*count) : NULL;
	if (!dft.sums)
	{
		// not fatal: without device sums the dump sums the fields on the host
		d->Free(lines);
		std::cerr << "GPU_Backend_HIP: not enough memory for a frequency domain dump on the device" << std::endl;
		return -1;
	}
	HIP_Check(hipMemsetAsync(dft.sums, 0, n*count*sizeof(float2), d->Stream()), "hipMemset");
	d->field_dfts.push_back(dft);
	return d->field_dfts.size()-1;
}

void GPU_Backend_HIP::AccumulateFieldDFT(int id, const std::vector<std::complex<float>>& weights)
{
	if ((id<0) || ((size_t)id>=d->field_dfts.size()))
		return;
	const Impl::FieldDFT& dft = d->field_dfts[id];
	const unsigned int count = std::min((size_t)dft.count, weights.size());
	for (unsigned int k0=0; k0<count; k0+=DFT_WEIGHTS)
	{
		HIP_DFTWeights W;
		const unsigned int nw = std::min(count-k0, (unsigned int)DFT_WEIGHTS);
		for (unsigned int k=0; k<nw; ++k)
			W.w[k] = make_float2(weights[k0+k].real(), weights[k0+k].imag());
		HIP_Launch(d, "field_dft", field_dft, dft.region.total, 1, 1, (const float*)(dft.currents ? d->curr : d->volt),
		            dft.region, d->dim, W, nw, dft.sums + k0*dft.region.total);
	}
}

bool GPU_Backend_HIP::ReadFieldDFT(int id, unsigned int k, std::vector<std::complex<float>>& sums)
{
	if ((id<0) || ((size_t)id>=d->field_dfts.size()) || (k>=d->field_dfts[id].count))
		return false;
	const Impl::FieldDFT& dft = d->field_dfts[id];
	sums.resize(dft.region.total);
	HIP_Check(hipMemcpyAsync(sums.data(), dft.sums + k*dft.region.total, sums.size()*sizeof(float2), hipMemcpyDeviceToHost, d->Stream()), "download");
	d->Flush();
	return true;
}

void GPU_Backend_HIP::WaitSnapshot()
{
	// called on the dump thread: report errors instead of throwing
	const hipError_t err = hipEventSynchronize(d->snap_done);
	if (err!=hipSuccess)
		std::cerr << "GPU_Backend_HIP: field snapshot failed: " << hipGetErrorString(err) << std::endl;
}

bool GPU_Backend_HIP::CalcFastEnergy(const unsigned int numNodes[3], double& E_energy, double& H_energy)
{
	E_energy = H_energy = 0;
	HIP_GridDim L = {numNodes[0], numNodes[1], numNodes[2]};
	if ((L.nx>d->dim.nx) || (L.ny>d->dim.ny) || (L.nz>d->dim.nz))
		return false;
	const size_t count = (size_t)L.ny*L.nz;
	if (L.nx==0 || count==0)
		return true;
	if (d->energy_count<count)
	{
		d->energy = d->Alloc<float2>(count);
		d->energy_count = count;
		d->energy_host.resize(count);   // small, not pinned
	}
	HIP_Launch(d, "energy_lines", energy_lines, L.nz, L.ny, 1, (const float*)d->volt, (const float*)d->curr, d->energy, d->dim, L);
	HIP_Check(hipMemcpyAsync(d->energy_host.data(), d->energy, count*sizeof(float2), hipMemcpyDeviceToHost, d->Stream()), "energy download");
	d->Flush();
	for (size_t k=0; k<count; ++k)
	{
		E_energy += d->energy_host[k].x;
		H_energy += d->energy_host[k].y;
	}
	return true;
}

// all HIP extensions, see hip_internal.h
static const HIP_ExtensionFactory HIP_EXTENSIONS[] = {
	HIP_CreateExt_Excitation,
	HIP_CreateExt_UPML,
	HIP_CreateExt_Mur_ABC,
	HIP_CreateExt_LorentzMaterial,
	HIP_CreateExt_DebyeMaterial,
	HIP_CreateExt_LumpedRLC,
	HIP_CreateExt_TFSF,
	HIP_CreateExt_Absorbing_BC,
	HIP_CreateExt_SteadyState,
	HIP_CreateExt_Cylinder,
};

GPU_Extension* GPU_Backend_HIP::CreateExtension(Engine_Extension* eng_ext, Engine* eng)
{
	// the fused step (see update_fused) knows the UPML, the excitation, lumped RLC elements and the steady-state detection
	if (!dynamic_cast<Engine_Ext_UPML*>(eng_ext) && !dynamic_cast<Engine_Ext_Excitation*>(eng_ext) && !dynamic_cast<Engine_Ext_LumpedRLC*>(eng_ext)
	    && !dynamic_cast<Engine_Ext_SteadyState*>(eng_ext))
		d->fused_blockers++;
	for (size_t n=0; n<sizeof(HIP_EXTENSIONS)/sizeof(HIP_EXTENSIONS[0]); ++n)
	{
		if (HIP_EXTENSIONS[n]==NULL)
			continue;
		GPU_Extension* gpu_ext = HIP_EXTENSIONS[n](d, eng_ext, eng);
		if (gpu_ext)
			return gpu_ext;
	}
	d->fused_blockers++;   // host fallback
	return NULL;
}
