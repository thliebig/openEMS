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
	const unsigned int sn = N.nx*N.ny*N.nz;
	const unsigned int i  = nijk(N, 0, x, y, z);
	const unsigned int xm = (x>0) ? N.ny*N.nz : 0;   // shift to the previous line, none on the first
	const unsigned int ym = (y>0) ? N.nz : 0;
	const unsigned int zm = (z>0) ? 1 : 0;
	const HIP_MainCoeff C = main_coeff(index, ca, cb, mode, sn, i, 0);
	const float* vv = C.a;
	const float* vi = C.b;
	const unsigned int s = C.s;
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
	const unsigned int sn = N.nx*N.ny*N.nz;
	const unsigned int i  = nijk(N, 0, x, y, z);
	const unsigned int xp = N.ny*N.nz;
	const unsigned int yp = N.nz;
	const HIP_MainCoeff C = main_coeff(index, ca, cb, mode, sn, i, 6);
	const float* ii = C.a;
	const float* iv = C.b;
	const unsigned int s = C.s;
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
__device__ __forceinline__ void fused_volt(const float* __restrict__ volt, const float* __restrict__ curr,
                                           const void* index, const float* va, const float* vb, unsigned int mode,
                                           const HIP_GridDim& N, unsigned int x, unsigned int y, unsigned int z, float v[3])
{
	const unsigned int sn = N.nx*N.ny*N.nz;
	const unsigned int i  = nijk(N, 0, x, y, z);
	const unsigned int xm = (x>0) ? N.ny*N.nz : 0;
	const unsigned int ym = (y>0) ? N.nz : 0;
	const unsigned int zm = (z>0) ? 1 : 0;
	const HIP_MainCoeff C = main_coeff(index, va, vb, mode, sn, i, 0);
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
template<unsigned int TZ, unsigned int TY>
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
	const unsigned int sn = N.nx*N.ny*N.nz;

	// the new voltages at x: computed on main nodes, else from the UPML regions (or unused)
	auto volt_at = [&](unsigned int x, float v[3])
	{
		if (main_yz && (x<E.nx))
		{
			const int r = F.count ? fused_region_of(F, x, y, z) : -1;
			if (r<0)
				fused_volt(volt_in, curr_in, index, va, vb, mode, N, x, y, z, v);
			else
				upml_volt_value(volt_in, curr_in, F.volt[r], N, F.R[r], x, y, z, owner && (x<xe), v);
			if (owner && (x<xe))
			{
				const unsigned int i = nijk(N, 0, x, y, z);
				volt_out[i] = v[0]; volt_out[sn+i] = v[1]; volt_out[2*sn+i] = v[2];
			}
		}
		else if (in_grid && (x<N.nx))
		{
			const unsigned int i = nijk(N, 0, x, y, z);
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
				upml_curr_value(curr_in, curr_out, F.curr[r], N, F.R[r], x, y, z, update, curl);
			else if (update)
			{
				const unsigned int i = nijk(N, 0, x, y, z);
				const HIP_MainCoeff C = main_coeff(index, ia, ib, mode, sn, i, 6);
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

template<unsigned int TZ, unsigned int TY>
static void launch_fused(hipStream_t stream, const HIP_FusedArgs& a)
{
	const dim3 block(TZ, TY, 1);
	const dim3 grid((a.E.nz-a.B.nz+TZ-2)/(TZ-1), (a.E.ny-a.B.ny+TY-2)/(TY-1), (a.E.nx-a.B.nx+a.xc-1)/a.xc);
	update_fused<TZ, TY><<<grid, block, 0, stream>>>(a.volt_in, a.curr_in, a.volt_out, a.curr_out, a.index,
	                                                 a.va, a.vb, a.ia, a.ib, a.mode, a.N, a.B, a.E, a.F, a.xc);
}

template<unsigned int TZ>
static void launch_fused_ty(hipStream_t stream, const HIP_FusedArgs& a, unsigned int ty)
{
	switch (ty)
	{
		case 4:  launch_fused<TZ, 4>(stream, a); break;
		case 16: launch_fused<TZ, 16>(stream, a); break;
		default: launch_fused<TZ, 8>(stream, a); break;
	}
}

// the dumped values of a snapshot: entries [0, nv) from the voltages, [nv, n) from the currents
// (see GPU_GatherEntry). Evaluated here so that only the dumped values go over PCIe, not the
// whole fields. Double precision without FMA contraction (--fmad=false), as on the host.
__device__ inline double gather_raw(const float* f, const GPU_GatherEntry& e, int k)
{
	if (e.delta[k])
		return (double)f[e.idx[k]] / e.delta[k];
	return 0.0;
}

__device__ inline float gather_value(const float* f, const GPU_GatherEntry& e)
{
	double v = 0;
	switch (e.form)
	{
	case GPU_GatherEntry::RAW:
		v = gather_raw(f, e, 0);
		break;
	case GPU_GatherEntry::LERP:
		v = gather_raw(f, e, 0)*(1.0-e.rel) + gather_raw(f, e, 1)*e.rel;
		break;
	case GPU_GatherEntry::AVG4:
		v = gather_raw(f, e, 0);
		v+= gather_raw(f, e, 1);
		v+= gather_raw(f, e, 2);
		v+= gather_raw(f, e, 3);
		v/= 4;
		break;
	default:
		break;
	}
	return v;
}

__global__ void gather_dumps(const float* __restrict__ volt, const float* __restrict__ curr, const GPU_GatherEntry* __restrict__ entries,
                             unsigned int nv, unsigned int n, float* __restrict__ out)
{
	const unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i>=n)
		return;
	out[i] = gather_value((i<nv) ? volt : curr, entries[i]);
}

// frequency domain dumps: sums[k*n + i] += value(i) * w[k] for the weights of one timestep,
// as ProcessFieldsFD::Process() (a float times a complex float, then a complex sum)
#define DFT_WEIGHTS 32   // weights per launch, passed by value in the kernel parameters
struct HIP_DFTWeights { float2 w[DFT_WEIGHTS]; };

__global__ void field_dft(const float* __restrict__ f, const GPU_GatherEntry* __restrict__ entries, unsigned int n,
                          HIP_DFTWeights W, unsigned int nw, float2* __restrict__ sums)
{
	const unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i>=n)
		return;
	const float v = gather_value(f, entries[i]);
	for (unsigned int k=0; k<nw; ++k)
	{
		float2 s = sums[(size_t)k*n + i];
		s.x += v*W.w[k].x;
		s.y += v*W.w[k].y;
		sums[(size_t)k*n + i] = s;
	}
}

// the currents of the given main nodes again (update_currents), after the voltage extensions
// changed voltages they read (see GPU_Backend_HIP::Impl::fixup)
__global__ void update_currents_nodes(const float* __restrict__ curr_in, float* __restrict__ curr_out, const float* __restrict__ volt,
                                      const unsigned int* nodes, unsigned int count,
                                      const void* index, const float* ia, const float* ib, unsigned int mode, HIP_GridDim N)
{
	const unsigned int k = blockIdx.x*blockDim.x + threadIdx.x;
	if (k>=count)
		return;
	const unsigned int sn = N.nx*N.ny*N.nz;
	const unsigned int i  = nodes[k];
	const unsigned int xp = N.ny*N.nz;
	const unsigned int yp = N.nz;
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
	const unsigned int sn = N.nx*N.ny*N.nz;
	float E = 0, H = 0;
	for (unsigned int x=0; x<L.nx; ++x)
	{
		const unsigned int i = nijk(N, 0, x, y, z);
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
	fixup_count = 0;
	snap_entries = NULL;
	snap_nv = snap_n = 0;
	snap_dev = NULL;
	snap_host[0] = snap_host[1] = NULL;
	snap_last = -1;
	HIP_Check(hipEventCreateWithFlags(&snap_evaluated, hipEventDisableTiming), "hipEventCreate");
	HIP_Check(hipEventCreateWithFlags(&snap_done[0], hipEventDisableTiming), "hipEventCreate");
	HIP_Check(hipEventCreateWithFlags(&snap_done[1], hipEventDisableTiming), "hipEventCreate");
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
	hipEventDestroy(snap_done[0]);
	hipEventDestroy(snap_done[1]);
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
	std::set<unsigned int> nodes;
	for (size_t k=0; k<volt_modified.size(); ++k)
	{
		const unsigned int node = volt_modified[k] % numCells;
		const unsigned int p[3] = {node/(dim.ny*dim.nz), (node/dim.nz)%dim.ny, node%dim.nz};
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
			nodes.insert((c[0]*dim.ny + c[1])*dim.nz + c[2]);
		}
	}
	std::vector<unsigned int> list;
	for (std::set<unsigned int>::const_iterator it=nodes.begin(); it!=nodes.end(); ++it)
	{
		const unsigned int c[3] = {*it/(dim.ny*dim.nz), (*it/dim.nz)%dim.ny, *it%dim.nz};
		if ((c[0]>=main_start.nx) && (c[0]<main_stop.nx) && (c[1]>=main_start.ny) && (c[1]<main_stop.ny) && (c[2]>=main_start.nz) && (c[2]<main_stop.nz))
			list.push_back(*it);   // main nodes only, the region kernels update the others after the extensions
		else if (fregions.count)
			list.push_back(*it);
	}
	fixup_count = list.size();
	fixup = Alloc<unsigned int>(fixup_count, list.data());

	// both buffers start with the current fields
	volt_next = Alloc<float>(3*numCells);
	curr_next = Alloc<float>(3*numCells);
	HIP_Check(hipMemcpyAsync(volt_next, volt, 3*numCells*sizeof(float), hipMemcpyDeviceToDevice, Stream()), "fused step buffers");
	HIP_Check(hipMemcpyAsync(curr_next, curr, 3*numCells*sizeof(float), hipMemcpyDeviceToDevice, Stream()), "fused step buffers");
	Flush();
	fused_step = 1;
	std::cout << "GPU_Backend_HIP: fused voltage and current updates";
	if (fregions.count)
		std::cout << " (" << fregions.count << " UPML regions in the kernel)";
	std::cout << std::endl;
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
	hipFree(snap_entries);
	hipFree(snap_dev);
	hipHostFree(snap_host[0]);
	hipHostFree(snap_host[1]);
	snap_entries = NULL;
	snap_dev = snap_host[0] = snap_host[1] = NULL;
	snap_nv = snap_n = 0;
	snap_last = -1;
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

	// the kernels index with 32 bit
	if (3*d->numCells > std::numeric_limits<unsigned int>::max())
	{
		std::cerr << "GPU_Backend_HIP::Init: Error: the mesh is too large for 32 bit indexing" << std::endl;
		return false;
	}

	const size_t count = 3*d->numCells;
	d->volt = d->Alloc<float>(count);
	d->curr = d->Alloc<float>(count);

	// upload the final operator coefficients, including all changes by operator extensions
	std::vector<float> vv(count), vi(count), ii(count), iv(count);
	unsigned int pos[3];
	size_t idx = 0;
	for (int n=0; n<3; ++n)
		for (pos[0]=0; pos[0]<numLines[0]; ++pos[0])
			for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
				for (pos[2]=0; pos[2]<numLines[2]; ++pos[2], ++idx)
				{
					vv[idx] = op->GetVV(n, pos[0], pos[1], pos[2]);
					vi[idx] = op->GetVI(n, pos[0], pos[1], pos[2]);
					ii[idx] = op->GetII(n, pos[0], pos[1], pos[2]);
					iv[idx] = op->GetIV(n, pos[0], pos[1], pos[2]);
				}

	// Most nodes share one of a few coefficient sets (same material and mesh spacing): store
	// every distinct set of the 12 coefficients of a node once and a set index per node
	const size_t sn = d->numCells;
	const std::vector<float>* src[4] = {&vv, &vi, &ii, &iv};
	GPU_CoeffSets sets;
	if (GPU_FindSets(sn, 12, [&](size_t i, float* values)
	    {
		    for (int c=0; c<4; ++c)
			    for (int n=0; n<3; ++n)
				    values[3*c+n] = (*src[c])[n*sn+i];
	    }, sets))
	{
		d->coeff_mode = sets.mode;
		d->index = d->AllocIndex(sets);
		d->coeff = d->Alloc<float>(sets.table.size(), sets.table.data());
		std::cout << "GPU_Backend_HIP: " << sets.count << " distinct coefficient sets, compressed update coefficients ("
		          << (sets.mode==1 ? 16 : 32) << " bit index)" << std::endl;
		return true;
	}

	d->coeff_mode = 0;
	d->vv = d->Alloc<float>(count, vv.data());
	d->vi = d->Alloc<float>(count, vi.data());
	d->ii = d->Alloc<float>(count, ii.data());
	d->iv = d->Alloc<float>(count, iv.data());
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
			launch_fused_ty<64>(d->Stream(), a, ty);
		else
#endif
			launch_fused_ty<32>(d->Stream(), a, ty);
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
			            (const unsigned int*)d->fixup, d->fixup_count, (const void*)d->index,
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

bool GPU_Backend_HIP::SetSnapshotGather(const std::vector<GPU_GatherEntry>& volt_entries, const std::vector<GPU_GatherEntry>& curr_entries)
{
	HIP_Check(hipStreamSynchronize(d->ctx->copy_stream), "stream synchronization");
	d->FreeSnapshots();
	const size_t n = volt_entries.size() + curr_entries.size();
	if ((n==0) || (n>UINT_MAX))
		return false;
	// not fatal: without snapshots the dumps read the fields from the host mirror
	if ((hipMalloc(&d->snap_entries, n*sizeof(GPU_GatherEntry))!=hipSuccess) || (hipMalloc(&d->snap_dev, n*sizeof(float))!=hipSuccess)
	    || (hipHostMalloc(&d->snap_host[0], n*sizeof(float))!=hipSuccess) || (hipHostMalloc(&d->snap_host[1], n*sizeof(float))!=hipSuccess))
	{
		hipGetLastError();
		d->FreeSnapshots();
		std::cerr << "GPU_Backend_HIP: not enough memory for field snapshots, the field dumps wait for the device" << std::endl;
		return false;
	}
	d->snap_nv = volt_entries.size();
	d->snap_n = n;
	HIP_Check(hipMemcpy(d->snap_entries, volt_entries.data(), volt_entries.size()*sizeof(GPU_GatherEntry), hipMemcpyHostToDevice), "hipMemcpy");
	HIP_Check(hipMemcpy(d->snap_entries + d->snap_nv, curr_entries.data(), curr_entries.size()*sizeof(GPU_GatherEntry), hipMemcpyHostToDevice), "hipMemcpy");
	return true;
}

bool GPU_Backend_HIP::SnapshotFields(unsigned int slot, const FDTD_FLOAT* &volt, const FDTD_FLOAT* &curr)
{
	if ((slot>1) || !d->snap_n)
		return false;
	// the dumped values on the work stream, once the last snapshot was downloaded from snap_dev
	if (d->snap_last>=0)
		HIP_Check(hipStreamWaitEvent(d->Stream(), d->snap_done[d->snap_last], 0), "hipStreamWaitEvent");
	HIP_Launch(d, "gather_dumps", gather_dumps, d->snap_n, 1, 1, (const float*)d->volt, (const float*)d->curr,
	            (const GPU_GatherEntry*)d->snap_entries, d->snap_nv, d->snap_n, d->snap_dev);
	HIP_Check(hipEventRecord(d->snap_evaluated, d->Stream()), "hipEventRecord");
	// the download on the copy stream, while the work stream continues with the next timesteps
	hipStream_t copy = d->ctx->copy_stream;
	HIP_Check(hipStreamWaitEvent(copy, d->snap_evaluated, 0), "hipStreamWaitEvent");
	HIP_Check(hipMemcpyAsync(d->snap_host[slot], d->snap_dev, d->snap_n*sizeof(float), hipMemcpyDeviceToHost, copy), "snapshot download");
	HIP_Check(hipEventRecord(d->snap_done[slot], copy), "hipEventRecord");
	d->snap_last = slot;
	volt = d->snap_host[slot];
	curr = d->snap_host[slot] + d->snap_nv;
	return true;
}

int GPU_Backend_HIP::AddFieldDFT(const std::vector<GPU_GatherEntry>& entries, bool currents, unsigned int count)
{
	const size_t n = entries.size();
	if ((n==0) || (count==0) || (n>UINT_MAX))
		return -1;
	Impl::FieldDFT dft;
	dft.n = n;
	dft.count = count;
	dft.currents = currents;
	try
	{
		dft.entries = d->Alloc<GPU_GatherEntry>(n, entries.data());
		dft.sums = d->Alloc<float2>(n*count);
	}
	catch (std::runtime_error&)
	{
		// not fatal: without device sums the dump accumulates the fields on the host
		hipGetLastError();
		std::cerr << "GPU_Backend_HIP: not enough memory for a frequency domain dump on the device" << std::endl;
		return -1;
	}
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
		HIP_Launch(d, "field_dft", field_dft, dft.n, 1, 1, (const float*)(dft.currents ? d->curr : d->volt),
		            (const GPU_GatherEntry*)dft.entries, dft.n, W, nw, dft.sums + (size_t)k0*dft.n);
	}
}

bool GPU_Backend_HIP::ReadFieldDFT(int id, std::vector<std::complex<float>>& sums)
{
	if ((id<0) || ((size_t)id>=d->field_dfts.size()))
		return false;
	const Impl::FieldDFT& dft = d->field_dfts[id];
	sums.resize((size_t)dft.n*dft.count);
	HIP_Check(hipMemcpyAsync(sums.data(), dft.sums, sums.size()*sizeof(float2), hipMemcpyDeviceToHost, d->Stream()), "download");
	d->Flush();
	return true;
}

void GPU_Backend_HIP::WaitSnapshot(unsigned int slot)
{
	// called on the dump thread: report errors instead of throwing
	if (slot<2)
	{
		const hipError_t err = hipEventSynchronize(d->snap_done[slot]);
		if (err!=hipSuccess)
			std::cerr << "GPU_Backend_HIP: field snapshot failed: " << hipGetErrorString(err) << std::endl;
	}
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
