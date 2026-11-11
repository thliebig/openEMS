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
#include <cstring>
#include <map>
#include <set>
#include <sstream>
#include <iostream>
#include <limits>
#include <stdexcept>

#include "hip_internal.h"
#include "gpu_build_id.h"
#include "FDTD/operator.h"
#include "FDTD/extensions/engine_ext_upml.h"
#include "FDTD/extensions/engine_ext_excitation.h"
#include "FDTD/extensions/engine_ext_lumpedRLC.h"
#include "FDTD/extensions/engine_ext_steadystate.h"
#include "FDTD/extensions/engine_ext_lorentzmaterial.h"
#include "FDTD/extensions/engine_ext_mur_abc.h"
#include "FDTD/extensions/engine_ext_tfsf.h"
#include "FDTD/extensions/engine_ext_debyematerial.h"
#include "FDTD/extensions/engine_ext_absorbing_bc.h"

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
	// an event pair per device: an event can only be recorded on a stream of its own device
	std::map<int, std::pair<hipEvent_t, hipEvent_t>> events;
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
	//! the event pair of the current device
	std::pair<hipEvent_t, hipEvent_t>& Events(int& device)
	{
		device = 0;
		(void)hipGetDevice(&device);
		std::pair<hipEvent_t, hipEvent_t>& e = events[device];
		if (!e.first)
		{
			(void)hipEventCreate(&e.first);
			(void)hipEventCreate(&e.second);
		}
		return e;
	}
};
KernelTimeTable& KernelTimes() {static KernelTimeTable table; return table;}
}

bool HIP_KernelTimes::Enabled()
{
	static const bool enabled = getenv("OPENEMS_HIP_KERNEL_TIMES") && (atoi(getenv("OPENEMS_HIP_KERNEL_TIMES"))!=0);
	return enabled;
}

// Begin() and End() take a stream of the current device
void HIP_KernelTimes::Begin(hipStream_t stream)
{
	int device;
	(void)hipEventRecord(KernelTimes().Events(device).first, stream);
}

void HIP_KernelTimes::End(hipStream_t stream, const char* name)
{
	KernelTimeTable& t = KernelTimes();
	int device;
	const std::pair<hipEvent_t, hipEvent_t>& e = t.Events(device);
	(void)hipEventRecord(e.second, stream);
	(void)hipEventSynchronize(e.second);
	float ms = 0;
	(void)hipEventElapsedTime(&ms, e.first, e.second);
	// the kernels of devices other than the first one apart (parts on several devices)
	std::pair<double, size_t>& entry = t.times[device ? std::string(name) + " (device " + std::to_string(device) + ")" : std::string(name)];
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

void HIP_Context::MakeCurrent() const
{
	HIP_Check(hipSetDevice(device), "hipSetDevice");
}

HIP_Context::~HIP_Context()
{
	(void)hipSetDevice(device);
	if (stream)
	{
		(void)hipStreamSynchronize(stream);
		(void)hipStreamDestroy(stream);
	}
	if (copy_stream)
	{
		(void)hipStreamSynchronize(copy_stream);
		(void)hipStreamDestroy(copy_stream);
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
	x0 = x1 = x_off = gnx = 0;
	pin_flags = hipHostRegisterDefault;
	// on the current device, that of the part (see New())
	HIP_Check(hipEventCreateWithFlags(&snap_evaluated, hipEventDisableTiming), "hipEventCreate");
	HIP_Check(hipEventCreateWithFlags(&snap_done, hipEventDisableTiming), "hipEventCreate");
	HIP_Check(hipEventCreateWithFlags(&halo_sent, hipEventDisableTiming), "hipEventCreate");
	HIP_Check(hipEventCreateWithFlags(&halo_received, hipEventDisableTiming), "hipEventCreate");
}

GPU_Backend_HIP::Impl::~Impl()
{
	if (ctx)
	{
		(void)hipSetDevice(ctx->device);
		(void)hipStreamSynchronize(ctx->stream);
		(void)hipStreamSynchronize(ctx->copy_stream);
	}
	FreeSnapshots();
	(void)hipEventDestroy(snap_evaluated);
	(void)hipEventDestroy(snap_done);
	(void)hipEventDestroy(halo_sent);
	(void)hipEventDestroy(halo_received);
	for (size_t n=0; n<m_Allocations.size(); ++n)
		(void)hipFree(m_Allocations.at(n));
	for (std::set<void*>::iterator it=m_Pinned.begin(); it!=m_Pinned.end(); ++it)
		(void)hipHostUnregister(*it);
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
	ctx->MakeCurrent();
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
	if (ctx)
		(void)hipSetDevice(ctx->device);
	(void)hipFree(snap_dev);
	(void)hipHostFree(snap_host);
	snap_dev = snap_host = NULL;
	snap_nv = snap_n = 0;
	snap_taken = false;
	snap_copies.clear();
}

void GPU_Backend_HIP::Impl::Flush()
{
	ctx->MakeCurrent();
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
	ctx->MakeCurrent();
	if (hipHostRegister(ptr, bytes, pin_flags)==hipSuccess)
		m_Pinned.insert(ptr);
	else
		(void)hipGetLastError();
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

int GPU_Backend_HIP::DeviceCount()
{
	int count = 0;
	if (hipGetDeviceCount(&count)!=hipSuccess)
	{
		(void)hipGetLastError();
		return 0;
	}
	return count;
}

bool GPU_Backend_HIP::DeviceInfo(int n, GPU_DeviceInfo* info)
{
	hipDeviceProp_t prop;
	if (hipGetDeviceProperties(&prop, n)!=hipSuccess)
	{
		(void)hipGetLastError();
		return false;
	}
	snprintf(info->name, sizeof(info->name), "%s", prop.name);
	info->memory = prop.totalGlobalMem;
	info->integrated = prop.integrated!=0;
	return true;
}

namespace
{
GPU_Backend* NewHIPBackend(int n)
{
	return GPU_Backend_HIP::New(n);
}

GPU_Backend* NewHIPBackendParts(const int* n, int count)
{
	return GPU_Backend_HIP::New(n, count);
}
}

const GPU_Plugin* GPU_Backend_HIP::Plugin()
{
	static const GPU_Plugin plugin = {OPENEMS_GPU_BUILD_ID, OPENEMS_GPU_PLUGIN_NAME, DeviceCount, DeviceInfo, NewHIPBackend, NewHIPBackendParts};
	return &plugin;
}

#ifdef OPENEMS_HIP_ARCHITECTURES
namespace
{
//! Whether this build has code for the AMD architecture \a arch (gcnArchName, e.g. "gfx1100:sramecc-:xnack-")
bool BuiltForArch(std::string arch)
{
	arch = arch.substr(0, arch.find(':'));
	// a generic target covers a family of architectures (see the AMDGPU docs of LLVM)
	static const char* const generic[][2] = {
		{"gfx9-generic", "gfx90"}, {"gfx9-4-generic", "gfx94"}, {"gfx9-4-generic", "gfx95"},
		{"gfx10-1-generic", "gfx101"}, {"gfx10-3-generic", "gfx103"}, {"gfx11-generic", "gfx11"}, {"gfx12-generic", "gfx12"}};
	const std::string list = OPENEMS_HIP_ARCHITECTURES;
	size_t start = 0;
	while (start<list.size())
	{
		size_t end = list.find(',', start);
		if (end==std::string::npos)
			end = list.size();
		const std::string target = list.substr(start, std::min(list.find(':', start), end)-start);
		start = end+1;
		if (target==arch)
			return true;
		for (const auto& g : generic)
			if ((target==g[0]) && (arch.compare(0, strlen(g[1]), g[1])==0) && (arch!="gfx908") && (arch!="gfx90a"))
				return true;
	}
	return false;
}
}
#endif

//! The state of a grid (part) on HIP device \a device, with its own streams, or NULL if it cannot run there
static GPU_Backend_HIP::Impl* NewPart(int device)
{
	typedef GPU_Backend_HIP::Impl Impl;
	if ((device<0) || (device>=GPU_Backend_HIP::DeviceCount()) || (hipSetDevice(device)!=hipSuccess))
	{
		(void)hipGetLastError();
		return NULL;
	}
	Impl* impl = new Impl();
	impl->ctx = std::make_shared<HIP_Context>();
	impl->ctx->device = device;
	hipDeviceProp_t prop;
	HIP_Check(hipGetDeviceProperties(&prop, impl->ctx->device), "hipGetDeviceProperties");
	impl->ctx->name = prop.name;
	impl->ctx->warp_size = prop.warpSize > 0 ? prop.warpSize : 32;
#ifdef OPENEMS_HIP_ARCHITECTURES
	// Some AMD runtimes (ROCm 7.2 with assertions) abort in the query below on a device
	// the build has no code for, instead of returning an error: compare with the list first.
	if (!BuiltForArch(prop.gcnArchName))
	{
		std::cerr << "GPU_Backend_HIP: no kernels for " << impl->ctx->name << " (" << prop.gcnArchName
		          << ") in this build, which is for " << OPENEMS_HIP_ARCHITECTURES << std::endl;
		delete impl;
		return NULL;
	}
#endif
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
		(void)hipGetLastError();
		delete impl;
		return NULL;
	}
	return impl;
}

GPU_Backend_HIP* GPU_Backend_HIP::New(int device)
{
	return New(&device, 1);
}

GPU_Backend_HIP* GPU_Backend_HIP::New(const int* devices, int count)
{
	std::vector<Impl*> parts;
	for (int n=0; n<count; ++n)
	{
		Impl* impl = NewPart(devices[n]);
		if (!impl)
		{
			for (Impl* p : parts)
				delete p;
			return NULL;
		}
		parts.push_back(impl);
	}
	if (parts.empty())
		return NULL;
	return new GPU_Backend_HIP(parts);
}

GPU_Backend* GPU_Backend_HIP::NewSubGridBackend()
{
	if (m_Parts.size()>1)
		throw std::runtime_error("GPU_Backend_HIP: a cylindrical multi-grid cannot be split over several GPU devices, use one (--gpu-devices)");
	// no fused step with sub-grids (cylindrical multi-grid)
	d->fused_blockers++;
	d->ctx->MakeCurrent();
	Impl* impl = new Impl();
	impl->fused_blockers++;
	impl->ctx = d->ctx;
	return new GPU_Backend_HIP(std::vector<Impl*>(1, impl));
}

GPU_Backend_HIP::GPU_Backend_HIP(const std::vector<Impl*>& parts)
{
	m_Parts = parts;
	d = m_Parts.front();
	m_SnapshotParts = false;
	m_SnapshotNV = 0;
	if (m_Parts.size()>1)
		d->pin_flags = hipHostRegisterPortable;   // the first part pins the host mirror for all
}

GPU_Backend_HIP::~GPU_Backend_HIP()
{
	for (Impl* p : m_Parts)
		delete p;
}

std::string GPU_Backend_HIP::GetName() const
{
	if (m_Parts.size()>1)
		return "HIP (" + d->ctx->name + ", " + std::to_string(m_Parts.size()) + " parts)";
	return "HIP (" + d->ctx->name + ")";
}

bool GPU_Backend_HIP::Init(const Operator* op)
{
	unsigned int numLines[3];
	for (int n=0; n<3; ++n)
		numLines[n] = op->GetNumberOfLines(n, true);
	if (m_Parts.size()==1)
		return InitPart(d, op, numLines, 0, numLines[0]);

	// Parts along x, at least 3 lines each: with the same number of lines (cells),
	// or between the lines of --gpu-domain-cuts (else OPENEMS_GPU_DOMAIN_CUTS)
	const unsigned int count = m_Parts.size();
	if (numLines[0]<3*count)
	{
		std::cerr << "GPU_Backend_HIP::Init: Error: " << numLines[0] << " mesh lines in x are too few for "
		          << count << " parts (at least 3 lines each)" << std::endl;
		return false;
	}
	std::vector<unsigned int> cuts;
	for (unsigned int n=0; n<=count; ++n)
		cuts.push_back(n*(numLines[0]/count) + std::min(n, numLines[0]%count));
	std::string option = "--gpu-domain-cuts=";
	std::string list = g_settings.GetGPUDomainCuts();
	if (list.empty() && getenv("OPENEMS_GPU_DOMAIN_CUTS"))
	{
		option = "OPENEMS_GPU_DOMAIN_CUTS=";
		list = getenv("OPENEMS_GPU_DOMAIN_CUTS");
	}
	if (!list.empty())
	{
		cuts.clear();
		std::istringstream items(list);
		std::string item;
		bool valid = true;
		while (std::getline(items, item, ','))
		{
			char* end = NULL;
			const long cut = strtol(item.c_str(), &end, 10);
			valid = valid && (end!=item.c_str()) && (*end=='\0') && (cut>=0) && (cut<=(long)numLines[0]);
			cuts.push_back(valid ? cut : 0);
		}
		valid = valid && (cuts.size()==count+1) && (cuts.front()==0) && (cuts.back()==numLines[0]);
		for (unsigned int n=0; valid && (n<count); ++n)
			valid = (cuts[n+1]>=cuts[n]+3);
		if (!valid)
		{
			std::cerr << "GPU_Backend_HIP::Init: Error: " << option << list << ": " << count << " parts need " << count+1
			          << " x lines from 0 to " << numLines[0] << " (the number of mesh lines in x), at least 3 lines apart" << std::endl;
			return false;
		}
	}
	for (unsigned int n=0; n<count; ++n)
	{
		const unsigned int x0 = cuts[n];
		const unsigned int x1 = cuts[n+1];
		Impl* p = m_Parts[n];
		p->fused_blockers++;   // no fused step with parts
		if (!InitPart(p, op, numLines, x0, x1))
			return false;
		if (g_settings.GetVerboseLevel()>0)
			std::cout << "GPU_Backend_HIP: part " << n << " on device " << p->ctx->device << " (" << p->ctx->name << "): x lines "
			          << x0 << " to " << x1-1 << ", " << (size_t)(x1-x0)*numLines[1]*numLines[2] << " cells" << std::endl;
	}
	return true;
}

bool GPU_Backend_HIP::InitPart(Impl* p, const Operator* op, const unsigned int numLines[3], unsigned int x0, unsigned int x1)
{
	// the owned lines and a halo line at each cut, see Impl::x0
	const unsigned int lo = (x0>0) ? 1 : 0;
	const unsigned int hi = (x1<numLines[0]) ? 1 : 0;
	p->x0 = x0;
	p->x1 = x1;
	p->x_off = x0 - lo;
	p->gnx = numLines[0];
	p->dim.nx = x1 - x0 + lo + hi;
	p->dim.ny = numLines[1];
	p->dim.nz = numLines[2];
	p->numCells = (size_t)p->dim.nx*numLines[1]*numLines[2];
	p->main_start.nx = p->main_start.ny = p->main_start.nz = 0;
	p->main_stop = p->dim;
	// the main updates cover the owned lines: at a cut x>0 and x+1<nx hold, as inside the grid
	p->main_start.nx = lo;
	p->main_stop.nx = lo + x1 - x0;
	p->upml_fused = -1;
	p->ResetZSlabs();

	// the main updates launch one block per x line along z, which the device limits (64k on AMD and NVIDIA);
	// from the properties: hipDeviceGetAttribute calls the CUDA driver API on NVIDIA, which the
	// plugin does not link (MSVC keeps that call, GCC folds it away)
	hipDeviceProp_t prop;
	HIP_Check(hipGetDeviceProperties(&prop, p->ctx->device), "hipGetDeviceProperties");
	int max_grid_z = prop.maxGridSize[2];
	if ((max_grid_z>0) && (p->dim.nx>(unsigned int)max_grid_z))
	{
		std::cerr << "GPU_Backend_HIP::Init: Error: " << p->dim.nx << " mesh lines in x, the device allows at most "
		          << max_grid_z << "; swap the axes, or use a CPU engine" << std::endl;
		return false;
	}

	// The fused step indexes with 32 bit on NVIDIA while the fields fit (faster there), with
	// 64 bit on AMD; OPENEMS_HIP_INDEX=32 or 64 overrides.
	const bool fits32 = 3*p->numCells <= UINT_MAX;
#ifdef __HIP_PLATFORM_AMD__
	p->index32 = false;
#else
	p->index32 = fits32;
#endif
	if (const char* env = getenv("OPENEMS_HIP_INDEX"))
		p->index32 = (atoi(env)==32);
	if (p->index32 && !fits32)
	{
		std::cerr << "GPU_Backend_HIP: the mesh is too large for 32 bit indices, the fused step uses 64 bit" << std::endl;
		p->index32 = false;
	}

	const size_t count = 3*p->numCells;
	p->volt = p->Alloc<float>(count);
	p->curr = p->Alloc<float>(count);

	// upload the final operator coefficients, including all changes by operator extensions.
	// Most nodes share one of a few coefficient sets (same material and mesh spacing): store
	// every distinct set of the 12 coefficients of a node once and a set index per node.
	// A part uploads those of its local nodes, the global node is offset by its first line.
	const size_t node_off = (size_t)p->x_off*numLines[1]*numLines[2];
	GPU_CoeffSets sets;
	if (GPU_FindSets(p->numCells, 12, [&](size_t i, float* values) {GPU_NodeCoefficients(op, numLines, node_off + i, values);}, sets))
	{
		p->coeff_mode = sets.mode;
		p->index = p->AllocIndex(sets);
		p->coeff = p->Alloc<float>(sets.table.size(), sets.table.data());
		if (g_settings.GetVerboseLevel()>1)
			std::cout << "GPU_Backend_HIP: " << sets.count << " distinct coefficient sets, compressed update coefficients ("
			          << (sets.mode==1 ? 16 : 32) << " bit index)" << std::endl;
		return true;
	}

	// the full arrays, one host copy at a time
	p->coeff_mode = 0;
	std::vector<float> values;
	auto upload = [&](int c)
	{
		if (m_Parts.size()==1)
			GPU_CoefficientArray(op, numLines, c, values);
		else
		{
			// the local nodes of the part, from the 12 coefficients of each (see GPU_NodeCoefficients())
			values.resize(count);
			float node[12];
			for (size_t i=0; i<p->numCells; ++i)
			{
				GPU_NodeCoefficients(op, numLines, node_off + i, node);
				for (int n=0; n<3; ++n)
					values[n*p->numCells + i] = node[3*c + n];
			}
		}
		return p->Alloc<float>(count, values.data());
	};
	p->vv = upload(0);
	p->vi = upload(1);
	p->ii = upload(2);
	p->iv = upload(3);
	return true;
}

// the main voltage update of one grid or part
static void UpdateVoltagesOf(GPU_Backend_HIP::Impl* d)
{
	d->ctx->MakeCurrent();
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

// the main current update of one grid or part
static void UpdateCurrentsOf(GPU_Backend_HIP::Impl* d)
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

void GPU_Backend_HIP::UpdateVoltages()
{
	if (m_Parts.size()>1)
		ExchangeHalos(true);
	for (Impl* p : m_Parts)
		UpdateVoltagesOf(p);
	if (m_Parts.size()>1)
		ReleaseHalos(true);
}

void GPU_Backend_HIP::UpdateCurrents()
{
	if (m_Parts.size()>1)
		ExchangeHalos(false);
	for (Impl* p : m_Parts)
		UpdateCurrentsOf(p);
	if (m_Parts.size()>1)
		ReleaseHalos(false);
}

// The halo lines are copied lazily, before the update that reads them: the voltage update of a part
// reads the currents of the line below its first one (the last line of the part below), the current
// update the voltages of the line above its last one (the first line of the part above). Only the y
// and z components are read across x (see update_voltages/update_currents), so only they are copied,
// on the stream of the receiving part once the sending part has computed the line (all extension
// hooks are done by then). This works across devices without peer access and between parts on one
// device.
void GPU_Backend_HIP::ExchangeHalos(bool currents)
{
	for (size_t n=0; n+1<m_Parts.size(); ++n)
	{
		Impl* src = m_Parts[currents ? n : n+1];
		Impl* dst = m_Parts[currents ? n+1 : n];
		const size_t line = (size_t)src->dim.ny*src->dim.nz;
		// currents: the last owned line of the lower part to the halo below the upper part;
		// voltages: the first owned line of the upper part to the halo above the lower part
		const size_t src_x = currents ? src->LocalX(src->x1)-1 : src->LocalX(src->x0);
		const size_t dst_x = currents ? 0 : dst->dim.nx-1;
		const float* sf = currents ? src->curr : src->volt;
		float* df = currents ? dst->curr : dst->volt;
		src->ctx->MakeCurrent();
		HIP_Check(hipEventRecord(src->halo_sent, src->Stream()), "hipEventRecord");
		dst->ctx->MakeCurrent();
		HIP_Check(hipStreamWaitEvent(dst->Stream(), src->halo_sent, 0), "hipStreamWaitEvent");
		const bool timed = HIP_KernelTimes::Enabled();
		if (timed)
			HIP_KernelTimes::Begin(dst->Stream());
		for (int c=1; c<3; ++c)
			HIP_Check(hipMemcpyAsync(df + c*dst->numCells + dst_x*line, sf + c*src->numCells + src_x*line, line*sizeof(float),
			                         hipMemcpyDeviceToDevice, dst->Stream()), "halo copy");
		if (timed)
			HIP_KernelTimes::End(dst->Stream(), currents ? "halo copy (currents)" : "halo copy (voltages)");
		HIP_Check(hipEventRecord(dst->halo_received, dst->Stream()), "hipEventRecord");
	}
}

// The sending part of each halo copy waits for it before it writes the sent line again. Called after
// the main updates are queued: they do not write the field that is sent (the voltage update reads the
// currents, the current update the voltages; no fused step with parts), so they run during the copy.
void GPU_Backend_HIP::ReleaseHalos(bool currents)
{
	for (size_t n=0; n+1<m_Parts.size(); ++n)
	{
		Impl* src = m_Parts[currents ? n : n+1];
		Impl* dst = m_Parts[currents ? n+1 : n];
		src->ctx->MakeCurrent();
		HIP_Check(hipStreamWaitEvent(src->Stream(), dst->halo_received, 0), "hipStreamWaitEvent");
	}
}

static void CheckSize(size_t host_count, size_t device_count)
{
	if (host_count!=device_count)
		throw std::runtime_error("GPU_Backend_HIP: host and device field size mismatch");
}

// The owned lines of every part to the host mirror \a host (\a size values), or the host mirror to
// the owned and halo lines of every part: one copy per component and part
void GPU_Backend_HIP::CopyParts(bool currents, bool upload, float* host, size_t size)
{
	const size_t line = (size_t)d->dim.ny*d->dim.nz;
	const size_t global = (size_t)d->gnx*line;
	CheckSize(size, 3*global);
	d->PinHostMemory(host, size*sizeof(float));
	d->ctx->MakeCurrent();
	for (Impl* p : m_Parts)
	{
		const size_t x = upload ? 0 : p->LocalX(p->x0);   // first local line
		const size_t lines = upload ? p->dim.nx : p->x1 - p->x0;
		float* field = currents ? p->curr : p->volt;
		p->ctx->MakeCurrent();
		for (int c=0; c<3; ++c)
		{
			float* dev = field + c*p->numCells + x*line;
			float* h = host + c*global + (p->x_off + x)*line;
			if (upload)
				HIP_Check(hipMemcpyAsync(dev, h, lines*line*sizeof(float), hipMemcpyHostToDevice, p->Stream()), "upload");
			else
				HIP_Check(hipMemcpyAsync(h, dev, lines*line*sizeof(float), hipMemcpyDeviceToHost, p->Stream()), "download");
		}
	}
	for (Impl* p : m_Parts)
		p->Flush();
}

void GPU_Backend_HIP::DownloadVoltages(ArrayLib::ArrayNIJK<FDTD_FLOAT>& volt)
{
	if (m_Parts.size()>1)
		return CopyParts(false, false, volt.data(), volt.size());
	CheckSize(volt.size(), 3*d->numCells);
	d->PinHostMemory(volt.data(), volt.size()*sizeof(FDTD_FLOAT));
	d->ctx->MakeCurrent();
	HIP_Check(hipMemcpyAsync(volt.data(), d->volt, volt.size()*sizeof(float), hipMemcpyDeviceToHost, d->Stream()), "download");
	d->Flush();
}

void GPU_Backend_HIP::DownloadCurrents(ArrayLib::ArrayNIJK<FDTD_FLOAT>& curr)
{
	if (m_Parts.size()>1)
		return CopyParts(true, false, curr.data(), curr.size());
	CheckSize(curr.size(), 3*d->numCells);
	d->PinHostMemory(curr.data(), curr.size()*sizeof(FDTD_FLOAT));
	d->ctx->MakeCurrent();
	HIP_Check(hipMemcpyAsync(curr.data(), d->curr, curr.size()*sizeof(float), hipMemcpyDeviceToHost, d->Stream()), "download");
	d->Flush();
}

void GPU_Backend_HIP::UploadVoltages(const ArrayLib::ArrayNIJK<FDTD_FLOAT>& volt)
{
	if (m_Parts.size()>1)
		return CopyParts(false, true, const_cast<FDTD_FLOAT*>(volt.data()), volt.size());
	CheckSize(volt.size(), 3*d->numCells);
	d->PinHostMemory(volt.data(), volt.size()*sizeof(FDTD_FLOAT));
	d->ctx->MakeCurrent();
	HIP_Check(hipMemcpyAsync(d->volt, volt.data(), volt.size()*sizeof(float), hipMemcpyHostToDevice, d->Stream()), "upload");
	d->Flush();
}

void GPU_Backend_HIP::UploadCurrents(const ArrayLib::ArrayNIJK<FDTD_FLOAT>& curr)
{
	if (m_Parts.size()>1)
		return CopyParts(true, true, const_cast<FDTD_FLOAT*>(curr.data()), curr.size());
	CheckSize(curr.size(), 3*d->numCells);
	d->PinHostMemory(curr.data(), curr.size()*sizeof(FDTD_FLOAT));
	d->ctx->MakeCurrent();
	HIP_Check(hipMemcpyAsync(d->curr, curr.data(), curr.size()*sizeof(float), hipMemcpyHostToDevice, d->Stream()), "upload");
	d->Flush();
}

bool GPU_Backend_HIP::DownloadRange(bool currents, size_t offset, size_t count, FDTD_FLOAT* dst)
{
	if (m_Parts.size()==1)
	{
		if (offset+count > 3*d->numCells)
			return false;
		const float* src = (currents ? d->curr : d->volt) + offset;
		d->ctx->MakeCurrent();
		HIP_Check(hipMemcpyAsync(dst, src, count*sizeof(float), hipMemcpyDeviceToHost, d->Stream()), "range download");
		d->Flush();
		return true;
	}
	// from the part that owns the range: one component and lines of one part (a z-line always is)
	if (count==0)
		return true;
	const size_t line = (size_t)d->dim.ny*d->dim.nz;
	const size_t global = (size_t)d->gnx*line;
	const size_t last = offset+count-1;
	if ((last>=3*global) || (offset/global!=last/global))
		return false;
	for (Impl* p : m_Parts)
	{
		if (!p->Owns((offset%global)/line) || !p->Owns((last%global)/line))
			continue;
		const size_t local = (offset/global)*p->numCells + (offset%global) - (size_t)p->x_off*line;
		const float* src = (currents ? p->curr : p->volt) + local;
		p->ctx->MakeCurrent();
		HIP_Check(hipMemcpyAsync(dst, src, count*sizeof(float), hipMemcpyDeviceToHost, p->Stream()), "range download");
		p->Flush();
		return true;
	}
	return false;
}

void GPU_Backend_HIP::Synchronize()
{
	for (Impl* p : m_Parts)
		p->Flush();
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
	ctx->MakeCurrent();
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

void GPU_Backend_HIP::Impl::OwnedRegion(const Engine_GatherRegion& region, size_t src, size_t dst, Engine_GatherRegion& owned, std::vector<RegionCopy>& copies) const
{
	for (int n=0; n<3; ++n)
	{
		// the lines are ascending: the owned x lines follow each other from the first one
		const std::vector<unsigned int>& xl = region.lines[n][0];
		owned.lines[n][0].clear();
		size_t first = 0;
		for (size_t a=0; a<xl.size(); ++a)
			if (Owns(xl[a]))
			{
				if (owned.lines[n][0].empty())
					first = a;
				owned.lines[n][0].push_back(LocalX(xl[a]));
			}
		owned.lines[n][1] = region.lines[n][1];
		owned.lines[n][2] = region.lines[n][2];
		const size_t count = owned.Size(n);
		const size_t to = dst + first*region.lines[n][1].size()*region.lines[n][2].size();
		if (count)
		{
			if (!copies.empty() && (copies.back().src+copies.back().count==src) && (copies.back().dst+copies.back().count==to))
				copies.back().count += count;   // continues the last copy
			else
				copies.push_back(RegionCopy{src, to, count});
		}
		src += count;
		dst += region.Size(n);
	}
}

bool GPU_Backend_HIP::SetSnapshotGather(const std::vector<Engine_GatherRegion>& volt_regions, const std::vector<Engine_GatherRegion>& curr_regions)
{
	if (m_Parts.size()>1)
		return SetSnapshotGatherParts(volt_regions, curr_regions);
	d->ctx->MakeCurrent();
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
	d->ctx->MakeCurrent();
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
		(void)hipGetLastError();
		d->FreeSnapshots();
		std::cerr << "GPU_Backend_HIP: not enough memory for field snapshots, the field dumps wait for the device" << std::endl;
		return false;
	}
	d->snap_nv = nv;
	d->snap_n = n;
	return true;
}

// Each part gathers the raw values of the regions on the x lines it owns into its own snap_dev, and
// downloads them to their places in one page-locked host buffer, the snapshot of all parts (see
// Impl::OwnedRegion()). The interpolation of the dumps reads neighbouring lines, but on the host from
// the snapshot, where every line comes from the part that owns and updates it: no halo is read.
bool GPU_Backend_HIP::SetSnapshotGatherParts(const std::vector<Engine_GatherRegion>& volt_regions, const std::vector<Engine_GatherRegion>& curr_regions)
{
	for (Impl* p : m_Parts)
	{
		p->ctx->MakeCurrent();
		HIP_Check(hipStreamSynchronize(p->ctx->copy_stream), "stream synchronization");
		p->FreeSnapshots();
	}
	m_SnapshotParts = false;
	size_t nv = 0, n = 0;
	for (size_t r=0; r<volt_regions.size(); ++r)
		nv += volt_regions[r].Size();
	n = nv;
	for (size_t r=0; r<curr_regions.size(); ++r)
		n += curr_regions[r].Size();
	if ((n==0) || (n>UINT_MAX))
		return false;
	// not fatal: without snapshots the dumps read the fields from the host mirror
	d->ctx->MakeCurrent();
	bool ok = (hipHostMalloc(&d->snap_host, n*sizeof(float), hipHostMallocPortable)==hipSuccess);   // written by all devices
	const std::vector<Engine_GatherRegion>* regions[2] = {&volt_regions, &curr_regions};
	for (Impl* p : m_Parts)
	{
		size_t src = 0, dst = 0;   // the first raw value of the region in snap_dev of the part, in the snapshot
		for (int f=0; f<2 && ok; ++f)
		{
			for (size_t r=0; r<regions[f]->size() && ok; ++r)
			{
				const Engine_GatherRegion& region = regions[f]->at(r);
				Engine_GatherRegion owned;
				p->OwnedRegion(region, src, dst, owned, p->snap_copies);
				dst += region.Size();
				if (owned.Size()==0)
					continue;
				HIP_GatherRegion dev;
				void* buf = p->UploadRegion(owned, dev);
				ok = (buf!=NULL);
				if (ok)
				{
					p->snap_lines.push_back(buf);
					p->snap_regions[f].push_back(dev);
					src += dev.total;
				}
			}
			if (f==0)
				p->snap_nv = src;
		}
		p->snap_n = src;
		p->ctx->MakeCurrent();
		ok = ok && (hipMalloc(&p->snap_dev, std::max(src, (size_t)1)*sizeof(float))==hipSuccess);
	}
	if (!ok)
	{
		(void)hipGetLastError();
		for (Impl* p : m_Parts)
			p->FreeSnapshots();
		std::cerr << "GPU_Backend_HIP: not enough memory for field snapshots, the field dumps wait for the device" << std::endl;
		return false;
	}
	m_SnapshotNV = nv;
	m_SnapshotParts = true;
	return true;
}

bool GPU_Backend_HIP::SnapshotFieldsParts(const FDTD_FLOAT* &volt, const FDTD_FLOAT* &curr)
{
	if (!m_SnapshotParts)
		return false;
	for (Impl* p : m_Parts)
	{
		// the raw values on the work stream, once the last snapshot was downloaded from snap_dev
		p->ctx->MakeCurrent();
		if (p->snap_taken)
			HIP_Check(hipStreamWaitEvent(p->Stream(), p->snap_done, 0), "hipStreamWaitEvent");
		size_t offset = 0;
		for (int f=0; f<2; ++f)
			for (size_t r=0; r<p->snap_regions[f].size(); ++r)
			{
				const HIP_GatherRegion& R = p->snap_regions[f][r];
				HIP_Launch(p, "gather_region", gather_region, R.total, 1, 1, (const float*)(f ? p->curr : p->volt), R, p->dim, p->snap_dev + offset);
				offset += R.total;
			}
		p->ctx->MakeCurrent();
		HIP_Check(hipEventRecord(p->snap_evaluated, p->Stream()), "hipEventRecord");
		// the download on the copy stream of the part, while its work stream continues
		hipStream_t copy = p->ctx->copy_stream;
		HIP_Check(hipStreamWaitEvent(copy, p->snap_evaluated, 0), "hipStreamWaitEvent");
		for (const Impl::RegionCopy& c : p->snap_copies)
			HIP_Check(hipMemcpyAsync(d->snap_host + c.dst, p->snap_dev + c.src, c.count*sizeof(float), hipMemcpyDeviceToHost, copy), "snapshot download");
		HIP_Check(hipEventRecord(p->snap_done, copy), "hipEventRecord");
		p->snap_taken = true;
	}
	volt = d->snap_host;
	curr = d->snap_host + m_SnapshotNV;
	return true;
}

bool GPU_Backend_HIP::SnapshotFields(const FDTD_FLOAT* &volt, const FDTD_FLOAT* &curr)
{
	if (m_Parts.size()>1)
		return SnapshotFieldsParts(volt, curr);
	if (!d->snap_n)
		return false;
	// the raw values on the work stream, once the last snapshot was downloaded from snap_dev
	d->ctx->MakeCurrent();
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
	d->ctx->MakeCurrent();
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
	if (m_Parts.size()>1)
	{
		// each part sums the raw values on the x lines it owns (see Impl::OwnedRegion()), all with the same id
		std::vector<Impl::FieldDFT> dfts(m_Parts.size(), Impl::FieldDFT());
		std::vector<void*> lines(m_Parts.size(), (void*)NULL);
		bool ok = true;
		for (size_t i=0; i<m_Parts.size() && ok; ++i)
		{
			Impl* p = m_Parts[i];
			Impl::FieldDFT& dft = dfts[i];
			dft.count = count;
			dft.currents = h_field;
			dft.all = n;
			Engine_GatherRegion owned;
			p->OwnedRegion(region, 0, 0, owned, dft.copies);
			const size_t pn = owned.Size();
			if (pn==0)
				continue;   // no raw values: no sums and no launches (region.total 0)
			lines[i] = p->UploadRegion(owned, dft.region);
			dft.sums = lines[i] ? p->TryAlloc<float2>(pn*count) : NULL;
			ok = (dft.sums!=NULL);
			if (ok)
			{
				p->ctx->MakeCurrent();
				HIP_Check(hipMemsetAsync(dft.sums, 0, pn*count*sizeof(float2), p->Stream()), "hipMemset");
			}
		}
		if (!ok)
		{
			// not fatal: without device sums the dump sums the fields on the host
			for (size_t i=0; i<m_Parts.size(); ++i)
			{
				m_Parts[i]->Free(dfts[i].sums);
				m_Parts[i]->Free(lines[i]);
			}
			std::cerr << "GPU_Backend_HIP: not enough memory for a frequency domain dump on the device" << std::endl;
			return -1;
		}
		for (size_t i=0; i<m_Parts.size(); ++i)
			m_Parts[i]->field_dfts.push_back(dfts[i]);
		return d->field_dfts.size()-1;
	}
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
	d->ctx->MakeCurrent();
	HIP_Check(hipMemsetAsync(dft.sums, 0, n*count*sizeof(float2), d->Stream()), "hipMemset");
	d->field_dfts.push_back(dft);
	return d->field_dfts.size()-1;
}

void GPU_Backend_HIP::AccumulateFieldDFT(int id, const std::vector<std::complex<float>>& weights)
{
	if ((id<0) || ((size_t)id>=d->field_dfts.size()))
		return;
	const unsigned int count = std::min((size_t)d->field_dfts[id].count, weights.size());
	for (unsigned int k0=0; k0<count; k0+=DFT_WEIGHTS)
	{
		HIP_DFTWeights W;
		const unsigned int nw = std::min(count-k0, (unsigned int)DFT_WEIGHTS);
		for (unsigned int k=0; k<nw; ++k)
			W.w[k] = make_float2(weights[k0+k].real(), weights[k0+k].imag());
		// several parts: each sums its own raw values, a part without any launches nothing
		for (Impl* p : m_Parts)
		{
			const Impl::FieldDFT& dft = p->field_dfts[id];
			HIP_Launch(p, "field_dft", field_dft, dft.region.total, 1, 1, (const float*)(dft.currents ? p->curr : p->volt),
			            dft.region, p->dim, W, nw, dft.sums + k0*dft.region.total);
		}
	}
}

bool GPU_Backend_HIP::ReadFieldDFT(int id, unsigned int k, std::vector<std::complex<float>>& sums)
{
	if ((id<0) || ((size_t)id>=d->field_dfts.size()) || (k>=d->field_dfts[id].count))
		return false;
	if (m_Parts.size()>1)
	{
		// the sums of each part to their places in those of the whole region
		sums.resize(d->field_dfts[id].all);
		for (Impl* p : m_Parts)
		{
			const Impl::FieldDFT& dft = p->field_dfts[id];
			p->ctx->MakeCurrent();
			for (const Impl::RegionCopy& c : dft.copies)
				HIP_Check(hipMemcpyAsync(sums.data() + c.dst, dft.sums + k*dft.region.total + c.src, c.count*sizeof(float2),
				                         hipMemcpyDeviceToHost, p->Stream()), "download");
		}
		for (Impl* p : m_Parts)
			p->Flush();
		return true;
	}
	const Impl::FieldDFT& dft = d->field_dfts[id];
	sums.resize(dft.region.total);
	d->ctx->MakeCurrent();
	HIP_Check(hipMemcpyAsync(sums.data(), dft.sums + k*dft.region.total, sums.size()*sizeof(float2), hipMemcpyDeviceToHost, d->Stream()), "download");
	d->Flush();
	return true;
}

void GPU_Backend_HIP::WaitSnapshot()
{
	// called on the dump thread: report errors instead of throwing
	for (Impl* p : m_Parts)
	{
		const hipError_t err = hipEventSynchronize(p->snap_done);
		if (err!=hipSuccess)
			std::cerr << "GPU_Backend_HIP: field snapshot failed: " << hipGetErrorString(err) << std::endl;
	}
}

bool GPU_Backend_HIP::CalcFastEnergy(const unsigned int numNodes[3], double& E_energy, double& H_energy)
{
	E_energy = H_energy = 0;
	const HIP_GridDim L = {numNodes[0], numNodes[1], numNodes[2]};
	if ((L.nx>d->gnx) || (L.ny>d->dim.ny) || (L.nz>d->dim.nz))
		return false;
	const size_t count = (size_t)L.ny*L.nz;
	if (L.nx==0 || count==0)
		return true;
	// each part sums its owned lines below L.nx, from its first owned line on
	const size_t line = (size_t)d->dim.ny*d->dim.nz;
	std::vector<Impl*> summed;
	for (Impl* p : m_Parts)
	{
		const HIP_GridDim Lp = {std::min(p->x1, L.nx) - std::min(p->x0, std::min(p->x1, L.nx)), L.ny, L.nz};
		if (Lp.nx==0)
			continue;
		if (p->energy_count<count)
		{
			p->energy = p->Alloc<float2>(count);
			p->energy_count = count;
			p->energy_host.resize(count);   // small, not pinned
		}
		const size_t first = (size_t)p->LocalX(p->x0)*line;
		HIP_Launch(p, "energy_lines", energy_lines, L.nz, L.ny, 1, (const float*)p->volt + first, (const float*)p->curr + first, p->energy, p->dim, Lp);
		HIP_Check(hipMemcpyAsync(p->energy_host.data(), p->energy, count*sizeof(float2), hipMemcpyDeviceToHost, p->Stream()), "energy download");
		summed.push_back(p);
	}
	for (Impl* p : summed)
	{
		p->Flush();
		for (size_t k=0; k<count; ++k)
		{
			E_energy += p->energy_host[k].x;
			H_energy += p->energy_host[k].y;
		}
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

namespace
{
//! An extension with a device implementation per part: the hooks of all of them, part by part
class HIP_PartsExtension : public GPU_Extension
{
public:
	HIP_PartsExtension(const std::vector<GPU_Extension*>& parts) : m_Parts(parts) {}
	virtual ~HIP_PartsExtension()
	{
		for (GPU_Extension* e : m_Parts)
			delete e;
	}

	virtual void DoPreVoltageUpdates()  {for (GPU_Extension* e : m_Parts) e->DoPreVoltageUpdates();}
	virtual void DoPostVoltageUpdates() {for (GPU_Extension* e : m_Parts) e->DoPostVoltageUpdates();}
	virtual void Apply2Voltages()       {for (GPU_Extension* e : m_Parts) e->Apply2Voltages();}
	virtual void DoPreCurrentUpdates()  {for (GPU_Extension* e : m_Parts) e->DoPreCurrentUpdates();}
	virtual void DoPostCurrentUpdates() {for (GPU_Extension* e : m_Parts) e->DoPostCurrentUpdates();}
	virtual void Apply2Current()        {for (GPU_Extension* e : m_Parts) e->Apply2Current();}
	virtual void Synchronize()          {for (GPU_Extension* e : m_Parts) e->Synchronize();}

protected:
	std::vector<GPU_Extension*> m_Parts;
};
}

GPU_Extension* GPU_Backend_HIP::CreateExtension(Engine_Extension* eng_ext, Engine* eng)
{
	// OPENEMS_HIP_EXTENSIONS=0: all extensions on the host, to compare with the host fallback of parts
	static const bool device_exts = !getenv("OPENEMS_HIP_EXTENSIONS") || (atoi(getenv("OPENEMS_HIP_EXTENSIONS"))!=0);
	if (!device_exts)
	{
		d->fused_blockers++;
		return NULL;
	}
	if (m_Parts.size()>1)
	{
		// The extensions keep only the nodes of their part (see Impl::Owns()), one per part, but the
		// steady-state detection. A factory returns NULL for a layout it does not split (e.g. a Mur
		// plane and the line it reads on two parts) or another extension: then all extensions run
		// on the host (host fallback).
		if (dynamic_cast<Engine_Ext_SteadyState*>(eng_ext))
			return HIP_CreateExt_SteadyStateParts(m_Parts, eng_ext, eng);   // one check over the probes of all parts
		HIP_ExtensionFactory factory = NULL;
		if (dynamic_cast<Engine_Ext_Excitation*>(eng_ext))
			factory = HIP_CreateExt_Excitation;
		else if (dynamic_cast<Engine_Ext_LorentzMaterial*>(eng_ext))
			factory = HIP_CreateExt_LorentzMaterial;
		else if (dynamic_cast<Engine_Ext_UPML*>(eng_ext))
			factory = HIP_CreateExt_UPML;   // the separate region kernels (see HIP_Ext_UPML::Fused())
		else if (dynamic_cast<Engine_Ext_Mur_ABC*>(eng_ext))
			factory = HIP_CreateExt_Mur_ABC;
		else if (dynamic_cast<Engine_Ext_LumpedRLC*>(eng_ext))
			factory = HIP_CreateExt_LumpedRLC;
		else if (dynamic_cast<Engine_Ext_TFSF*>(eng_ext))
			factory = HIP_CreateExt_TFSF;
		else if (dynamic_cast<Engine_Ext_DebyeMaterial*>(eng_ext))
			factory = HIP_CreateExt_DebyeMaterial;
		else if (dynamic_cast<Engine_Ext_Absorbing_BC*>(eng_ext))
			factory = HIP_CreateExt_Absorbing_BC;
		else
			factory = HIP_CreateExt_Cylinder;   // detected by its operator extension, NULL for any other
		std::vector<GPU_Extension*> exts;
		for (size_t n=0; factory && (n<m_Parts.size()); ++n)
		{
			GPU_Extension* e = factory(m_Parts[n], eng_ext, eng);
			if (!e)
				break;
			exts.push_back(e);
		}
		if (factory && (exts.size()==m_Parts.size()))
			return new HIP_PartsExtension(exts);
		for (GPU_Extension* e : exts)
			delete e;
		return NULL;
	}

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
