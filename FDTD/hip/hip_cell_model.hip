/*
*	Copyright (C) 2026 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

// The main updates on the cell model (GPU_Backend::Init()): the per-node body of
// FDTD/kernels/cell_update.h, one thread per node.

#include "hip/hip_runtime.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

#include "hip_internal.h"
#include "FDTD/gpu_cell_model.h"
#include "tools/global.h"
#include "FDTD/engine_plan.h"
#include <algorithm>

template <typename Index>
__global__ void cm_voltages(float* __restrict__ V, const float* __restrict__ I, const Index* __restrict__ index, CM_Field F, HIP_GridDim N, unsigned int x0)
{
	const unsigned int z = blockIdx.x*blockDim.x + threadIdx.x;
	const unsigned int y = blockIdx.y*blockDim.y + threadIdx.y;
	const unsigned int x = x0 + blockIdx.z;
	if ((z>=N.nz) || (y>=N.ny))
		return;
	cm_update_voltages<Index>(V, I, index, F, N.nx, N.ny, N.nz, x, y, z);
}

template <typename Index>
__global__ void cm_currents(float* __restrict__ I, const float* __restrict__ V, const Index* __restrict__ index, CM_Field F, HIP_GridDim N, unsigned int x0)
{
	const unsigned int z = blockIdx.x*blockDim.x + threadIdx.x;
	const unsigned int y = blockIdx.y*blockDim.y + threadIdx.y;
	const unsigned int x = x0 + blockIdx.z;
	if ((z>=N.nz) || (y>=N.ny))
		return;
	cm_update_currents<Index>(I, V, index, F, N.nx, N.ny, N.nz, x, y, z);
}

// The fused step (GPU_Backend::PrepareFused()): the voltages and currents of a timestep in one pass, from the fields
// into the second buffers, swapped after it (the step reads only old values, so the tiles need no order). A block of
// TZ x TY threads walks through the XC lines of its chunk along x: the new voltages of line x+1, then the currents of
// line x from the new voltages of x (shared memory) and x+1 (registers). It owns (TZ-1) x (TY-1) nodes; its last
// thread along z and y computes the voltages next to them, and every thread those of the line after the chunk, which
// the neighbouring tiles own: the same values, without storing the state (cm_apply_component()). The new state of the
// voltages goes to a second buffer too: the neighbouring tiles read the old one. Only the owner computes a current, its
// state stays in place. The late currents are left for SeamUpdate(), after the hooks between the two updates.

//! Whether the current of node \a i waits for SeamUpdate()
__device__ __forceinline__ bool cm_is_late(const uint32_t* late, uint64_t i)
{
	return late && ((late[i>>5]>>(i&31)) & 1u);
}

//! The currents of node \a i at (x, y, z) from \a Iin into \a Iout, as cm_update_currents(), from the new voltages in shared
//! memory at \a s (component stride \a sc, the next node along y at +TZ, along z at +1) and of the next node along x, \a vx1 and
//! \a vx2 (the curls read no x component there).
//! The components written out, a constant n each: no array indexed at run time (local memory on a device)
template <typename Index, typename A>
__device__ __forceinline__ void cm_fused_currents(const float* __restrict__ Iin, float* __restrict__ Iout, const Index* __restrict__ index, const CM_Field& F,
                                                  const HIP_GridDim& N, unsigned int x, unsigned int y, unsigned int z, A i, A sN,
                                                  const float* s, unsigned int sc, unsigned int TZ, float vx1, float vx2)
{
	const float o0 = s[0], o1 = s[sc], o2 = s[2*sc];
	const float c0 = o2 - s[2*sc + TZ] - o1 + s[sc + 1];
	const float c1 = o0 - s[1] - o2 + vx2;
	const float c2 = o1 - vx1 - o0 + s[TZ];
	const CM_Combo& cb = F.combos[index[i]];
	const A line = (A)x*N.ny + y;
	if (cb.state & CM_LOSSLESS3)
	{
		const float a0 = CM_FACTOR(cb, F, 0, line), a1 = CM_FACTOR(cb, F, 1, line), a2 = CM_FACTOR(cb, F, 2, line);
		Iout[i]      = Iin[i]      + a0*F.gz[0][z]*c0;
		Iout[sN+i]   = Iin[sN+i]   + a1*F.gz[1][z]*c1;
		Iout[2*sN+i] = Iin[2*sN+i] + a2*F.gz[2][z]*c2;
		return;
	}
	const uint32_t slot = (cb.state & CM_STATE) ? F.slot[i] : 0;
	Iout[i]      = cm_component_value(cb, F, 0, line, z, Iin[i], c0, slot, true, F.state);
	Iout[sN+i]   = cm_component_value(cb, F, 1, line, z, Iin[sN+i], c1, slot, true, F.state);
	Iout[2*sN+i] = cm_component_value(cb, F, 2, line, z, Iin[2*sN+i], c2, slot, true, F.state);
}

template <typename Index, typename A>
__global__ void __launch_bounds__(512) cm_fused(const float* __restrict__ Vin, const float* __restrict__ Iin, float* __restrict__ Vout, float* __restrict__ Iout,
                                                const Index* __restrict__ iv, const Index* __restrict__ ii, CM_Field FV, CM_Field FI, HIP_GridDim N,
                                                unsigned int XC, const uint32_t* __restrict__ late, float* __restrict__ vstate)
{
	extern __shared__ float sV[];   // [3][TY][TZ]: the new voltages of line x
	const unsigned int TZ = blockDim.x, TY = blockDim.y, tz = threadIdx.x, ty = threadIdx.y;
	const unsigned int z = blockIdx.x*(TZ-1) + tz, y = blockIdx.y*(TY-1) + ty;
	const unsigned int xs = blockIdx.z*XC, xe = min(xs+XC, N.nx);
	const bool in_grid = (z<N.nz) && (y<N.ny);
	const bool owner = in_grid && (tz+1<TZ) && (ty+1<TY);
	const A sy = N.nz, sx = (A)N.ny*N.nz, sN = (A)N.nx*sx;
	const unsigned int sc = TY*TZ;
	float* s = sV + ty*TZ + tz;
	// the new voltages of line x, stored if owned
	auto volt = [&](unsigned int x, float& v0, float& v1, float& v2) {
		v0 = v1 = v2 = 0;
		if (in_grid && (x<N.nx))
		{
			const bool own = owner && (x<xe);
			cm_voltage_values<Index, A>(Vin, Iin, iv, FV, N.nx, N.ny, N.nz, x, y, z, v0, v1, v2, own, vstate);
			if (own)
			{
				const A i = x*sx + y*sy + z;
				Vout[i] = v0;
				Vout[sN+i] = v1;
				Vout[2*sN+i] = v2;
			}
		}
	};
	float v0, v1, v2;
	volt(xs, v0, v1, v2);
	s[0] = v0;
	s[sc] = v1;
	s[2*sc] = v2;
	__syncthreads();
	for (unsigned int x=xs; x<xe; ++x)
	{
		volt(x+1, v0, v1, v2);
		if (owner)
		{
			const A i = x*sx + y*sy + z;
			if (!cm_is_late(late, i))
				cm_fused_currents<Index, A>(Iin, Iout, ii, FI, N, x, y, z, i, sN, s, sc, TZ, v1, v2);
		}
		__syncthreads();
		s[0] = v0;
		s[sc] = v1;
		s[2*sc] = v2;
		__syncthreads();
	}
}

//! The late currents, after the hooks: in place in the old currents \a I, from the new voltages, then copied to \a Iout
template <typename Index>
__global__ void cm_late_currents(float* __restrict__ I, float* __restrict__ Iout, const float* __restrict__ V, const Index* __restrict__ index, CM_Field F,
                                 HIP_GridDim N, const uint64_t* __restrict__ list, uint64_t count)
{
	const uint64_t n = (uint64_t)blockIdx.x*blockDim.x + threadIdx.x;
	if (n>=count)
		return;
	const uint64_t i = list[n], sN = (uint64_t)N.nx*N.ny*N.nz;
	const unsigned int z = i%N.nz, y = (i/N.nz)%N.ny, x = i/((uint64_t)N.nz*N.ny);
	cm_update_currents<Index>(I, V, index, F, N.nx, N.ny, N.nz, x, y, z);
	Iout[i] = I[i];
	Iout[sN+i] = I[sN+i];
	Iout[2*sN+i] = I[2*sN+i];
}

bool GPU_Backend_HIP::Init(const Operator* op, const GPU_CellModel& m)
{
	if (!InitFields(op))
		return false;
	if ((m.numLines[0]!=d->dim.nx) || (m.numLines[1]!=d->dim.ny) || (m.numLines[2]!=d->dim.nz))
		return false;
	d->cm_index16 = m.Index16();
	size_t bytes = 0;
	for (int f=0; f<2; ++f)
	{
		const GPU_CellModel::Field& src = m.field[f];
		CM_Field& F = d->cm_field[f];
		if (d->cm_index16)
		{
			// converted in chunks, no host copy of the whole index
			const size_t n = src.index.size();
			uint16_t* dev = d->Alloc<uint16_t>(n);
			const size_t chunk = (size_t)1<<22;
			std::vector<uint16_t> part(std::min(n, chunk));
			for (size_t i=0; i<n; i+=chunk)
			{
				const size_t m = std::min(chunk, n-i);
				std::copy(src.index.begin()+i, src.index.begin()+i+m, part.begin());
				HIP_Check(hipMemcpy(dev+i, part.data(), m*sizeof(uint16_t), hipMemcpyHostToDevice), "hipMemcpy");
			}
			d->cm_index[f] = dev;
			bytes += n*2;
		}
		else
		{
			d->cm_index[f] = d->Alloc<uint32_t>(src.index.size(), src.index.data());
			bytes += src.index.size()*4;
		}
		F.combos = d->Alloc<CM_Combo>(src.combos.size(), src.combos.data());
		F.slot = d->Alloc<uint32_t>(src.slot.size(), src.slot.data());
		F.numSlots = src.numSlots;
		F.state = d->Alloc<float>(src.numSlots*src.slotFloats);   // zero
		d->cm_state_size[f] = src.numSlots*src.slotFloats;
		F.poles = d->Alloc<CM_Pole>(src.poles.size(), src.poles.data());
		for (int q=0; q<5; ++q)
			F.profile[q] = d->Alloc<float>(src.profile[q].size(), src.profile[q].data());
		for (int n=0; n<3; ++n)
		{
			F.gxy[n] = d->Alloc<float>(src.gxy[n].size(), src.gxy[n].data());
			F.gz[n] = d->Alloc<float>(src.gz[n].size(), src.gz[n].data());
		}
		bytes += src.combos.size()*sizeof(CM_Combo) + src.slot.size()*4 + src.numSlots*src.slotFloats*4;
		if (g_settings.GetVerboseLevel()>0)
			std::cout << "GPU_Backend_HIP: cell model, " << (f ? "currents" : "voltages") << ": " << src.combos.size() << " combinations, "
			          << src.numSlots << " state slots of " << src.slotFloats << " floats, " << src.poles.size() << " poles" << std::endl;
	}
	if (g_settings.GetVerboseLevel()>0)
		std::cout << "GPU_Backend_HIP: cell model on the device, " << bytes/1048576 << " MiB, " << (d->cm_index16 ? 16 : 32) << " bit node index" << std::endl;
	return true;
}

//! The main voltage or current update of the x lines [x0, x1)
static void HIP_CM_Update(GPU_Backend_HIP::Impl* d, bool voltages, unsigned int x0, unsigned int x1)
{
	const HIP_GridDim N = d->dim;
	x1 = std::min(x1, N.nx);
	// a block per x line: in launches of at most max_grid_z lines
	for (unsigned int xa=x0; xa<x1; xa=(unsigned int)std::min<uint64_t>(x1, (uint64_t)xa + d->max_grid_z))
	{
		const unsigned int n = (unsigned int)std::min<uint64_t>(x1-xa, d->max_grid_z);
		if (voltages)
		{
			if (d->cm_index16)
				HIP_Launch(d, "cm_voltages", cm_voltages<uint16_t>, N.nz, N.ny, n, d->volt, (const float*)d->curr, (const uint16_t*)d->cm_index[0], d->cm_field[0], N, xa);
			else
				HIP_Launch(d, "cm_voltages", cm_voltages<uint32_t>, N.nz, N.ny, n, d->volt, (const float*)d->curr, (const uint32_t*)d->cm_index[0], d->cm_field[0], N, xa);
		}
		else
		{
			if (d->cm_index16)
				HIP_Launch(d, "cm_currents", cm_currents<uint16_t>, N.nz, N.ny, n, d->curr, (const float*)d->volt, (const uint16_t*)d->cm_index[1], d->cm_field[1], N, xa);
			else
				HIP_Launch(d, "cm_currents", cm_currents<uint32_t>, N.nz, N.ny, n, d->curr, (const float*)d->volt, (const uint32_t*)d->cm_index[1], d->cm_field[1], N, xa);
		}
	}
}

void GPU_Backend_HIP::UpdateVoltages(const Engine_Box& box)
{
	HIP_CM_Update(d, true, box.start[0], box.stop[0]);
}

void GPU_Backend_HIP::UpdateCurrents(const Engine_Box& box)
{
	HIP_CM_Update(d, false, box.start[0], box.stop[0]);
}

bool GPU_Backend_HIP::PrepareFused(const std::vector<uint64_t>& late)
{
	Impl::Fused& f = d->fused;
	ReleaseFused();
	// as the fused step before the cell model had it: a wavefront along z, 8 rows, short chunks (8 lines: +3 %
	// over 4 on an RTX 3060)
	f.TZ = 32;
	f.TY = 8;
	f.XC = 8;
	if (const char* env = getenv("OPENEMS_HIP_FUSED"))   // TZ:TY:XC, for measurements
	{
		unsigned int v[3] = {0, 0, 0};
		if ((sscanf(env, "%u:%u:%u", &v[0], &v[1], &v[2])!=3) || (v[0]<2) || (v[1]<2) || (v[0]*v[1]>512) || (v[2]<1))
			throw std::invalid_argument("GPU_Backend_HIP: OPENEMS_HIP_FUSED must be TZ:TY:XC, TZ and TY at least 2, TZ*TY at most 512");
		f.TZ = v[0];
		f.TY = v[1];
		f.XC = v[2];
	}
	// one launch: a block per chunk of XC lines along x
	if ((d->dim.nx + f.XC - 1)/f.XC > d->max_grid_z)
	{
		if (g_settings.GetVerboseLevel()>0)
			std::cout << "GPU_Backend_HIP: no fused step: more chunks of " << f.XC << " x lines than a launch allows" << std::endl;
		return false;
	}
	// the second buffers: only a speed-up, without the memory the main updates stay separate. They
	// leave some memory free (256 MiB or 5 %), for what is allocated later (energy sums, snapshots)
	const size_t need = (2*3*d->numCells + d->cm_state_size[0])*sizeof(float) + (late.empty() ? 0 : d->numCells/8 + 8*late.size());
	size_t free_bytes = 0, total_bytes = 0;
	if ((hipMemGetInfo(&free_bytes, &total_bytes)==hipSuccess)
	    && (free_bytes < need + std::max<size_t>((size_t)256 << 20, total_bytes/20)))
	{
		if (g_settings.GetVerboseLevel()>0)
			std::cout << "GPU_Backend_HIP: no fused step: " << need/1048576 << " MiB for the second field buffers, "
			          << free_bytes/1048576 << " MiB free" << std::endl;
		return false;
	}
	(void)hipGetLastError();
	f.volt = d->TryAlloc<float>(3*d->numCells);
	f.curr = f.volt ? d->TryAlloc<float>(3*d->numCells) : NULL;
	f.vstate = f.curr ? d->TryAlloc<float>(std::max<size_t>(d->cm_state_size[0], 1)) : NULL;
	if (f.vstate==NULL)
	{
		ReleaseFused();
		if (g_settings.GetVerboseLevel()>0)
			std::cout << "GPU_Backend_HIP: no fused step: not enough device memory for the second field buffers" << std::endl;
		return false;
	}
	if (!late.empty())
	{
		std::vector<uint32_t> bits((d->numCells + 31)/32, 0);
		for (uint64_t n : late)
			bits[n>>5] |= 1u<<(n&31);
		f.late = d->Alloc<uint32_t>(bits.size(), bits.data());
		f.list = d->Alloc<uint64_t>(late.size(), late.data());
		f.count = late.size();
	}
	if (g_settings.GetVerboseLevel()>0)
		std::cout << "GPU_Backend_HIP: fused step, tiles of " << f.TZ-1 << "x" << f.TY-1 << " nodes (z, y), chunks of " << f.XC << " lines along x, "
		          << (2*3*d->numCells + d->cm_state_size[0])*sizeof(float)/1048576 << " MiB for the second field buffers" << std::endl;
	return true;
}

void GPU_Backend_HIP::ReleaseFused()
{
	Impl::Fused& f = d->fused;
	d->Free(f.volt);
	d->Free(f.curr);
	d->Free(f.vstate);
	d->Free(f.late);
	d->Free(f.list);
	f.volt = f.curr = f.vstate = NULL;
	f.late = NULL;
	f.list = NULL;
	f.count = 0;
}

template <typename Index, typename A>
static void HIP_LaunchFused(GPU_Backend_HIP::Impl* d)
{
	const GPU_Backend_HIP::Impl::Fused& f = d->fused;
	const HIP_GridDim N = d->dim;
	const dim3 block(f.TZ, f.TY, 1);
	const dim3 grid((N.nz + f.TZ - 2)/(f.TZ - 1), (N.ny + f.TY - 2)/(f.TY - 1), (N.nx + f.XC - 1)/f.XC);
	cm_fused<Index, A><<<grid, block, 3*f.TY*f.TZ*sizeof(float), d->Stream()>>>(d->volt, d->curr, f.volt, f.curr, (const Index*)d->cm_index[0],
	                                                                      (const Index*)d->cm_index[1], d->cm_field[0], d->cm_field[1], N, f.XC, f.late, f.vstate);
}

void GPU_Backend_HIP::FusedUpdate()
{
	const bool timed = HIP_KernelTimes::Enabled();
	if (timed)
		HIP_KernelTimes::Begin(d->Stream());
	// the index arithmetic in 32 bit where all values of a field have such an index (fewer instructions)
	const bool a32 = 3*d->numCells <= UINT32_MAX;
	if (d->cm_index16)
		a32 ? HIP_LaunchFused<uint16_t, uint32_t>(d) : HIP_LaunchFused<uint16_t, uint64_t>(d);
	else
		a32 ? HIP_LaunchFused<uint32_t, uint32_t>(d) : HIP_LaunchFused<uint32_t, uint64_t>(d);
	d->CheckLaunch("cm_fused");
	if (timed)
		HIP_KernelTimes::End(d->Stream(), "cm_fused");
	// the new voltages for the hooks; the old currents stay until SeamUpdate(), the hooks between read them
	std::swap(d->volt, d->fused.volt);
	std::swap(d->cm_field[0].state, d->fused.vstate);
}

void GPU_Backend_HIP::SeamUpdate()
{
	Impl::Fused& f = d->fused;
	if (f.count)
	{
		if (d->cm_index16)
			HIP_Launch(d, "cm_late_currents", cm_late_currents<uint16_t>, f.count, 1, 1, d->curr, f.curr, (const float*)d->volt,
			           (const uint16_t*)d->cm_index[1], d->cm_field[1], d->dim, (const uint64_t*)f.list, (uint64_t)f.count);
		else
			HIP_Launch(d, "cm_late_currents", cm_late_currents<uint32_t>, f.count, 1, 1, d->curr, f.curr, (const float*)d->volt,
			           (const uint32_t*)d->cm_index[1], d->cm_field[1], d->dim, (const uint64_t*)f.list, (uint64_t)f.count);
	}
	std::swap(d->curr, f.curr);
}
