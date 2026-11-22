/*
*	Copyright (C) 2026 Thorsten Liebig (Thorsten.Liebig@gmx.de)
*	Copyright (C) 2026 Alexandru Iercosan-Lucaci (alex.lucaci@alacrity.ro)
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
// FDTD/kernels/cell_update.h, one work-item per node, as in hip/hip_cell_model.hip.

#include <algorithm>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

#include "sycl_internal.h"
#include "FDTD/gpu_cell_model.h"
#include "FDTD/engine_plan.h"
#include "tools/global.h"

namespace
{
//! The most x lines of one launch of a main update and the most x chunks of the fused step: the
//! work-groups along the slowest dimension, which devices limit (65535 on some)
const unsigned int MAX_GROUPS_X = 65535;
}

bool GPU_Backend_SYCL::Init(const Operator* op, const GPU_CellModel& m)
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
				const size_t c = std::min(chunk, n-i);
				std::copy(src.index.begin()+i, src.index.begin()+i+c, part.begin());
				SYCL_Check("upload", [&] {d->Queue().memcpy(dev+i, part.data(), c*sizeof(uint16_t));});
				d->Flush();
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
			std::cout << "GPU_Backend_SYCL: cell model, " << (f ? "currents" : "voltages") << ": " << src.combos.size() << " combinations, "
			          << src.numSlots << " state slots of " << src.slotFloats << " floats, " << src.poles.size() << " poles" << std::endl;
	}
	if (g_settings.GetVerboseLevel()>0)
		std::cout << "GPU_Backend_SYCL: cell model on the device, " << bytes/1048576 << " MiB, " << (d->cm_index16 ? 16 : 32) << " bit node index" << std::endl;
	return true;
}

//! The main voltage or current update of the nodes of \a box
template <typename Index>
static void SYCL_CM_Update(GPU_Backend_SYCL::Impl* d, bool voltages, const Engine_Box& box)
{
	const SYCL_GridDim N = d->dim;
	CM_Box all;
	if (!GPU_UpdateBox(box, N.nx, N.ny, N.nz, all))
		return;
	const unsigned int nz = all.z1-all.z0, ny = all.y1-all.y0;
	const CM_Field F = d->cm_field[voltages ? 0 : 1];
	const Index* index = (const Index*)d->cm_index[voltages ? 0 : 1];
	float* out = voltages ? d->volt : d->curr;
	const float* in = voltages ? d->curr : d->volt;
	// a work-group row per x line: in launches of at most MAX_GROUPS_X lines
	for (unsigned int xa=all.x0; xa<all.x1; xa=(unsigned int)std::min<uint64_t>(all.x1, (uint64_t)xa + MAX_GROUPS_X))
	{
		const unsigned int n = (unsigned int)std::min<uint64_t>(all.x1-xa, MAX_GROUPS_X);
		const unsigned int x0 = xa, y0 = all.y0, z0 = all.z0;
		if (voltages)
			SYCL_Launch(d, "cm_voltages", nz, ny, n, [=](size_t i, size_t j, size_t k) {
				cm_update_voltages<Index>(out, in, index, F, N.nx, N.ny, N.nz, x0+k, y0+j, z0+i);
			});
		else
			SYCL_Launch(d, "cm_currents", nz, ny, n, [=](size_t i, size_t j, size_t k) {
				cm_update_currents<Index>(out, in, index, F, N.nx, N.ny, N.nz, x0+k, y0+j, z0+i);
			});
	}
}

void GPU_Backend_SYCL::UpdateVoltages(const Engine_Box& box)
{
	d->cm_index16 ? SYCL_CM_Update<uint16_t>(d, true, box) : SYCL_CM_Update<uint32_t>(d, true, box);
}

void GPU_Backend_SYCL::UpdateCurrents(const Engine_Box& box)
{
	d->cm_index16 ? SYCL_CM_Update<uint16_t>(d, false, box) : SYCL_CM_Update<uint32_t>(d, false, box);
}

// The fused step (GPU_Backend::PrepareFused()), as cm_fused of the HIP backend: the voltages and currents of a timestep
// in one pass, from the fields into the second buffers, swapped after it (the step reads only old values, so the tiles
// need no order). A work-group of TZ x TY items walks through the XC lines of its chunk along x: the new voltages of
// line x+1, then the currents of line x from the new voltages of x (local memory) and x+1 (registers). It owns
// (TZ-1) x (TY-1) nodes; its last item along z and y computes the voltages next to them, and every item those of the
// line after the chunk, which the neighbouring tiles own: the same values, without storing the state
// (cm_apply_component()). The new state of the voltages goes to a second buffer too: the neighbouring tiles read the
// old one. Only the owner computes a current, its state stays in place. The late currents are left for SeamUpdate(),
// after the hooks between the two updates.

namespace
{
//! Whether the current of node \a i waits for SeamUpdate()
inline bool cm_is_late(const uint32_t* late, uint64_t i)
{
	return late && ((late[i>>5]>>(i&31)) & 1u);
}

//! The currents of node \a i at (x, y, z) from \a Iin into \a Iout, as cm_update_currents(), from the new voltages in local
//! memory at \a s (component stride \a sc, the next node along y at +TZ, along z at +1) and of the next node along x, \a vx1 and
//! \a vx2 (the curls read no x component there).
//! The components written out, a constant n each: no array indexed at run time (private memory on a device)
template <typename Index, typename A>
inline void cm_fused_currents(const float* Iin, float* Iout, const Index* index, const CM_Field& F,
                              const SYCL_GridDim& N, unsigned int x, unsigned int y, unsigned int z, A i, A sN,
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
}

bool GPU_Backend_SYCL::PrepareFused(const std::vector<uint64_t>& late)
{
	Impl::Fused& f = d->fused;
	ReleaseFused();
	// as the HIP backend: a wavefront along z, 8 rows, chunks of 8 lines along x
	f.TZ = 32;
	f.TY = 8;
	f.XC = 8;
	if (const char* env = getenv("OPENEMS_SYCL_FUSED"))   // TZ:TY:XC, for measurements
	{
		unsigned int v[3] = {0, 0, 0};
		if ((sscanf(env, "%u:%u:%u", &v[0], &v[1], &v[2])!=3) || (v[0]<2) || (v[1]<2) || (v[2]<1))
			throw std::invalid_argument("GPU_Backend_SYCL: OPENEMS_SYCL_FUSED must be TZ:TY:XC, TZ and TY at least 2");
		f.TZ = v[0];
		f.TY = v[1];
		f.XC = v[2];
	}
	if (f.TZ*f.TY > d->ctx->max_group)
	{
		if (g_settings.GetVerboseLevel()>0)
			std::cout << "GPU_Backend_SYCL: no fused step: work-groups of " << f.TZ*f.TY << " items, the device allows " << d->ctx->max_group << std::endl;
		return false;
	}
	if ((d->dim.nx + f.XC - 1)/f.XC > MAX_GROUPS_X)
	{
		if (g_settings.GetVerboseLevel()>0)
			std::cout << "GPU_Backend_SYCL: no fused step: more chunks of " << f.XC << " x lines than a launch allows" << std::endl;
		return false;
	}
	// the second buffers: only a speed-up, without the memory the main updates stay separate. They
	// leave some memory free (256 MiB or 5 %), for what is allocated later (energy sums, snapshots).
	// An Intel GPU tells its free memory, and is asked: its driver gives out more than the card
	// holds and spills to host memory, which slows the fused step instead of failing the
	// allocation (a B580 spilling 2 % of the fields ran it 24 % slower). A device that does not
	// tell it holds that much while they are allocated
	const size_t need = (2*3*d->numCells + d->cm_state_size[0])*sizeof(float) + (late.empty() ? 0 : d->numCells/8 + 8*late.size());
	const size_t total = d->ctx->device.get_info<sycl::info::device::global_mem_size>();
	const size_t spare = std::max<size_t>((size_t)256 << 20, total/20);
	void* reserve = NULL;
	bool room = true;
	if (d->ctx->device.has(sycl::aspect::ext_intel_free_memory))
	{
		const size_t free_bytes = d->ctx->device.get_info<sycl::ext::intel::info::device::free_memory>();
		if (free_bytes < need + spare)
		{
			if (g_settings.GetVerboseLevel()>0)
				std::cout << "GPU_Backend_SYCL: no fused step: " << need/1048576 << " MiB for the second field buffers, "
				          << free_bytes/1048576 << " MiB free" << std::endl;
			return false;
		}
	}
	else
	{
		reserve = d->TryAlloc<char>(std::min<size_t>(spare, d->ctx->max_alloc));
		room = (reserve!=NULL);
	}
	f.volt = room ? d->TryAlloc<float>(3*d->numCells) : NULL;
	f.curr = f.volt ? d->TryAlloc<float>(3*d->numCells) : NULL;
	f.vstate = f.curr ? d->TryAlloc<float>(std::max<size_t>(d->cm_state_size[0], 1)) : NULL;
	d->Free(reserve);
	if (f.vstate==NULL)
	{
		ReleaseFused();
		if (g_settings.GetVerboseLevel()>0)
			std::cout << "GPU_Backend_SYCL: no fused step: not enough device memory for the second field buffers" << std::endl;
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
		std::cout << "GPU_Backend_SYCL: fused step, tiles of " << f.TZ-1 << "x" << f.TY-1 << " nodes (z, y), chunks of " << f.XC << " lines along x, "
		          << (2*3*d->numCells + d->cm_state_size[0])*sizeof(float)/1048576 << " MiB for the second field buffers" << std::endl;
	return true;
}

void GPU_Backend_SYCL::ReleaseFused()
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
static void SYCL_LaunchFused(GPU_Backend_SYCL::Impl* d)
{
	const GPU_Backend_SYCL::Impl::Fused& f = d->fused;
	const SYCL_GridDim N = d->dim;
	const unsigned int TZ = f.TZ, TY = f.TY, XC = f.XC;
	const sycl::range<3> groups(((size_t)N.nx + XC - 1)/XC, ((size_t)N.ny + TY - 2)/(TY - 1), ((size_t)N.nz + TZ - 2)/(TZ - 1));
	const float* Vin = d->volt;
	const float* Iin = d->curr;
	float* Vout = f.volt;
	float* Iout = f.curr;
	const Index* iv = (const Index*)d->cm_index[0];
	const Index* ii = (const Index*)d->cm_index[1];
	const CM_Field FV = d->cm_field[0];
	const CM_Field FI = d->cm_field[1];
	const uint32_t* late = f.late;
	float* vstate = f.vstate;
	SYCL_Check("cm_fused", [&] {
		d->Queue().submit([&](sycl::handler& h) {
			sycl::local_accessor<float, 1> sV(sycl::range<1>(3*TY*TZ), h);   // [3][TY][TZ]: the new voltages of line x
			h.parallel_for(sycl::nd_range<3>(groups*sycl::range<3>(1, TY, TZ), sycl::range<3>(1, TY, TZ)), [=](sycl::nd_item<3> it) {
				const unsigned int tz = it.get_local_id(2), ty = it.get_local_id(1);
				const unsigned int z = it.get_group(2)*(TZ-1) + tz, y = it.get_group(1)*(TY-1) + ty;
				const unsigned int xs = it.get_group(0)*XC, xe = sycl::min(xs+XC, N.nx);
				const bool in_grid = (z<N.nz) && (y<N.ny);
				const bool owner = in_grid && (tz+1<TZ) && (ty+1<TY);
				const A sy = N.nz, sx = (A)N.ny*N.nz, sN = (A)N.nx*sx;
				const unsigned int sc = TY*TZ;
				float* s = sV.template get_multi_ptr<sycl::access::decorated::no>().get() + ty*TZ + tz;
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
				sycl::group_barrier(it.get_group());
				for (unsigned int x=xs; x<xe; ++x)
				{
					volt(x+1, v0, v1, v2);
					if (owner)
					{
						const A i = x*sx + y*sy + z;
						if (!cm_is_late(late, i))
							cm_fused_currents<Index, A>(Iin, Iout, ii, FI, N, x, y, z, i, sN, s, sc, TZ, v1, v2);
					}
					sycl::group_barrier(it.get_group());
					s[0] = v0;
					s[sc] = v1;
					s[2*sc] = v2;
					sycl::group_barrier(it.get_group());
				}
			});
		});
	});
}

void GPU_Backend_SYCL::FusedUpdate()
{
	// the index arithmetic in 32 bit where all values of a field have such an index (fewer instructions)
	const bool a32 = 3*d->numCells <= UINT32_MAX;
	if (d->cm_index16)
		a32 ? SYCL_LaunchFused<uint16_t, uint32_t>(d) : SYCL_LaunchFused<uint16_t, uint64_t>(d);
	else
		a32 ? SYCL_LaunchFused<uint32_t, uint32_t>(d) : SYCL_LaunchFused<uint32_t, uint64_t>(d);
	// the new voltages for the hooks; the old currents stay until SeamUpdate(), the hooks between read them
	std::swap(d->volt, d->fused.volt);
	std::swap(d->cm_field[0].state, d->fused.vstate);
}

//! The late currents, after the hooks: in place in the old currents, from the new voltages, then copied to the new ones
template <typename Index>
static void SYCL_LateCurrents(GPU_Backend_SYCL::Impl* d)
{
	const GPU_Backend_SYCL::Impl::Fused& f = d->fused;
	const SYCL_GridDim N = d->dim;
	float* I = d->curr;
	float* Iout = f.curr;
	const float* V = d->volt;
	const Index* index = (const Index*)d->cm_index[1];
	const CM_Field F = d->cm_field[1];
	const uint64_t* list = f.list;
	SYCL_Launch(d, "cm_late_currents", f.count, 1, 1, [=](size_t n, size_t, size_t) {
		const uint64_t i = list[n], sN = (uint64_t)N.nx*N.ny*N.nz;
		const unsigned int z = i%N.nz, y = (i/N.nz)%N.ny, x = i/((uint64_t)N.nz*N.ny);
		cm_update_currents<Index>(I, V, index, F, N.nx, N.ny, N.nz, x, y, z);
		Iout[i] = I[i];
		Iout[sN+i] = I[sN+i];
		Iout[2*sN+i] = I[2*sN+i];
	});
}

void GPU_Backend_SYCL::SeamUpdate()
{
	Impl::Fused& f = d->fused;
	if (f.count)
		d->cm_index16 ? SYCL_LateCurrents<uint16_t>(d) : SYCL_LateCurrents<uint32_t>(d);
	std::swap(d->curr, f.curr);
}
