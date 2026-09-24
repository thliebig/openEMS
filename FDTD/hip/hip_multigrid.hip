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
#include "hip_internal.h"

// Coupling of a cylindrical multi-grid level (base) and its sub-grid, see
// Engine_CylinderMultiGrid and Engine_Ext_CylinderMultiGrid. Both grids share the
// HIP stream, so the kernels run in order with the updates of both grids.

// base grid voltages at r line r to the sub-grid (every other alpha line)
__global__ void mg_sync_volt(float* sub, const float* base, HIP_GridDim NB, HIP_GridDim NS, unsigned int r)
{
	const unsigned int z = blockIdx.x*blockDim.x + threadIdx.x;
	const unsigned int h = blockIdx.y*blockDim.y + threadIdx.y;
	const unsigned int a = 2*h;
	if (z>=NB.nz || a+1>=NB.ny)
		return;
	sub[nijk(NS, 0, r, h, z)] = 0;
	sub[nijk(NS, 2, r, h, z)] = base[nijk(NB, 2, r, a, z)];
	sub[nijk(NS, 1, r, h, z)] = base[nijk(NB, 1, r, a, z)] + base[nijk(NB, 1, r, a+1, z)];
}

//! interpolation positions and weights, [n][a] with n=0 for the r and z components, n=1 for the alpha component
struct MGWeights
{
	const unsigned int* pos_2p; const unsigned int* pos_2pp;
	const float* w_2p; const float* w_2pp;
};

// sub-grid fields to the base grid at the r lines r0 ... r0+count-1
__global__ void mg_interpolate(float* base, const float* sub, MGWeights W, HIP_GridDim NB, HIP_GridDim NS, unsigned int r0, unsigned int count)
{
	const unsigned int z = blockIdx.x*blockDim.x + threadIdx.x;
	const unsigned int a = blockIdx.y*blockDim.y + threadIdx.y;
	const unsigned int dr = blockIdx.z;
	if (z>=NB.nz || a>=NB.ny || dr>=count)
		return;
	const unsigned int r = r0 + dr;
	const unsigned int k1 = NB.ny + a;
	base[nijk(NB, 0, r, a, z)] = W.w_2p[a]*sub[nijk(NS, 0, r, W.pos_2p[a], z)] + W.w_2pp[a]*sub[nijk(NS, 0, r, W.pos_2pp[a], z)];
	base[nijk(NB, 2, r, a, z)] = W.w_2p[a]*sub[nijk(NS, 2, r, W.pos_2p[a], z)] + W.w_2pp[a]*sub[nijk(NS, 2, r, W.pos_2pp[a], z)];
	base[nijk(NB, 1, r, a, z)] = W.w_2p[k1]*sub[nijk(NS, 1, r, W.pos_2p[k1], z)] + W.w_2pp[k1]*sub[nijk(NS, 1, r, W.pos_2pp[k1], z)];
}

class HIP_MultiGridLink : public GPU_MultiGridLink
{
public:
	HIP_MultiGridLink(GPU_Backend_HIP::Impl* base, GPU_Backend_HIP::Impl* sub, const GPU_MultiGridInterpolation& interpol);

	virtual void SyncVoltages();
	virtual void SyncCurrents();
	virtual void InterpolateToBase();

protected:
	MGWeights Upload(const std::vector<unsigned int> pos_2p[2], const std::vector<unsigned int> pos_2pp[2],
	                 const std::vector<float> w_2p[2], const std::vector<float> w_2pp[2]);
	void Interpolate(float* base_field, const float* sub_field, const MGWeights& w, unsigned int r0, unsigned int count);

	GPU_Backend_HIP::Impl* b;
	GPU_Backend_HIP::Impl* s;
	unsigned int m_SplitPos;
	MGWeights m_Volt, m_Curr;
};

HIP_MultiGridLink::HIP_MultiGridLink(GPU_Backend_HIP::Impl* base, GPU_Backend_HIP::Impl* sub, const GPU_MultiGridInterpolation& interpol)
{
	b = base;
	s = sub;
	m_SplitPos = interpol.split_pos;
	m_Volt = Upload(interpol.pos_v_2p, interpol.pos_v_2pp, interpol.w_v_2p, interpol.w_v_2pp);
	m_Curr = Upload(interpol.pos_i_2p, interpol.pos_i_2pp, interpol.w_i_2p, interpol.w_i_2pp);
}

MGWeights HIP_MultiGridLink::Upload(const std::vector<unsigned int> pos_2p[2], const std::vector<unsigned int> pos_2pp[2],
                                     const std::vector<float> w_2p[2], const std::vector<float> w_2pp[2])
{
	// [n][a] with n=0,1 and a < number of base alpha lines
	std::vector<unsigned int> p2p, p2pp;
	std::vector<float> f2p, f2pp;
	for (int n=0; n<2; ++n)
	{
		p2p.insert(p2p.end(), pos_2p[n].begin(), pos_2p[n].end());
		p2pp.insert(p2pp.end(), pos_2pp[n].begin(), pos_2pp[n].end());
		f2p.insert(f2p.end(), w_2p[n].begin(), w_2p[n].end());
		f2pp.insert(f2pp.end(), w_2pp[n].begin(), w_2pp[n].end());
	}
	// the buffers belong to the base grid, which is deleted after the link
	MGWeights w;
	w.pos_2p = b->Alloc<unsigned int>(p2p.size(), p2p.data());
	w.pos_2pp = b->Alloc<unsigned int>(p2pp.size(), p2pp.data());
	w.w_2p = b->Alloc<float>(f2p.size(), f2p.data());
	w.w_2pp = b->Alloc<float>(f2pp.size(), f2pp.data());
	return w;
}

void HIP_MultiGridLink::SyncVoltages()
{
	if (m_SplitPos<1)
		return;
	HIP_Launch(b, "mg_sync_volt", mg_sync_volt, b->dim.nz, b->dim.ny/2, 1, s->volt, (const float*)b->volt, b->dim, s->dim, m_SplitPos-1);
}

void HIP_MultiGridLink::Interpolate(float* base_field, const float* sub_field, const MGWeights& w, unsigned int r0, unsigned int count)
{
	HIP_Launch(b, "mg_interpolate", mg_interpolate, b->dim.nz, b->dim.ny, count, base_field, sub_field, w, b->dim, s->dim, r0, count);
}

void HIP_MultiGridLink::SyncCurrents()
{
	if (m_SplitPos>=2)
		Interpolate(b->curr, s->curr, m_Curr, m_SplitPos-2, 1);
}

void HIP_MultiGridLink::InterpolateToBase()
{
	if (m_SplitPos>=1)
		Interpolate(b->volt, s->volt, m_Volt, 0, m_SplitPos-1);
	if (m_SplitPos>=2)
		Interpolate(b->curr, s->curr, m_Curr, 0, m_SplitPos-2);
}

GPU_MultiGridLink* GPU_Backend_HIP::CreateMultiGridLink(GPU_Backend* sub_grid, const GPU_MultiGridInterpolation& interpol)
{
	GPU_Backend_HIP* sub = dynamic_cast<GPU_Backend_HIP*>(sub_grid);
	// only for a sub-grid on the same context: the coupling kernels are ordered by the shared stream
	if (!sub || sub->d->ctx!=d->ctx)
		return NULL;
	return new HIP_MultiGridLink(d, sub->d, interpol);
}
