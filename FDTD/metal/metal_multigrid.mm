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

#include "metal_internal.h"

// Coupling of a cylindrical multi-grid level (base) and its sub-grid, see
// Engine_CylinderMultiGrid and Engine_Ext_CylinderMultiGrid. Both grids share the
// Metal context, so the kernels run in order with the updates of both grids.
static const char* MULTIGRID_SOURCE = R"MSL(
// base grid voltages at r line r to the sub-grid (every other alpha line); gid = (z, a/2)
kernel void mg_sync_volt(device float* sub           [[buffer(0)]],
                         const device float* base    [[buffer(1)]],
                         constant GridDim& NB        [[buffer(4)]],
                         constant GridDim& NS        [[buffer(5)]],
                         constant uint& r            [[buffer(6)]],
                         uint2 gid [[thread_position_in_grid]])
{
	const uint z = gid.x, h = gid.y, a = 2*h;
	if (a+1 >= NB.ny)
		return;
	sub[nijk(NS, 0, r, h, z)] = 0;
	sub[nijk(NS, 2, r, h, z)] = base[nijk(NB, 2, r, a, z)];
	sub[nijk(NS, 1, r, h, z)] = base[nijk(NB, 1, r, a, z)] + base[nijk(NB, 1, r, a+1, z)];
}

// sub-grid fields to the base grid at the r lines r0 ... r0+count-1; gid = (z, a, r-r0)
// pos/w: [n][a], n=0 for the r and z components, n=1 for the alpha component
kernel void mg_interpolate(device float* base          [[buffer(0)]],
                           const device float* sub     [[buffer(1)]],
                           const device uint* pos_2p   [[buffer(2)]],
                           const device uint* pos_2pp  [[buffer(3)]],
                           constant GridDim& NB        [[buffer(4)]],
                           constant GridDim& NS        [[buffer(5)]],
                           constant uint& r0           [[buffer(6)]],
                           const device float* w_2p    [[buffer(7)]],
                           const device float* w_2pp   [[buffer(8)]],
                           uint3 gid [[thread_position_in_grid]])
{
	const uint z = gid.x, a = gid.y, r = r0 + gid.z;
	const uint k1 = NB.ny + a;
	base[nijk(NB, 0, r, a, z)] = w_2p[a]*sub[nijk(NS, 0, r, pos_2p[a], z)] + w_2pp[a]*sub[nijk(NS, 0, r, pos_2pp[a], z)];
	base[nijk(NB, 2, r, a, z)] = w_2p[a]*sub[nijk(NS, 2, r, pos_2p[a], z)] + w_2pp[a]*sub[nijk(NS, 2, r, pos_2pp[a], z)];
	base[nijk(NB, 1, r, a, z)] = w_2p[k1]*sub[nijk(NS, 1, r, pos_2p[k1], z)] + w_2pp[k1]*sub[nijk(NS, 1, r, pos_2pp[k1], z)];
}
)MSL";

class Metal_MultiGridLink : public GPU_MultiGridLink
{
public:
	Metal_MultiGridLink(GPU_Backend_Metal::Impl* base, GPU_Backend_Metal::Impl* sub, const GPU_MultiGridInterpolation& interpol);

	virtual void SyncVoltages();
	virtual void SyncCurrents();
	virtual void InterpolateToBase();

protected:
	struct Weights { id<MTLBuffer> pos_2p, pos_2pp, w_2p, w_2pp; };
	Weights Upload(const std::vector<unsigned int> pos_2p[2], const std::vector<unsigned int> pos_2pp[2],
	               const std::vector<float> w_2p[2], const std::vector<float> w_2pp[2]);
	void Interpolate(id<MTLBuffer> base_field, id<MTLBuffer> sub_field, const Weights& w, unsigned int r0, unsigned int count);

	GPU_Backend_Metal::Impl* b;   //!< base grid, owned by its backend
	GPU_Backend_Metal::Impl* s;   //!< sub-grid, owned by its backend
	unsigned int m_SplitPos;
	Weights m_Volt, m_Curr;
};

Metal_MultiGridLink::Metal_MultiGridLink(GPU_Backend_Metal::Impl* base, GPU_Backend_Metal::Impl* sub, const GPU_MultiGridInterpolation& interpol)
{
	b = base;
	s = sub;
	m_SplitPos = interpol.split_pos;
	m_Volt = Upload(interpol.pos_v_2p, interpol.pos_v_2pp, interpol.w_v_2p, interpol.w_v_2pp);
	m_Curr = Upload(interpol.pos_i_2p, interpol.pos_i_2pp, interpol.w_i_2p, interpol.w_i_2pp);
	b->Pipeline(MULTIGRID_SOURCE, "mg_sync_volt");
	b->Pipeline(MULTIGRID_SOURCE, "mg_interpolate");
}

Metal_MultiGridLink::Weights Metal_MultiGridLink::Upload(const std::vector<unsigned int> pos_2p[2], const std::vector<unsigned int> pos_2pp[2],
                                                         const std::vector<float> w_2p[2], const std::vector<float> w_2pp[2])
{
	// [n][a] with n=0,1 and a < number of base alpha lines
	std::vector<uint32_t> p2p, p2pp;
	std::vector<float> f2p, f2pp;
	for (int n=0; n<2; ++n)
	{
		p2p.insert(p2p.end(), pos_2p[n].begin(), pos_2p[n].end());
		p2pp.insert(p2pp.end(), pos_2pp[n].begin(), pos_2pp[n].end());
		f2p.insert(f2p.end(), w_2p[n].begin(), w_2p[n].end());
		f2pp.insert(f2pp.end(), w_2pp[n].begin(), w_2pp[n].end());
	}
	Weights w;
	w.pos_2p = b->NewBuffer(p2p.size()*sizeof(uint32_t), p2p.data());
	w.pos_2pp = b->NewBuffer(p2pp.size()*sizeof(uint32_t), p2pp.data());
	w.w_2p = b->NewBuffer(f2p.size()*sizeof(float), f2p.data());
	w.w_2pp = b->NewBuffer(f2pp.size()*sizeof(float), f2pp.data());
	return w;
}

void Metal_MultiGridLink::SyncVoltages()
{
	if (m_SplitPos<1)
		return;
	const uint32_t r = m_SplitPos-1;
	id<MTLComputePipelineState> pso = b->Pipeline(MULTIGRID_SOURCE, "mg_sync_volt");
	id<MTLComputeCommandEncoder> enc = b->Encoder();
	[enc setComputePipelineState:pso];
	[enc setBuffer:s->volt offset:0 atIndex:0];
	[enc setBuffer:b->volt offset:0 atIndex:1];
	[enc setBytes:&b->dim length:sizeof(b->dim) atIndex:4];
	[enc setBytes:&s->dim length:sizeof(s->dim) atIndex:5];
	[enc setBytes:&r length:sizeof(r) atIndex:6];
	b->Dispatch(pso, b->dim.nz, b->dim.ny/2);
}

void Metal_MultiGridLink::Interpolate(id<MTLBuffer> base_field, id<MTLBuffer> sub_field, const Weights& w, unsigned int r0, unsigned int count)
{
	if (count==0)
		return;
	const uint32_t r = r0;
	id<MTLComputePipelineState> pso = b->Pipeline(MULTIGRID_SOURCE, "mg_interpolate");
	id<MTLComputeCommandEncoder> enc = b->Encoder();
	[enc setComputePipelineState:pso];
	[enc setBuffer:base_field offset:0 atIndex:0];
	[enc setBuffer:sub_field offset:0 atIndex:1];
	[enc setBuffer:w.pos_2p offset:0 atIndex:2];
	[enc setBuffer:w.pos_2pp offset:0 atIndex:3];
	[enc setBytes:&b->dim length:sizeof(b->dim) atIndex:4];
	[enc setBytes:&s->dim length:sizeof(s->dim) atIndex:5];
	[enc setBytes:&r length:sizeof(r) atIndex:6];
	[enc setBuffer:w.w_2p offset:0 atIndex:7];
	[enc setBuffer:w.w_2pp offset:0 atIndex:8];
	b->Dispatch(pso, b->dim.nz, b->dim.ny, count);
}

void Metal_MultiGridLink::SyncCurrents()
{
	if (m_SplitPos>=2)
		Interpolate(b->curr, s->curr, m_Curr, m_SplitPos-2, 1);
}

void Metal_MultiGridLink::InterpolateToBase()
{
	if (m_SplitPos>=1)
		Interpolate(b->volt, s->volt, m_Volt, 0, m_SplitPos-1);
	if (m_SplitPos>=2)
		Interpolate(b->curr, s->curr, m_Curr, 0, m_SplitPos-2);
}

GPU_MultiGridLink* GPU_Backend_Metal::CreateMultiGridLink(GPU_Backend* sub_grid, const GPU_MultiGridInterpolation& interpol)
{
	GPU_Backend_Metal* sub = dynamic_cast<GPU_Backend_Metal*>(sub_grid);
	if (!sub || sub->d->ctx!=d->ctx)
		return NULL;
	return new Metal_MultiGridLink(d, sub->d, interpol);
}
