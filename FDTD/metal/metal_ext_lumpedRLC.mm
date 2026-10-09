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

#include <set>

#include "metal_internal.h"
#include "FDTD/engine.h"
#include "FDTD/extensions/operator_ext_lumpedRLC.h"
#include "FDTD/extensions/engine_ext_lumpedRLC.h"

// Series/parallel lumped RLC elements, see Engine_Ext_LumpedRLC.
// The histories Vd[n], Vd[n-1], Vd[n-2] and J[n], J[n-1], J[n-2] are three slots
// of count values each; the host rotates the slot indices instead of pointers.
static const char* RLC_SOURCE = R"MSL(
struct RLCParam { uint count; uint vd0, vd1, vd2; uint j0, j1, j2; };

struct RLCCoeff
{
	const device float* ilv; const device float* i2v;
	const device float* vv2; const device float* vj1; const device float* vj2; const device float* vvd;
	const device float* ib0; const device float* b1;  const device float* b2;
};

// parallel inductor current
kernel void rlc_pre(device float* il           [[buffer(0)]],
                    const device float* vd     [[buffer(1)]],
                    const device float* i2v    [[buffer(2)]],
                    const device float* ilv    [[buffer(3)]],
                    constant RLCParam& P       [[buffer(4)]],
                    uint p [[thread_position_in_grid]])
{
	if (p<P.count)
		il[p] = il[p] + (i2v[p]*ilv[p])*vd[P.vd1*P.count + p];
}

inline float rlc_vd(uint p, float v, const device float* il, const device float* vd, const device float* j,
                    const device float* vvd, const device float* vv2, const device float* vj1, const device float* vj2,
                    constant RLCParam& P)
{
	return vvd[p]*(v - il[p] + vv2[p]*vd[P.vd2*P.count + p] + vj1[p]*j[P.j1*P.count + p] + vj2[p]*j[P.j2*P.count + p]);
}

inline float rlc_j(uint p, const device float* vd, const device float* j,
                   const device float* ib0, const device float* b1, const device float* b2, constant RLCParam& P)
{
	return ib0[p]*(vd[P.vd0*P.count + p] - vd[P.vd2*P.count + p])
	       - (b1[p]*ib0[p])*j[P.j1*P.count + p]
	       - (b2[p]*ib0[p])*j[P.j2*P.count + p];
}

// node voltage update, one thread per element (all elements on different edges)
kernel void rlc_apply(device float* volt          [[buffer(0)]],
                      const device float* il      [[buffer(1)]],
                      device float* vd            [[buffer(2)]],
                      device float* j             [[buffer(3)]],
                      const device uint* index    [[buffer(4)]],
                      const device float* vvd     [[buffer(5)]],
                      const device float* vv2     [[buffer(6)]],
                      const device float* vj1     [[buffer(7)]],
                      const device float* vj2     [[buffer(8)]],
                      const device float* ib0     [[buffer(9)]],
                      const device float* b1      [[buffer(10)]],
                      const device float* b2      [[buffer(11)]],
                      constant RLCParam& P        [[buffer(12)]],
                      uint p [[thread_position_in_grid]])
{
	if (p>=P.count)
		return;
	vd[P.vd0*P.count + p] = rlc_vd(p, volt[index[p]], il, vd, j, vvd, vv2, vj1, vj2, P);
	j[P.j0*P.count + p] = rlc_j(p, vd, j, ib0, b1, b2, P);
	volt[index[p]] = vd[P.vd0*P.count + p];
}

// same in the three passes of the host (read all, update all, write all), for elements sharing an edge
kernel void rlc_apply_serial(device float* volt          [[buffer(0)]],
                             const device float* il      [[buffer(1)]],
                             device float* vd            [[buffer(2)]],
                             device float* j             [[buffer(3)]],
                             const device uint* index    [[buffer(4)]],
                             const device float* vvd     [[buffer(5)]],
                             const device float* vv2     [[buffer(6)]],
                             const device float* vj1     [[buffer(7)]],
                             const device float* vj2     [[buffer(8)]],
                             const device float* ib0     [[buffer(9)]],
                             const device float* b1      [[buffer(10)]],
                             const device float* b2      [[buffer(11)]],
                             constant RLCParam& P        [[buffer(12)]])
{
	for (uint p=0; p<P.count; ++p)
		vd[P.vd0*P.count + p] = volt[index[p]];
	for (uint p=0; p<P.count; ++p)
	{
		vd[P.vd0*P.count + p] = rlc_vd(p, vd[P.vd0*P.count + p], il, vd, j, vvd, vv2, vj1, vj2, P);
		j[P.j0*P.count + p] = rlc_j(p, vd, j, ib0, b1, b2, P);
	}
	for (uint p=0; p<P.count; ++p)
		volt[index[p]] = vd[P.vd0*P.count + p];
}
)MSL";

class Metal_Ext_LumpedRLC : public GPU_Extension
{
public:
	Metal_Ext_LumpedRLC(GPU_Backend_Metal::Impl* impl, Operator_Ext_LumpedRLC* op_ext);

	virtual void DoPreVoltageUpdates();
	virtual void Apply2Voltages();

protected:
	GPU_Backend_Metal::Impl* d;
	bool m_Serial;   //!< elements share an edge
	struct {uint32_t count, vd0, vd1, vd2, j0, j1, j2;} m_Param;

	id<MTLBuffer> m_Il, m_Vd, m_J, m_Index;
	id<MTLBuffer> m_ilv, m_i2v, m_vv2, m_vj1, m_vj2, m_vvd, m_ib0, m_b1, m_b2;
};

Metal_Ext_LumpedRLC::Metal_Ext_LumpedRLC(GPU_Backend_Metal::Impl* impl, Operator_Ext_LumpedRLC* op_ext)
{
	d = impl;
	const unsigned int count = op_ext->RLC_count;
	m_Param.count = count;
	// slot order as initialized by the host: v_Vdn[k] and v_Jn[k] in slot k
	m_Param.vd0 = 0; m_Param.vd1 = 1; m_Param.vd2 = 2;
	m_Param.j0 = 0;  m_Param.j1 = 1;  m_Param.j2 = 2;
	m_Serial = false;
	if (count==0)
		return;

	std::vector<uint32_t> index(count);
	std::set<uint32_t> edges;
	for (unsigned int p=0; p<count; ++p)
	{
		index[p] = ((op_ext->v_RLC_dir[p]*d->dim.nx + op_ext->v_RLC_pos[0][p])*d->dim.ny + op_ext->v_RLC_pos[1][p])*d->dim.nz + op_ext->v_RLC_pos[2][p];
		if (!edges.insert(index[p]).second)
			m_Serial = true;
	}
	const size_t bytes = count*sizeof(float);
	m_Index = d->NewBuffer(count*sizeof(uint32_t), index.data());
	m_Il = d->NewBuffer(bytes);
	m_Vd = d->NewBuffer(3*bytes);
	m_J  = d->NewBuffer(3*bytes);
	m_ilv = d->NewBuffer(bytes, op_ext->v_RLC_ilv);
	m_i2v = d->NewBuffer(bytes, op_ext->v_RLC_i2v);
	m_vv2 = d->NewBuffer(bytes, op_ext->v_RLC_vv2);
	m_vj1 = d->NewBuffer(bytes, op_ext->v_RLC_vj1);
	m_vj2 = d->NewBuffer(bytes, op_ext->v_RLC_vj2);
	m_vvd = d->NewBuffer(bytes, op_ext->v_RLC_vvd);
	m_ib0 = d->NewBuffer(bytes, op_ext->v_RLC_ib0);
	m_b1  = d->NewBuffer(bytes, op_ext->v_RLC_b1);
	m_b2  = d->NewBuffer(bytes, op_ext->v_RLC_b2);

	d->Pipeline(RLC_SOURCE, "rlc_pre");
	d->Pipeline(RLC_SOURCE, m_Serial ? "rlc_apply_serial" : "rlc_apply");
}

void Metal_Ext_LumpedRLC::DoPreVoltageUpdates()
{
	if (m_Param.count==0)
		return;
	// rotate the Vd history: [0,1,2] <- [2,0,1]
	const uint32_t vd2 = m_Param.vd2;
	m_Param.vd2 = m_Param.vd1;
	m_Param.vd1 = m_Param.vd0;
	m_Param.vd0 = vd2;

	id<MTLComputePipelineState> pso = d->Pipeline(RLC_SOURCE, "rlc_pre");
	id<MTLComputeCommandEncoder> enc = d->Encoder();
	[enc setComputePipelineState:pso];
	[enc setBuffer:m_Il offset:0 atIndex:0];
	[enc setBuffer:m_Vd offset:0 atIndex:1];
	[enc setBuffer:m_i2v offset:0 atIndex:2];
	[enc setBuffer:m_ilv offset:0 atIndex:3];
	[enc setBytes:&m_Param length:sizeof(m_Param) atIndex:4];
	d->Dispatch(pso, m_Param.count);
}

void Metal_Ext_LumpedRLC::Apply2Voltages()
{
	if (m_Param.count==0)
		return;
	// rotate the J history: [0,1,2] <- [2,0,1]
	const uint32_t j2 = m_Param.j2;
	m_Param.j2 = m_Param.j1;
	m_Param.j1 = m_Param.j0;
	m_Param.j0 = j2;

	id<MTLComputePipelineState> pso = d->Pipeline(RLC_SOURCE, m_Serial ? "rlc_apply_serial" : "rlc_apply");
	id<MTLComputeCommandEncoder> enc = d->Encoder();
	[enc setComputePipelineState:pso];
	[enc setBuffer:d->volt offset:0 atIndex:0];
	[enc setBuffer:m_Il offset:0 atIndex:1];
	[enc setBuffer:m_Vd offset:0 atIndex:2];
	[enc setBuffer:m_J offset:0 atIndex:3];
	[enc setBuffer:m_Index offset:0 atIndex:4];
	[enc setBuffer:m_vvd offset:0 atIndex:5];
	[enc setBuffer:m_vv2 offset:0 atIndex:6];
	[enc setBuffer:m_vj1 offset:0 atIndex:7];
	[enc setBuffer:m_vj2 offset:0 atIndex:8];
	[enc setBuffer:m_ib0 offset:0 atIndex:9];
	[enc setBuffer:m_b1 offset:0 atIndex:10];
	[enc setBuffer:m_b2 offset:0 atIndex:11];
	[enc setBytes:&m_Param length:sizeof(m_Param) atIndex:12];
	d->Dispatch(pso, m_Serial ? 1 : m_Param.count);
}

GPU_Extension* Metal_CreateExt_LumpedRLC(GPU_Backend_Metal::Impl* d, Engine_Extension* eng_ext, Engine* eng)
{
	UNUSED(eng);
	if (!dynamic_cast<Engine_Ext_LumpedRLC*>(eng_ext))
		return NULL;
	Operator_Ext_LumpedRLC* op_ext = dynamic_cast<Operator_Ext_LumpedRLC*>(eng_ext->GetOperatorExtension());
	if (!op_ext)
		return NULL;
	return new Metal_Ext_LumpedRLC(d, op_ext);
}
