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
#include "FDTD/engine.h"
#include "FDTD/extensions/operator_ext_absorbing_bc.h"
#include "FDTD/extensions/engine_ext_absorbing_bc.h"

// Local absorbing sheet (first order Mur, optionally with super-absorption of the
// currents), see Engine_Ext_Absorbing_BC. One thread per (i,j) of the sheet,
// i along m_nyP, j along m_nyPP; gid = (j, i). The state and coefficient arrays
// have the full sheet size ni x nj, the currents use ni-1 x nj-1 of it.
// The host divides in double precision, this runs in float.
static const char* ABSORBING_BC_SOURCE = R"MSL(
struct ABCParam { int ny, nyP, nyPP; uint line, line_shift; uint startP, startPP; uint nj; };

inline uint abc_index(constant GridDim& N, constant ABCParam& P, uint n, uint line, uint i, uint j)
{
	uint pos[3];
	pos[P.ny] = line;
	pos[P.nyP] = P.startP + i;
	pos[P.nyPP] = P.startPP + j;
	return nijk(N, n, pos[0], pos[1], pos[2]);
}

// F(i+s,n) - K1*F(i,n)
kernel void abc_pre(const device float* field  [[buffer(0)]],
                    device float* s_nyP        [[buffer(1)]],
                    device float* s_nyPP       [[buffer(2)]],
                    const device float* k1_nyP [[buffer(3)]],
                    const device float* k1_nyPP[[buffer(4)]],
                    constant GridDim& N        [[buffer(7)]],
                    constant ABCParam& P       [[buffer(8)]],
                    uint2 gid [[thread_position_in_grid]])
{
	const uint i = gid.y, j = gid.x, ij = i*P.nj + j;
	s_nyP[ij]  = field[abc_index(N, P, P.nyP,  P.line_shift, i, j)] - k1_nyP[ij]  * field[abc_index(N, P, P.nyP,  P.line, i, j)];
	s_nyPP[ij] = field[abc_index(N, P, P.nyPP, P.line_shift, i, j)] - k1_nyPP[ij] * field[abc_index(N, P, P.nyPP, P.line, i, j)];
}

// + K1*F(i+s,n+1)
kernel void abc_post(const device float* field  [[buffer(0)]],
                     device float* s_nyP        [[buffer(1)]],
                     device float* s_nyPP       [[buffer(2)]],
                     const device float* k1_nyP [[buffer(3)]],
                     const device float* k1_nyPP[[buffer(4)]],
                     constant GridDim& N        [[buffer(7)]],
                     constant ABCParam& P       [[buffer(8)]],
                     uint2 gid [[thread_position_in_grid]])
{
	const uint i = gid.y, j = gid.x, ij = i*P.nj + j;
	s_nyP[ij]  = s_nyP[ij]  + k1_nyP[ij]  * field[abc_index(N, P, P.nyP,  P.line_shift, i, j)];
	s_nyPP[ij] = s_nyPP[ij] + k1_nyPP[ij] * field[abc_index(N, P, P.nyPP, P.line_shift, i, j)];
}

// voltages: replace by the Mur value
kernel void abc_apply_volt(device float* field        [[buffer(0)]],
                           const device float* s_nyP  [[buffer(1)]],
                           const device float* s_nyPP [[buffer(2)]],
                           constant GridDim& N        [[buffer(7)]],
                           constant ABCParam& P       [[buffer(8)]],
                           uint2 gid [[thread_position_in_grid]])
{
	const uint i = gid.y, j = gid.x, ij = i*P.nj + j;
	field[abc_index(N, P, P.nyP,  P.line, i, j)] = s_nyP[ij];
	field[abc_index(N, P, P.nyPP, P.line, i, j)] = s_nyPP[ij];
}

// currents, super-absorption: (Hsa*K2 + Hc)/(K2 + 1)
kernel void abc_apply_curr(device float* field        [[buffer(0)]],
                           const device float* s_nyP  [[buffer(1)]],
                           const device float* s_nyPP [[buffer(2)]],
                           const device float* k2_nyP [[buffer(5)]],
                           const device float* k2_nyPP[[buffer(6)]],
                           constant GridDim& N        [[buffer(7)]],
                           constant ABCParam& P       [[buffer(8)]],
                           uint2 gid [[thread_position_in_grid]])
{
	const uint i = gid.y, j = gid.x, ij = i*P.nj + j;
	const uint gP  = abc_index(N, P, P.nyP,  P.line, i, j);
	const uint gPP = abc_index(N, P, P.nyPP, P.line, i, j);
	field[gP]  = (s_nyP[ij] *k2_nyP[ij]  + field[gP]) /(k2_nyP[ij]  + 1.0f);
	field[gPP] = (s_nyPP[ij]*k2_nyPP[ij] + field[gPP])/(k2_nyPP[ij] + 1.0f);
}
)MSL";

class Metal_Ext_Absorbing_BC : public GPU_Extension
{
public:
	Metal_Ext_Absorbing_BC(GPU_Backend_Metal::Impl* impl, Operator_Ext_Absorbing_BC* op_ext);

	virtual void DoPreVoltageUpdates()  {Run("abc_pre", d->volt, m_Volt_nyP, m_Volt_nyPP, m_VoltParam, m_ni, m_nj);}
	virtual void DoPostVoltageUpdates() {Run("abc_post", d->volt, m_Volt_nyP, m_Volt_nyPP, m_VoltParam, m_ni, m_nj);}
	virtual void Apply2Voltages()       {Run("abc_apply_volt", d->volt, m_Volt_nyP, m_Volt_nyPP, m_VoltParam, m_ni, m_nj);}

	// super-absorption of the currents, on the dual grid one line less
	virtual void DoPreCurrentUpdates()  {if (m_SuperAbsorption) Run("abc_pre", d->curr, m_Curr_nyP, m_Curr_nyPP, m_CurrParam, m_ni-1, m_nj-1);}
	virtual void DoPostCurrentUpdates() {if (m_SuperAbsorption) Run("abc_post", d->curr, m_Curr_nyP, m_Curr_nyPP, m_CurrParam, m_ni-1, m_nj-1);}
	virtual void Apply2Current()        {if (m_SuperAbsorption) Run("abc_apply_curr", d->curr, m_Curr_nyP, m_Curr_nyPP, m_CurrParam, m_ni-1, m_nj-1);}

protected:
	struct Param {int32_t ny, nyP, nyPP; uint32_t line, line_shift, startP, startPP, nj;};

	void Run(const char* kernel, id<MTLBuffer> field, id<MTLBuffer> s_nyP, id<MTLBuffer> s_nyPP,
	         const Param& param, unsigned int ni, unsigned int nj);

	GPU_Backend_Metal::Impl* d;
	bool m_SuperAbsorption;
	unsigned int m_ni, m_nj;
	Param m_VoltParam, m_CurrParam;

	id<MTLBuffer> m_Volt_nyP, m_Volt_nyPP, m_Curr_nyP, m_Curr_nyPP;
	id<MTLBuffer> m_K1_nyP, m_K1_nyPP, m_K2_nyP, m_K2_nyPP;
};

Metal_Ext_Absorbing_BC::Metal_Ext_Absorbing_BC(GPU_Backend_Metal::Impl* impl, Operator_Ext_Absorbing_BC* op_ext)
{
	d = impl;
	m_SuperAbsorption = (op_ext->m_ABCtype==Operator_Ext_Absorbing_BC::MUR_1ST_SA);
	m_ni = op_ext->m_numLines[0];
	m_nj = op_ext->m_numLines[1];

	// positions as in Engine_Ext_Absorbing_BC::Engine_Ext_Absorbing_BC()
	const int ny = op_ext->m_ny;
	const unsigned int start = op_ext->m_sheetX0[ny];
	const bool positive = op_ext->m_normalSignPositive;
	Param base;
	base.ny = ny;
	base.nyP = op_ext->m_nyP;
	base.nyPP = op_ext->m_nyPP;
	base.startP = op_ext->m_sheetX0[op_ext->m_nyP];
	base.startPP = op_ext->m_sheetX0[op_ext->m_nyPP];
	base.nj = m_nj;

	m_VoltParam = base;
	m_VoltParam.line = start;
	m_VoltParam.line_shift = start + (positive ? 1 : -1);

	m_CurrParam = base;
	m_CurrParam.line = start + (positive ? 0 : -1);
	m_CurrParam.line_shift = start + (positive ? 1 : -2);

	const size_t bytes = (size_t)m_ni*m_nj*sizeof(FDTD_FLOAT);
	m_Volt_nyP  = d->NewBuffer(bytes);
	m_Volt_nyPP = d->NewBuffer(bytes);
	m_Curr_nyP  = d->NewBuffer(bytes);
	m_Curr_nyPP = d->NewBuffer(bytes);
	m_K1_nyP  = d->NewBuffer(bytes, op_ext->m_K1_nyP.data());
	m_K1_nyPP = d->NewBuffer(bytes, op_ext->m_K1_nyPP.data());
	m_K2_nyP  = d->NewBuffer(bytes, op_ext->m_K2_nyP.valid() ? op_ext->m_K2_nyP.data() : NULL);
	m_K2_nyPP = d->NewBuffer(bytes, op_ext->m_K2_nyPP.valid() ? op_ext->m_K2_nyPP.data() : NULL);

	d->Pipeline(ABSORBING_BC_SOURCE, "abc_pre");
	d->Pipeline(ABSORBING_BC_SOURCE, "abc_post");
	d->Pipeline(ABSORBING_BC_SOURCE, "abc_apply_volt");
	d->Pipeline(ABSORBING_BC_SOURCE, "abc_apply_curr");
}

void Metal_Ext_Absorbing_BC::Run(const char* kernel, id<MTLBuffer> field, id<MTLBuffer> s_nyP, id<MTLBuffer> s_nyPP,
                                 const Param& param, unsigned int ni, unsigned int nj)
{
	id<MTLComputePipelineState> pso = d->Pipeline(ABSORBING_BC_SOURCE, kernel);
	id<MTLComputeCommandEncoder> enc = d->Encoder();
	[enc setComputePipelineState:pso];
	[enc setBuffer:field offset:0 atIndex:0];
	[enc setBuffer:s_nyP offset:0 atIndex:1];
	[enc setBuffer:s_nyPP offset:0 atIndex:2];
	[enc setBuffer:m_K1_nyP offset:0 atIndex:3];
	[enc setBuffer:m_K1_nyPP offset:0 atIndex:4];
	[enc setBuffer:m_K2_nyP offset:0 atIndex:5];
	[enc setBuffer:m_K2_nyPP offset:0 atIndex:6];
	d->SetGridDim(7);
	[enc setBytes:&param length:sizeof(param) atIndex:8];
	d->Dispatch(pso, nj, ni);
}

GPU_Extension* Metal_CreateExt_Absorbing_BC(GPU_Backend_Metal::Impl* d, Engine_Extension* eng_ext, Engine* eng)
{
	UNUSED(eng);
	if (!dynamic_cast<Engine_Ext_Absorbing_BC*>(eng_ext))
		return NULL;
	Operator_Ext_Absorbing_BC* op_ext = dynamic_cast<Operator_Ext_Absorbing_BC*>(eng_ext->GetOperatorExtension());
	if (!op_ext)
		return NULL;
	return new Metal_Ext_Absorbing_BC(d, op_ext);
}
