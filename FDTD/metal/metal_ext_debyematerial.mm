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
#include "FDTD/extensions/operator_ext_debyematerial.h"
#include "FDTD/extensions/engine_ext_debyematerial.h"

// Debye materials, see Engine_Ext_DebyeMaterial. One thread per mesh position, all
// poles in the same thread, since the implicit solve couples them. The pole states
// and coefficients are stored as [pole][direction][position], the rest as
// [direction][position].
static const char* DEBYE_SOURCE = R"MSL(
struct DebyeParam { uint count; uint poles; uint sn; };

// before the main update, the only place V^n exists: the part of the new pole
// states V^n fixes, and its sum
kernel void debye_pre(const device float* field   [[buffer(0)]],
                      device float* pole          [[buffer(1)]],
                      device float* pre           [[buffer(2)]],
                      const device float* relax   [[buffer(3)]],
                      const device float* drive   [[buffer(4)]],
                      const device uint* pos      [[buffer(5)]],
                      constant DebyeParam& P      [[buffer(6)]],
                      uint i [[thread_position_in_grid]])
{
	if (i>=P.count)
		return;
	for (uint n=0; n<3; ++n)
	{
		const float volt = field[n*P.sn + pos[i]];
		float sum = 0;
		for (uint o=0; o<P.poles; ++o)
		{
			const uint k = (o*3 + n)*P.count + i;
			const float W = relax[k]*pole[k] + drive[k]*volt;
			sum = sum + (W - pole[k]);
			pole[k] = W;
		}
		pre[n*P.count + i] = sum;
	}
}

// after the main update: solve for V^(n+1), keep the correction for
// dispersive_apply and finish the pole states
kernel void debye_post(const device float* field   [[buffer(0)]],
                       device float* pole          [[buffer(1)]],
                       const device float* pre     [[buffer(2)]],
                       device float* ade           [[buffer(3)]],
                       const device float* drive   [[buffer(4)]],
                       const device float* solve   [[buffer(5)]],
                       const device uint* pos      [[buffer(6)]],
                       constant DebyeParam& P      [[buffer(7)]],
                       uint i [[thread_position_in_grid]])
{
	if (i>=P.count)
		return;
	for (uint n=0; n<3; ++n)
	{
		const uint j = n*P.count + i;
		const float v_raw = field[n*P.sn + pos[i]];
		const float volt = (v_raw - pre[j]) * solve[j];
		ade[j] = v_raw - volt;
		for (uint o=0; o<P.poles; ++o)
		{
			const uint k = (o*3 + n)*P.count + i;
			pole[k] = pole[k] + drive[k]*volt;
		}
	}
}

// subtract the correction from the field
kernel void debye_apply(device float* field      [[buffer(0)]],
                        const device float* ade  [[buffer(1)]],
                        const device uint* pos   [[buffer(2)]],
                        constant DebyeParam& P   [[buffer(3)]],
                        uint i [[thread_position_in_grid]])
{
	if (i>=P.count)
		return;
	for (uint n=0; n<3; ++n)
	{
		const uint g = n*P.sn + pos[i];
		field[g] = field[g] - ade[n*P.count + i];
	}
}
)MSL";

class Metal_Ext_DebyeMaterial : public GPU_Extension
{
public:
	Metal_Ext_DebyeMaterial(GPU_Backend_Metal::Impl* impl, Operator_Ext_DebyeMaterial* op_ext);

	virtual void DoPreVoltageUpdates();
	virtual void DoPostVoltageUpdates();
	virtual void Apply2Voltages();

protected:
	id<MTLBuffer> Coefficients(int poles, FDTD_FLOAT*** c);

	GPU_Backend_Metal::Impl* d;
	struct {uint32_t count, poles, sn;} m_Param;
	id<MTLBuffer> m_Pos;              //!< flat NIJK index of direction 0
	id<MTLBuffer> m_Pole, m_Pre, m_ADE;   //!< state
	id<MTLBuffer> m_Relax, m_Drive, m_Solve;
};

Metal_Ext_DebyeMaterial::Metal_Ext_DebyeMaterial(GPU_Backend_Metal::Impl* impl, Operator_Ext_DebyeMaterial* op_ext)
{
	d = impl;
	const unsigned int count = op_ext->m_LM_Count.at(0);
	m_Param.count = count;
	m_Param.poles = op_ext->m_PoleCount;
	m_Param.sn = d->numCells;
	if (count==0)
		return;

	unsigned int** pos = op_ext->m_LM_pos[0];
	std::vector<uint32_t> flat(count);
	for (unsigned int i=0; i<count; ++i)
		flat[i] = (pos[0][i]*d->dim.ny + pos[1][i])*d->dim.nz + pos[2][i];
	m_Pos = d->NewBuffer(count*sizeof(uint32_t), flat.data());

	m_Pole  = d->NewBuffer((size_t)m_Param.poles*3*count*sizeof(float));
	m_Pre   = d->NewBuffer(3*(size_t)count*sizeof(float));
	m_ADE   = d->NewBuffer(3*(size_t)count*sizeof(float));
	m_Relax = Coefficients(m_Param.poles, op_ext->v_relax_ADE);
	m_Drive = Coefficients(m_Param.poles, op_ext->v_drive_ADE);
	m_Solve = Coefficients(1, &op_ext->v_solve_ADE);

	d->Pipeline(DEBYE_SOURCE, "debye_pre");
	d->Pipeline(DEBYE_SOURCE, "debye_post");
	d->Pipeline(DEBYE_SOURCE, "debye_apply");
}

id<MTLBuffer> Metal_Ext_DebyeMaterial::Coefficients(int poles, FDTD_FLOAT*** c)
{
	const unsigned int count = m_Param.count;
	std::vector<float> data((size_t)poles*3*count);
	for (int o=0; o<poles; ++o)
		for (int n=0; n<3; ++n)
			for (unsigned int i=0; i<count; ++i)
				data[((size_t)o*3 + n)*count + i] = c[o][n][i];
	return d->NewBuffer(data.size()*sizeof(float), data.data());
}

void Metal_Ext_DebyeMaterial::DoPreVoltageUpdates()
{
	if (m_Param.count==0)
		return;
	id<MTLComputePipelineState> pso = d->Pipeline(DEBYE_SOURCE, "debye_pre");
	id<MTLComputeCommandEncoder> enc = d->Encoder();
	[enc setComputePipelineState:pso];
	[enc setBuffer:d->volt offset:0 atIndex:0];
	[enc setBuffer:m_Pole offset:0 atIndex:1];
	[enc setBuffer:m_Pre offset:0 atIndex:2];
	[enc setBuffer:m_Relax offset:0 atIndex:3];
	[enc setBuffer:m_Drive offset:0 atIndex:4];
	[enc setBuffer:m_Pos offset:0 atIndex:5];
	[enc setBytes:&m_Param length:sizeof(m_Param) atIndex:6];
	d->Dispatch(pso, m_Param.count);
}

void Metal_Ext_DebyeMaterial::DoPostVoltageUpdates()
{
	if (m_Param.count==0)
		return;
	id<MTLComputePipelineState> pso = d->Pipeline(DEBYE_SOURCE, "debye_post");
	id<MTLComputeCommandEncoder> enc = d->Encoder();
	[enc setComputePipelineState:pso];
	[enc setBuffer:d->volt offset:0 atIndex:0];
	[enc setBuffer:m_Pole offset:0 atIndex:1];
	[enc setBuffer:m_Pre offset:0 atIndex:2];
	[enc setBuffer:m_ADE offset:0 atIndex:3];
	[enc setBuffer:m_Drive offset:0 atIndex:4];
	[enc setBuffer:m_Solve offset:0 atIndex:5];
	[enc setBuffer:m_Pos offset:0 atIndex:6];
	[enc setBytes:&m_Param length:sizeof(m_Param) atIndex:7];
	d->Dispatch(pso, m_Param.count);
}

void Metal_Ext_DebyeMaterial::Apply2Voltages()
{
	if (m_Param.count==0)
		return;
	id<MTLComputePipelineState> pso = d->Pipeline(DEBYE_SOURCE, "debye_apply");
	id<MTLComputeCommandEncoder> enc = d->Encoder();
	[enc setComputePipelineState:pso];
	[enc setBuffer:d->volt offset:0 atIndex:0];
	[enc setBuffer:m_ADE offset:0 atIndex:1];
	[enc setBuffer:m_Pos offset:0 atIndex:2];
	[enc setBytes:&m_Param length:sizeof(m_Param) atIndex:3];
	d->Dispatch(pso, m_Param.count);
}

GPU_Extension* Metal_CreateExt_DebyeMaterial(GPU_Backend_Metal::Impl* d, Engine_Extension* eng_ext, Engine* eng)
{
	UNUSED(eng);
	if (!dynamic_cast<Engine_Ext_DebyeMaterial*>(eng_ext))
		return NULL;
	Operator_Ext_DebyeMaterial* op_ext = dynamic_cast<Operator_Ext_DebyeMaterial*>(eng_ext->GetOperatorExtension());
	if (!op_ext)
		return NULL;
	return new Metal_Ext_DebyeMaterial(d, op_ext);
}
