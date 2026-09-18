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
#include "FDTD/extensions/operator_ext_lorentzmaterial.h"
#include "FDTD/extensions/engine_ext_lorentzmaterial.h"

// Drude/Lorentz/Debye materials and conducting sheets, see Engine_Ext_LorentzMaterial
// and Engine_Ext_Dispersive. One dispatch per dispersion order, one thread per mesh
// position of that order (unique within an order). The ADE state and coefficients
// are stored per order as [direction][position].
static const char* LORENTZ_SOURCE = R"MSL(
struct ADEParam { uint count; uint lorentz; uint sn; };

// ADE update before the main update, reads the field
kernel void lorentz_pre(const device float* field  [[buffer(0)]],
                        device float* ade          [[buffer(1)]],
                        device float* lor          [[buffer(2)]],
                        const device float* c_int  [[buffer(3)]],
                        const device float* c_ext  [[buffer(4)]],
                        const device float* c_lor  [[buffer(5)]],
                        const device uint* pos     [[buffer(6)]],
                        constant ADEParam& P       [[buffer(7)]],
                        uint i [[thread_position_in_grid]])
{
	if (i>=P.count)
		return;
	for (uint n=0; n<3; ++n)
	{
		const uint k = n*P.count + i;
		const float f = field[n*P.sn + pos[i]];
		if (P.lorentz)
		{
			lor[k] = lor[k] + c_lor[k]*ade[k];
			ade[k] = ade[k] * c_int[k];
			ade[k] = ade[k] + c_ext[k]*(f - lor[k]);
		}
		else
		{
			ade[k] = ade[k] * c_int[k];
			ade[k] = ade[k] + c_ext[k]*f;
		}
	}
}

// subtract the ADE current from the field
kernel void dispersive_apply(device float* field      [[buffer(0)]],
                             const device float* ade  [[buffer(1)]],
                             const device uint* pos   [[buffer(6)]],
                             constant ADEParam& P     [[buffer(7)]],
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

class Metal_Ext_LorentzMaterial : public GPU_Extension
{
public:
	Metal_Ext_LorentzMaterial(GPU_Backend_Metal::Impl* impl, Operator_Ext_LorentzMaterial* op_ext);

	virtual void DoPreVoltageUpdates() {for (size_t o=0; o<m_Volt.size(); ++o) Pre(m_Volt[o], d->volt);}
	virtual void Apply2Voltages()      {for (size_t o=0; o<m_Volt.size(); ++o) Apply(m_Volt[o], d->volt);}
	virtual void DoPreCurrentUpdates() {for (size_t o=0; o<m_Curr.size(); ++o) Pre(m_Curr[o], d->curr);}
	virtual void Apply2Current()       {for (size_t o=0; o<m_Curr.size(); ++o) Apply(m_Curr[o], d->curr);}

protected:
	//! one dispersion order of the voltages or currents
	struct Order
	{
		struct {uint32_t count, lorentz, sn;} param;
		id<MTLBuffer> pos;          //!< flat NIJK index of direction 0
		id<MTLBuffer> ade, lor;     //!< state
		id<MTLBuffer> c_int, c_ext, c_lor;
	};

	void Setup(std::vector<Order>& orders, int o, unsigned int count, unsigned int** pos,
	           bool lorentz, FDTD_FLOAT** c_int, FDTD_FLOAT** c_ext, FDTD_FLOAT** c_lor);
	id<MTLBuffer> Coefficients(unsigned int count, FDTD_FLOAT** c);
	void Pre(Order& o, id<MTLBuffer> field);
	void Apply(Order& o, id<MTLBuffer> field);

	GPU_Backend_Metal::Impl* d;
	std::vector<Order> m_Volt;
	std::vector<Order> m_Curr;
};

Metal_Ext_LorentzMaterial::Metal_Ext_LorentzMaterial(GPU_Backend_Metal::Impl* impl, Operator_Ext_LorentzMaterial* op_ext)
{
	d = impl;
	for (int o=0; o<op_ext->m_Order; ++o)
	{
		const unsigned int count = op_ext->m_LM_Count.at(o);
		if (op_ext->m_volt_ADE_On[o])
			Setup(m_Volt, o, count, op_ext->m_LM_pos[o], op_ext->m_volt_Lor_ADE_On[o],
			      op_ext->v_int_ADE[o], op_ext->v_ext_ADE[o], op_ext->m_volt_Lor_ADE_On[o] ? op_ext->v_Lor_ADE[o] : NULL);
		if (op_ext->m_curr_ADE_On[o])
			Setup(m_Curr, o, count, op_ext->m_LM_pos[o], op_ext->m_curr_Lor_ADE_On[o],
			      op_ext->i_int_ADE[o], op_ext->i_ext_ADE[o], op_ext->m_curr_Lor_ADE_On[o] ? op_ext->i_Lor_ADE[o] : NULL);
	}
	d->Pipeline(LORENTZ_SOURCE, "lorentz_pre");
	d->Pipeline(LORENTZ_SOURCE, "dispersive_apply");
}

id<MTLBuffer> Metal_Ext_LorentzMaterial::Coefficients(unsigned int count, FDTD_FLOAT** c)
{
	std::vector<float> data(3*(size_t)count, 0);
	if (c)
		for (int n=0; n<3; ++n)
			for (unsigned int i=0; i<count; ++i)
				data[n*count + i] = c[n][i];
	return d->NewBuffer(data.size()*sizeof(float), data.data());
}

void Metal_Ext_LorentzMaterial::Setup(std::vector<Order>& orders, int o, unsigned int count, unsigned int** pos,
                                      bool lorentz, FDTD_FLOAT** c_int, FDTD_FLOAT** c_ext, FDTD_FLOAT** c_lor)
{
	UNUSED(o);
	if (count==0)
		return;
	Order ord;
	ord.param.count = count;
	ord.param.lorentz = lorentz;
	ord.param.sn = d->numCells;

	std::vector<uint32_t> flat(count);
	for (unsigned int i=0; i<count; ++i)
		flat[i] = (pos[0][i]*d->dim.ny + pos[1][i])*d->dim.nz + pos[2][i];
	ord.pos = d->NewBuffer(count*sizeof(uint32_t), flat.data());

	ord.ade = d->NewBuffer(3*(size_t)count*sizeof(float));
	ord.lor = d->NewBuffer(3*(size_t)count*sizeof(float));
	ord.c_int = Coefficients(count, c_int);
	ord.c_ext = Coefficients(count, c_ext);
	ord.c_lor = Coefficients(count, c_lor);
	orders.push_back(ord);
}

void Metal_Ext_LorentzMaterial::Pre(Order& o, id<MTLBuffer> field)
{
	id<MTLComputePipelineState> pso = d->Pipeline(LORENTZ_SOURCE, "lorentz_pre");
	id<MTLComputeCommandEncoder> enc = d->Encoder();
	[enc setComputePipelineState:pso];
	[enc setBuffer:field offset:0 atIndex:0];
	[enc setBuffer:o.ade offset:0 atIndex:1];
	[enc setBuffer:o.lor offset:0 atIndex:2];
	[enc setBuffer:o.c_int offset:0 atIndex:3];
	[enc setBuffer:o.c_ext offset:0 atIndex:4];
	[enc setBuffer:o.c_lor offset:0 atIndex:5];
	[enc setBuffer:o.pos offset:0 atIndex:6];
	[enc setBytes:&o.param length:sizeof(o.param) atIndex:7];
	d->Dispatch(pso, o.param.count);
}

void Metal_Ext_LorentzMaterial::Apply(Order& o, id<MTLBuffer> field)
{
	id<MTLComputePipelineState> pso = d->Pipeline(LORENTZ_SOURCE, "dispersive_apply");
	id<MTLComputeCommandEncoder> enc = d->Encoder();
	[enc setComputePipelineState:pso];
	[enc setBuffer:field offset:0 atIndex:0];
	[enc setBuffer:o.ade offset:0 atIndex:1];
	[enc setBuffer:o.pos offset:0 atIndex:6];
	[enc setBytes:&o.param length:sizeof(o.param) atIndex:7];
	d->Dispatch(pso, o.param.count);
}

GPU_Extension* Metal_CreateExt_LorentzMaterial(GPU_Backend_Metal::Impl* d, Engine_Extension* eng_ext, Engine* eng)
{
	UNUSED(eng);
	if (!dynamic_cast<Engine_Ext_LorentzMaterial*>(eng_ext))
		return NULL;
	Operator_Ext_LorentzMaterial* op_ext = dynamic_cast<Operator_Ext_LorentzMaterial*>(eng_ext->GetOperatorExtension());
	if (!op_ext)
		return NULL;
	return new Metal_Ext_LorentzMaterial(d, op_ext);
}
