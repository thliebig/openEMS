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
#include "FDTD/engine.h"
#include "FDTD/extensions/operator_ext_lorentzmaterial.h"
#include "FDTD/extensions/engine_ext_lorentzmaterial.h"

// Drude/Lorentz/Debye materials and conducting sheets, see Engine_Ext_LorentzMaterial
// and Engine_Ext_Dispersive. One launch per dispersion order, one thread per mesh
// position of that order (unique within an order). The ADE state and coefficients
// are stored per order as [direction][position].
struct ADEParam { unsigned int count; unsigned int lorentz; unsigned int sn; };

// ADE update before the main update, reads the field
__global__ void lorentz_pre(const float* field, float* ade, float* lor, const float* c_int, const float* c_ext,
                            const float* c_lor, const unsigned int* pos, ADEParam P)
{
	const unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i>=P.count)
		return;
	for (unsigned int n=0; n<3; ++n)
	{
		const unsigned int k = n*P.count + i;
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
__global__ void dispersive_apply(float* field, const float* ade, const unsigned int* pos, ADEParam P)
{
	const unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i>=P.count)
		return;
	for (unsigned int n=0; n<3; ++n)
	{
		const unsigned int g = n*P.sn + pos[i];
		field[g] = field[g] - ade[n*P.count + i];
	}
}

class HIP_Ext_LorentzMaterial : public GPU_Extension
{
public:
	HIP_Ext_LorentzMaterial(GPU_Backend_HIP::Impl* impl, Operator_Ext_LorentzMaterial* op_ext);

	virtual void DoPreVoltageUpdates() {for (size_t o=0; o<m_Volt.size(); ++o) Pre(m_Volt[o], d->volt);}
	virtual void Apply2Voltages()      {for (size_t o=0; o<m_Volt.size(); ++o) Apply(m_Volt[o], d->volt);}
	virtual void DoPreCurrentUpdates() {for (size_t o=0; o<m_Curr.size(); ++o) Pre(m_Curr[o], d->curr);}
	virtual void Apply2Current()       {for (size_t o=0; o<m_Curr.size(); ++o) Apply(m_Curr[o], d->curr);}

protected:
	//! one dispersion order of the voltages or currents
	struct Order
	{
		ADEParam param;
		unsigned int* pos;          //!< flat NIJK index of direction 0
		float *ade, *lor;           //!< state
		float *c_int, *c_ext, *c_lor;
	};

	void Setup(std::vector<Order>& orders, unsigned int count, unsigned int** pos,
	           bool lorentz, FDTD_FLOAT** c_int, FDTD_FLOAT** c_ext, FDTD_FLOAT** c_lor);
	float* Coefficients(unsigned int count, FDTD_FLOAT** c);
	void Pre(Order& o, float* field)
	{
		HIP_Launch(d, "lorentz_pre", lorentz_pre, o.param.count, 1, 1, (const float*)field, o.ade, o.lor,
		            (const float*)o.c_int, (const float*)o.c_ext, (const float*)o.c_lor, (const unsigned int*)o.pos, o.param);
	}
	void Apply(Order& o, float* field)
	{
		HIP_Launch(d, "dispersive_apply", dispersive_apply, o.param.count, 1, 1, field, (const float*)o.ade, (const unsigned int*)o.pos, o.param);
	}

	GPU_Backend_HIP::Impl* d;
	std::vector<Order> m_Volt;
	std::vector<Order> m_Curr;
};

HIP_Ext_LorentzMaterial::HIP_Ext_LorentzMaterial(GPU_Backend_HIP::Impl* impl, Operator_Ext_LorentzMaterial* op_ext)
{
	d = impl;
	for (int o=0; o<op_ext->m_Order; ++o)
	{
		const unsigned int count = op_ext->m_LM_Count.at(o);
		if (op_ext->m_volt_ADE_On[o])
			Setup(m_Volt, count, op_ext->m_LM_pos[o], op_ext->m_volt_Lor_ADE_On[o],
			      op_ext->v_int_ADE[o], op_ext->v_ext_ADE[o], op_ext->m_volt_Lor_ADE_On[o] ? op_ext->v_Lor_ADE[o] : NULL);
		if (op_ext->m_curr_ADE_On[o])
			Setup(m_Curr, count, op_ext->m_LM_pos[o], op_ext->m_curr_Lor_ADE_On[o],
			      op_ext->i_int_ADE[o], op_ext->i_ext_ADE[o], op_ext->m_curr_Lor_ADE_On[o] ? op_ext->i_Lor_ADE[o] : NULL);
	}
}

float* HIP_Ext_LorentzMaterial::Coefficients(unsigned int count, FDTD_FLOAT** c)
{
	std::vector<float> data(3*(size_t)count, 0);
	if (c)
		for (int n=0; n<3; ++n)
			for (unsigned int i=0; i<count; ++i)
				data[n*count + i] = c[n][i];
	return d->Alloc<float>(data.size(), data.data());
}

void HIP_Ext_LorentzMaterial::Setup(std::vector<Order>& orders, unsigned int count, unsigned int** pos,
                                     bool lorentz, FDTD_FLOAT** c_int, FDTD_FLOAT** c_ext, FDTD_FLOAT** c_lor)
{
	if (count==0)
		return;
	Order ord;
	ord.param.count = count;
	ord.param.lorentz = lorentz;
	ord.param.sn = d->numCells;

	std::vector<unsigned int> flat(count);
	for (unsigned int i=0; i<count; ++i)
		flat[i] = (pos[0][i]*d->dim.ny + pos[1][i])*d->dim.nz + pos[2][i];
	ord.pos = d->Alloc<unsigned int>(count, flat.data());

	ord.ade = d->Alloc<float>(3*(size_t)count);
	ord.lor = d->Alloc<float>(3*(size_t)count);
	ord.c_int = Coefficients(count, c_int);
	ord.c_ext = Coefficients(count, c_ext);
	ord.c_lor = Coefficients(count, c_lor);
	orders.push_back(ord);
}

GPU_Extension* HIP_CreateExt_LorentzMaterial(GPU_Backend_HIP::Impl* d, Engine_Extension* eng_ext, Engine* eng)
{
	UNUSED(eng);
	if (!dynamic_cast<Engine_Ext_LorentzMaterial*>(eng_ext))
		return NULL;
	Operator_Ext_LorentzMaterial* op_ext = dynamic_cast<Operator_Ext_LorentzMaterial*>(eng_ext->GetOperatorExtension());
	if (!op_ext)
		return NULL;
	return new HIP_Ext_LorentzMaterial(d, op_ext);
}
