/*
*	Copyright (C) 2026 Sean Mollet (sean@malmoset.com)
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

#include "sycl_internal.h"
#include "FDTD/engine.h"
#include "FDTD/extensions/operator_ext_lorentzmaterial.h"
#include "FDTD/extensions/engine_ext_lorentzmaterial.h"

// Drude/Lorentz materials and conducting sheets, see Engine_Ext_LorentzMaterial
// and Engine_Ext_Dispersive. The ADE state and coefficients are stored per order as
// [direction][position]. The orders of one side (voltages or currents) on the same
// positions (e.g. the two orders of a conducting sheet) run in one launch, one thread
// per position: the field is read once for all of them, and the corrections are
// subtracted one after the other in the order of the orders, as separate launches
// would do. A direction no order of a launch has a coefficient for is skipped: its
// state stays zero, so its correction would be a subtraction of zero (e.g. the
// normal direction of a conducting sheet).
#define ADE_MAX_ORDERS 4

struct ADEOrder
{
	float *ade, *lor;
	const float *c_int, *c_ext, *c_lor;
	unsigned int lorentz;
	unsigned int dirs;     //!< bit n: direction n has a nonzero coefficient
};

struct ADEGroup
{
	size_t count;
	size_t sn;
	const size_t* pos;     //!< flat NIJK index of direction 0, in memory order
	unsigned int orders;
	unsigned int dirs;     //!< union of the directions of the orders
	ADEOrder o[ADE_MAX_ORDERS];
};

// ADE update before the main update, reads the field
static inline void lorentz_pre(const sycl::nd_item<3>& it, const float* field, const ADEGroup& G)
{
	const size_t i = (size_t)it.get_group(2)*it.get_local_range(2) + it.get_local_id(2);
	if (i>=G.count)
		return;
	const size_t p = G.pos[i];
	for (unsigned int n=0; n<3; ++n)
	{
		if (!((G.dirs>>n) & 1))
			continue;
		const size_t k = n*G.count + i;
		const float f = field[n*G.sn + p];
		#pragma unroll
		for (unsigned int q=0; q<ADE_MAX_ORDERS; ++q)   // fixed indices: the arguments stay in registers
		{
			const ADEOrder& O = G.o[q];
			if ((q>=G.orders) || !((O.dirs>>n) & 1))
				continue;
			if (O.lorentz)
			{
				O.lor[k] = O.lor[k] + O.c_lor[k]*O.ade[k];
				O.ade[k] = O.ade[k] * O.c_int[k];
				O.ade[k] = O.ade[k] + O.c_ext[k]*(f - O.lor[k]);
			}
			else
			{
				O.ade[k] = O.ade[k] * O.c_int[k];
				O.ade[k] = O.ade[k] + O.c_ext[k]*f;
			}
		}
	}
}

// subtract the ADE currents of the orders from the field
static inline void dispersive_apply(const sycl::nd_item<3>& it, float* field, const ADEGroup& G)
{
	const size_t i = (size_t)it.get_group(2)*it.get_local_range(2) + it.get_local_id(2);
	if (i>=G.count)
		return;
	const size_t p = G.pos[i];
	for (unsigned int n=0; n<3; ++n)
	{
		if (!((G.dirs>>n) & 1))
			continue;
		const size_t g = n*G.sn + p;
		float f = field[g];
		#pragma unroll
		for (unsigned int q=0; q<ADE_MAX_ORDERS; ++q)
			if ((q<G.orders) && ((G.o[q].dirs>>n) & 1))
				f = f - G.o[q].ade[n*G.count + i];
		field[g] = f;
	}
}

class SYCL_Ext_LorentzMaterial : public GPU_Extension
{
public:
	SYCL_Ext_LorentzMaterial(GPU_Backend_SYCL::Impl* impl, Operator_Ext_LorentzMaterial* op_ext);

	virtual void DoPreVoltageUpdates() {for (size_t g=0; g<m_Volt.size(); ++g) Pre(m_Volt[g], d->volt);}
	virtual void Apply2Voltages()      {for (size_t g=0; g<m_Volt.size(); ++g) Apply(m_Volt[g], d->volt);}
	virtual void DoPreCurrentUpdates() {for (size_t g=0; g<m_Curr.size(); ++g) Pre(m_Curr[g], d->curr);}
	virtual void Apply2Current()       {for (size_t g=0; g<m_Curr.size(); ++g) Apply(m_Curr[g], d->curr);}

protected:
	void Setup(std::vector<ADEGroup>& groups, std::vector<std::vector<size_t>>& group_pos, size_t count, unsigned int** pos,
	           bool lorentz, FDTD_FLOAT** c_int, FDTD_FLOAT** c_ext, FDTD_FLOAT** c_lor);
	float* Coefficients(size_t count, FDTD_FLOAT** c);
	void Pre(const ADEGroup& group, float* field)
	{
		if (group.dirs==0)
			return;
		const float* f = field;
		const ADEGroup G = group;
		SYCL_Launch(d, "lorentz_pre", G.count, 1, 1, [=](sycl::nd_item<3> it)
		{
			lorentz_pre(it, f, G);
		});
	}
	void Apply(const ADEGroup& group, float* field)
	{
		if (group.dirs==0)
			return;
		float* f = field;
		const ADEGroup G = group;
		SYCL_Launch(d, "dispersive_apply", G.count, 1, 1, [=](sycl::nd_item<3> it)
		{
			dispersive_apply(it, f, G);
		});
	}

	GPU_Backend_SYCL::Impl* d;
	std::vector<ADEGroup> m_Volt;
	std::vector<ADEGroup> m_Curr;
};

SYCL_Ext_LorentzMaterial::SYCL_Ext_LorentzMaterial(GPU_Backend_SYCL::Impl* impl, Operator_Ext_LorentzMaterial* op_ext)
{
	d = impl;
	std::vector<std::vector<size_t>> volt_pos, curr_pos;
	for (int o=0; o<op_ext->m_Order; ++o)
	{
		const size_t count = op_ext->m_LM_Count.at(o);
		if (op_ext->m_volt_ADE_On[o])
			Setup(m_Volt, volt_pos, count, op_ext->m_LM_pos[o], op_ext->m_volt_Lor_ADE_On[o],
			      op_ext->v_int_ADE[o], op_ext->v_ext_ADE[o], op_ext->m_volt_Lor_ADE_On[o] ? op_ext->v_Lor_ADE[o] : NULL);
		if (op_ext->m_curr_ADE_On[o])
			Setup(m_Curr, curr_pos, count, op_ext->m_LM_pos[o], op_ext->m_curr_Lor_ADE_On[o],
			      op_ext->i_int_ADE[o], op_ext->i_ext_ADE[o], op_ext->m_curr_Lor_ADE_On[o] ? op_ext->i_Lor_ADE[o] : NULL);
	}
	// fused step (see update_fused): Apply2Voltages() changes the voltages of these nodes after the
	// voltage update, the currents next to them are recomputed; Apply2Current() changes their
	// currents after the current update
	for (size_t g=0; g<m_Volt.size(); ++g)
		for (unsigned int n=0; n<3; ++n)
			if ((m_Volt[g].dirs>>n) & 1)
				for (size_t k=0; k<volt_pos[g].size(); ++k)
					d->volt_modified.push_back(n*d->numCells + volt_pos[g][k]);
	for (size_t g=0; g<m_Curr.size(); ++g)
		if (m_Curr[g].dirs)
			for (unsigned int n=0; n<3; ++n)
				if ((m_Curr[g].dirs>>n) & 1)
					for (size_t k=0; k<curr_pos[g].size(); ++k)
						d->curr_modified.push_back(n*d->numCells + curr_pos[g][k]);
}

float* SYCL_Ext_LorentzMaterial::Coefficients(size_t count, FDTD_FLOAT** c)
{
	std::vector<float> data(3*(size_t)count, 0);
	if (c)
		for (int n=0; n<3; ++n)
			for (size_t i=0; i<count; ++i)
				data[n*count + i] = c[n][i];
	return d->Alloc<float>(data.size(), data.data());
}

void SYCL_Ext_LorentzMaterial::Setup(std::vector<ADEGroup>& groups, std::vector<std::vector<size_t>>& group_pos, size_t count, unsigned int** pos,
                                     bool lorentz, FDTD_FLOAT** c_int, FDTD_FLOAT** c_ext, FDTD_FLOAT** c_lor)
{
	if (count==0)
		return;
	std::vector<size_t> flat(count);
	for (size_t i=0; i<count; ++i)
		flat[i] = ((size_t)pos[0][i]*d->dim.ny + pos[1][i])*d->dim.nz + pos[2][i];

	ADEOrder ord;
	ord.lorentz = lorentz;
	// a direction without a coefficient keeps a zero state (ade = ade*c_int + 0, lor = lor + c_lor*0)
	ord.dirs = 0;
	for (int n=0; n<3; ++n)
		for (size_t i=0; c_ext && (i<count); ++i)
			if (c_ext[n][i]!=0)
			{
				ord.dirs |= 1u<<n;
				break;
			}
	ord.ade = d->Alloc<float>(3*(size_t)count);
	ord.lor = d->Alloc<float>(3*(size_t)count);
	ord.c_int = Coefficients(count, c_int);
	ord.c_ext = Coefficients(count, c_ext);
	ord.c_lor = Coefficients(count, c_lor);

	// the same positions as the previous order of this side: one launch for both
	if (!groups.empty() && (group_pos.back()==flat) && (groups.back().orders<ADE_MAX_ORDERS))
	{
		ADEGroup& G = groups.back();
		G.o[G.orders++] = ord;
		G.dirs |= ord.dirs;
		return;
	}
	ADEGroup G;
	G.count = count;
	G.sn = d->numCells;
	G.pos = d->Alloc<size_t>(count, flat.data());
	G.orders = 1;
	G.dirs = ord.dirs;
	G.o[0] = ord;
	groups.push_back(G);
	group_pos.push_back(flat);
}

GPU_Extension* SYCL_CreateExt_LorentzMaterial(GPU_Backend_SYCL::Impl* d, Engine_Extension* eng_ext, Engine* eng)
{
	(void)eng;
	if (!dynamic_cast<Engine_Ext_LorentzMaterial*>(eng_ext))
		return NULL;
	Operator_Ext_LorentzMaterial* op_ext = dynamic_cast<Operator_Ext_LorentzMaterial*>(eng_ext->GetOperatorExtension());
	if (!op_ext)
		return NULL;
	return new SYCL_Ext_LorentzMaterial(d, op_ext);
}
