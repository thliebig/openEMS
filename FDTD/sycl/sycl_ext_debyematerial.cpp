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
#include "FDTD/extensions/operator_ext_debyematerial.h"
#include "FDTD/extensions/engine_ext_debyematerial.h"

// Debye materials, see Engine_Ext_DebyeMaterial. One thread per mesh position, all
// poles in the same thread, since the implicit solve couples them. The pole states
// and coefficients are stored as [pole][direction][position], the rest as
// [direction][position].
struct DebyeParam { size_t count; unsigned int poles; size_t sn; };

// before the main update, the only place V^n exists: the part of the new pole
// states V^n fixes, and its sum
static inline void debye_pre(const sycl::nd_item<3>& it, const float* field, float* pole, float* pre, const float* relax,
                             const float* drive, const size_t* pos, DebyeParam P)
{
	const size_t i = (size_t)it.get_group(2)*it.get_local_range(2) + it.get_local_id(2);
	if (i>=P.count)
		return;
	for (unsigned int n=0; n<3; ++n)
	{
		const float volt = field[n*P.sn + pos[i]];
		float sum = 0;
		for (unsigned int o=0; o<P.poles; ++o)
		{
			const size_t k = (o*3 + n)*P.count + i;
			const float W = relax[k]*pole[k] + drive[k]*volt;
			sum = sum + (W - pole[k]);
			pole[k] = W;
		}
		pre[n*P.count + i] = sum;
	}
}

// after the main update: solve for V^(n+1), keep the correction for
// debye_apply and finish the pole states
static inline void debye_post(const sycl::nd_item<3>& it, const float* field, float* pole, const float* pre, float* ade,
                              const float* drive, const float* solve, const size_t* pos, DebyeParam P)
{
	const size_t i = (size_t)it.get_group(2)*it.get_local_range(2) + it.get_local_id(2);
	if (i>=P.count)
		return;
	for (unsigned int n=0; n<3; ++n)
	{
		const size_t j = n*P.count + i;
		const float v_raw = field[n*P.sn + pos[i]];
		const float volt = (v_raw - pre[j]) * solve[j];
		ade[j] = v_raw - volt;
		for (unsigned int o=0; o<P.poles; ++o)
		{
			const size_t k = (o*3 + n)*P.count + i;
			pole[k] = pole[k] + drive[k]*volt;
		}
	}
}

// subtract the correction from the field
static inline void debye_apply(const sycl::nd_item<3>& it, float* field, const float* ade, const size_t* pos, DebyeParam P)
{
	const size_t i = (size_t)it.get_group(2)*it.get_local_range(2) + it.get_local_id(2);
	if (i>=P.count)
		return;
	for (unsigned int n=0; n<3; ++n)
	{
		const size_t g = n*P.sn + pos[i];
		field[g] = field[g] - ade[n*P.count + i];
	}
}

class SYCL_Ext_DebyeMaterial : public GPU_Extension
{
public:
	SYCL_Ext_DebyeMaterial(GPU_Backend_SYCL::Impl* impl, Operator_Ext_DebyeMaterial* op_ext);

	virtual void DoPreVoltageUpdates()
	{
		if (m_Param.count)
		{
			const float* f = d->volt;
			float* pole = m_Pole;
			float* pre = m_Pre;
			const float* relax = m_Relax;
			const float* drive = m_Drive;
			const size_t* pos = m_Pos;
			const DebyeParam param = m_Param;
			SYCL_Launch(d, "debye_pre", param.count, 1, 1, [=](sycl::nd_item<3> it)
			{
				debye_pre(it, f, pole, pre, relax, drive, pos, param);
			});
		}
	}
	virtual void DoPostVoltageUpdates()
	{
		if (m_Param.count)
		{
			const float* f = d->volt;
			float* pole = m_Pole;
			const float* pre = m_Pre;
			float* ade = m_ADE;
			const float* drive = m_Drive;
			const float* solve = m_Solve;
			const size_t* pos = m_Pos;
			const DebyeParam param = m_Param;
			SYCL_Launch(d, "debye_post", param.count, 1, 1, [=](sycl::nd_item<3> it)
			{
				debye_post(it, f, pole, pre, ade, drive, solve, pos, param);
			});
		}
	}
	virtual void Apply2Voltages()
	{
		if (m_Param.count)
		{
			float* f = d->volt;
			const float* ade = m_ADE;
			const size_t* pos = m_Pos;
			const DebyeParam param = m_Param;
			SYCL_Launch(d, "debye_apply", param.count, 1, 1, [=](sycl::nd_item<3> it)
			{
				debye_apply(it, f, ade, pos, param);
			});
		}
	}

protected:
	float* Coefficients(int poles, FDTD_FLOAT*** c);

	GPU_Backend_SYCL::Impl* d;
	DebyeParam m_Param;
	size_t* m_Pos;                     //!< flat NIJK index of direction 0
	float *m_Pole, *m_Pre, *m_ADE;     //!< state
	float *m_Relax, *m_Drive, *m_Solve;
};

SYCL_Ext_DebyeMaterial::SYCL_Ext_DebyeMaterial(GPU_Backend_SYCL::Impl* impl, Operator_Ext_DebyeMaterial* op_ext)
{
	d = impl;
	const size_t count = op_ext->m_LM_Count.at(0);
	m_Param.count = count;
	m_Param.poles = op_ext->m_PoleCount;
	m_Param.sn = d->numCells;
	if (count==0)
		return;

	unsigned int** pos = op_ext->m_LM_pos[0];
	std::vector<size_t> flat(count);
	for (size_t i=0; i<count; ++i)
		flat[i] = ((size_t)pos[0][i]*d->dim.ny + pos[1][i])*d->dim.nz + pos[2][i];
	m_Pos = d->Alloc<size_t>(count, flat.data());
	// fused step (see update_fused): Apply2Voltages() changes the voltages of these nodes after the
	// voltage update, the currents next to them are recomputed
	for (unsigned int n=0; n<3; ++n)
		for (size_t i=0; i<count; ++i)
			d->volt_modified.push_back(n*d->numCells + flat[i]);

	m_Pole  = d->Alloc<float>((size_t)m_Param.poles*3*count);
	m_Pre   = d->Alloc<float>(3*(size_t)count);
	m_ADE   = d->Alloc<float>(3*(size_t)count);
	m_Relax = Coefficients(m_Param.poles, op_ext->v_relax_ADE);
	m_Drive = Coefficients(m_Param.poles, op_ext->v_drive_ADE);
	m_Solve = Coefficients(1, &op_ext->v_solve_ADE);
}

float* SYCL_Ext_DebyeMaterial::Coefficients(int poles, FDTD_FLOAT*** c)
{
	const size_t count = m_Param.count;
	std::vector<float> data((size_t)poles*3*count);
	for (int o=0; o<poles; ++o)
		for (int n=0; n<3; ++n)
			for (size_t i=0; i<count; ++i)
				data[((size_t)o*3 + n)*count + i] = c[o][n][i];
	return d->Alloc<float>(data.size(), data.data());
}

GPU_Extension* SYCL_CreateExt_DebyeMaterial(GPU_Backend_SYCL::Impl* d, Engine_Extension* eng_ext, Engine* eng)
{
	(void)eng;
	if (!dynamic_cast<Engine_Ext_DebyeMaterial*>(eng_ext))
		return NULL;
	Operator_Ext_DebyeMaterial* op_ext = dynamic_cast<Operator_Ext_DebyeMaterial*>(eng_ext->GetOperatorExtension());
	if (!op_ext)
		return NULL;
	return new SYCL_Ext_DebyeMaterial(d, op_ext);
}
