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
#include "FDTD/extensions/operator_ext_absorbing_bc.h"
#include "FDTD/extensions/engine_ext_absorbing_bc.h"

// Local absorbing sheet (first order Mur, optionally with super-absorption of the
// currents), see Engine_Ext_Absorbing_BC. One thread per (i,j) of the sheet,
// i along m_nyP, j along m_nyPP. The state and coefficient arrays have the full
// sheet size ni x nj, the currents use ni-1 x nj-1 of it.
struct ABCParam { int ny, nyP, nyPP; unsigned int line, line_shift; unsigned int startP, startPP; unsigned int ni, nj, stride; unsigned int swap; };

//! state and coefficients of one field
struct ABCData
{
	float* s_nyP; float* s_nyPP;
	const float* k1_nyP; const float* k1_nyPP;
	const float* k2_nyP; const float* k2_nyPP;
};

inline size_t abc_index(const SYCL_GridDim& N, const ABCParam& P, unsigned int n, unsigned int line, unsigned int i, unsigned int j)
{
	unsigned int pos[3];
	pos[P.ny] = line;
	pos[P.nyP] = P.startP + i;
	pos[P.nyPP] = P.startPP + j;
	return nijk(N, n, pos[0], pos[1], pos[2]);
}

// thread index: j along dimension 2, i along dimension 1, or the other way round (swap: the
// contiguous dimension along the axis with the smaller memory stride); the arrays have a row
// length of stride
#define ABC_THREAD \
	const unsigned int t2 = it.get_group(2)*it.get_local_range(2) + it.get_local_id(2); \
	const unsigned int t1 = it.get_group(1)*it.get_local_range(1) + it.get_local_id(1); \
	const unsigned int j = P.swap ? t1 : t2; \
	const unsigned int i = P.swap ? t2 : t1; \
	if (i>=P.ni || j>=P.nj) \
		return; \
	const unsigned int ij = i*P.stride + j;

// F(i+s,n) - K1*F(i,n)
static inline void abc_pre(const sycl::nd_item<3>& it, const float* field, ABCData D, SYCL_GridDim N, ABCParam P)
{
	ABC_THREAD
	D.s_nyP[ij]  = field[abc_index(N, P, P.nyP,  P.line_shift, i, j)] - D.k1_nyP[ij]  * field[abc_index(N, P, P.nyP,  P.line, i, j)];
	D.s_nyPP[ij] = field[abc_index(N, P, P.nyPP, P.line_shift, i, j)] - D.k1_nyPP[ij] * field[abc_index(N, P, P.nyPP, P.line, i, j)];
}

// + K1*F(i+s,n+1)
static inline void abc_post(const sycl::nd_item<3>& it, const float* field, ABCData D, SYCL_GridDim N, ABCParam P)
{
	ABC_THREAD
	D.s_nyP[ij]  = D.s_nyP[ij]  + D.k1_nyP[ij]  * field[abc_index(N, P, P.nyP,  P.line_shift, i, j)];
	D.s_nyPP[ij] = D.s_nyPP[ij] + D.k1_nyPP[ij] * field[abc_index(N, P, P.nyPP, P.line_shift, i, j)];
}

// voltages: replace by the Mur value
static inline void abc_apply_volt(const sycl::nd_item<3>& it, float* field, ABCData D, SYCL_GridDim N, ABCParam P)
{
	ABC_THREAD
	field[abc_index(N, P, P.nyP,  P.line, i, j)] = D.s_nyP[ij];
	field[abc_index(N, P, P.nyPP, P.line, i, j)] = D.s_nyPP[ij];
}

// currents, super-absorption: (Hsa*K2 + Hc)/(K2 + 1)
static inline void abc_apply_curr(const sycl::nd_item<3>& it, float* field, ABCData D, SYCL_GridDim N, ABCParam P)
{
	ABC_THREAD
	const size_t gP  = abc_index(N, P, P.nyP,  P.line, i, j);
	const size_t gPP = abc_index(N, P, P.nyPP, P.line, i, j);
	field[gP]  = (D.s_nyP[ij] *D.k2_nyP[ij]  + field[gP]) /(D.k2_nyP[ij]  + 1.0f);
	field[gPP] = (D.s_nyPP[ij]*D.k2_nyPP[ij] + field[gPP])/(D.k2_nyPP[ij] + 1.0f);
}

class SYCL_Ext_Absorbing_BC : public GPU_Extension
{
public:
	SYCL_Ext_Absorbing_BC(GPU_Backend_SYCL::Impl* impl, Operator_Ext_Absorbing_BC* op_ext);

	virtual void DoPreVoltageUpdates()  {LaunchRead("abc_pre", 0, d->volt, m_VoltData, m_VoltParam);}
	virtual void DoPostVoltageUpdates() {LaunchRead("abc_post", 1, d->volt, m_VoltData, m_VoltParam);}
	virtual void Apply2Voltages()       {LaunchApply("abc_apply_volt", d->volt, m_VoltData, m_VoltParam, false);}

	// super-absorption of the currents, on the dual grid one line less
	virtual void DoPreCurrentUpdates()  {if (m_SuperAbsorption) LaunchRead("abc_pre", 0, d->curr, m_CurrData, m_CurrParam);}
	virtual void DoPostCurrentUpdates() {if (m_SuperAbsorption) LaunchRead("abc_post", 1, d->curr, m_CurrData, m_CurrParam);}
	virtual void Apply2Current()        {if (m_SuperAbsorption) LaunchApply("abc_apply_curr", d->curr, m_CurrData, m_CurrParam, true);}

protected:
	//! abc_pre (post==0) or abc_post (post==1): the field is only read
	void LaunchRead(const char* name, int post, const float* field_in, const ABCData& data, const ABCParam& param)
	{
		const float* field = field_in;
		const ABCData D = data;
		const SYCL_GridDim N = d->dim;
		const ABCParam P = param;
		if (post)
			SYCL_Launch(d, name, P.swap ? P.ni : P.nj, P.swap ? P.nj : P.ni, 1, [=](sycl::nd_item<3> it) {abc_post(it, field, D, N, P);});
		else
			SYCL_Launch(d, name, P.swap ? P.ni : P.nj, P.swap ? P.nj : P.ni, 1, [=](sycl::nd_item<3> it) {abc_pre(it, field, D, N, P);});
	}

	//! abc_apply_volt or abc_apply_curr (curr==true)
	void LaunchApply(const char* name, float* field_in, const ABCData& data, const ABCParam& param, bool curr)
	{
		float* field = field_in;
		const ABCData D = data;
		const SYCL_GridDim N = d->dim;
		const ABCParam P = param;
		if (curr)
			SYCL_Launch(d, name, P.swap ? P.ni : P.nj, P.swap ? P.nj : P.ni, 1, [=](sycl::nd_item<3> it) {abc_apply_curr(it, field, D, N, P);});
		else
			SYCL_Launch(d, name, P.swap ? P.ni : P.nj, P.swap ? P.nj : P.ni, 1, [=](sycl::nd_item<3> it) {abc_apply_volt(it, field, D, N, P);});
	}

	GPU_Backend_SYCL::Impl* d;
	bool m_SuperAbsorption;
	ABCParam m_VoltParam, m_CurrParam;
	ABCData m_VoltData, m_CurrData;
};

SYCL_Ext_Absorbing_BC::SYCL_Ext_Absorbing_BC(GPU_Backend_SYCL::Impl* impl, Operator_Ext_Absorbing_BC* op_ext)
{
	d = impl;
	m_SuperAbsorption = (op_ext->m_ABCtype==Operator_Ext_Absorbing_BC::MUR_1ST_SA);
	const unsigned int ni = op_ext->m_numLines[0];
	const unsigned int nj = op_ext->m_numLines[1];

	// positions as in Engine_Ext_Absorbing_BC::Engine_Ext_Absorbing_BC()
	const int ny = op_ext->m_ny;
	const unsigned int start = op_ext->m_sheetX0[ny];
	const bool positive = op_ext->m_normalSignPositive;
	ABCParam base;
	base.ny = ny;
	base.nyP = op_ext->m_nyP;
	base.nyPP = op_ext->m_nyPP;
	base.startP = op_ext->m_sheetX0[op_ext->m_nyP];
	base.startPP = op_ext->m_sheetX0[op_ext->m_nyPP];
	base.stride = nj;
	base.swap = base.nyP > base.nyPP;   // the higher axis has the smaller stride

	m_VoltParam = base;
	m_VoltParam.line = start;
	m_VoltParam.line_shift = start + (positive ? 1 : -1);
	m_VoltParam.ni = ni;
	m_VoltParam.nj = nj;

	m_CurrParam = base;
	m_CurrParam.line = start + (positive ? 0 : -1);
	m_CurrParam.line_shift = start + (positive ? 1 : -2);
	m_CurrParam.ni = ni>0 ? ni-1 : 0;
	m_CurrParam.nj = nj>0 ? nj-1 : 0;

	const size_t count = (size_t)ni*nj;
	const float* k1_nyP  = d->Alloc<float>(count, op_ext->m_K1_nyP.data());
	const float* k1_nyPP = d->Alloc<float>(count, op_ext->m_K1_nyPP.data());
	const float* k2_nyP  = d->Alloc<float>(count, op_ext->m_K2_nyP.valid() ? op_ext->m_K2_nyP.data() : NULL);
	const float* k2_nyPP = d->Alloc<float>(count, op_ext->m_K2_nyPP.valid() ? op_ext->m_K2_nyPP.data() : NULL);

	m_VoltData.s_nyP  = d->Alloc<float>(count);
	m_VoltData.s_nyPP = d->Alloc<float>(count);
	m_CurrData.s_nyP  = d->Alloc<float>(count);
	m_CurrData.s_nyPP = d->Alloc<float>(count);
	m_VoltData.k1_nyP = m_CurrData.k1_nyP = k1_nyP;
	m_VoltData.k1_nyPP = m_CurrData.k1_nyPP = k1_nyPP;
	m_VoltData.k2_nyP = m_CurrData.k2_nyP = k2_nyP;
	m_VoltData.k2_nyPP = m_CurrData.k2_nyPP = k2_nyPP;

	// fused step (see update_fused): Apply2Voltages() replaces the voltages of the sheet after the
	// voltage update, the currents next to them are recomputed; Apply2Current() changes the
	// currents of its plane after the current update
	for (unsigned int i=0; i<m_VoltParam.ni; ++i)
		for (unsigned int j=0; j<m_VoltParam.nj; ++j)
		{
			d->volt_modified.push_back(abc_index(d->dim, m_VoltParam, m_VoltParam.nyP, m_VoltParam.line, i, j));
			d->volt_modified.push_back(abc_index(d->dim, m_VoltParam, m_VoltParam.nyPP, m_VoltParam.line, i, j));
		}
	if (m_SuperAbsorption)
		for (unsigned int i=0; i<m_CurrParam.ni; ++i)
			for (unsigned int j=0; j<m_CurrParam.nj; ++j)
				d->curr_modified.push_back(abc_index(d->dim, m_CurrParam, 0, m_CurrParam.line, i, j));
}

GPU_Extension* SYCL_CreateExt_Absorbing_BC(GPU_Backend_SYCL::Impl* d, Engine_Extension* eng_ext, Engine* eng)
{
	UNUSED(eng);
	if (!dynamic_cast<Engine_Ext_Absorbing_BC*>(eng_ext))
		return NULL;
	Operator_Ext_Absorbing_BC* op_ext = dynamic_cast<Operator_Ext_Absorbing_BC*>(eng_ext->GetOperatorExtension());
	if (!op_ext)
		return NULL;
	return new SYCL_Ext_Absorbing_BC(d, op_ext);
}
