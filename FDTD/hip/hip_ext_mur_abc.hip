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
#include "FDTD/extensions/operator_ext_mur_abc.h"
#include "FDTD/extensions/engine_ext_mur_abc.h"

// First order Mur ABC on one boundary plane, see Engine_Ext_Mur_ABC.
// One thread per (i,j) of the plane, i along m_nyP, j along m_nyPP.
struct MurParam { int ny, nyP, nyPP; unsigned int line, line_shift; unsigned int ni, nj; };

__device__ __forceinline__ unsigned int mur_index(const HIP_GridDim& N, const MurParam& P, unsigned int n, unsigned int line, unsigned int i, unsigned int j)
{
	unsigned int pos[3];
	pos[P.ny] = line;
	pos[P.nyP] = i;
	pos[P.nyPP] = j;
	return nijk(N, n, pos[0], pos[1], pos[2]);
}

// thread index: j along x-threads, i along y-threads
#define MUR_THREAD \
	const unsigned int j = blockIdx.x*blockDim.x + threadIdx.x; \
	const unsigned int i = blockIdx.y*blockDim.y + threadIdx.y; \
	if (i>=P.ni || j>=P.nj) \
		return; \
	const unsigned int ij = i*P.nj + j;

__global__ void mur_pre(const float* volt, float* v_nyP, float* v_nyPP, const float* c_nyP, const float* c_nyPP, HIP_GridDim N, MurParam P)
{
	MUR_THREAD
	v_nyP[ij]  = volt[mur_index(N, P, P.nyP,  P.line_shift, i, j)] - c_nyP[ij]  * volt[mur_index(N, P, P.nyP,  P.line, i, j)];
	v_nyPP[ij] = volt[mur_index(N, P, P.nyPP, P.line_shift, i, j)] - c_nyPP[ij] * volt[mur_index(N, P, P.nyPP, P.line, i, j)];
}

__global__ void mur_post(const float* volt, float* v_nyP, float* v_nyPP, const float* c_nyP, const float* c_nyPP, HIP_GridDim N, MurParam P)
{
	MUR_THREAD
	v_nyP[ij]  = v_nyP[ij]  + c_nyP[ij]  * volt[mur_index(N, P, P.nyP,  P.line_shift, i, j)];
	v_nyPP[ij] = v_nyPP[ij] + c_nyPP[ij] * volt[mur_index(N, P, P.nyPP, P.line_shift, i, j)];
}

__global__ void mur_apply(float* volt, const float* v_nyP, const float* v_nyPP, HIP_GridDim N, MurParam P)
{
	MUR_THREAD
	volt[mur_index(N, P, P.nyP,  P.line, i, j)] = v_nyP[ij];
	volt[mur_index(N, P, P.nyPP, P.line, i, j)] = v_nyPP[ij];
}

class HIP_Ext_Mur_ABC : public GPU_Extension
{
public:
	HIP_Ext_Mur_ABC(GPU_Backend_HIP::Impl* impl, Operator_Ext_Mur_ABC* op_ext, Engine_Ext_Mur_ABC* eng_ext, Engine* eng);

	virtual void DoPreVoltageUpdates()
	{
		if (IsActive())
			HIP_Launch(d, "mur_pre", mur_pre, m_Param.nj, m_Param.ni, 1, (const float*)d->volt, m_Volt_nyP, m_Volt_nyPP, (const float*)m_Coeff_nyP, (const float*)m_Coeff_nyPP, d->dim, m_Param);
	}
	virtual void DoPostVoltageUpdates()
	{
		if (IsActive())
			HIP_Launch(d, "mur_post", mur_post, m_Param.nj, m_Param.ni, 1, (const float*)d->volt, m_Volt_nyP, m_Volt_nyPP, (const float*)m_Coeff_nyP, (const float*)m_Coeff_nyPP, d->dim, m_Param);
	}
	virtual void Apply2Voltages()
	{
		if (IsActive())
			HIP_Launch(d, "mur_apply", mur_apply, m_Param.nj, m_Param.ni, 1, d->volt, (const float*)m_Volt_nyP, (const float*)m_Volt_nyPP, d->dim, m_Param);
	}

protected:
	//! the ABC is off until an excitation on its plane is done
	bool IsActive() {return m_Eng->GetNumberOfTimesteps()>=m_StartTS;}

	GPU_Backend_HIP::Impl* d;
	Engine* m_Eng;
	unsigned int m_StartTS;
	MurParam m_Param;

	float *m_Volt_nyP, *m_Volt_nyPP;
	float *m_Coeff_nyP, *m_Coeff_nyPP;
};

HIP_Ext_Mur_ABC::HIP_Ext_Mur_ABC(GPU_Backend_HIP::Impl* impl, Operator_Ext_Mur_ABC* op_ext, Engine_Ext_Mur_ABC* eng_ext, Engine* eng)
{
	d = impl;
	m_Eng = eng;
	m_StartTS = eng_ext->GetStartTimestep();
	m_Param.ny = op_ext->m_ny;
	m_Param.nyP = op_ext->m_nyP;
	m_Param.nyPP = op_ext->m_nyPP;
	m_Param.line = op_ext->m_LineNr;
	m_Param.line_shift = op_ext->m_LineNr_Shift;
	m_Param.ni = op_ext->m_numLines[0];
	m_Param.nj = op_ext->m_numLines[1];

	const size_t count = (size_t)m_Param.ni*m_Param.nj;
	m_Volt_nyP  = d->Alloc<float>(count);
	m_Volt_nyPP = d->Alloc<float>(count);
	m_Coeff_nyP  = d->Alloc<float>(count, op_ext->m_Mur_Coeff_nyP.data());
	m_Coeff_nyPP = d->Alloc<float>(count, op_ext->m_Mur_Coeff_nyPP.data());
}

GPU_Extension* HIP_CreateExt_Mur_ABC(GPU_Backend_HIP::Impl* d, Engine_Extension* eng_ext, Engine* eng)
{
	Engine_Ext_Mur_ABC* mur_ext = dynamic_cast<Engine_Ext_Mur_ABC*>(eng_ext);
	if (!mur_ext)
		return NULL;
	Operator_Ext_Mur_ABC* op_ext = dynamic_cast<Operator_Ext_Mur_ABC*>(eng_ext->GetOperatorExtension());
	if (!op_ext)
		return NULL;
	return new HIP_Ext_Mur_ABC(d, op_ext, mur_ext, eng);
}
