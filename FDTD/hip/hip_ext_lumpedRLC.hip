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
#include <set>

#include "hip_internal.h"
#include "FDTD/engine.h"
#include "FDTD/extensions/operator_ext_lumpedRLC.h"
#include "FDTD/extensions/engine_ext_lumpedRLC.h"

// Series/parallel lumped RLC elements, see Engine_Ext_LumpedRLC.
// The histories Vd[n], Vd[n-1], Vd[n-2] and J[n], J[n-1], J[n-2] are three slots
// of count values each; the host rotates the slot indices instead of pointers.
struct RLCParam { unsigned int count; unsigned int vd0, vd1, vd2; unsigned int j0, j1, j2; };

//! state and coefficients of all elements
struct RLCData
{
	float* il; float* vd; float* j;
	const unsigned int* index;
	const float* ilv; const float* i2v;
	const float* vv2; const float* vj1; const float* vj2; const float* vvd;
	const float* ib0; const float* b1;  const float* b2;
};

// parallel inductor current
__global__ void rlc_pre(RLCData D, RLCParam P)
{
	const unsigned int p = blockIdx.x*blockDim.x + threadIdx.x;
	if (p<P.count)
		D.il[p] = D.il[p] + (D.i2v[p]*D.ilv[p])*D.vd[P.vd1*P.count + p];
}

__device__ __forceinline__ float rlc_vd(unsigned int p, float v, const RLCData& D, const RLCParam& P)
{
	return D.vvd[p]*(v - D.il[p] + D.vv2[p]*D.vd[P.vd2*P.count + p] + D.vj1[p]*D.j[P.j1*P.count + p] + D.vj2[p]*D.j[P.j2*P.count + p]);
}

__device__ __forceinline__ float rlc_j(unsigned int p, const RLCData& D, const RLCParam& P)
{
	return D.ib0[p]*(D.vd[P.vd0*P.count + p] - D.vd[P.vd2*P.count + p])
	       - (D.b1[p]*D.ib0[p])*D.j[P.j1*P.count + p]
	       - (D.b2[p]*D.ib0[p])*D.j[P.j2*P.count + p];
}

// node voltage update, one thread per element (all elements on different edges)
__global__ void rlc_apply(float* volt, RLCData D, RLCParam P)
{
	const unsigned int p = blockIdx.x*blockDim.x + threadIdx.x;
	if (p>=P.count)
		return;
	D.vd[P.vd0*P.count + p] = rlc_vd(p, volt[D.index[p]], D, P);
	D.j[P.j0*P.count + p] = rlc_j(p, D, P);
	volt[D.index[p]] = D.vd[P.vd0*P.count + p];
}

// same in the three passes of the host (read all, update all, write all), for elements sharing an edge
__global__ void rlc_apply_serial(float* volt, RLCData D, RLCParam P)
{
	if (blockIdx.x*blockDim.x + threadIdx.x != 0)
		return;
	for (unsigned int p=0; p<P.count; ++p)
		D.vd[P.vd0*P.count + p] = volt[D.index[p]];
	for (unsigned int p=0; p<P.count; ++p)
	{
		D.vd[P.vd0*P.count + p] = rlc_vd(p, D.vd[P.vd0*P.count + p], D, P);
		D.j[P.j0*P.count + p] = rlc_j(p, D, P);
	}
	for (unsigned int p=0; p<P.count; ++p)
		volt[D.index[p]] = D.vd[P.vd0*P.count + p];
}

class HIP_Ext_LumpedRLC : public GPU_Extension
{
public:
	HIP_Ext_LumpedRLC(GPU_Backend_HIP::Impl* impl, Operator_Ext_LumpedRLC* op_ext);

	virtual void DoPreVoltageUpdates();
	virtual void Apply2Voltages();

protected:
	GPU_Backend_HIP::Impl* d;
	bool m_Serial;   //!< elements share an edge
	RLCParam m_Param;
	RLCData m_Data;
};

HIP_Ext_LumpedRLC::HIP_Ext_LumpedRLC(GPU_Backend_HIP::Impl* impl, Operator_Ext_LumpedRLC* op_ext)
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

	std::vector<unsigned int> index(count);
	std::set<unsigned int> edges;
	for (unsigned int p=0; p<count; ++p)
	{
		index[p] = ((op_ext->v_RLC_dir[p]*d->dim.nx + op_ext->v_RLC_pos[0][p])*d->dim.ny + op_ext->v_RLC_pos[1][p])*d->dim.nz + op_ext->v_RLC_pos[2][p];
		if (!edges.insert(index[p]).second)
			m_Serial = true;
	}
	// fused step (see update_fused): Apply2Voltages() changes these voltages after the voltage
	// update, the currents next to them are recomputed
	d->volt_modified.insert(d->volt_modified.end(), index.begin(), index.end());
	m_Data.index = d->Alloc<unsigned int>(count, index.data());
	m_Data.il = d->Alloc<float>(count);
	m_Data.vd = d->Alloc<float>(3*(size_t)count);
	m_Data.j  = d->Alloc<float>(3*(size_t)count);
	m_Data.ilv = d->Alloc<float>(count, op_ext->v_RLC_ilv);
	m_Data.i2v = d->Alloc<float>(count, op_ext->v_RLC_i2v);
	m_Data.vv2 = d->Alloc<float>(count, op_ext->v_RLC_vv2);
	m_Data.vj1 = d->Alloc<float>(count, op_ext->v_RLC_vj1);
	m_Data.vj2 = d->Alloc<float>(count, op_ext->v_RLC_vj2);
	m_Data.vvd = d->Alloc<float>(count, op_ext->v_RLC_vvd);
	m_Data.ib0 = d->Alloc<float>(count, op_ext->v_RLC_ib0);
	m_Data.b1  = d->Alloc<float>(count, op_ext->v_RLC_b1);
	m_Data.b2  = d->Alloc<float>(count, op_ext->v_RLC_b2);
}

void HIP_Ext_LumpedRLC::DoPreVoltageUpdates()
{
	if (m_Param.count==0)
		return;
	// rotate the Vd history: [0,1,2] <- [2,0,1]
	const unsigned int vd2 = m_Param.vd2;
	m_Param.vd2 = m_Param.vd1;
	m_Param.vd1 = m_Param.vd0;
	m_Param.vd0 = vd2;
	HIP_Launch(d, "rlc_pre", rlc_pre, m_Param.count, 1, 1, m_Data, m_Param);
}

void HIP_Ext_LumpedRLC::Apply2Voltages()
{
	if (m_Param.count==0)
		return;
	// rotate the J history: [0,1,2] <- [2,0,1]
	const unsigned int j2 = m_Param.j2;
	m_Param.j2 = m_Param.j1;
	m_Param.j1 = m_Param.j0;
	m_Param.j0 = j2;
	if (m_Serial)
		HIP_Launch(d, "rlc_apply_serial", rlc_apply_serial, 1, 1, 1, d->volt, m_Data, m_Param);
	else
		HIP_Launch(d, "rlc_apply", rlc_apply, m_Param.count, 1, 1, d->volt, m_Data, m_Param);
}

GPU_Extension* HIP_CreateExt_LumpedRLC(GPU_Backend_HIP::Impl* d, Engine_Extension* eng_ext, Engine* eng)
{
	UNUSED(eng);
	if (!dynamic_cast<Engine_Ext_LumpedRLC*>(eng_ext))
		return NULL;
	Operator_Ext_LumpedRLC* op_ext = dynamic_cast<Operator_Ext_LumpedRLC*>(eng_ext->GetOperatorExtension());
	if (!op_ext)
		return NULL;
	return new HIP_Ext_LumpedRLC(d, op_ext);
}
