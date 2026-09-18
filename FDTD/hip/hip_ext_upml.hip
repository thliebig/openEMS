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

#include "hip_internal.h"
#include "FDTD/engine.h"
#include "FDTD/extensions/operator_ext_upml.h"
#include "FDTD/extensions/engine_ext_upml.h"

// Uniaxial PML, see Engine_Ext_UPML. One thread per cell of the PML region.
// The flux and coefficient arrays use the local region size.
struct UPMLRegion { unsigned int sx, sy, sz; unsigned int lx, ly, lz; };

__device__ __forceinline__ unsigned int upml_local(const UPMLRegion& R, unsigned int n, unsigned int x, unsigned int y, unsigned int z)
{
	return ((n*R.lx + x)*R.ly + y)*R.lz + z;
}

// before the main update: store the scaled old field in the flux, hand the old flux to the engine
__global__ void upml_pre(float* field, float* flux, const float* c_old, const float* c_fo, HIP_GridDim N, UPMLRegion R)
{
	const unsigned int z = blockIdx.x*blockDim.x + threadIdx.x;
	const unsigned int y = blockIdx.y*blockDim.y + threadIdx.y;
	const unsigned int x = blockIdx.z;
	if (z>=R.lz || y>=R.ly || x>=R.lx)
		return;
	for (unsigned int n=0; n<3; ++n)
	{
		const unsigned int l = upml_local(R, n, x, y, z);
		const unsigned int g = nijk(N, n, R.sx + x, R.sy + y, R.sz + z);
		const float f_help = c_old[l]*field[g] - c_fo[l]*flux[l];
		field[g] = flux[l];
		flux[l] = f_help;
	}
}

// after the main update: keep the updated flux, combine it with the stored field
__global__ void upml_post(float* field, float* flux, const float* c_fn, HIP_GridDim N, UPMLRegion R)
{
	const unsigned int z = blockIdx.x*blockDim.x + threadIdx.x;
	const unsigned int y = blockIdx.y*blockDim.y + threadIdx.y;
	const unsigned int x = blockIdx.z;
	if (z>=R.lz || y>=R.ly || x>=R.lx)
		return;
	for (unsigned int n=0; n<3; ++n)
	{
		const unsigned int l = upml_local(R, n, x, y, z);
		const unsigned int g = nijk(N, n, R.sx + x, R.sy + y, R.sz + z);
		const float f_help = flux[l];
		flux[l] = field[g];
		field[g] = f_help + c_fn[l]*flux[l];
	}
}

class HIP_Ext_UPML : public GPU_Extension
{
public:
	HIP_Ext_UPML(GPU_Backend_HIP::Impl* impl, Operator_Ext_UPML* op_ext);

	virtual void DoPreVoltageUpdates()  {Pre(d->volt, m_VoltFlux, m_VV, m_VVFO);}
	virtual void DoPostVoltageUpdates() {Post(d->volt, m_VoltFlux, m_VVFN);}
	virtual void DoPreCurrentUpdates()  {Pre(d->curr, m_CurrFlux, m_II, m_IIFO);}
	virtual void DoPostCurrentUpdates() {Post(d->curr, m_CurrFlux, m_IIFN);}

protected:
	void Pre(float* field, float* flux, const float* c_old, const float* c_fo)
	{
		HIP_Launch(d, "upml_pre", upml_pre, m_Region.lz, m_Region.ly, m_Region.lx, field, flux, c_old, c_fo, d->dim, m_Region);
	}
	void Post(float* field, float* flux, const float* c_fn)
	{
		HIP_Launch(d, "upml_post", upml_post, m_Region.lz, m_Region.ly, m_Region.lx, field, flux, c_fn, d->dim, m_Region);
	}
	float* Upload(const ArrayLib::ArrayNIJK<FDTD_FLOAT>& arr) {return d->Alloc<float>(arr.size(), arr.data());}

	GPU_Backend_HIP::Impl* d;
	UPMLRegion m_Region;

	float *m_VoltFlux, *m_CurrFlux;
	float *m_VV, *m_VVFO, *m_VVFN;
	float *m_II, *m_IIFO, *m_IIFN;
};

HIP_Ext_UPML::HIP_Ext_UPML(GPU_Backend_HIP::Impl* impl, Operator_Ext_UPML* op_ext)
{
	d = impl;
	m_Region.sx = op_ext->m_StartPos[0];
	m_Region.sy = op_ext->m_StartPos[1];
	m_Region.sz = op_ext->m_StartPos[2];
	m_Region.lx = op_ext->m_numLines[0];
	m_Region.ly = op_ext->m_numLines[1];
	m_Region.lz = op_ext->m_numLines[2];

	const size_t count = 3*(size_t)m_Region.lx*m_Region.ly*m_Region.lz;
	m_VoltFlux = d->Alloc<float>(count);
	m_CurrFlux = d->Alloc<float>(count);
	m_VV   = Upload(op_ext->vv);
	m_VVFO = Upload(op_ext->vvfo);
	m_VVFN = Upload(op_ext->vvfn);
	m_II   = Upload(op_ext->ii);
	m_IIFO = Upload(op_ext->iifo);
	m_IIFN = Upload(op_ext->iifn);
}

GPU_Extension* HIP_CreateExt_UPML(GPU_Backend_HIP::Impl* d, Engine_Extension* eng_ext, Engine* eng)
{
	UNUSED(eng);
	if (!dynamic_cast<Engine_Ext_UPML*>(eng_ext))
		return NULL;
	Operator_Ext_UPML* op_ext = dynamic_cast<Operator_Ext_UPML*>(eng_ext->GetOperatorExtension());
	if (!op_ext)
		return NULL;
	return new HIP_Ext_UPML(d, op_ext);
}
