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
#include "FDTD/extensions/operator_ext_upml.h"
#include "FDTD/extensions/engine_ext_upml.h"

// Uniaxial PML, see Engine_Ext_UPML. One thread per cell of the PML region,
// gid = local (z, y, x). The flux and coefficient arrays use the local region size.
static const char* UPML_SOURCE = R"MSL(
struct UPMLParam { uint sx, sy, sz; uint lx, ly, lz; };

inline uint upml_local(constant UPMLParam& P, uint n, uint3 gid)
{
	return ((n*P.lx + gid.z)*P.ly + gid.y)*P.lz + gid.x;
}

inline uint upml_global(constant GridDim& N, constant UPMLParam& P, uint n, uint3 gid)
{
	return nijk(N, n, P.sx + gid.z, P.sy + gid.y, P.sz + gid.x);
}

// before the main update: store the scaled old field in the flux, hand the old flux to the engine
kernel void upml_pre(device float* field       [[buffer(0)]],
                     device float* flux        [[buffer(1)]],
                     const device float* c_old [[buffer(2)]],
                     const device float* c_fo  [[buffer(3)]],
                     constant GridDim& N       [[buffer(4)]],
                     constant UPMLParam& P     [[buffer(5)]],
                     uint3 gid [[thread_position_in_grid]])
{
	for (uint n=0; n<3; ++n)
	{
		const uint l = upml_local(P, n, gid);
		const uint g = upml_global(N, P, n, gid);
		const float f_help = c_old[l]*field[g] - c_fo[l]*flux[l];
		field[g] = flux[l];
		flux[l] = f_help;
	}
}

// after the main update: keep the updated flux, combine it with the stored field
kernel void upml_post(device float* field       [[buffer(0)]],
                      device float* flux        [[buffer(1)]],
                      const device float* c_fn  [[buffer(2)]],
                      constant GridDim& N       [[buffer(4)]],
                      constant UPMLParam& P     [[buffer(5)]],
                      uint3 gid [[thread_position_in_grid]])
{
	for (uint n=0; n<3; ++n)
	{
		const uint l = upml_local(P, n, gid);
		const uint g = upml_global(N, P, n, gid);
		const float f_help = flux[l];
		flux[l] = field[g];
		field[g] = f_help + c_fn[l]*flux[l];
	}
}
)MSL";

class Metal_Ext_UPML : public GPU_Extension
{
public:
	Metal_Ext_UPML(GPU_Backend_Metal::Impl* impl, Operator_Ext_UPML* op_ext);

	virtual void DoPreVoltageUpdates()  {Pre(d->volt, m_VoltFlux, m_VV, m_VVFO);}
	virtual void DoPostVoltageUpdates() {Post(d->volt, m_VoltFlux, m_VVFN);}
	virtual void DoPreCurrentUpdates()  {Pre(d->curr, m_CurrFlux, m_II, m_IIFO);}
	virtual void DoPostCurrentUpdates() {Post(d->curr, m_CurrFlux, m_IIFN);}

protected:
	void Pre(id<MTLBuffer> field, id<MTLBuffer> flux, id<MTLBuffer> c_old, id<MTLBuffer> c_fo);
	void Post(id<MTLBuffer> field, id<MTLBuffer> flux, id<MTLBuffer> c_fn);
	void SetRegion(id<MTLComputeCommandEncoder> enc);

	GPU_Backend_Metal::Impl* d;
	struct {uint32_t sx, sy, sz, lx, ly, lz;} m_Region;

	id<MTLBuffer> m_VoltFlux, m_CurrFlux;
	id<MTLBuffer> m_VV, m_VVFO, m_VVFN;
	id<MTLBuffer> m_II, m_IIFO, m_IIFN;
};

static id<MTLBuffer> UploadArray(GPU_Backend_Metal::Impl* d, const ArrayLib::ArrayNIJK<FDTD_FLOAT>& arr)
{
	return d->NewBuffer(arr.size()*sizeof(FDTD_FLOAT), arr.data());
}

Metal_Ext_UPML::Metal_Ext_UPML(GPU_Backend_Metal::Impl* impl, Operator_Ext_UPML* op_ext)
{
	d = impl;
	m_Region.sx = op_ext->m_StartPos[0];
	m_Region.sy = op_ext->m_StartPos[1];
	m_Region.sz = op_ext->m_StartPos[2];
	m_Region.lx = op_ext->m_numLines[0];
	m_Region.ly = op_ext->m_numLines[1];
	m_Region.lz = op_ext->m_numLines[2];

	const size_t bytes = 3*(size_t)m_Region.lx*m_Region.ly*m_Region.lz*sizeof(FDTD_FLOAT);
	m_VoltFlux = d->NewBuffer(bytes);
	m_CurrFlux = d->NewBuffer(bytes);
	m_VV   = UploadArray(d, op_ext->vv);
	m_VVFO = UploadArray(d, op_ext->vvfo);
	m_VVFN = UploadArray(d, op_ext->vvfn);
	m_II   = UploadArray(d, op_ext->ii);
	m_IIFO = UploadArray(d, op_ext->iifo);
	m_IIFN = UploadArray(d, op_ext->iifn);

	d->Pipeline(UPML_SOURCE, "upml_pre");
	d->Pipeline(UPML_SOURCE, "upml_post");
}

void Metal_Ext_UPML::SetRegion(id<MTLComputeCommandEncoder> enc)
{
	d->SetGridDim(4);
	[enc setBytes:&m_Region length:sizeof(m_Region) atIndex:5];
}

void Metal_Ext_UPML::Pre(id<MTLBuffer> field, id<MTLBuffer> flux, id<MTLBuffer> c_old, id<MTLBuffer> c_fo)
{
	id<MTLComputePipelineState> pso = d->Pipeline(UPML_SOURCE, "upml_pre");
	id<MTLComputeCommandEncoder> enc = d->Encoder();
	[enc setComputePipelineState:pso];
	[enc setBuffer:field offset:0 atIndex:0];
	[enc setBuffer:flux offset:0 atIndex:1];
	[enc setBuffer:c_old offset:0 atIndex:2];
	[enc setBuffer:c_fo offset:0 atIndex:3];
	SetRegion(enc);
	d->Dispatch(pso, m_Region.lz, m_Region.ly, m_Region.lx);
}

void Metal_Ext_UPML::Post(id<MTLBuffer> field, id<MTLBuffer> flux, id<MTLBuffer> c_fn)
{
	id<MTLComputePipelineState> pso = d->Pipeline(UPML_SOURCE, "upml_post");
	id<MTLComputeCommandEncoder> enc = d->Encoder();
	[enc setComputePipelineState:pso];
	[enc setBuffer:field offset:0 atIndex:0];
	[enc setBuffer:flux offset:0 atIndex:1];
	[enc setBuffer:c_fn offset:0 atIndex:2];
	SetRegion(enc);
	d->Dispatch(pso, m_Region.lz, m_Region.ly, m_Region.lx);
}

GPU_Extension* Metal_CreateExt_UPML(GPU_Backend_Metal::Impl* d, Engine_Extension* eng_ext, Engine* eng)
{
	UNUSED(eng);
	if (!dynamic_cast<Engine_Ext_UPML*>(eng_ext))
		return NULL;
	Operator_Ext_UPML* op_ext = dynamic_cast<Operator_Ext_UPML*>(eng_ext->GetOperatorExtension());
	if (!op_ext)
		return NULL;
	return new Metal_Ext_UPML(d, op_ext);
}
