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

#include <algorithm>
#include <iostream>

// Uniaxial PML, see Engine_Ext_UPML. One thread per cell of the PML region,
// gid = local (z, y, x). The flux and coefficient arrays use the local region size.
static const char* UPML_SOURCE = R"MSL(
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

// the fused kernels of a whole region, see upml_fused_volt_node()
kernel void upml_fused_volt(device float* volt         [[buffer(0)]],
                            const device float* curr   [[buffer(1)]],
                            device float* flux         [[buffer(2)]],
                            const device float* c_old  [[buffer(3)]],
                            const device float* c_fo   [[buffer(4)]],
                            const device float* c_fn   [[buffer(5)]],
                            const device void* index   [[buffer(6)]],
                            const device float* ca     [[buffer(7)]],
                            const device float* cb     [[buffer(10)]],
                            constant uint& mode        [[buffer(11)]],
                            const device void* u_index [[buffer(12)]],
                            const device float* u_sets [[buffer(13)]],
                            constant uint& u_mode      [[buffer(14)]],
                            constant GridDim& N        [[buffer(8)]],
                            constant UPMLParam& P      [[buffer(9)]],
                            uint3 gid [[thread_position_in_grid]])
{
	upml_fused_volt_node(volt, curr, flux, c_old, c_fo, c_fn, index, ca, cb, mode, u_index, u_sets, u_mode, N, P, gid);
}

kernel void upml_fused_curr(device float* curr         [[buffer(0)]],
                            const device float* volt   [[buffer(1)]],
                            device float* flux         [[buffer(2)]],
                            const device float* c_old  [[buffer(3)]],
                            const device float* c_fo   [[buffer(4)]],
                            const device float* c_fn   [[buffer(5)]],
                            const device void* index   [[buffer(6)]],
                            const device float* ca     [[buffer(7)]],
                            const device float* cb     [[buffer(10)]],
                            constant uint& mode        [[buffer(11)]],
                            const device void* u_index [[buffer(12)]],
                            const device float* u_sets [[buffer(13)]],
                            constant uint& u_mode      [[buffer(14)]],
                            constant GridDim& N        [[buffer(8)]],
                            constant UPMLParam& P      [[buffer(9)]],
                            uint3 gid [[thread_position_in_grid]])
{
	upml_fused_curr_node(curr, volt, flux, c_old, c_fo, c_fn, index, ca, cb, mode, u_index, u_sets, u_mode, N, P, gid);
}
)MSL";

// The regions of all UPML extensions of a grid are disjoint: the dispatches of one
// hook form a group, which can run concurrently (see GPU_Backend_Metal::Impl::Dispatch()).
// The UPML extensions have the same priority, so they are called one after the other.
static const char UPML_GROUP[4] = {0, 0, 0, 0};

class Metal_Ext_UPML : public GPU_Extension
{
public:
	Metal_Ext_UPML(GPU_Backend_Metal::Impl* impl, Operator_Ext_UPML* op_ext, Engine* eng);
	virtual ~Metal_Ext_UPML();

	// fused: nothing before the main update, the region update after it
	virtual void DoPreVoltageUpdates()  {if (!Fused()) Pre(d->volt, m_VoltFlux, m_VV, m_VVFO, &UPML_GROUP[0]);}
	virtual void DoPostVoltageUpdates() {if (Fused()) {if (!m_InMain) Fuse("upml_fused_volt", d->volt, d->curr, m_VoltFlux, m_VV, m_VVFO, m_VVFN, d->vv, d->vi, &UPML_GROUP[1]);} else Post(d->volt, m_VoltFlux, m_VVFN, &UPML_GROUP[1]);}
	virtual void DoPreCurrentUpdates()  {if (!Fused()) Pre(d->curr, m_CurrFlux, m_II, m_IIFO, &UPML_GROUP[2]);}
	virtual void DoPostCurrentUpdates() {if (Fused()) {if (!m_InMain) Fuse("upml_fused_curr", d->curr, d->volt, m_CurrFlux, m_II, m_IIFO, m_IIFN, d->ii, d->iv, &UPML_GROUP[3]);} else Post(d->curr, m_CurrFlux, m_IIFN, &UPML_GROUP[3]);}

protected:
	void Pre(id<MTLBuffer> field, id<MTLBuffer> flux, id<MTLBuffer> c_old, id<MTLBuffer> c_fo, const void* group);
	void Post(id<MTLBuffer> field, id<MTLBuffer> flux, id<MTLBuffer> c_fn, const void* group);
	void Fuse(const char* kernel, id<MTLBuffer> field, id<MTLBuffer> other, id<MTLBuffer> flux,
	          id<MTLBuffer> c_old, id<MTLBuffer> c_fo, id<MTLBuffer> c_fn,
	          id<MTLBuffer> full_a, id<MTLBuffer> full_b, const void* group);
	void SetRegion(id<MTLComputeCommandEncoder> enc);

	//! Whether the UPML of this grid runs fused with the main updates, decided at the first call
	bool Fused();
	bool CanFuse();

	GPU_Backend_Metal::Impl* d;
	Engine* m_Eng;
	bool m_InMain;   //!< fused: updated by the main kernels (see Metal_ZSlab)
	struct {uint32_t sx, sy, sz, lx, ly, lz;} m_Region;

	id<MTLBuffer> m_VoltFlux, m_CurrFlux;
	id<MTLBuffer> m_VV, m_VVFO, m_VVFN;
	id<MTLBuffer> m_II, m_IIFO, m_IIFN;

	//! Coefficients for the fused kernels: 0: the full arrays, 1/2: a 16/32 bit set index per cell and the sets (see upml_coeff())
	uint32_t m_CoeffMode;
	id<MTLBuffer> m_SetIndex, m_Sets;
};

static id<MTLBuffer> UploadArray(GPU_Backend_Metal::Impl* d, const ArrayLib::ArrayNIJK<FDTD_FLOAT>& arr)
{
	return d->NewBuffer(arr.size()*sizeof(FDTD_FLOAT), arr.data());
}

Metal_Ext_UPML::Metal_Ext_UPML(GPU_Backend_Metal::Impl* impl, Operator_Ext_UPML* op_ext, Engine* eng)
{
	d = impl;
	m_Eng = eng;
	m_InMain = false;
	d->upml.push_back(this);
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

	// the 18 coefficients of a cell for the fused kernels: vv[3], vvfo[3], vvfn[3], ii[3], iifo[3], iifn[3]
	const size_t cells = (size_t)m_Region.lx*m_Region.ly*m_Region.lz;
	const FDTD_FLOAT* src[6] = {op_ext->vv.data(), op_ext->vvfo.data(), op_ext->vvfn.data(),
	                            op_ext->ii.data(), op_ext->iifo.data(), op_ext->iifn.data()};
	GPU_CoeffSets sets;
	m_CoeffMode = 0;
	if (GPU_FindSets(cells, 18, [&](size_t l0, float* values)
	    {
		    for (int c=0; c<6; ++c)
			    for (int n=0; n<3; ++n)
				    values[3*c+n] = src[c][n*cells + l0];
	    }, sets))
	{
		m_CoeffMode = sets.mode;
		m_SetIndex = d->NewIndexBuffer(sets);
		m_Sets = d->NewBuffer(sets.table.size()*sizeof(float), sets.table.data());
	}

	d->Pipeline(UPML_SOURCE, "upml_pre");
	d->Pipeline(UPML_SOURCE, "upml_post");
	d->Pipeline(UPML_SOURCE, "upml_fused_volt");
	d->Pipeline(UPML_SOURCE, "upml_fused_curr");
}

Metal_Ext_UPML::~Metal_Ext_UPML()
{
	d->upml.erase(std::remove(d->upml.begin(), d->upml.end(), this), d->upml.end());
	// the main updates cover the whole grid again
	d->upml_fused = -1;
	d->main_start.nx = d->main_start.ny = d->main_start.nz = 0;
	d->main_stop = d->dim;
	d->ResetZSlabs();
	for (size_t r=0; r<d->upml.size(); ++r)
		static_cast<Metal_Ext_UPML*>(d->upml[r])->m_InMain = false;
}

bool Metal_Ext_UPML::Fused()
{
	if (d->upml_fused<0)
	{
		d->upml_fused = CanFuse();
		if (d->upml_fused)
			std::cout << "GPU_Backend_Metal: " << d->upml.size() << " UPML regions fused with the main updates ("
			          << d->zslab[0].active+d->zslab[1].active << " along z in the main kernels)" << std::endl;
	}
	return d->upml_fused>0;
}

// The fused kernels replace the UPML hooks and the main update of the UPML regions,
// see GPU_UPMLFusionBox() for the conditions.
bool Metal_Ext_UPML::CanFuse()
{
	std::vector<GPU_UPMLRegion> regions;
	for (size_t r=0; r<d->upml.size(); ++r)
	{
		const Metal_Ext_UPML* u = static_cast<const Metal_Ext_UPML*>(d->upml[r]);
		const GPU_UPMLRegion R = {{u->m_Region.sx, u->m_Region.sy, u->m_Region.sz}, {u->m_Region.lx, u->m_Region.ly, u->m_Region.lz}};
		regions.push_back(R);
	}
	const unsigned int numLines[3] = {d->dim.nx, d->dim.ny, d->dim.nz};
	unsigned int start[3], stop[3];
	if (!GPU_UPMLFusionBox(m_Eng, regions, numLines, start, stop))
		return false;
	d->main_start.nx = start[0]; d->main_start.ny = start[1]; d->main_start.nz = start[2];
	d->main_stop.nx = stop[0]; d->main_stop.ny = stop[1]; d->main_stop.nz = stop[2];
	d->ResetZSlabs();

	// The regions along z (which cover the x/y range of the main updates) are updated by the
	// main kernels, along full z lines: a region of their own would only have short rows
	// (e.g. 9 nodes), which use the memory transactions poorly.
	for (size_t r=0; r<d->upml.size(); ++r)
	{
		Metal_Ext_UPML* u = static_cast<Metal_Ext_UPML*>(d->upml[r]);
		const auto& R = u->m_Region;
		if ((R.sx!=start[0]) || (R.lx!=stop[0]-start[0]) || (R.sy!=start[1]) || (R.ly!=stop[1]-start[1]))
			continue;
		int s = -1;
		if (!d->zslab[0].active && (R.sz+R.lz==start[2]))
			s = 0;
		else if (!d->zslab[1].active && (R.sz==stop[2]))
			s = 1;
		if (s<0)
			continue;
		Metal_ZSlab& z = d->zslab[s];
		z.active = true;
		z.region.sx = R.sx; z.region.sy = R.sy; z.region.sz = R.sz;
		z.region.lx = R.lx; z.region.ly = R.ly; z.region.lz = R.lz;
		z.flux[0] = u->m_VoltFlux; z.c_old[0] = u->m_VV; z.c_fo[0] = u->m_VVFO; z.c_fn[0] = u->m_VVFN;
		z.flux[1] = u->m_CurrFlux; z.c_old[1] = u->m_II; z.c_fo[1] = u->m_IIFO; z.c_fn[1] = u->m_IIFN;
		z.u_mode = u->m_CoeffMode;
		z.u_index = u->m_SetIndex;
		z.u_sets = u->m_Sets;
		u->m_InMain = true;
	}
	return true;
}

void Metal_Ext_UPML::Fuse(const char* kernel, id<MTLBuffer> field, id<MTLBuffer> other, id<MTLBuffer> flux,
                          id<MTLBuffer> c_old, id<MTLBuffer> c_fo, id<MTLBuffer> c_fn,
                          id<MTLBuffer> full_a, id<MTLBuffer> full_b, const void* group)
{
	id<MTLComputePipelineState> pso = d->Pipeline(UPML_SOURCE, kernel);
	id<MTLComputeCommandEncoder> enc = d->Encoder();
	[enc setComputePipelineState:pso];
	[enc setBuffer:field offset:0 atIndex:0];
	[enc setBuffer:other offset:0 atIndex:1];
	[enc setBuffer:flux offset:0 atIndex:2];
	[enc setBuffer:c_old offset:0 atIndex:3];
	[enc setBuffer:c_fo offset:0 atIndex:4];
	[enc setBuffer:c_fn offset:0 atIndex:5];
	d->SetCoefficients(6, 7, 10, 11, full_a, full_b);
	// unused arguments get a valid buffer
	[enc setBuffer:(m_CoeffMode ? m_SetIndex : c_old) offset:0 atIndex:12];
	[enc setBuffer:(m_CoeffMode ? m_Sets : c_old) offset:0 atIndex:13];
	[enc setBytes:&m_CoeffMode length:sizeof(m_CoeffMode) atIndex:14];
	d->SetGridDim(8);
	[enc setBytes:&m_Region length:sizeof(m_Region) atIndex:9];
	d->Dispatch(pso, m_Region.lz, m_Region.ly, m_Region.lx, group);
}

void Metal_Ext_UPML::SetRegion(id<MTLComputeCommandEncoder> enc)
{
	d->SetGridDim(4);
	[enc setBytes:&m_Region length:sizeof(m_Region) atIndex:5];
}

void Metal_Ext_UPML::Pre(id<MTLBuffer> field, id<MTLBuffer> flux, id<MTLBuffer> c_old, id<MTLBuffer> c_fo, const void* group)
{
	id<MTLComputePipelineState> pso = d->Pipeline(UPML_SOURCE, "upml_pre");
	id<MTLComputeCommandEncoder> enc = d->Encoder();
	[enc setComputePipelineState:pso];
	[enc setBuffer:field offset:0 atIndex:0];
	[enc setBuffer:flux offset:0 atIndex:1];
	[enc setBuffer:c_old offset:0 atIndex:2];
	[enc setBuffer:c_fo offset:0 atIndex:3];
	SetRegion(enc);
	d->Dispatch(pso, m_Region.lz, m_Region.ly, m_Region.lx, group);
}

void Metal_Ext_UPML::Post(id<MTLBuffer> field, id<MTLBuffer> flux, id<MTLBuffer> c_fn, const void* group)
{
	id<MTLComputePipelineState> pso = d->Pipeline(UPML_SOURCE, "upml_post");
	id<MTLComputeCommandEncoder> enc = d->Encoder();
	[enc setComputePipelineState:pso];
	[enc setBuffer:field offset:0 atIndex:0];
	[enc setBuffer:flux offset:0 atIndex:1];
	[enc setBuffer:c_fn offset:0 atIndex:2];
	SetRegion(enc);
	d->Dispatch(pso, m_Region.lz, m_Region.ly, m_Region.lx, group);
}

GPU_Extension* Metal_CreateExt_UPML(GPU_Backend_Metal::Impl* d, Engine_Extension* eng_ext, Engine* eng)
{
	if (!dynamic_cast<Engine_Ext_UPML*>(eng_ext))
		return NULL;
	Operator_Ext_UPML* op_ext = dynamic_cast<Operator_Ext_UPML*>(eng_ext->GetOperatorExtension());
	if (!op_ext)
		return NULL;
	return new Metal_Ext_UPML(d, op_ext, eng);
}
