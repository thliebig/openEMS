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
#include "FDTD/extensions/engine_ext_steadystate.h"

#include <algorithm>
#include <iostream>

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

// Fused: upml_pre, the main update (update_voltages) and upml_post of the nodes of
// the region, operation by operation the same. Valid because the UPML hooks run
// directly before and after the main update (see Metal_Ext_UPML::CanFuse()), and the
// voltage update only reads the voltage of its own node.
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
                            constant GridDim& N        [[buffer(8)]],
                            constant UPMLParam& P      [[buffer(9)]],
                            uint3 gid [[thread_position_in_grid]])
{
	const uint x = P.sx + gid.z, y = P.sy + gid.y, z = P.sz + gid.x;
	const uint sn = N.nx*N.ny*N.nz;
	const uint i  = nijk(N, 0, x, y, z);
	const uint xm = (x>0) ? N.ny*N.nz : 0;
	const uint ym = (y>0) ? N.nz : 0;
	const uint zm = (z>0) ? 1 : 0;
	const MainCoeff C = main_coeff(index, ca, cb, mode, sn, i, 0);
	float curl[3];
	curl[0] = (curr[2*sn+i] - curr[2*sn+i-ym] - curr[sn+i] + curr[sn+i-zm]);
	curl[1] = (curr[i] - curr[i-zm] - curr[2*sn+i] + curr[2*sn+i-xm]);
	curl[2] = (curr[sn+i] - curr[sn+i-xm] - curr[i] + curr[i-ym]);
	for (uint n=0; n<3; ++n)
	{
		const uint l = upml_local(P, n, gid);
		const uint g = n*sn + i;
		const float f_help = c_old[l]*volt[g] - c_fo[l]*flux[l];
		float v;
		v  = flux[l] * C.a[n*C.s];
		v += C.b[n*C.s] * curl[n];
		flux[l] = v;
		volt[g] = f_help + c_fn[l]*v;
	}
}

// same for the currents (update_currents), which are not updated on the last mesh lines
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
                            constant GridDim& N        [[buffer(8)]],
                            constant UPMLParam& P      [[buffer(9)]],
                            uint3 gid [[thread_position_in_grid]])
{
	const uint x = P.sx + gid.z, y = P.sy + gid.y, z = P.sz + gid.x;
	const uint sn = N.nx*N.ny*N.nz;
	const uint i  = nijk(N, 0, x, y, z);
	const bool update = (x+1<N.nx) && (y+1<N.ny) && (z+1<N.nz);
	float curl[3] = {0, 0, 0};
	MainCoeff C = {ca, cb, 0};
	if (update)
	{
		const uint xp = N.ny*N.nz;
		const uint yp = N.nz;
		C = main_coeff(index, ca, cb, mode, sn, i, 6);
		curl[0] = (volt[2*sn+i] - volt[2*sn+i+yp] - volt[sn+i] + volt[sn+i+1]);
		curl[1] = (volt[i] - volt[i+1] - volt[2*sn+i] + volt[2*sn+i+xp]);
		curl[2] = (volt[sn+i] - volt[sn+i+xp] - volt[i] + volt[i+yp]);
	}
	for (uint n=0; n<3; ++n)
	{
		const uint l = upml_local(P, n, gid);
		const uint g = n*sn + i;
		const float f_help = c_old[l]*curr[g] - c_fo[l]*flux[l];
		float c = flux[l];
		if (update)
		{
			c  = c * C.a[n*C.s];
			c += C.b[n*C.s] * curl[n];
		}
		flux[l] = c;
		curr[g] = f_help + c_fn[l]*c;
	}
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
	virtual void DoPostVoltageUpdates() {if (Fused()) Fuse("upml_fused_volt", d->volt, d->curr, m_VoltFlux, m_VV, m_VVFO, m_VVFN, d->vv, d->vi, &UPML_GROUP[1]); else Post(d->volt, m_VoltFlux, m_VVFN, &UPML_GROUP[1]);}
	virtual void DoPreCurrentUpdates()  {if (!Fused()) Pre(d->curr, m_CurrFlux, m_II, m_IIFO, &UPML_GROUP[2]);}
	virtual void DoPostCurrentUpdates() {if (Fused()) Fuse("upml_fused_curr", d->curr, d->volt, m_CurrFlux, m_II, m_IIFO, m_IIFN, d->ii, d->iv, &UPML_GROUP[3]); else Post(d->curr, m_CurrFlux, m_IIFN, &UPML_GROUP[3]);}

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
	struct {uint32_t sx, sy, sz, lx, ly, lz;} m_Region;

	id<MTLBuffer> m_VoltFlux, m_CurrFlux;
	id<MTLBuffer> m_VV, m_VVFO, m_VVFN;
	id<MTLBuffer> m_II, m_IIFO, m_IIFN;
};

static id<MTLBuffer> UploadArray(GPU_Backend_Metal::Impl* d, const ArrayLib::ArrayNIJK<FDTD_FLOAT>& arr)
{
	return d->NewBuffer(arr.size()*sizeof(FDTD_FLOAT), arr.data());
}

Metal_Ext_UPML::Metal_Ext_UPML(GPU_Backend_Metal::Impl* impl, Operator_Ext_UPML* op_ext, Engine* eng)
{
	d = impl;
	m_Eng = eng;
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
}

bool Metal_Ext_UPML::Fused()
{
	if (d->upml_fused<0)
	{
		d->upml_fused = CanFuse();
		if (d->upml_fused)
			std::cout << "GPU_Backend_Metal: " << d->upml.size() << " UPML regions fused with the main updates" << std::endl;
	}
	return d->upml_fused>0;
}

// The fused kernels replace the UPML hooks and the main update of the UPML regions.
// This requires:
// - the UPML hooks run directly before and after the main update: only extensions
//   without device pre/post hooks (steady-state) come before the UPML extensions,
// - all UPML extensions of the grid run on the device,
// - the regions cover exactly the nodes outside a box, which the main updates then cover.
bool Metal_Ext_UPML::CanFuse()
{
	size_t num_upml = 0;
	size_t last_upml = 0;
	for (size_t n=0; n<m_Eng->GetExtensionCount(); ++n)
		if (dynamic_cast<Engine_Ext_UPML*>(m_Eng->GetExtension(n)))
		{
			++num_upml;
			last_upml = n;
		}
	if (num_upml!=d->upml.size())
		return false;
	for (size_t n=0; n<last_upml; ++n)
	{
		Engine_Extension* ext = m_Eng->GetExtension(n);
		if (!dynamic_cast<Engine_Ext_UPML*>(ext) && !dynamic_cast<Engine_Ext_SteadyState*>(ext))
			return false;
	}

	// mark the region nodes, the others must form a box
	const Metal_GridDim& N = d->dim;
	std::vector<unsigned char> mask(d->numCells, 0);
	for (size_t r=0; r<d->upml.size(); ++r)
	{
		const Metal_Ext_UPML* u = static_cast<const Metal_Ext_UPML*>(d->upml[r]);
		const unsigned int sx = u->m_Region.sx, sy = u->m_Region.sy, sz = u->m_Region.sz;
		if ((sx+u->m_Region.lx>N.nx) || (sy+u->m_Region.ly>N.ny) || (sz+u->m_Region.lz>N.nz))
			return false;
		for (unsigned int x=sx; x<sx+u->m_Region.lx; ++x)
			for (unsigned int y=sy; y<sy+u->m_Region.ly; ++y)
				for (unsigned int z=sz; z<sz+u->m_Region.lz; ++z)
				{
					unsigned char& m = mask[((size_t)x*N.ny + y)*N.nz + z];
					if (m)
						return false;   // overlapping regions
					m = 1;
				}
	}
	Metal_GridDim lo = N, hi = {0, 0, 0};
	size_t free_nodes = 0;
	for (unsigned int x=0; x<N.nx; ++x)
		for (unsigned int y=0; y<N.ny; ++y)
			for (unsigned int z=0; z<N.nz; ++z)
				if (!mask[((size_t)x*N.ny + y)*N.nz + z])
				{
					++free_nodes;
					lo.nx = std::min(lo.nx, x); hi.nx = std::max(hi.nx, x+1);
					lo.ny = std::min(lo.ny, y); hi.ny = std::max(hi.ny, y+1);
					lo.nz = std::min(lo.nz, z); hi.nz = std::max(hi.nz, z+1);
				}
	if (free_nodes==0)
		return false;
	if (free_nodes!=(size_t)(hi.nx-lo.nx)*(hi.ny-lo.ny)*(hi.nz-lo.nz))
		return false;

	d->main_start = lo;
	d->main_stop = hi;
	return true;
}

// full_a/full_b: the full main coefficient arrays, used if the coefficients are not compressed
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
