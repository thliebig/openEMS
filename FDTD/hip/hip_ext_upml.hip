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

#include <algorithm>
#include <iostream>

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

// thread index: local (z, y, x) of the region
#define UPML_THREAD \
	const unsigned int z = blockIdx.x*blockDim.x + threadIdx.x; \
	const unsigned int y = blockIdx.y*blockDim.y + threadIdx.y; \
	const unsigned int x = blockIdx.z; \
	if (z>=R.lz || y>=R.ly || x>=R.lx) \
		return;

// before the main update: store the scaled old field in the flux, hand the old flux to the engine
__global__ void upml_pre(float* field, float* flux, const float* c_old, const float* c_fo, HIP_GridDim N, UPMLRegion R)
{
	UPML_THREAD
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
	UPML_THREAD
	for (unsigned int n=0; n<3; ++n)
	{
		const unsigned int l = upml_local(R, n, x, y, z);
		const unsigned int g = nijk(N, n, R.sx + x, R.sy + y, R.sz + z);
		const float f_help = flux[l];
		flux[l] = field[g];
		field[g] = f_help + c_fn[l]*flux[l];
	}
}

//! Arguments of the fused kernels: the flux and UPML coefficients of the region, the main coefficients
struct UPMLFusedArgs
{
	float* flux;
	const float* c_old; const float* c_fo; const float* c_fn;   //!< full UPML coefficient arrays (mode 0)
	const void* u_index; const float* u_sets; unsigned int u_mode;
	const void* index; const float* ca; const float* cb; unsigned int mode;   //!< main coefficients, see main_coeff()
};

//! UPML coefficients of the region cell l0: old[n*s], fo[n*s] and fn[n*s] are vv/vvfo/vvfn (set_offset 0)
//! or ii/iifo/iifn (set_offset 9) of direction n. Mode 0: the full arrays, else a 16/32 bit set index and the sets of 18.
struct UPMLCoeff { const float* old; const float* fo; const float* fn; unsigned int s; };

__device__ __forceinline__ UPMLCoeff upml_coeff(const UPMLFusedArgs& A, unsigned int l0, unsigned int cells, unsigned int set_offset)
{
	UPMLCoeff c;
	if (A.u_mode==0)
	{
		c.old = A.c_old + l0; c.fo = A.c_fo + l0; c.fn = A.c_fn + l0; c.s = cells;
		return c;
	}
	const unsigned int set = (A.u_mode==1) ? (unsigned int)((const unsigned short*)A.u_index)[l0] : ((const unsigned int*)A.u_index)[l0];
	c.old = A.u_sets + 18*set + set_offset; c.fo = c.old + 3; c.fn = c.old + 6; c.s = 1;
	return c;
}

// Fused: upml_pre, the main update (update_voltages) and upml_post of the nodes of the
// region, operation by operation the same (see GPU_UPMLFusionBox() for the conditions).
__global__ void upml_fused_volt(float* volt, const float* curr, UPMLFusedArgs A, HIP_GridDim N, UPMLRegion R)
{
	UPML_THREAD
	const unsigned int gx = R.sx + x, gy = R.sy + y, gz = R.sz + z;
	const unsigned int sn = N.nx*N.ny*N.nz;
	const unsigned int i  = nijk(N, 0, gx, gy, gz);
	const unsigned int xm = (gx>0) ? N.ny*N.nz : 0;
	const unsigned int ym = (gy>0) ? N.nz : 0;
	const unsigned int zm = (gz>0) ? 1 : 0;
	const HIP_MainCoeff C = main_coeff(A.index, A.ca, A.cb, A.mode, sn, i, 0);
	const unsigned int cells = R.lx*R.ly*R.lz;
	const UPMLCoeff U = upml_coeff(A, (x*R.ly + y)*R.lz + z, cells, 0);
	float curl[3];
	curl[0] = (curr[2*sn+i] - curr[2*sn+i-ym] - curr[sn+i] + curr[sn+i-zm]);
	curl[1] = (curr[i] - curr[i-zm] - curr[2*sn+i] + curr[2*sn+i-xm]);
	curl[2] = (curr[sn+i] - curr[sn+i-xm] - curr[i] + curr[i-ym]);
	for (unsigned int n=0; n<3; ++n)
	{
		const unsigned int l = upml_local(R, n, x, y, z);
		const unsigned int g = n*sn + i;
		const float f_help = U.old[n*U.s]*volt[g] - U.fo[n*U.s]*A.flux[l];
		float v;
		v  = A.flux[l] * C.a[n*C.s];
		v += C.b[n*C.s] * curl[n];
		A.flux[l] = v;
		volt[g] = f_help + U.fn[n*U.s]*v;
	}
}

// same for the currents (update_currents), which are not updated on the last mesh lines
__global__ void upml_fused_curr(float* curr, const float* volt, UPMLFusedArgs A, HIP_GridDim N, UPMLRegion R)
{
	UPML_THREAD
	const unsigned int gx = R.sx + x, gy = R.sy + y, gz = R.sz + z;
	const unsigned int sn = N.nx*N.ny*N.nz;
	const unsigned int i  = nijk(N, 0, gx, gy, gz);
	const bool update = (gx+1<N.nx) && (gy+1<N.ny) && (gz+1<N.nz);
	const unsigned int cells = R.lx*R.ly*R.lz;
	const UPMLCoeff U = upml_coeff(A, (x*R.ly + y)*R.lz + z, cells, 9);
	float curl[3] = {0, 0, 0};
	HIP_MainCoeff C = {A.ca, A.cb, 0};
	if (update)
	{
		const unsigned int xp = N.ny*N.nz;
		const unsigned int yp = N.nz;
		C = main_coeff(A.index, A.ca, A.cb, A.mode, sn, i, 6);
		curl[0] = (volt[2*sn+i] - volt[2*sn+i+yp] - volt[sn+i] + volt[sn+i+1]);
		curl[1] = (volt[i] - volt[i+1] - volt[2*sn+i] + volt[2*sn+i+xp]);
		curl[2] = (volt[sn+i] - volt[sn+i+xp] - volt[i] + volt[i+yp]);
	}
	for (unsigned int n=0; n<3; ++n)
	{
		const unsigned int l = upml_local(R, n, x, y, z);
		const unsigned int g = n*sn + i;
		const float f_help = U.old[n*U.s]*curr[g] - U.fo[n*U.s]*A.flux[l];
		float c = A.flux[l];
		if (update)
		{
			c  = c * C.a[n*C.s];
			c += C.b[n*C.s] * curl[n];
		}
		A.flux[l] = c;
		curr[g] = f_help + U.fn[n*U.s]*c;
	}
}

class HIP_Ext_UPML : public GPU_Extension
{
public:
	HIP_Ext_UPML(GPU_Backend_HIP::Impl* impl, Operator_Ext_UPML* op_ext, Engine* eng);
	virtual ~HIP_Ext_UPML();

	// fused: nothing before the main update, the region update after it
	virtual void DoPreVoltageUpdates()  {if (!Fused()) Pre(d->volt, m_VoltFlux, m_VV, m_VVFO);}
	virtual void DoPostVoltageUpdates() {if (Fused()) FuseVolt(); else Post(d->volt, m_VoltFlux, m_VVFN);}
	virtual void DoPreCurrentUpdates()  {if (!Fused()) Pre(d->curr, m_CurrFlux, m_II, m_IIFO);}
	virtual void DoPostCurrentUpdates() {if (Fused()) FuseCurr(); else Post(d->curr, m_CurrFlux, m_IIFN);}

protected:
	void Pre(float* field, float* flux, const float* c_old, const float* c_fo)
	{
		HIP_Launch(d, "upml_pre", upml_pre, m_Region.lz, m_Region.ly, m_Region.lx, field, flux, c_old, c_fo, d->dim, m_Region);
	}
	void Post(float* field, float* flux, const float* c_fn)
	{
		HIP_Launch(d, "upml_post", upml_post, m_Region.lz, m_Region.ly, m_Region.lx, field, flux, c_fn, d->dim, m_Region);
	}
	void FuseVolt()
	{
		UPMLFusedArgs A = Args(m_VoltFlux, m_VV, m_VVFO, m_VVFN, d->vv, d->vi);
		HIP_Launch(d, "upml_fused_volt", upml_fused_volt, m_Region.lz, m_Region.ly, m_Region.lx, d->volt, (const float*)d->curr, A, d->dim, m_Region);
	}
	void FuseCurr()
	{
		UPMLFusedArgs A = Args(m_CurrFlux, m_II, m_IIFO, m_IIFN, d->ii, d->iv);
		HIP_Launch(d, "upml_fused_curr", upml_fused_curr, m_Region.lz, m_Region.ly, m_Region.lx, d->curr, (const float*)d->volt, A, d->dim, m_Region);
	}
	UPMLFusedArgs Args(float* flux, const float* c_old, const float* c_fo, const float* c_fn, const float* full_a, const float* full_b)
	{
		UPMLFusedArgs A;
		A.flux = flux;
		A.c_old = c_old; A.c_fo = c_fo; A.c_fn = c_fn;
		A.u_index = m_SetIndex; A.u_sets = m_Sets; A.u_mode = m_CoeffMode;
		A.index = d->index; A.mode = d->coeff_mode;
		A.ca = d->coeff_mode ? d->coeff : full_a;
		A.cb = d->coeff_mode ? d->coeff : full_b;
		return A;
	}
	float* Upload(const ArrayLib::ArrayNIJK<FDTD_FLOAT>& arr) {return d->Alloc<float>(arr.size(), arr.data());}

	//! Whether the UPML of this grid runs fused with the main updates, decided at the first call
	bool Fused();
	bool CanFuse();

	GPU_Backend_HIP::Impl* d;
	Engine* m_Eng;
	UPMLRegion m_Region;

	float *m_VoltFlux, *m_CurrFlux;
	float *m_VV, *m_VVFO, *m_VVFN;
	float *m_II, *m_IIFO, *m_IIFN;

	//! Coefficients for the fused kernels: 0: the full arrays, 1/2: a 16/32 bit set index per cell and the sets (see upml_coeff())
	unsigned int m_CoeffMode;
	void* m_SetIndex;
	float* m_Sets;
};

HIP_Ext_UPML::HIP_Ext_UPML(GPU_Backend_HIP::Impl* impl, Operator_Ext_UPML* op_ext, Engine* eng)
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

	const size_t cells = (size_t)m_Region.lx*m_Region.ly*m_Region.lz;
	m_VoltFlux = d->Alloc<float>(3*cells);
	m_CurrFlux = d->Alloc<float>(3*cells);
	m_VV   = Upload(op_ext->vv);
	m_VVFO = Upload(op_ext->vvfo);
	m_VVFN = Upload(op_ext->vvfn);
	m_II   = Upload(op_ext->ii);
	m_IIFO = Upload(op_ext->iifo);
	m_IIFN = Upload(op_ext->iifn);

	// the 18 coefficients of a cell for the fused kernels: vv[3], vvfo[3], vvfn[3], ii[3], iifo[3], iifn[3]
	const FDTD_FLOAT* src[6] = {op_ext->vv.data(), op_ext->vvfo.data(), op_ext->vvfn.data(),
	                            op_ext->ii.data(), op_ext->iifo.data(), op_ext->iifn.data()};
	GPU_CoeffSets sets;
	m_CoeffMode = 0;
	m_SetIndex = NULL;
	m_Sets = NULL;
	if (GPU_FindSets(cells, 18, [&](size_t l0, float* values)
	    {
		    for (int c=0; c<6; ++c)
			    for (int n=0; n<3; ++n)
				    values[3*c+n] = src[c][n*cells + l0];
	    }, sets))
	{
		m_CoeffMode = sets.mode;
		m_SetIndex = d->AllocIndex(sets);
		m_Sets = d->Alloc<float>(sets.table.size(), sets.table.data());
	}
}

HIP_Ext_UPML::~HIP_Ext_UPML()
{
	d->upml.erase(std::remove(d->upml.begin(), d->upml.end(), this), d->upml.end());
	// the main updates cover the whole grid again
	d->upml_fused = -1;
	d->main_start.nx = d->main_start.ny = d->main_start.nz = 0;
	d->main_stop = d->dim;
}

bool HIP_Ext_UPML::Fused()
{
	if (d->upml_fused<0)
	{
		d->upml_fused = CanFuse();
		if (d->upml_fused)
			std::cout << "GPU_Backend_HIP: " << d->upml.size() << " UPML regions fused with the main updates" << std::endl;
	}
	return d->upml_fused>0;
}

// The fused kernels replace the UPML hooks and the main update of the UPML regions,
// see GPU_UPMLFusionBox() for the conditions.
bool HIP_Ext_UPML::CanFuse()
{
	std::vector<GPU_UPMLRegion> regions;
	for (size_t r=0; r<d->upml.size(); ++r)
	{
		const HIP_Ext_UPML* u = static_cast<const HIP_Ext_UPML*>(d->upml[r]);
		const GPU_UPMLRegion R = {{u->m_Region.sx, u->m_Region.sy, u->m_Region.sz}, {u->m_Region.lx, u->m_Region.ly, u->m_Region.lz}};
		regions.push_back(R);
	}
	const unsigned int numLines[3] = {d->dim.nx, d->dim.ny, d->dim.nz};
	unsigned int start[3], stop[3];
	if (!GPU_UPMLFusionBox(m_Eng, regions, numLines, start, stop))
		return false;
	d->main_start.nx = start[0]; d->main_start.ny = start[1]; d->main_start.nz = start[2];
	d->main_stop.nx = stop[0]; d->main_stop.ny = stop[1]; d->main_stop.nz = stop[2];
	return true;
}

GPU_Extension* HIP_CreateExt_UPML(GPU_Backend_HIP::Impl* d, Engine_Extension* eng_ext, Engine* eng)
{
	if (!dynamic_cast<Engine_Ext_UPML*>(eng_ext))
		return NULL;
	Operator_Ext_UPML* op_ext = dynamic_cast<Operator_Ext_UPML*>(eng_ext->GetOperatorExtension());
	if (!op_ext)
		return NULL;
	return new HIP_Ext_UPML(d, op_ext, eng);
}
