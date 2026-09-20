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
#include <algorithm>
#include <iostream>

#include "hip_internal.h"
#include "FDTD/engine.h"
#include "FDTD/extensions/operator_ext_upml.h"
#include "FDTD/extensions/engine_ext_upml.h"

// thread index: local (z, y) flattened along the threads, x along the grid z-dimension;
// keeps the warps full for thin regions (e.g. 9 nodes along z)
#define UPML_THREAD \
	const unsigned int t = blockIdx.x*blockDim.x + threadIdx.x; \
	const unsigned int x = blockIdx.z; \
	if (t>=R.lz*R.ly || x>=R.lx) \
		return; \
	const unsigned int z = t % R.lz; \
	const unsigned int y = t / R.lz;

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

// the fused kernels of a whole region, in place or from one buffer of the fused step into the other
__global__ void upml_fused_volt(float* volt_out, const float* volt_in, const float* curr, UPMLFusedArgs A, HIP_GridDim N, UPMLRegion R)
{
	UPML_THREAD
	upml_fused_volt_node(volt_out, volt_in, curr, A, N, R, x, y, z);
}

__global__ void upml_fused_curr(float* curr_out, const float* curr_in, const float* volt, UPMLFusedArgs A, HIP_GridDim N, UPMLRegion R)
{
	UPML_THREAD
	upml_fused_curr_node(curr_out, curr_in, volt, A, N, R, x, y, z);
}

class HIP_Ext_UPML : public GPU_Extension
{
public:
	HIP_Ext_UPML(GPU_Backend_HIP::Impl* impl, Operator_Ext_UPML* op_ext, Engine* eng);
	virtual ~HIP_Ext_UPML();

	// fused: nothing before the main update, the region update after it (unless done by the main kernels);
	// fused step (see update_fused): the voltages of the region before it (the main nodes next to the
	// region need them), the currents after the voltage extensions, into the next buffers
	virtual void DoPreVoltageUpdates()
	{
		if (!Fused())
			Pre(d->volt, m_VoltFlux, m_VV, m_VVFO);
		else if ((d->fused_step>0) && d->fregions.count)
		{
			// in the kernel: the current voltage flux in, the new one into the second buffer,
			// as blocks that only need a border voltage compute the node again on the old flux
			d->fregions.volt[m_Side].flux = m_VoltFlux;
			d->fregions.volt[m_Side].flux_out = m_VoltFlux2;
			std::swap(m_VoltFlux, m_VoltFlux2);
		}
		else if (d->fused_step>0)
			FuseVolt(d->volt_next);
	}
	virtual void DoPostVoltageUpdates()
	{
		if (Fused() && (d->fused_step<=0) && !m_InMain)
			FuseVolt(d->volt);
		else if (!Fused())
			Post(d->volt, m_VoltFlux, m_VVFN);
	}
	virtual void DoPreCurrentUpdates()
	{
		if (!Fused())
			Pre(d->curr, m_CurrFlux, m_II, m_IIFO);
		else if ((d->fused_step>0) && !d->fregions.count)
			FuseCurr(d->curr_next);
	}
	virtual void DoPostCurrentUpdates()
	{
		if (Fused() && (d->fused_step<=0) && !m_InMain)
			FuseCurr(d->curr);
		else if (!Fused())
			Post(d->curr, m_CurrFlux, m_IIFN);
	}

protected:
	void Pre(float* field, float* flux, const float* c_old, const float* c_fo)
	{
		HIP_Launch(d, "upml_pre", upml_pre, m_Region.lz*m_Region.ly, 1, m_Region.lx, field, flux, c_old, c_fo, d->dim, m_Region);
	}
	void Post(float* field, float* flux, const float* c_fn)
	{
		HIP_Launch(d, "upml_post", upml_post, m_Region.lz*m_Region.ly, 1, m_Region.lx, field, flux, c_fn, d->dim, m_Region);
	}
	UPMLFusedArgs VoltArgs() {return Args(m_VoltFlux, m_VV, m_VVFO, m_VVFN, d->vv, d->vi);}
	UPMLFusedArgs CurrArgs() {return Args(m_CurrFlux, m_II, m_IIFO, m_IIFN, d->ii, d->iv);}
	// in place, or in the fused step from the current buffers into the next ones (see update_fused)
	void FuseVolt(float* volt_out)
	{
		HIP_Launch(d, "upml_fused_volt", upml_fused_volt, m_Region.lz*m_Region.ly, 1, m_Region.lx, volt_out, (const float*)d->volt, (const float*)d->curr, VoltArgs(), d->dim, m_Region);
	}
	void FuseCurr(float* curr_out)
	{
		HIP_Launch(d, "upml_fused_curr", upml_fused_curr, m_Region.lz*m_Region.ly, 1, m_Region.lx, curr_out, (const float*)d->curr, (const float*)d->volt, CurrArgs(), d->dim, m_Region);
	}
	UPMLFusedArgs Args(float* flux, const float* c_old, const float* c_fo, const float* c_fn, const float* full_a, const float* full_b)
	{
		UPMLFusedArgs A;
		A.flux = flux;
		A.flux_out = flux;
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
	bool m_InMain;   //!< fused: updated by the main kernels (see HIP_ZSlabs)
	int m_Side;      //!< fused step with the regions in the kernel: the side of the region (see HIP_FusedRegions)
	float* m_VoltFlux2;   //!< ... and the second buffer of the voltage flux, written by the block that owns the node

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
	m_InMain = false;
	m_Side = -1;
	m_VoltFlux2 = NULL;
	// the UPML extensions of a grid decide together whether they fuse (see CanFuse())
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
	d->ResetZSlabs();
	for (size_t r=0; r<d->upml.size(); ++r)
		static_cast<HIP_Ext_UPML*>(d->upml[r])->m_InMain = false;
}

bool HIP_Ext_UPML::Fused()
{
	if (d->upml_fused<0)
	{
		d->upml_fused = CanFuse();
		if (d->upml_fused)
			std::cout << "GPU_Backend_HIP: " << d->upml.size() << " UPML regions fused with the main updates ("
			          << d->zslab_lo+d->zslab_hi << " along z in the main kernels)" << std::endl;
		else
			d->fused_blockers++;   // the fused step needs the fused UPML
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
	d->ResetZSlabs();

	// The fused step (see update_fused) can update the regions in its kernel, if they are on the
	// sides of the main nodes in the order x, y, z (as Operator_Ext_UPML::Create_UPML() makes them),
	// else the region kernels update them.
	HIP_FusedRegions& F = d->fregions;
	F.count = 0;
	F.lo = d->main_start;
	F.hi = d->main_stop;
	bool sides_ok = true;
	for (int s=0; s<6; ++s)
		F.side[s] = -1;
	for (size_t r=0; r<d->upml.size(); ++r)
	{
		HIP_Ext_UPML* u = static_cast<HIP_Ext_UPML*>(d->upml[r]);
		const UPMLRegion& R = u->m_Region;
		const unsigned int sR[3] = {R.sx, R.sy, R.sz};
		const unsigned int lR[3] = {R.lx, R.ly, R.lz};
		const unsigned int lo[3] = {F.lo.nx, F.lo.ny, F.lo.nz};
		const unsigned int hi[3] = {F.hi.nx, F.hi.ny, F.hi.nz};
		const unsigned int N[3] = {d->dim.nx, d->dim.ny, d->dim.nz};
		int side = -1;
		for (int n=0; n<3 && side<0; ++n)
			for (int h=0; h<2 && side<0; ++h)
			{
				// along n below lo or from hi on; along the directions before n the main nodes, after it all
				bool match = h ? ((sR[n]==hi[n]) && (sR[n]+lR[n]==N[n])) : ((sR[n]==0) && (lR[n]==lo[n]));
				for (int m=0; m<3; ++m)
					if (m<n)
						match &= (sR[m]==lo[m]) && (lR[m]==hi[m]-lo[m]);
					else if (m>n)
						match &= (sR[m]==0) && (lR[m]==N[m]);
				if (match && (F.side[2*n+h]<0))
					side = 2*n+h;
			}
		if (side<0)
		{
			sides_ok = false;
			break;
		}
		F.side[side] = F.count;
		F.R[F.count] = R;
		F.volt[F.count] = u->VoltArgs();
		F.curr[F.count] = u->CurrArgs();
		u->m_Side = F.count;
		++F.count;
	}
	if (!sides_ok)
		F.count = 0;
	if (F.count)
		for (size_t r=0; r<d->upml.size(); ++r)
		{
			HIP_Ext_UPML* u = static_cast<HIP_Ext_UPML*>(d->upml[r]);
			const size_t n = 3*(size_t)u->m_Region.lx*u->m_Region.ly*u->m_Region.lz;
			u->m_VoltFlux2 = d->Alloc<float>(n);
			HIP_Check(hipMemcpyAsync(u->m_VoltFlux2, u->m_VoltFlux, n*sizeof(float), hipMemcpyDeviceToDevice, d->Stream()), "UPML flux");
		}
	if (d->DecideFusedStep())
		return true;
	F.count = 0;

	// The regions along z (which cover the x/y range of the main updates) are updated by the
	// main kernels, along full z lines: a region of their own would only have short rows
	// (e.g. 9 nodes), which use the memory transactions poorly.
	for (size_t r=0; r<d->upml.size(); ++r)
	{
		HIP_Ext_UPML* u = static_cast<HIP_Ext_UPML*>(d->upml[r]);
		const UPMLRegion& R = u->m_Region;
		if ((R.sx!=start[0]) || (R.lx!=stop[0]-start[0]) || (R.sy!=start[1]) || (R.ly!=stop[1]-start[1]))
			continue;
		if (!d->zslab_lo && (R.sz+R.lz==start[2]))
		{
			d->zslab_lo = true;
			d->zslabs.R_lo = R;
			d->zslabs.volt_lo = u->VoltArgs();
			d->zslabs.curr_lo = u->CurrArgs();
			u->m_InMain = true;
		}
		else if (!d->zslab_hi && (R.sz==stop[2]))
		{
			d->zslab_hi = true;
			d->zslabs.R_hi = R;
			d->zslabs.volt_hi = u->VoltArgs();
			d->zslabs.curr_hi = u->CurrArgs();
			u->m_InMain = true;
		}
	}
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
