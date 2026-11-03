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

#include <algorithm>
#include <cstdlib>
#include <iostream>

#include "sycl_internal.h"
#include "FDTD/engine.h"
#include "FDTD/extensions/operator_ext_upml.h"
#include "FDTD/extensions/engine_ext_upml.h"

// thread index: local (z, y) flattened along the work-items, x along the dimension 0 of the groups;
// keeps the sub-groups full for thin regions (e.g. 9 nodes along z)
#define UPML_THREAD \
	const unsigned int t = it.get_group(2)*it.get_local_range(2) + it.get_local_id(2); \
	const unsigned int x = it.get_group(0); \
	if (t>=R.lz*R.ly || x>=R.lx) \
		return; \
	const unsigned int z = t % R.lz; \
	const unsigned int y = t / R.lz;

// before the main update: store the scaled old field in the flux, hand the old flux to the engine
static inline void upml_pre(const sycl::nd_item<3>& it, float* field, UPMLFusedArgs A, unsigned int set_offset, SYCL_GridDim N, UPMLRegion R)
{
	UPML_THREAD
	const UPMLCoeff U = upml_coeff(A, (x*R.ly + y)*R.lz + z, (size_t)R.lx*R.ly*R.lz, set_offset);
	for (unsigned int n=0; n<3; ++n)
	{
		const size_t l = upml_local(R, n, x, y, z);
		const size_t g = nijk(N, n, R.sx + x, R.sy + y, R.sz + z);
		const float f_help = U.old[n*U.s]*field[g] - U.fo[n*U.s]*A.flux[l];
		field[g] = A.flux[l];
		A.flux[l] = f_help;
	}
}

// after the main update: keep the updated flux, combine it with the stored field
static inline void upml_post(const sycl::nd_item<3>& it, float* field, UPMLFusedArgs A, unsigned int set_offset, SYCL_GridDim N, UPMLRegion R)
{
	UPML_THREAD
	const UPMLCoeff U = upml_coeff(A, (x*R.ly + y)*R.lz + z, (size_t)R.lx*R.ly*R.lz, set_offset);
	for (unsigned int n=0; n<3; ++n)
	{
		const size_t l = upml_local(R, n, x, y, z);
		const size_t g = nijk(N, n, R.sx + x, R.sy + y, R.sz + z);
		const float f_help = A.flux[l];
		A.flux[l] = field[g];
		field[g] = f_help + U.fn[n*U.s]*A.flux[l];
	}
}

// the fused kernels of a whole region, in place or from one buffer of the fused step into the other
// (FACTOR: coefficient modes 3/4, see main_coeff())
template<bool FACTOR>
static inline void upml_fused_volt(const sycl::nd_item<3>& it, float* volt_out, const float* volt_in, const float* curr, UPMLFusedArgs A, SYCL_GridDim N, UPMLRegion R)
{
	UPML_THREAD
	upml_fused_volt_node<FACTOR>(volt_out, volt_in, curr, A, N, R, x, y, z);
}

template<bool FACTOR>
static inline void upml_fused_curr(const sycl::nd_item<3>& it, float* curr_out, const float* curr_in, const float* volt, UPMLFusedArgs A, SYCL_GridDim N, UPMLRegion R)
{
	UPML_THREAD
	upml_fused_curr_node<FACTOR>(curr_out, curr_in, volt, A, N, R, x, y, z);
}

// The UPML regions the region kernels update in the fused step (those along x and y), all in one
// launch per half-step: the regions are small (e.g. 8 nodes thick), as separate kernels the launch
// and the ramp of each one cost more than their nodes on small grids. A table in device memory
// (constant, as these regions update their flux in place), one region per work-group.
#define UPML_GROUP_SIZE 256   // work-items per work-group: cells of an x plane of the region
struct SYCL_UPMLGroup
{
	UPMLRegion R;
	UPMLFusedArgs volt, curr;
	unsigned int first;    //!< the first work-group of the region
	unsigned int blocks;   //!< work-groups per x plane of the region
};

// the region (the same for the work-group) and the cell of the work-item; returns past the region's cells
#define UPML_GROUP_THREAD \
	const unsigned int g = it.get_group(2); \
	unsigned int r = 0; \
	for (unsigned int k=1; k<count; ++k) \
		if (g>=G[k].first) \
			r = k; \
	const UPMLRegion R = G[r].R; \
	const unsigned int l = g - G[r].first; \
	const unsigned int x = l / G[r].blocks; \
	const unsigned int t = (l % G[r].blocks)*UPML_GROUP_SIZE + it.get_local_id(2); \
	if (t>=R.lz*R.ly) \
		return; \
	const unsigned int z = t % R.lz; \
	const unsigned int y = t / R.lz;

template<bool FACTOR>
static inline void upml_group_volt(const sycl::nd_item<3>& it, float* volt_out, const float* volt_in, const float* curr,
                                   const SYCL_UPMLGroup* G, unsigned int count, SYCL_GridDim N)
{
	UPML_GROUP_THREAD
	const UPMLFusedArgs A = G[r].volt;
	upml_fused_volt_node<FACTOR>(volt_out, volt_in, curr, A, N, R, x, y, z);
}

template<bool FACTOR>
static inline void upml_group_curr(const sycl::nd_item<3>& it, float* curr_out, const float* curr_in, const float* volt,
                                   const SYCL_UPMLGroup* G, unsigned int count, SYCL_GridDim N)
{
	UPML_GROUP_THREAD
	const UPMLFusedArgs A = G[r].curr;
	upml_fused_curr_node<FACTOR>(curr_out, curr_in, volt, A, N, R, x, y, z);
}

class SYCL_Ext_UPML : public GPU_Extension
{
public:
	SYCL_Ext_UPML(GPU_Backend_SYCL::Impl* impl, Operator_Ext_UPML* op_ext, Engine* eng);
	virtual ~SYCL_Ext_UPML();

	// fused: nothing before the main update, the region update after it (unless done by the main kernels);
	// fused step (see update_fused): the voltages of the region before it (the main nodes next to the
	// region need them), the currents after the voltage extensions, into the next buffers; a z slab
	// in the kernel only alternates its voltage flux buffers
	virtual void DoPreVoltageUpdates()
	{
		if (!Fused())
			Pre(d->volt, VoltArgs(), 0);
		else if ((d->fused_step>0) && m_InMain)
		{
			// in the kernel: the current voltage flux in, the new one into the second buffer,
			// as work-groups that only need a border voltage compute the node again on the old flux
			UPMLFusedArgs& A = (m_ZSlab<0) ? d->zslabs.volt_lo : d->zslabs.volt_hi;
			A.flux = m_VoltFlux;
			A.flux_out = m_VoltFlux2;
			std::swap(m_VoltFlux, m_VoltFlux2);
		}
		else if ((d->fused_step>0) && m_GroupCount)
			FuseGroup(false);
		else if ((d->fused_step>0) && !m_Grouped)
			FuseVolt(d->volt_next);
	}
	virtual void DoPostVoltageUpdates()
	{
		if (Fused() && (d->fused_step<=0) && !m_InMain)
			FuseVolt(d->volt);
		else if (!Fused())
			Post(d->volt, VoltArgs(), 0);
	}
	virtual void DoPreCurrentUpdates()
	{
		if (!Fused())
			Pre(d->curr, CurrArgs(), 9);
		else if ((d->fused_step>0) && m_GroupCount)
			FuseGroup(true);
		else if ((d->fused_step>0) && !m_InMain && !m_Grouped)
			FuseCurr(d->curr_next);
	}
	virtual void DoPostCurrentUpdates()
	{
		if (Fused() && (d->fused_step<=0) && !m_InMain)
			FuseCurr(d->curr);
		else if (!Fused())
			Post(d->curr, CurrArgs(), 9);
	}

protected:
	// set_offset: the voltage (0) or current (9) coefficients, see upml_coeff()
	void Pre(float* field, UPMLFusedArgs A, unsigned int set_offset)
	{
		const SYCL_GridDim N = d->dim;
		const UPMLRegion R = m_Region;
		SYCL_Launch(d, "upml_pre", R.lz*R.ly, 1, R.lx, [=](sycl::nd_item<3> it)
		{
			upml_pre(it, field, A, set_offset, N, R);
		});
	}
	void Post(float* field, UPMLFusedArgs A, unsigned int set_offset)
	{
		const SYCL_GridDim N = d->dim;
		const UPMLRegion R = m_Region;
		SYCL_Launch(d, "upml_post", R.lz*R.ly, 1, R.lx, [=](sycl::nd_item<3> it)
		{
			upml_post(it, field, A, set_offset, N, R);
		});
	}
	UPMLFusedArgs VoltArgs() {return Args(m_VoltFlux, m_VV, m_VVFO, m_VVFN, d->vv, d->vi);}
	UPMLFusedArgs CurrArgs() {return Args(m_CurrFlux, m_II, m_IIFO, m_IIFN, d->ii, d->iv);}
	// in place, or in the fused step from the current buffers into the next ones (see update_fused)
	void FuseVolt(float* volt_out)
	{
		const float* volt_in = d->volt;
		const float* curr = d->curr;
		const UPMLFusedArgs A = VoltArgs();
		const SYCL_GridDim N = d->dim;
		const UPMLRegion R = m_Region;
		if (d->coeff_mode>=3)
			SYCL_Launch(d, "upml_fused_volt", R.lz*R.ly, 1, R.lx, [=](sycl::nd_item<3> it)
			{
				upml_fused_volt<true>(it, volt_out, volt_in, curr, A, N, R);
			});
		else
			SYCL_Launch(d, "upml_fused_volt", R.lz*R.ly, 1, R.lx, [=](sycl::nd_item<3> it)
			{
				upml_fused_volt<false>(it, volt_out, volt_in, curr, A, N, R);
			});
	}
	void FuseCurr(float* curr_out)
	{
		const float* curr_in = d->curr;
		const float* volt = d->volt;
		const UPMLFusedArgs A = CurrArgs();
		const SYCL_GridDim N = d->dim;
		const UPMLRegion R = m_Region;
		if (d->coeff_mode>=3)
			SYCL_Launch(d, "upml_fused_curr", R.lz*R.ly, 1, R.lx, [=](sycl::nd_item<3> it)
			{
				upml_fused_curr<true>(it, curr_out, curr_in, volt, A, N, R);
			});
		else
			SYCL_Launch(d, "upml_fused_curr", R.lz*R.ly, 1, R.lx, [=](sycl::nd_item<3> it)
			{
				upml_fused_curr<false>(it, curr_out, curr_in, volt, A, N, R);
			});
	}
	// the fused step: the regions of the group (see SYCL_UPMLGroup), into the next buffers
	void FuseGroup(bool currents)
	{
		const float* volt = d->volt;
		const float* curr = d->curr;
		float* out = currents ? d->curr_next : d->volt_next;
		const SYCL_UPMLGroup* G = m_Group;
		const unsigned int count = m_GroupCount;
		const SYCL_GridDim N = d->dim;
		const bool factor = d->coeff_mode>=3;
		if (currents && factor)
			SYCL_Launch(d, "upml_group_curr", (size_t)m_GroupBlocks*UPML_GROUP_SIZE, 1, 1, [=](sycl::nd_item<3> it)
			{
				upml_group_curr<true>(it, out, curr, volt, G, count, N);
			});
		else if (currents)
			SYCL_Launch(d, "upml_group_curr", (size_t)m_GroupBlocks*UPML_GROUP_SIZE, 1, 1, [=](sycl::nd_item<3> it)
			{
				upml_group_curr<false>(it, out, curr, volt, G, count, N);
			});
		else if (factor)
			SYCL_Launch(d, "upml_group_volt", (size_t)m_GroupBlocks*UPML_GROUP_SIZE, 1, 1, [=](sycl::nd_item<3> it)
			{
				upml_group_volt<true>(it, out, volt, curr, G, count, N);
			});
		else
			SYCL_Launch(d, "upml_group_volt", (size_t)m_GroupBlocks*UPML_GROUP_SIZE, 1, 1, [=](sycl::nd_item<3> it)
			{
				upml_group_volt<false>(it, out, volt, curr, G, count, N);
			});
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

	GPU_Backend_SYCL::Impl* d;
	Engine* m_Eng;
	UPMLRegion m_Region;
	bool m_InMain;   //!< fused: updated by the main kernels (see SYCL_ZSlabs)
	int m_ZSlab;     //!< -1: slab lo, 1: slab hi, 0: not a z slab
	float* m_VoltFlux2;   //!< fused step: the second buffer of the voltage flux of a z slab, written by the work-group that owns the node
	bool m_Grouped;              //!< fused step: updated by the group launch (see SYCL_UPMLGroup) ...
	SYCL_UPMLGroup* m_Group;     //!< ... of this extension, the first of the group: the device table,
	unsigned int m_GroupCount;   //!< its number of regions (0: not the first)
	unsigned int m_GroupBlocks;  //!< and of work-groups

	float *m_VoltFlux, *m_CurrFlux;
	float *m_VV, *m_VVFO, *m_VVFN;   //!< mode 0 only, else NULL
	float *m_II, *m_IIFO, *m_IIFN;
	//! UPML coefficients, 0: the full arrays, 1/2: a 16/32 bit set index per cell and the sets (see upml_coeff())
	unsigned int m_CoeffMode;
	void* m_SetIndex;
	float* m_Sets;
};

SYCL_Ext_UPML::SYCL_Ext_UPML(GPU_Backend_SYCL::Impl* impl, Operator_Ext_UPML* op_ext, Engine* eng)
{
	d = impl;
	m_Eng = eng;
	m_InMain = false;
	m_ZSlab = 0;
	m_VoltFlux2 = NULL;
	m_Grouped = false;
	m_Group = NULL;
	m_GroupCount = m_GroupBlocks = 0;
	m_Region.sx = op_ext->m_StartPos[0];
	m_Region.sy = op_ext->m_StartPos[1];
	m_Region.sz = op_ext->m_StartPos[2];
	m_Region.lx = op_ext->m_numLines[0];
	m_Region.ly = op_ext->m_numLines[1];
	m_Region.lz = op_ext->m_numLines[2];

	const size_t cells = (size_t)m_Region.lx*m_Region.ly*m_Region.lz;
	m_VoltFlux = d->Alloc<float>(3*cells);
	m_CurrFlux = d->Alloc<float>(3*cells);

	// the 18 coefficients of a cell: vv[3], vvfo[3], vvfn[3], ii[3], iifo[3], iifn[3]; as sets,
	// or as the full arrays if nearly every cell has its own set (e.g. a graded material)
	const FDTD_FLOAT* src[6] = {op_ext->vv.data(), op_ext->vvfo.data(), op_ext->vvfn.data(),
	                            op_ext->ii.data(), op_ext->iifo.data(), op_ext->iifn.data()};
	GPU_CoeffSets sets;
	m_CoeffMode = 0;
	m_SetIndex = NULL;
	m_Sets = NULL;
	m_VV = m_VVFO = m_VVFN = NULL;
	m_II = m_IIFO = m_IIFN = NULL;
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
	else
	{
		m_VV   = Upload(op_ext->vv);
		m_VVFO = Upload(op_ext->vvfo);
		m_VVFN = Upload(op_ext->vvfn);
		m_II   = Upload(op_ext->ii);
		m_IIFO = Upload(op_ext->iifo);
		m_IIFN = Upload(op_ext->iifn);
		if (g_settings.GetVerboseLevel()>1)
			std::cout << "GPU_Backend_SYCL: UPML region of " << cells << " cells: full coefficient arrays" << std::endl;
	}

	// the UPML extensions of a grid decide together whether they fuse (see CanFuse()); only once
	// complete: a constructor that throws (out of device memory) must not stay in the list
	d->upml.push_back(this);
}

SYCL_Ext_UPML::~SYCL_Ext_UPML()
{
	d->upml.erase(std::remove(d->upml.begin(), d->upml.end(), this), d->upml.end());
	// the main updates cover the whole grid again
	d->upml_fused = -1;
	d->main_start.nx = d->main_start.ny = d->main_start.nz = 0;
	d->main_stop = d->dim;
	d->ResetZSlabs();
	for (size_t r=0; r<d->upml.size(); ++r)
	{
		SYCL_Ext_UPML* u = static_cast<SYCL_Ext_UPML*>(d->upml[r]);
		u->m_InMain = false;
		u->m_ZSlab = 0;
		u->m_Grouped = false;
		u->m_GroupCount = 0;
	}
}

bool SYCL_Ext_UPML::Fused()
{
	if (d->upml_fused<0)
	{
		// OPENEMS_SYCL_UPML_FUSION=0: the separate UPML kernels, which the usual region layout never needs
		const char* env = getenv("OPENEMS_SYCL_UPML_FUSION");
		d->upml_fused = (env && (atoi(env)==0)) ? 0 : CanFuse();
		if (d->upml_fused && (g_settings.GetVerboseLevel()>1))
			std::cout << "GPU_Backend_SYCL: " << d->upml.size() << " UPML regions fused with the main updates ("
			          << d->zslab_lo+d->zslab_hi << " along z in the main " << ((d->fused_step>0) ? "kernel" : "kernels") << ")" << std::endl;
		if (!d->upml_fused)
			d->fused_blockers++;   // the fused step needs the fused UPML
	}
	return d->upml_fused>0;
}

// The fused kernels replace the UPML hooks and the main update of the UPML regions,
// see GPU_UPMLFusionBox() for the conditions.
bool SYCL_Ext_UPML::CanFuse()
{
	std::vector<GPU_UPMLRegion> regions;
	for (size_t r=0; r<d->upml.size(); ++r)
	{
		const SYCL_Ext_UPML* u = static_cast<const SYCL_Ext_UPML*>(d->upml[r]);
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
	// (e.g. 9 nodes), which use the memory transactions poorly. The fused step can leave them
	// to the region kernels (see DecideFusedStep()).
	for (size_t r=0; r<d->upml.size(); ++r)
	{
		SYCL_Ext_UPML* u = static_cast<SYCL_Ext_UPML*>(d->upml[r]);
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
			u->m_ZSlab = -1;
		}
		else if (!d->zslab_hi && (R.sz==stop[2]))
		{
			d->zslab_hi = true;
			d->zslabs.R_hi = R;
			d->zslabs.volt_hi = u->VoltArgs();
			d->zslabs.curr_hi = u->CurrArgs();
			u->m_InMain = true;
			u->m_ZSlab = 1;
		}
	}

	// In the fused step a z slab node is computed by several work-groups: the voltage flux goes
	// into a second buffer (optional, as the fused step: without it, its region kernels update it)
	for (size_t r=0; r<d->upml.size(); ++r)
	{
		SYCL_Ext_UPML* u = static_cast<SYCL_Ext_UPML*>(d->upml[r]);
		if (!u->m_ZSlab)
			continue;
		const size_t n = 3*(size_t)u->m_Region.lx*u->m_Region.ly*u->m_Region.lz;
		u->m_VoltFlux2 = d->TryAlloc<float>(n);
		if (u->m_VoltFlux2)
			SYCL_Check("UPML flux", [&] {d->Queue().memcpy(u->m_VoltFlux2, u->m_VoltFlux, n*sizeof(float));});
		else if (u->m_ZSlab<0)
			d->zslab_lo = false;
		else
			d->zslab_hi = false;
	}
	const bool fused_step = d->DecideFusedStep();
	for (size_t r=0; r<d->upml.size(); ++r)
	{
		SYCL_Ext_UPML* u = static_cast<SYCL_Ext_UPML*>(d->upml[r]);
		// the decision may leave a slab to its region kernels
		if (u->m_ZSlab && !((u->m_ZSlab<0) ? d->zslab_lo : d->zslab_hi))
		{
			u->m_InMain = false;
			u->m_ZSlab = 0;
		}
		if (!fused_step || !u->m_ZSlab)
		{
			d->Free(u->m_VoltFlux2);
			u->m_VoltFlux2 = NULL;
		}
	}

	// the fused step: the other regions in one launch per half-step (optional: else one each)
	std::vector<SYCL_UPMLGroup> group;
	unsigned int blocks = 0;
	for (size_t r=0; fused_step && (r<d->upml.size()); ++r)
	{
		SYCL_Ext_UPML* u = static_cast<SYCL_Ext_UPML*>(d->upml[r]);
		if (u->m_InMain)
			continue;
		SYCL_UPMLGroup e;
		e.R = u->m_Region;
		e.volt = u->VoltArgs();
		e.curr = u->CurrArgs();
		e.first = blocks;
		e.blocks = (e.R.lz*e.R.ly + UPML_GROUP_SIZE-1)/UPML_GROUP_SIZE;
		blocks += e.blocks*e.R.lx;
		group.push_back(e);
	}
	if (group.size()>1)
	{
		SYCL_UPMLGroup* table = d->TryAlloc<SYCL_UPMLGroup>(group.size());
		if (table)
		{
			SYCL_Check("UPML group", [&] {d->Queue().memcpy(table, group.data(), group.size()*sizeof(SYCL_UPMLGroup)).wait();});
			bool first = true;   // the first one launches the group
			for (size_t r=0; r<d->upml.size(); ++r)
			{
				SYCL_Ext_UPML* u = static_cast<SYCL_Ext_UPML*>(d->upml[r]);
				u->m_Grouped = !u->m_InMain;
				if (u->m_Grouped && first)
				{
					u->m_Group = table;
					u->m_GroupCount = group.size();
					u->m_GroupBlocks = blocks;
					first = false;
				}
			}
		}
	}
	return true;
}

GPU_Extension* SYCL_CreateExt_UPML(GPU_Backend_SYCL::Impl* d, Engine_Extension* eng_ext, Engine* eng)
{
	if (!dynamic_cast<Engine_Ext_UPML*>(eng_ext))
		return NULL;
	Operator_Ext_UPML* op_ext = dynamic_cast<Operator_Ext_UPML*>(eng_ext->GetOperatorExtension());
	if (!op_ext)
		return NULL;
	return new SYCL_Ext_UPML(d, op_ext, eng);
}
