/*
*	Copyright (C) 2026 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#include "gpu_cell_model.h"
#include "cell_model.h"

#include <unordered_map>
#include <cstring>
#include <stdexcept>
#include <algorithm>

#include <cstddef>

// the kernels read the classes of the model as they are (kernels/cell_update.h)
static_assert(CM_ZERO==(int)CellModel_Set::ZERO && CM_LOSSLESS==(int)CellModel_Set::LOSSLESS && CM_LOSSY==(int)CellModel_Set::LOSSY
              && CM_UPML==(int)CellModel_Set::UPML && CM_DISP==(int)CellModel_Set::DISP
              && CM_UPML_PROFILE==(int)CellModel_Set::UPML_PROFILE,
              "CM_* and CellModel_Set::Class differ");
// the layout of the shared structs on 64 bit hosts, which the Metal kernels share (8 byte pointers
// there): a change must hold for MSL too, so it shows here
constexpr bool CM_HOST64 = (sizeof(void*)==8);
static_assert(!CM_HOST64 || (sizeof(CM_Comp)==24 && offsetof(CM_Comp, p0)==16 && offsetof(CM_Comp, state)==20
                             && offsetof(CM_Comp, nd)==23), "CM_Comp: a layout the MSL side does not share");
static_assert(!CM_HOST64 || (sizeof(CM_Combo)==100 && offsetof(CM_Combo, vv)==4 && offsetof(CM_Combo, m)==16
                             && offsetof(CM_Combo, c)==28), "CM_Combo: a layout the MSL side does not share");
static_assert(!CM_HOST64 || (sizeof(CM_Field)==128 && offsetof(CM_Field, numSlots)==24 && offsetof(CM_Field, gxy)==80
                             && offsetof(CM_Field, gz)==104),
              "CM_Field: a layout the MSL side does not share");

namespace
{
//! the floats of state a component of class \a s needs per node
unsigned int StateFloats(const CellModel_Set& s, const std::vector<CellModel_Pole>& poles)
{
	switch (s.cls)
	{
	case CellModel_Set::UPML:
	case CellModel_Set::UPML_PROFILE:
		return 1;
	case CellModel_Set::DISP:
	{
		unsigned int f = s.nd;
		for (unsigned int p=0; p<s.nl; ++p)
			f += (poles[s.offset+p].k_lor!=0) ? 2 : 1;
		return f;
	}
	default:
		return 0;
	}
}

//! a combination by its set ids and run offsets: equal keys give equal combinations
struct ComboKey
{
	uint32_t set[3];
	int32_t p0[3];
	bool operator==(const ComboKey& o) const {return std::memcmp(this, &o, sizeof(ComboKey))==0;}
};
struct ComboHash
{
	size_t operator()(const ComboKey& k) const
	{
		uint64_t h = 1469598103934665603ull;
		const uint32_t* w = (const uint32_t*)&k;
		for (size_t n=0; n<sizeof(ComboKey)/4; ++n)
			h = (h ^ w[n])*1099511628211ull;
		return (size_t)h;
	}
};
}

GPU_CellModel::GPU_CellModel(const CellModel& model, const unsigned int lines[3])
{
	for (int a=0; a<3; ++a)
		numLines[a] = lines[a];
	const size_t nodes = (size_t)numLines[0]*numLines[1]*numLines[2];
	const unsigned int nx = numLines[0], ny = numLines[1], nz = numLines[2];

	for (int f=0; f<2; ++f)
	{
		Field& F = field[f];
		// the global pole and profile offsets of each component
		const CellModel_Component* comp[3];
		uint32_t pole_base[3], prof_base[3];
		for (int n=0; n<3; ++n)
		{
			const CellModel_Component& c = model.Component(f, n);
			if (c.first.size()!=(size_t)nx*ny+1)
				throw std::runtime_error("GPU_CellModel: the cell model is not built along z (split x, mid y)");
			comp[n] = &c;
			pole_base[n] = (uint32_t)F.poles.size();
			for (const CellModel_Pole& p : c.poles)
				F.poles.push_back({p.k_int, p.k_ext, p.k_lor});
			prof_base[n] = (uint32_t)F.profile[0].size();
			for (int q=0; q<5; ++q)
				F.profile[q].insert(F.profile[q].end(), c.profile[q].begin(), c.profile[q].end());
			const std::vector<double>& gx = model.G(f, n, 0);
			const std::vector<double>& gy = model.G(f, n, 1);
			F.gxy[n].resize((size_t)nx*ny);
			for (unsigned int x=0; x<nx; ++x)
				for (unsigned int y=0; y<ny; ++y)
					F.gxy[n][(size_t)x*ny + y] = (float)(gx[x]*gy[y]);
			const std::vector<double>& gz = model.G(f, n, 2);
			F.gz[n].assign(gz.begin(), gz.end());
		}
		// the offsets of the combinations (CM_Comp::p0) are 32 bit with sign
		if ((F.profile[0].size() > INT32_MAX) || (F.poles.size() > INT32_MAX))
			throw std::length_error("GPU_CellModel: more than 2^31 cells of UPML profiles or dispersive poles in a field");

		// per line, the segments along z on which none of the three components changes its run: one
		// combination each (no per node copy of the runs). The slots are numbered in node order.
		std::unordered_map<ComboKey, uint32_t, ComboHash> ids;
		std::vector<unsigned int> slot_floats;   // per combination
		F.index.resize(nodes);
		for (size_t line=0; line<(size_t)nx*ny; ++line)
		{
			uint64_t r[3];
			for (int n=0; n<3; ++n)
				r[n] = comp[n]->first[line];
			for (uint32_t z=0; z<nz; )
			{
				// the run of each component at z (or none), and where the segment ends
				ComboKey key;
				uint32_t z1 = nz;
				for (int n=0; n<3; ++n)
				{
					const CellModel_Component& c = *comp[n];
					while ((r[n]<c.first[line+1]) && (c.runs[r[n]].stop<=z))
						++r[n];
					key.set[n] = UINT32_MAX;
					key.p0[n] = 0;
					if (r[n]==c.first[line+1])
						continue;
					const CellModel_Run& run = c.runs[r[n]];
					if (run.start>z)   // a gap before the next run: no update (CM_SKIP)
					{
						z1 = std::min(z1, run.start);
						continue;
					}
					z1 = std::min(z1, run.stop);
					const CellModel_Set& set = c.sets[run.set];
					key.set[n] = run.set;
					key.p0[n] = (set.cls==CellModel_Set::UPML_PROFILE) ? (int32_t)(prof_base[n] + set.offset) - (int32_t)run.start
					          : ((set.cls==CellModel_Set::DISP) ? (int32_t)(pole_base[n] + set.offset) : 0);
				}

				uint32_t id;
				auto it = ids.find(key);
				if (it!=ids.end())
					id = it->second;
				else
				{
					CM_Combo cb;
					std::memset(&cb, 0, sizeof(cb));
					unsigned int floats = 0;
					for (int n=0; n<3; ++n)
					{
						CM_Comp& cc = cb.c[n];
						if (key.set[n]==UINT32_MAX)
						{
							cb.cls[n] = CM_SKIP;
							continue;
						}
						const CellModel_Set& set = comp[n]->sets[key.set[n]];
						cb.cls[n] = (uint8_t)set.cls;
						cb.m[n] = (float)set.m;
						cb.vv[n] = set.vv;
						cc.uvv = set.uvv;
						cc.ufo = set.ufo;
						cc.ufn = set.ufn;
						cc.solve = set.solve;
						cc.nl = set.nl;
						cc.nd = set.nd;
						cc.p0 = key.p0[n];
						cc.state = (uint16_t)floats;
						floats += StateFloats(set, comp[n]->poles);
					}
					cb.state = (floats>0) ? CM_STATE : 0;
					if ((cb.cls[0]==CM_LOSSLESS) && (cb.cls[1]==CM_LOSSLESS) && (cb.cls[2]==CM_LOSSLESS))
						cb.state |= CM_LOSSLESS3;
					id = (uint32_t)F.combos.size();
					ids.emplace(key, id);
					F.combos.push_back(cb);
					slot_floats.push_back(floats);
					F.slotFloats = std::max(F.slotFloats, floats);
				}

				const size_t i0 = line*nz;
				std::fill(F.index.begin() + i0 + z, F.index.begin() + i0 + z1, id);
				if (slot_floats[id]>0)
				{
					// the per node slots, allocated with the first node that has state (none: no slot array)
					if (F.slot.empty())
						F.slot.assign(nodes, UINT32_MAX);
					// a slot of 32 bit per node (64 would cost 4 bytes more per node of the grid): UINT32_MAX is none
					if (F.numSlots + (z1-z) > UINT32_MAX)
						throw std::length_error("GPU_CellModel: more than 2^32-1 nodes with state (UPML, dispersive materials) in a field");
					for (uint32_t k=z; k<z1; ++k)
						F.slot[i0+k] = (uint32_t)F.numSlots++;
				}
				z = z1;
			}
		}
	}
}

double GPU_CellModel::NodeBytes() const
{
	const double nodes = (double)numLines[0]*numLines[1]*numLines[2];
	double bytes = 6*sizeof(float);
	for (int f=0; f<2; ++f)
	{
		const Field& F = field[f];
		bytes += Index16() ? sizeof(uint16_t) : sizeof(uint32_t);
		if (!F.slot.empty())
			bytes += sizeof(uint32_t) + (double)F.numSlots*F.slotFloats*sizeof(float)/nodes;
	}
	return bytes;
}
