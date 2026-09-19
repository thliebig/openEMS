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
#include <cstring>
#include <limits>
#include <unordered_map>

#include "gpu_coeff_sets.h"

namespace
{
//! bit pattern of the values of a set, see GPU_FindSets()
struct SetKey
{
	uint32_t v[GPU_MAX_SET_WIDTH];
	unsigned int width;
	bool operator==(const SetKey& o) const {return std::memcmp(v, o.v, width*sizeof(uint32_t))==0;}
};

struct SetKeyHash
{
	size_t operator()(const SetKey& s) const
	{
		uint64_t h = 1469598103934665603ULL;   // FNV-1a
		for (unsigned int k=0; k<s.width; ++k)
			h = (h ^ s.v[k]) * 1099511628211ULL;
		return h;
	}
};
}

bool GPU_FindSets(size_t count, unsigned int width, const std::function<void(size_t, float*)>& get, GPU_CoeffSets& sets)
{
	if ((width==0) || (width>GPU_MAX_SET_WIDTH))
		return false;
	// at most half the size of the full arrays (width floats per item)
	const size_t max_sets16 = (size_t)std::numeric_limits<uint16_t>::max()+1;
	const size_t max_sets = std::max(max_sets16, (2*width*count - 4*count)/(4*width));
	sets.index.resize(count);
	sets.table.clear();
	std::unordered_map<SetKey, uint32_t, SetKeyHash> found;

	SetKey prev;
	bool have_prev = false;
	float values[GPU_MAX_SET_WIDTH];
	for (size_t i=0; i<count; ++i)
	{
		SetKey key;
		key.width = width;
		get(i, values);
		std::memcpy(key.v, values, width*sizeof(float));
		// neighbours mostly share the set
		if (have_prev && (key==prev))
		{
			sets.index[i] = sets.index[i-1];
			continue;
		}
		std::unordered_map<SetKey, uint32_t, SetKeyHash>::const_iterator it = found.find(key);
		if (it==found.end())
		{
			if (found.size()>=max_sets)
				return false;
			it = found.insert(std::make_pair(key, (uint32_t)found.size())).first;
			sets.table.insert(sets.table.end(), values, values+width);
		}
		sets.index[i] = it->second;
		prev = key;
		have_prev = true;
	}
	sets.count = found.size();
	sets.mode = (sets.count<=max_sets16) ? 1 : 2;
	return true;
}
