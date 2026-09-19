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

#ifndef GPU_COEFF_SETS_H
#define GPU_COEFF_SETS_H

#include <cstddef>
#include <cstdint>
#include <functional>
#include <vector>

//! Coefficients stored as distinct sets and a set index per item (node or cell), see GPU_FindSets()
#define GPU_MAX_SET_WIDTH 18
struct GPU_CoeffSets
{
	uint32_t mode;                //!< 1: 16 bit index, 2: 32 bit index
	size_t count;                 //!< number of sets
	std::vector<uint32_t> index;  //!< set of each item
	std::vector<float> table;     //!< the sets, \a width values each
};

//! Find the distinct sets of \a width values of \a count items, \a get(i, values) returns the values of item i
/*!
  Returns false if a set index plus the sets would be larger than half the full
  arrays (\a width floats per item). 16 bit indices for up to 65536 sets, else 32 bit.
  The sets hold the bit patterns of the values, so reading them is exact.
  */
bool GPU_FindSets(size_t count, unsigned int width, const std::function<void(size_t, float*)>& get, GPU_CoeffSets& sets);

#endif // GPU_COEFF_SETS_H
