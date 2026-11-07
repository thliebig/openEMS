/*
 * Copyright (C) 2024 Yifeng Li <tomli@tomli.me>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef ARRAYLIB_ARRAY_NIJK_H
#define ARRAYLIB_ARRAY_NIJK_H

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <array>
#include <string>
#include <stdexcept>

#include "impl/array_base.h"

namespace ArrayLib
{
	template <typename T, typename IndexType=uint32_t, size_t extentN=3>
	class ArrayNIJK;
};

template <typename T, typename IndexType, size_t extentN>
class ArrayLib::ArrayNIJK :
	public ArrayBase<ArrayNIJK<T, IndexType, extentN>, T, 4, IndexType>
{
public:
	using Base = ArrayBase<ArrayNIJK<T, IndexType>, T, 4, IndexType>;
	using Base::operator();

	// 2-phase initialization: user declares a dummy array object first,
	// and calls Init() later.
	ArrayNIJK() {}

	virtual void Init(std::string name, std::array<IndexType, 3> extent)
	{
		InitView(name, extent, NULL);
	}

	// Allocate with the memory layout \a layout, see layout().
	void Init(std::string name, std::array<IndexType, 3> extent, const std::string& layout)
	{
		InitView(name, extent, NULL, layout);
	}

	// Use the external memory at ptr (at least extentN*i*j*k elements) instead
	// of allocating, e.g. memory shared with a GPU. The array does not free it.
	// With ptr==NULL, allocate and own the memory like Init().
	void InitView(std::string name, std::array<IndexType, 3> extent, T* ptr,
	              const std::string& layout="NIJK")
	{
		std::array<int, 4> order = ParseLayout(layout);   // throws before allocating
		this->Reset();

		this->m_name = name;
		this->m_size = (size_t)extent[0] * extent[1] * extent[2] * extentN;
		this->m_bytes = sizeof(T) * this->m_size;
		this->m_owner = (ptr==NULL);
		this->m_ptr = ptr ? ptr : Base::AllocatorType::alloc(this->m_size);

		this->m_extent = {extentN, extent[0], extent[1], extent[2]};
		size_t stride = 1;
		for (int m=3; m>=0; --m)
		{
			this->m_stride[order[m]] = stride;
			stride *= this->m_extent[order[m]];
		}
		m_layout = layout;
	}

	// The order of the four axes in memory, slowest first: a permutation of
	// "NIJK", e.g. "NIJK" (the default, one block per component), "IJNK" (the
	// component rows of a line next to each other) or "IJKN" (the components
	// of a node interleaved). Only the strides change; access via operator()
	// is the same for every layout. Code that walks data() linearly must
	// check layout().
	std::string layout() const {return m_layout;}

	ArrayNIJK(std::string name, std::array<IndexType, 3> extent, const std::string& layout)
	{
		Init(name, extent, layout);
	}

	void Init(std::string name, IndexType extent[3])
	{
		Init(name, {extent[0], extent[1], extent[2]});
	}

	// standard RAII initialization
	ArrayNIJK(std::string name, std::array<IndexType, 3> extent)
	{
		Init(name, extent);
	}

	ArrayNIJK(std::string name, IndexType extent[3])
	{
		Init(name, {extent[0], extent[1], extent[2]});
	}

	size_t linearIndex(std::array<IndexType, 4> tupleIndex) const
	{
		return  this->m_stride[0] * tupleIndex[0] +
		        this->m_stride[1] * tupleIndex[1] +
		        this->m_stride[2] * tupleIndex[2] +
		        this->m_stride[3] * tupleIndex[3];
	}

	T& operator() (IndexType n, IndexType i, IndexType j, IndexType k) const
	{
		return this->m_ptr[linearIndex({n, i, j, k})];
	}

private:
	std::string m_layout = "NIJK";

	// the axis (0=N .. 3=K) at each position of layout, slowest first
	static std::array<int, 4> ParseLayout(const std::string& layout)
	{
		const std::string axes = "NIJK";
		bool seen[4] = {false, false, false, false};
		std::array<int, 4> order = {0, 1, 2, 3};
		for (size_t m=0; m<layout.size() && m<4; ++m)
		{
			size_t a = axes.find(layout[m]);
			if (a==std::string::npos || seen[a])
				break;
			seen[a] = true;
			order[m] = (int)a;
		}
		if (layout.size()!=4 || !(seen[0] && seen[1] && seen[2] && seen[3]))
			throw std::invalid_argument("ArrayNIJK: layout \"" + layout + "\" is not a permutation of NIJK");
		return order;
	}
};

#endif // ARRAYLIB_ARRAY_NIJK_H
