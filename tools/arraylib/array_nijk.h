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

#include "impl/array_base.h"

namespace ArrayLib
{
	template <typename T, typename IndexType=uint32_t, size_t extentN=3>
	class ArrayNIJK;
};

template <typename T, typename IndexType, size_t extentN>
class ArrayLib::ArrayNIJK :
	public ArrayBase<ArrayNIJK<T, IndexType>, T, 4, IndexType>
{
public:
	using Base = ArrayBase<ArrayNIJK<T, IndexType>, T, 4, IndexType>;
	using Base::operator();

	// 2-phase initialization: user declares a dummy array object first,
	// and calls Init() later.
	ArrayNIJK() {}

	virtual void Init(std::string name, std::array<IndexType, 3> extent)
	{
		if (this->m_ptr != NULL)
			Base::AllocatorType::free(this->m_ptr, this->m_size);

		this->m_name = name;
		this->m_size = extent[0] * extent[1] * extent[2] * extentN;
		this->m_bytes = sizeof(T) * this->m_size;
		this->m_ptr = Base::AllocatorType::alloc(this->m_size);

		this->m_extent = {extentN, extent[0], extent[1], extent[2]};

		// N-I-J-K ordering
		this->m_stride[0] = extent[0] * extent[1] * extent[2];
		this->m_stride[1] = extent[1] * extent[2];
		this->m_stride[2] = extent[2];
		this->m_stride[3] = 1;
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

	IndexType linearIndex(std::array<IndexType, 4> tupleIndex) const
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
};

#endif // ARRAYLIB_ARRAY_NIJK_H
