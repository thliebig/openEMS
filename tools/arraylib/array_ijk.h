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

#ifndef ARRAYLIB_ARRAY_IJK_H
#define ARRAYLIB_ARRAY_IJK_H

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <array>

#include "impl/array_base.h"

namespace ArrayLib
{
	template <typename T, typename IndexType=uint32_t>
	class ArrayIJK;
};

template <typename T, typename IndexType>
class ArrayLib::ArrayIJK :
	public ArrayBase<ArrayIJK<T, IndexType>, T, 3, IndexType>
{
private:
	std::array<IndexType, 3> m_stride;

public:
	using Base = ArrayBase<ArrayIJK<T, IndexType>, T, 3, IndexType>;
	using Base::operator();

	// 2-phase initialization: user declares a dummy array object first,
	// and calls Init() later. Ugly but needed because openEMS's FDTD
	// engine itself uses the 2-phase initialization pattern, we can't
	// determine simulation domain size in the constructor.
	ArrayIJK() {}

	void Init(std::string name, std::array<IndexType, 3> extent)
	{
		if (this->m_ptr != NULL)
			Base::AllocatorType::free(this->m_ptr, this->m_size);

		this->m_name = name;
		this->m_size = extent[0] * extent[1] * extent[2];
		this->m_bytes = sizeof(T) * this->m_size;
		this->m_ptr = Base::AllocatorType::alloc(this->m_size);

		this->m_extent = extent;
		this->m_stride[0] = extent[1] * extent[2];
		this->m_stride[1] = extent[2];
	}

	void Init(std::string name, IndexType extent[3])
	{
		Init(name, {extent[0], extent[1], extent[2]});
	}

	// standard RAII initialization
	ArrayIJK(std::string name, std::array<IndexType, 3> extent)
	{
		Init(name, extent);
	}

	ArrayIJK(std::string name, IndexType extent[3])
	{
		Init(name, {extent[0], extent[1], extent[2]});
	}

	IndexType linearIndex(std::array<IndexType, 3> tupleIndex) const
	{
		return m_stride[0] * tupleIndex[0] +
		       m_stride[1] * tupleIndex[1] +
		                     tupleIndex[2];
	}

	T& operator() (IndexType i, IndexType j, IndexType k) const
	{
		return this->m_ptr[linearIndex({i, j, k})];
	}

	~ArrayIJK()
	{
		Base::AllocatorType::free(this->m_ptr, this->m_size);
	}
};

#endif // ARRAYLIB_ARRAY_IJK_H
