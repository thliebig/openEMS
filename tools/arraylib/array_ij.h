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

#ifndef ARRAYLIB_ARRAY_IJ_H
#define ARRAYLIB_ARRAY_IJ_H

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <array>

#include "impl/array_base.h"

namespace ArrayLib
{
	template <typename T, typename IndexType=uint32_t>
	class ArrayIJ;
};

template <typename T, typename IndexType>
class ArrayLib::ArrayIJ :
	public ArrayBase<ArrayIJ<T, IndexType>, T, 2, IndexType>
{
private:
	IndexType m_stride;

public:
	using Base = ArrayBase<ArrayIJ<T, IndexType>, T, 2, IndexType>;
	using Base::operator();

	// 2-phase initialization: user declares a dummy array object first,
	// and calls Init() later. Ugly but needed because openEMS's FDTD
	// engine itself uses the 2-phase initialization pattern, we can't
	// determine simulation domain size in the constructor.
	ArrayIJ() {}

	void Init(std::string name, std::array<IndexType, 2> extent)
	{
		if (this->m_ptr != NULL)
			Base::AllocatorType::free(this->m_ptr, this->m_size);

		this->m_name = name;
		this->m_size = extent[0] * extent[1];
		this->m_bytes = sizeof(T) * this->m_size;
		this->m_ptr = Base::AllocatorType::alloc(this->m_size);

		this->m_extent = extent;
		this->m_stride = extent[1];
	}

	void Init(std::string name, IndexType extent[2])
	{
		Init(name, {extent[0], extent[1]});
	}

	// standard RAII initialization
	ArrayIJ(std::string name, std::array<IndexType, 2> extent)
	{
		Init(name, extent);
	}

	ArrayIJ(std::string name, IndexType extent[2])
	{
		Init(name, {extent[0], extent[1]});
	}

	IndexType linearIndex(std::array<IndexType, 2> tupleIndex) const
	{
		return m_stride * tupleIndex[0] + tupleIndex[1];
	}

	T& operator() (IndexType i, IndexType j) const
	{
		return this->m_ptr[linearIndex({i, j})];
	}

	~ArrayIJ()
	{
		Base::AllocatorType::free(this->m_ptr, this->m_size);
	}
};

#endif // ARRAYLIB_ARRAY_IJ_H
