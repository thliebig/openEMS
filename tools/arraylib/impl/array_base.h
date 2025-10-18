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

#ifndef ARRAYLIB_ARRAY_BASE_H
#define ARRAYLIB_ARRAY_BASE_H

#include <string>

#include "subscript.h"
#include "allocator.h"

// ArrayBase is a base class of which 2D, 3D and 4D arrays are derived from.
// this reduces code duplication and allows minimum code in derived classes.
//
// Each derived class should implement:
//
// * call base constructor
// * m_ptr, m_size, and m_bytes (by assignment in derived constructor)
// * IndexType linearIndex(std::array<IndexType, dim> tupleIndex) const
// * T& operator() (IndexType i, IndexType j, IndexType k)
// * using Base::operator(); (so both base and derived operator() are visible)

namespace ArrayLib
{
	template <
		typename Derived,
		typename T,
		size_t rank,
		typename Index,
		typename Allocator=AlignedAllocator<T>
	>
	class ArrayBase;
};

template <
	typename Derived,
	typename T,
	size_t rank,
	typename Index,
	typename Allocator
>
class ArrayLib::ArrayBase
{
protected:
	using IndexType = Index;
	using AllocatorType = Allocator;

	std::string m_name="";
	std::array<IndexType, rank> m_extent = {};
	std::array<IndexType, rank> m_stride = {};
	IndexType m_size=0, m_bytes=0;

	T* __restrict m_ptr = NULL;

	// 2-phase initialization
	ArrayBase() {}

public:
	// Access array via arr({i, j, k}) syntax with one array of indices.
	// Each derived class should also implement operator(i, j, k), which
	// is the recommended syntax.
	T&
	operator() (std::array<IndexType, rank> tupleIndex) const
	{
		auto derived = static_cast<const Derived*>(this);
		return m_ptr[derived->linearIndex(tupleIndex)];
	}

	void Reset()
	{
		if (this->m_ptr==NULL) return;
		AllocatorType::free(this->m_ptr, this->m_size);
		this->m_name = "";
		this->m_ptr  = NULL;
		this->m_size = 0;
		for (int n=0;n<rank;++n)
		{
			this->m_extent[n] = 0;
			this->m_stride[n] = 0;
		}
	}

	// This array must always be passed via reference, not value, because
	// it's prohibitively expensive for GiB-sized simulation data. If the
	// rule is ignored, different copies would hold the same underlying
	// pointer. When RAII frees a single array within a scope, it affects
	// the entire program.
	ArrayBase (const ArrayBase&) = delete;
	ArrayBase& operator= (const ArrayBase&) = delete;
	~ArrayBase() {Reset();}

	std::string                 name()              const { return m_name;     }
	std::array<IndexType, rank> extent()            const { return m_extent;   }
	IndexType                   extent(IndexType n) const { return m_extent[n];}
	IndexType                   stride(IndexType n) const { return m_stride[n];}
	bool                        valid()             const { return m_ptr!=NULL;}

	// return the number of array elements
	IndexType                   size()              const { return m_size;     }

	// return occupied memory
	IndexType                   bytes()             const { return m_bytes;    }

	// return raw pointer to underlying data
	T*                          data()              const { return m_ptr;      }
	T&                          data(IndexType n)   const { return m_ptr[n];   }
};

#endif // ARRAYLIB_ARRAY_BASE_H
