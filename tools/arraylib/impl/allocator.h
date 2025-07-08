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

#ifndef ARRAYLIB_ALLOCATOR_H
#define ARRAYLIB_ALLOCATOR_H

#include <cstddef>
#include <cmath>
#include <cstring>

namespace ArrayLib
{
	template <typename T>
	class SimpleAllocator;

	template <typename T>
	class AlignedAllocator;
};

template <typename T>
class ArrayLib::SimpleAllocator
{
public:
	static T* alloc(size_t numelem)
	{
		T* ptr = new T[numelem]();
	}

	static void free(T* ptr, size_t numelem)
	{
		delete[] ptr;
	}
};

#ifdef WIN32
#include <malloc.h>
#endif

template <typename T>
class ArrayLib::AlignedAllocator
{
public:
	static T* alloc(size_t numelem)
	{
		size_t alignment;
		// User-defined SIMD variables can be larger than the built-in C
		// types, so they need custom alignments to their own size.
		if (sizeof(T) <= sizeof(void*))
			// POSIX says alignment >= sizeof(void*)
			alignment = sizeof(void *);
		else
			// POSIX says alignment must be a power of 2
			alignment = 1 << (size_t) ceil(log2(sizeof(T)));

		T* buf;
#ifdef WIN32
		buf = (T*) _mm_malloc(numelem * sizeof(T), alignment);
		if (buf == NULL)
		{
			std::cerr << "Failed to allocate aligned memory" << std::endl;
			throw std::bad_alloc();
		}
#else
		int retval = posix_memalign((void**) &buf, alignment, numelem * sizeof(T));
		if (retval != 0)
		{
			std::cerr << "Failed to allocate aligned memory" << std::endl;
			throw std::bad_alloc();
		}
#endif
		memset(buf, 0, numelem * sizeof(T));
		for (size_t i = 0; i < numelem; i++)
			new (buf + i) T();

		return buf;
	}

	static void free(T* ptr, size_t numelem)
	{
		if (ptr)
		{
			for (size_t i = 0; i < numelem; i++)
				(&ptr[i])->~T();
#ifdef WIN32
			_mm_free(ptr);
#else
			std::free(ptr);
#endif
		}

	}
};

#endif // ARRAYLIB_ALLOCATOR_H
