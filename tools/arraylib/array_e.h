/*
 * Copyright (C) 2025 Thorsten Liebig <thorsten.liebig@gmx.de>
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

#ifndef ARRAYLIB_ARRAY_E_H
#define ARRAYLIB_ARRAY_E_H

#include "impl/array_base.h"
#include "array_nijk.h"

namespace ArrayLib
{
	template <typename T>
	class ArrayENG;
};

// This is an 4D array specific for the FDTD engine & operator.
// Do not rely on a specific data ordering inside the buffer.
// Use ArrayNIJK for propper N-I-J-K data ordering
template <typename T>
class ArrayLib::ArrayENG :
	public ArrayNIJK<T>
{
public:
	using Base = ArrayNIJK<T>;
	using Base::operator();
	ArrayENG() {}

	// standard RAII initialization
	ArrayENG(std::string name, std::array<uint32_t, 3> extent)
	{
		Init(name, extent);
	}

	ArrayENG(std::string name, uint32_t extent[3])
	{
		Init(name, {extent[0], extent[1], extent[2]});
	}

	virtual void Init(std::string name, std::array<uint32_t, 3> extent)
	{
		Base::Init(name, extent);
		// I-J-K-N ordering
		this->m_stride[0] = 1;
		this->m_stride[1] = this->m_extent[2] * this->m_extent[3] * this->m_extent[0];
		this->m_stride[2] = this->m_extent[3] * this->m_extent[0];
		this->m_stride[3] = this->m_extent[0];

		// // N-I-J-K ordering
#ifdef OPENEMS_ENGINE_NIJK
		this->m_stride[0] = this->m_extent[1] * this->m_extent[2] * this->m_extent[3];
		this->m_stride[1] = this->m_extent[2] * this->m_extent[3];
		this->m_stride[2] = this->m_extent[3];
		this->m_stride[3] = 1;
#endif
	}
};

#endif // ARRAYLIB_ARRAY_E_H
