/*
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

// Row-pointer access to the UPML coefficient and flux arrays. Every row base
// and the z increment come from the linearIndex() and stride() of that array.
#ifndef ENGINE_EXT_UPML_SSE_ROWS_H
#define ENGINE_EXT_UPML_SSE_ROWS_H

#include "tools/arraylib/array_nijk.h"

#include <cstddef>
#include <cstdint>

template <typename T, typename IndexType, std::size_t ExtentN>
inline bool UpmlSSERowsHaveBounds(
	const ArrayLib::ArrayNIJK<T, IndexType, ExtentN>& array,
	std::size_t xStart,
	std::size_t xCount,
	std::size_t yCount,
	std::size_t zCount
)
{
	const std::size_t xExtent = static_cast<std::size_t>(array.extent(1));
	return array.extent(0) >= 3 &&
		xStart <= xExtent && xCount <= xExtent - xStart &&
		yCount <= static_cast<std::size_t>(array.extent(2)) &&
		zCount <= static_cast<std::size_t>(array.extent(3)) &&
		((xCount == 0 || yCount == 0 || zCount == 0) || array.data() != NULL);
}

template <typename T, typename IndexType = std::uint32_t, std::size_t ExtentN = 3>
class UpmlSSENIJKRow3
{
public:
	void BeginRow(
		ArrayLib::ArrayNIJK<T, IndexType, ExtentN>& array,
		IndexType x,
		IndexType y
	)
	{
		T* data = array.data();
		m_component[0] = data + array.linearIndex({0, x, y, 0});
		m_component[1] = data + array.linearIndex({1, x, y, 0});
		m_component[2] = data + array.linearIndex({2, x, y, 0});
		m_zStride = static_cast<std::size_t>(array.stride(3));
	}

	T& Component(std::size_t component) const
	{
		return *m_component[component];
	}

	void Advance()
	{
		m_component[0] += m_zStride;
		m_component[1] += m_zStride;
		m_component[2] += m_zStride;
	}

private:
	T* m_component[3] = {NULL, NULL, NULL};
	std::size_t m_zStride = 0;
};

#endif // ENGINE_EXT_UPML_SSE_ROWS_H
