/*
*	Copyright (C) 2026 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#ifndef ENGINE_EXT_CELLINDEX_H
#define ENGINE_EXT_CELLINDEX_H

#include <cstddef>
#include <cstdint>
#include <vector>

#include "FDTD/engine_plan.h"

//! The elements of an extension's cell list in a box, without walking the whole list
/*!
  Per axis the element indices ordered by their line on that axis (stable,
  so elements on the same cell keep their order), with the first element of
  each line. A box query walks the candidates of the axis that has fewest in
  the box. Built once, before the threads run (Engine_Extension::PrepareBoxes()),
  read-only after that.
  */
class Engine_Ext_CellIndex
{
public:
	//! The list: \a count elements at (pos[0][i], pos[1][i], pos[2][i]); ForEach() walks it in order until Build()
	void SetList(size_t count, unsigned int* const pos[3])
	{
		m_Count = (pos==NULL) ? 0 : count;   // an extension without elements may have no list at all
		for (int a=0; a<3; ++a)
			m_Pos[a] = (pos==NULL) ? NULL : pos[a];
		m_Ready = false;
	}

	//! The index by box, on a mesh of numLines
	void Build(const unsigned int numLines[3])
	{
		const size_t count = m_Count;
		unsigned int* const* pos = m_Pos;
		for (int a=0; a<3; ++a)
		{
			m_First[a].assign(numLines[a]+2, 0);
			for (size_t i=0; i<count; ++i)
				++m_First[a][pos[a][i]+1];
			for (unsigned int l=0; l<=numLines[a]; ++l)
				m_First[a][l+1] += m_First[a][l];
			m_Order[a].resize(count);
			std::vector<uint32_t> next(m_First[a].begin(), m_First[a].end()-1);
			for (size_t i=0; i<count; ++i)
				m_Order[a][next[pos[a][i]]++] = (uint32_t)i;
		}
		m_Ready = true;
	}

	bool IsReady() const {return m_Ready;}

	//! f(i) for the elements in \a box, all in list order if \a box is NULL
	template <typename F>
	void ForEach(const Engine_Box* box, F f) const
	{
		if ((box==NULL) || !m_Ready)
		{
			for (size_t i=0; i<m_Count; ++i)
				if (Engine_Box_Contains(box, m_Pos[0][i], m_Pos[1][i], m_Pos[2][i]))
					f(i);
			return;
		}
		int best = 0;
		size_t fewest = (size_t)-1;
		for (int a=0; a<3; ++a)
		{
			const size_t n = Range(a, box->start[a], box->stop[a], NULL);
			if (n<fewest)
			{
				fewest = n;
				best = a;
			}
		}
		size_t k0;
		const size_t n = Range(best, box->start[best], box->stop[best], &k0);
		for (size_t k=k0; k<k0+n; ++k)
		{
			const uint32_t i = m_Order[best][k];
			if (box->Contains(m_Pos[0][i], m_Pos[1][i], m_Pos[2][i]))
				f(i);
		}
	}

protected:
	//! the number of elements on the lines [lo, hi) of axis a, and the first of them
	size_t Range(int a, unsigned int lo, unsigned int hi, size_t* first) const
	{
		const size_t L = m_First[a].size()-2;
		if (lo>L) lo = L;
		if (hi>L) hi = L;
		if (first)
			*first = m_First[a][lo];
		return (hi>lo) ? m_First[a][hi]-m_First[a][lo] : 0;
	}

	size_t m_Count = 0;
	unsigned int* m_Pos[3] = {NULL, NULL, NULL};
	std::vector<uint32_t> m_First[3];
	std::vector<uint32_t> m_Order[3];
	bool m_Ready = false;
};

#endif // ENGINE_EXT_CELLINDEX_H
