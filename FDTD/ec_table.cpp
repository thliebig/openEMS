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

#include "ec_table.h"

#include <cstring>
#include <stdexcept>

// file-local, so that Intern() has them inline (a member would be called through the PLT)
namespace
{

// field by field, not the bytes of the struct: no assumption about its padding
bool Same(const EC_Table::Entry& a, const EC_Table::Entry& b)
{
	return (a.abs_mask==b.abs_mask) && (a.keep_mask==b.keep_mask) && (a.base==b.base) && (std::memcmp(a.v, b.v, sizeof(a.v))==0);
}

size_t Hash(const EC_Table::Entry& e)
{
	uint64_t h = 0x9e3779b97f4a7c15ull ^ e.abs_mask ^ ((uint64_t)e.keep_mask<<16) ^ ((uint64_t)e.base<<32);
	for (int k=0; k<12; ++k)
	{
		uint64_t b;
		std::memcpy(&b, &e.v[k], sizeof(b));
		h ^= b + 0x9e3779b97f4a7c15ull + (h<<6) + (h>>2);
		h *= 0xbf58476d1ce4e5b9ull;
		h ^= h>>31;
	}
	// the slot is the low bits: the high bits of the last word down to them as well
	h ^= h>>33;
	h *= 0xff51afd7ed558ccdull;
	h ^= h>>33;
	return (size_t)h;
}

}

EC_Table::EC_Table() : m_Size(0)
{
	m_Slots.assign(64, 0);
	Entry zero;
	std::memset(&zero, 0, sizeof(zero));
	Intern(zero);
}

uint32_t EC_Table::Intern(const Entry& e)
{
	// the operator goes dense long before (Operator::GetECTableLimit()): never a wrapped index
	if (m_Size>=UINT32_MAX-1)
		throw std::length_error("EC_Table::Intern: more entries than a 32 bit index can address");
	if (2*(m_Size+1)>m_Slots.size())
		Grow();
	const size_t mask = m_Slots.size()-1;
	for (size_t i=Hash(e)&mask; ; i=(i+1)&mask)
	{
		const uint32_t slot = m_Slots[i];
		if (slot==0)
		{
			if ((m_Size & (BLOCK-1))==0)
				m_Blocks.emplace_back(new Entry[BLOCK]);
			Get(m_Size) = e;
			m_Slots[i] = (uint32_t)(++m_Size);
			return (uint32_t)(m_Size-1);
		}
		if (Same(Get(slot-1), e))
			return slot-1;
	}
}

void EC_Table::Grow()
{
	std::vector<uint32_t> slots(2*m_Slots.size(), 0);
	const size_t mask = slots.size()-1;
	for (size_t id=0; id<m_Size; ++id)
	{
		size_t i = Hash(Get(id))&mask;
		while (slots[i])
			i = (i+1)&mask;
		slots[i] = (uint32_t)(id+1);
	}
	m_Slots.swap(slots);
}
