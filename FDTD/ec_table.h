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

#ifndef EC_TABLE_H
#define EC_TABLE_H

#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <unordered_set>
#include <vector>

//! The distinct effective materials of the nodes of an operator, each stored once (Operator::GetEC())
/*!
  An entry holds the effective material of a node for all three directions, the
  circuit value of an edge is its material times the geometry of the edge. An
  extension that sets a circuit value itself (a lumped element, a conducting
  sheet) gives the node an entry with that value and its bit in abs_mask set.

  Entry 0 is all zero: no material, e.g. the nodes of a multi-grid operator left to
  its sub-grid. Intern() returns the index of an entry, adding it if it is new.
  */
class EC_Table
{
public:
	struct Entry
	{
		double v[12];        //!< [n*4 + k]: epsilon, kappa, mue, sigma (k: Operator::ECValue) of direction n, or the circuit value itself (abs_mask)
		uint16_t abs_mask;   //!< bit n*4 + k: v[n*4 + k] is the circuit value itself, not to be multiplied with the geometry
	};

	EC_Table() : m_Lookup(64, Hash(&m_Entries), Equal(&m_Entries))
	{
		Entry zero;
		std::memset(&zero, 0, sizeof(zero));
		Intern(zero);
	}
	// the lookup refers to m_Entries of this object
	EC_Table(const EC_Table&) = delete;
	EC_Table& operator=(const EC_Table&) = delete;

	//! The index of \a e, added if it is new
	uint32_t Intern(const Entry& e)
	{
		// the operator goes dense long before (Operator::GetECTableLimit()): never a wrapped index
		if (m_Entries.size()>=UINT32_MAX)
			throw std::length_error("EC_Table::Intern: more entries than a 32 bit index can address");
		m_Entries.push_back(e);
		const uint32_t id = (uint32_t)(m_Entries.size()-1);
		auto r = m_Lookup.insert(id);
		if (!r.second)
		{
			m_Entries.pop_back();
			return *r.first;
		}
		return id;
	}

	size_t Size() const {return m_Entries.size();}
	const Entry& operator[](size_t i) const {return m_Entries[i];}
	//! The memory held, in bytes, about
	size_t Bytes() const {return m_Entries.capacity()*sizeof(Entry) + m_Lookup.size()*(sizeof(void*)+sizeof(uint32_t)+sizeof(size_t)) + m_Lookup.bucket_count()*sizeof(void*);}

private:
	// the values of an entry, field by field: its padding is not initialised
	static bool Same(const Entry& a, const Entry& b) {return (a.abs_mask==b.abs_mask) && (std::memcmp(a.v, b.v, sizeof(a.v))==0);}
	struct Hash
	{
		const std::vector<Entry>* entries;
		explicit Hash(const std::vector<Entry>* e) : entries(e) {}
		size_t operator()(uint32_t id) const
		{
			const Entry& e = (*entries)[id];
			uint64_t h = 0x9e3779b97f4a7c15ull ^ e.abs_mask;
			for (int k=0; k<12; ++k)
			{
				uint64_t b;
				std::memcpy(&b, &e.v[k], sizeof(b));
				h ^= b + 0x9e3779b97f4a7c15ull + (h<<6) + (h>>2);
				h *= 0xbf58476d1ce4e5b9ull;
				h ^= h>>31;
			}
			return (size_t)h;
		}
	};
	struct Equal
	{
		const std::vector<Entry>* entries;
		explicit Equal(const std::vector<Entry>* e) : entries(e) {}
		bool operator()(uint32_t a, uint32_t b) const {return Same((*entries)[a], (*entries)[b]);}
	};

	std::vector<Entry> m_Entries;
	std::unordered_set<uint32_t, Hash, Equal> m_Lookup;
};

#endif // EC_TABLE_H
