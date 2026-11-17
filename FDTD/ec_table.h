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
#include <memory>
#include <vector>

//! The distinct effective materials of the nodes of an operator, each stored once (Operator::GetEC())
/*!
  An entry holds the effective material of a node for all three directions, the
  circuit value of an edge is its material times the geometry of the edge. An
  extension that sets a circuit value itself (a lumped element, a conducting
  sheet) gives the node an entry with that value and its bit in abs_mask set.
  The material dumps see the effective material of such a value, unless its bit
  in keep_mask is set too: then they see the material of the entry \a base, the
  node's entry before the extension changed it (the cylinder axis, whose
  circuit values are those of all the sectors around it).

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
		uint16_t keep_mask;  //!< bit n*4 + k (with abs_mask): the dumps see the material of entry \a base for it, not the effective material of v[n*4 + k]
		uint32_t base;       //!< with keep_mask: an entry without circuit values of its own, the material of the node; else 0
	};

	EC_Table();

	//! The index of \a e, added if it is new
	uint32_t Intern(const Entry& e);

	size_t Size() const {return m_Size;}
	const Entry& operator[](size_t i) const {return m_Blocks[i/BLOCK][i%BLOCK];}
	//! The memory held, in bytes
	size_t Bytes() const {return m_Blocks.size()*BLOCK*sizeof(Entry) + m_Slots.size()*sizeof(uint32_t);}

private:
	// the entries in blocks: no copies when it grows, and no small allocations
	static const size_t BLOCK = 4096;
	std::vector<std::unique_ptr<Entry[]> > m_Blocks;
	size_t m_Size;
	//! open addressing, linear probing: the index of an entry plus one, 0 for a free slot; at most half of them used
	std::vector<uint32_t> m_Slots;

	Entry& Get(size_t i) {return m_Blocks[i/BLOCK][i%BLOCK];}
	//! Twice the slots, the entries hashed into them again
	void Grow();
};

#endif // EC_TABLE_H
