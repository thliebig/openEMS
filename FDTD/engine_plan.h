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

#ifndef ENGINE_PLAN_H
#define ENGINE_PLAN_H

#include <vector>

//! A box of mesh lines, [start, stop) per axis
struct Engine_Box
{
	unsigned int start[3];
	unsigned int stop[3];

	bool Empty() const {return (start[0]>=stop[0]) || (start[1]>=stop[1]) || (start[2]>=stop[2]);}
};

//! The work of the threads of Engine_ThreadPool in a half step
/*!
  Per thread the boxes of its voltage and of its current update. The plan is
  in mesh coordinates; how a box maps to memory is the engine's business
  (Engine::UpdateVoltages(const Engine_Box&)). A generator builds the plan from
  what the engine can update; the pool only walks it.

  The voltages are updated on the lines [0, numLines) of each axis, the
  currents on [0, numLines-1): the last line of an axis has no current.
  */
struct Engine_Plan
{
	std::vector< std::vector<Engine_Box> > volt;   //!< per thread its voltage boxes
	std::vector< std::vector<Engine_Box> > curr;   //!< per thread its current boxes

	unsigned int GetNumThreads() const {return volt.size();}

	//! The flat sweep: the domain split into one box per thread along \a splitAxis; fewer threads than requested if the axis has fewer lines
	static Engine_Plan Flat(const unsigned int numLines[3], unsigned int numThreads, int splitAxis);
};

#endif // ENGINE_PLAN_H
