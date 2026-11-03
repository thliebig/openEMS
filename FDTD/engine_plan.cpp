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

#include "engine_plan.h"

#include <algorithm>

#include "tools/useful.h"

Engine_Plan Engine_Plan::Flat(const unsigned int numLines[3], unsigned int numThreads, int splitAxis)
{
	Engine_Plan plan;
	std::vector<unsigned int> jpt = AssignJobs2Threads(numLines[splitAxis], numThreads, true);
	plan.volt.resize(jpt.size());
	plan.curr.resize(jpt.size());

	unsigned int start = 0;
	for (size_t t=0; t<jpt.size(); ++t)
	{
		Engine_Box v, c;
		for (int n=0; n<3; ++n)
		{
			v.start[n] = c.start[n] = 0;
			v.stop[n] = numLines[n];
			c.stop[n] = numLines[n]-1;
		}
		v.start[splitAxis] = c.start[splitAxis] = start;
		v.stop[splitAxis] = start + jpt[t];
		c.stop[splitAxis] = std::min(start + jpt[t], numLines[splitAxis]-1);
		start += jpt[t];

		plan.volt[t].push_back(v);
		if (!c.Empty())
			plan.curr[t].push_back(c);
	}
	return plan;
}
