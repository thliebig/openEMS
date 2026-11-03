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

/*
  Tests for the plans of Engine_ThreadPool (Engine_Plan): in every half step
  each cell the engine updates is in exactly one box of exactly one thread,
  the voltages on [0, numLines), the currents on [0, numLines-1) per axis,
  and no thread gets more boxes than the plan has threads for.

  Build with -DOPENEMS_BUILD_TESTS=ON and run it through ctest:

    cmake -DOPENEMS_BUILD_TESTS=ON <the usual options> ..
    make && ctest --output-on-failure

  Exits non-zero and prints "FAIL: ..." per failed check.
*/

#include "FDTD/engine_plan.h"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <vector>

static int fails = 0;
#define CHECK(cond, msg) do { if (!(cond)) { std::cout << "FAIL: " << msg << "\n"; ++fails; } } while (0)

// count per cell how often the boxes of all threads cover it, and check each
// cell inside [0, numLines - last) once and none outside
static void CheckCover(const std::vector< std::vector<Engine_Box> >& boxes, const unsigned int numLines[3],
                       unsigned int last, const std::string& what)
{
	const unsigned int n0 = numLines[0], n1 = numLines[1], n2 = numLines[2];
	std::vector<int> count((size_t)n0*n1*n2, 0);
	bool inside = true;
	for (const std::vector<Engine_Box>& thread : boxes)
		for (const Engine_Box& b : thread)
		{
			for (int n=0; n<3; ++n)
				inside = inside && b.stop[n]<=numLines[n];
			if (!inside)
				break;
			for (unsigned int i=b.start[0]; i<b.stop[0]; ++i)
				for (unsigned int j=b.start[1]; j<b.stop[1]; ++j)
					for (unsigned int k=b.start[2]; k<b.stop[2]; ++k)
						++count[((size_t)i*n1 + j)*n2 + k];
		}
	CHECK(inside, what << ": a box reaches outside the domain");
	int bad = 0;
	for (unsigned int i=0; i<n0; ++i)
		for (unsigned int j=0; j<n1; ++j)
			for (unsigned int k=0; k<n2; ++k)
			{
				const bool updated = (i+last<n0) && (j+last<n1) && (k+last<n2);
				bad += count[((size_t)i*n1 + j)*n2 + k] != (updated ? 1 : 0);
			}
	CHECK(bad==0, what << ": " << bad << " cells not covered exactly once");
}

int main()
{
	// odd sizes, a split axis with fewer lines than threads, and the minimum of 2 lines
	const unsigned int sizes[][3] = {{7, 5, 3}, {2, 9, 4}, {13, 2, 6}, {5, 6, 2}};
	for (const auto& numLines : sizes)
		for (int axis=0; axis<3; ++axis)
			for (unsigned int threads=1; threads<=9; ++threads)
			{
				std::ostringstream what;
				what << numLines[0] << "x" << numLines[1] << "x" << numLines[2]
				     << ", split " << "xyz"[axis] << ", " << threads << " threads";
				Engine_Plan plan = Engine_Plan::Flat(numLines, threads, axis);
				const unsigned int expect = std::min(threads, numLines[axis]);
				CHECK(plan.GetNumThreads()==expect && plan.curr.size()==expect,
				      what.str() << ": " << plan.GetNumThreads() << " threads in the plan, not " << expect);
				CheckCover(plan.volt, numLines, 0, what.str() + ", voltages");
				CheckCover(plan.curr, numLines, 1, what.str() + ", currents");
			}

	if (fails==0)
		std::cout << "all Engine_Plan checks passed\n";
	return fails ? 1 : 0;
}
