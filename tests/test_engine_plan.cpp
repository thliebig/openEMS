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
  Tests for the plans of Engine_ThreadPool. The flat plan (Engine_Plan): in
  every half step each cell the engine updates is in exactly one box of
  exactly one thread, the voltages on [0, numLines), the currents on
  [0, numLines-1) per axis, and no thread gets more boxes than the plan has
  threads for. The blocked plans (Engine_BlockedPlan) by a stamp check: a
  time level per cell, each update checked to find the cell and its
  neighbours at the levels the Yee scheme needs, no cell touched by two
  threads in one phase, and every cell at k after the block.

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

// the time levels of a block: a voltage at p (level t -> t+1) reads the currents at p and
// p-e_a (level t), a current at p (t -> t+1) the voltages at p and p+e_a (level t+1); along
// all three axes, which covers the curl of every component
static void CheckBlocked(const Engine_BlockedPlan& plan, const unsigned int numLines[3], const std::string& what)
{
	const unsigned int n0 = numLines[0], n1 = numLines[1], n2 = numLines[2];
	const size_t cells = (size_t)n0*n1*n2;
	std::vector<int> vlev(cells, 0), ilev(cells, 0);
	auto idx = [&](const unsigned int p[3]) {return ((size_t)p[0]*n1 + p[1])*n2 + p[2];};
	int bad = 0, race = 0;
	for (const auto& phase : plan.phases)
	{
		// who writes and who reads each cell in this phase, per field
		std::vector<int> writer[2] = {std::vector<int>(cells, -1), std::vector<int>(cells, -1)};
		std::vector<int> reader[2] = {std::vector<int>(cells, -1), std::vector<int>(cells, -1)};
		auto touch = [&](std::vector<int>& who, size_t i, int th) {
			if ((who[i]>=0) && (who[i]!=th)) return false;
			who[i] = th; return true;};
		for (size_t th=0; th<phase.size(); ++th)
			for (const Engine_Step& st : phase[th])
			{
				const unsigned int last = st.volt ? 0 : 1;
				for (int a=0; a<3; ++a)
					if (st.box.stop[a]+last>numLines[a])
						++bad;
				unsigned int p[3];
				for (p[0]=st.box.start[0]; p[0]<st.box.stop[0]; ++p[0])
					for (p[1]=st.box.start[1]; p[1]<st.box.stop[1]; ++p[1])
						for (p[2]=st.box.start[2]; p[2]<st.box.stop[2]; ++p[2])
						{
							const size_t i = idx(p);
							const int t = st.t;
							std::vector<int>& own = st.volt ? vlev : ilev;
							std::vector<int>& other = st.volt ? ilev : vlev;
							const int f = st.volt ? 0 : 1;
							if (own[i]!=t)
								++bad;
							if (!touch(writer[f], i, th) || ((reader[f][i]>=0) && (reader[f][i]!=(int)th)))
								++race;
							const int need = st.volt ? t : t+1;
							// the cell itself and its neighbours of the other field
							for (int a=-1; a<3; ++a)
							{
								unsigned int q[3] = {p[0], p[1], p[2]};
								if (a>=0)
								{
									if (st.volt) { if (q[a]==0) continue; --q[a]; }
									else         { ++q[a]; }
								}
								const bool exists = st.volt ? (q[0]+1<n0 && q[1]+1<n1 && q[2]+1<n2) : true;
								if (!exists)
									continue;   // a voltage on the last line of an axis has no current there
								const size_t j = idx(q);
								if (other[j]!=need)
									++bad;
								if ((writer[1-f][j]>=0) && (writer[1-f][j]!=(int)th))
									++race;
								reader[1-f][j] = (reader[1-f][j]<0 || reader[1-f][j]==(int)th) ? th : -2;
							}
							own[i] = t+1;
						}
			}
	}
	int unfinished = 0;
	for (unsigned int i=0; i<n0; ++i)
		for (unsigned int j=0; j<n1; ++j)
			for (unsigned int k=0; k<n2; ++k)
			{
				const unsigned int p[3] = {i, j, k};
				unfinished += vlev[idx(p)] != (int)plan.k;
				if ((i+1<n0) && (j+1<n1) && (k+1<n2))
					unfinished += ilev[idx(p)] != (int)plan.k;
			}
	CHECK(bad==0, what << ": " << bad << " updates out of order or out of the domain");
	CHECK(race==0, what << ": " << race << " cells touched by two threads in one phase");
	CHECK(unfinished==0, what << ": " << unfinished << " cells not at level k after the block");
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

	// blocked plans: tiles that divide the axis and ones that do not, a short tail, a single tile
	const unsigned int bsizes[][3] = {{23, 4, 5}, {16, 5, 3}, {9, 3, 4}, {41, 3, 3}};
	const unsigned int kw[][2] = {{1, 3}, {2, 5}, {3, 7}, {2, 8}, {4, 9}};
	for (const auto& numLines : bsizes)
		for (int S=0; S<3; ++S)
			for (int T=0; T<3; ++T)
			{
				if (T==S)
					continue;
				for (const auto& p : kw)
					for (unsigned int threads : {1u, 2u, 3u})
					{
						std::ostringstream what;
						what << numLines[0] << "x" << numLines[1] << "x" << numLines[2] << ", split " << "xyz"[S]
						     << ", threads along " << "xyz"[T] << ", k " << p[0] << ", width " << p[1] << ", " << threads << " threads";
						Engine_BlockedPlan plan = Engine_BlockedPlan::Trapezoids(numLines, threads, S, T, p[1], p[0]);
						CHECK(plan.k==p[0], what.str() << ": k " << plan.k);
						CheckBlocked(plan, numLines, what.str());
					}
			}

	if (fails==0)
		std::cout << "all Engine_Plan checks passed\n";
	return fails ? 1 : 0;
}
