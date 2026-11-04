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

// the phase of one half step: the range [lo, hi) along the split axis, split over the threads along the thread axis
static void AddPhase(Engine_BlockedPlan& plan, const unsigned int numLines[3], unsigned int numThreads, int S, int T,
                     bool volt, unsigned int t, unsigned int lo, unsigned int hi)
{
	if (lo>=hi)
		return;
	const unsigned int last = volt ? 0 : 1;   // the currents end a line earlier on every axis
	std::vector< std::vector<Engine_Step> > phase(numThreads);
	std::vector<unsigned int> jpt = AssignJobs2Threads(numLines[T]-last, numThreads, false);
	unsigned int start = 0;
	for (unsigned int n=0; n<numThreads; ++n)
	{
		Engine_Step step;
		step.volt = volt;
		step.t = t;
		for (int a=0; a<3; ++a)
		{
			step.box.start[a] = 0;
			step.box.stop[a] = numLines[a]-last;
		}
		step.box.start[S] = lo;
		step.box.stop[S] = hi;
		step.box.start[T] = start;
		step.box.stop[T] = start + jpt[n];
		start += jpt[n];
		if (jpt[n]>0)
			phase[n].push_back(step);
	}
	plan.phases.push_back(phase);
}

Engine_BlockedPlan Engine_BlockedPlan::Trapezoids(const unsigned int numLines[3], unsigned int numThreads, int S, int T,
                                                  unsigned int width, unsigned int k)
{
	Engine_BlockedPlan plan;
	plan.k = k;
	const unsigned int N = numLines[S];
	if (width < 2*k+1)
		width = 2*k+1;
	// the tile borders, a short last tile folded into the one before
	std::vector<unsigned int> border(1, 0);
	while (N - border.back() >= 2*width)
		border.push_back(border.back() + width);
	border.push_back(N);

	// the cores: V_t on [a+t, b-t), I_t on [a+t, b-t-1), not shrinking at the domain edges
	for (size_t i=0; i+1<border.size(); ++i)
	{
		const unsigned int a = border[i], b = border[i+1];
		for (unsigned int t=0; t<k; ++t)
		{
			const unsigned int lo = (a>0) ? a+t : 0;
			AddPhase(plan, numLines, numThreads, S, T, true,  t, lo, (b<N) ? b-t   : N);
			AddPhase(plan, numLines, numThreads, S, T, false, t, lo, (b<N) ? b-t-1 : N-1);
		}
	}
	// the wedges at the inner borders: V_t on [a-t, a+t), I_t on [a-t-1, a+t)
	for (size_t i=1; i+1<border.size(); ++i)
	{
		const unsigned int a = border[i];
		for (unsigned int t=0; t<k; ++t)
		{
			AddPhase(plan, numLines, numThreads, S, T, true,  t, a-t,   a+t);
			AddPhase(plan, numLines, numThreads, S, T, false, t, a-t-1, a+t);
		}
	}
	return plan;
}
