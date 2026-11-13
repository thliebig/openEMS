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
#include <climits>
#include <cstdlib>
#include <stdexcept>
#include <string>

#include "tools/useful.h"

Engine_Plan Engine_Plan::Flat(const unsigned int numLines[3], unsigned int numThreads, int splitAxis, const unsigned int* first)
{
	Engine_Plan plan;
	const unsigned int f[3] = {first ? first[0] : 0, first ? first[1] : 0, first ? first[2] : 0};
	std::vector<unsigned int> jpt = AssignJobs2Threads(numLines[splitAxis]-f[splitAxis], numThreads, true);
	plan.volt.resize(jpt.size());
	plan.curr.resize(jpt.size());

	unsigned int start = f[splitAxis];
	for (size_t t=0; t<jpt.size(); ++t)
	{
		Engine_Box v, c;
		for (int n=0; n<3; ++n)
		{
			v.start[n] = c.start[n] = f[n];
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
	// at least two lines per thread: an extension on a plane normal to the thread axis
	// (Mur) needs the line inward of it in the same box
	const unsigned int lines = numLines[T]-last;
	const unsigned int used = std::max(1u, std::min(numThreads, lines/2));
	std::vector<unsigned int> jpt = AssignJobs2Threads(lines, used, false);
	jpt.resize(numThreads, 0);
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

//! The tile borders along an axis of \a N lines, tiles of \a width lines: 0, width, ..., N, a short last tile folded into the one before
static std::vector<unsigned int> TileBorders(unsigned int N, unsigned int width)
{
	std::vector<unsigned int> border(1, 0);
	while (N - border.back() >= 2*width)
		border.push_back(border.back() + width);
	border.push_back(N);
	return border;
}

Engine_BlockedPlan Engine_BlockedPlan::Trapezoids(const unsigned int numLines[3], unsigned int numThreads, int S, int T,
                                                  unsigned int width, unsigned int k)
{
	Engine_BlockedPlan plan;
	plan.k = k;
	const unsigned int N = numLines[S];
	const std::vector<unsigned int> border = TileBorders(N, std::max(width, 2*k+1));
	plan.tiles = border.size()-1;

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

Engine_BlockedPlan Engine_BlockedPlan::Parallelograms(const unsigned int numLines[3], unsigned int numThreads, int S, int T,
                                                      unsigned int width, unsigned int k)
{
	Engine_BlockedPlan plan;
	plan.k = k;
	const unsigned int N = numLines[S];
	const std::vector<unsigned int> border = TileBorders(N, std::max(width, 2*k+1));
	plan.tiles = border.size()-1;
	for (size_t i=0; i+1<border.size(); ++i)
	{
		const unsigned int a = border[i], b = border[i+1];
		for (unsigned int t=0; t<k; ++t)
		{
			// the first tile starts at line 0, the last ends at the end of the axis, at every step
			AddPhase(plan, numLines, numThreads, S, T, true,  t, (a>0) ? a-2*t   : 0, (b<N) ? b-2*t   : N);
			AddPhase(plan, numLines, numThreads, S, T, false, t, (a>0) ? a-2*t-1 : 0, (b<N) ? b-2*t-1 : N-1);
		}
	}
	return plan;
}

bool Engine_TilesFromEnv(const char* name, std::vector< std::pair<unsigned int, unsigned int> >& tiles)
{
	tiles.clear();
	const char* env = std::getenv(name);
	if ((env==NULL) || (env[0]=='\0'))
		return false;
	if (std::string(env)=="0")
		return true;
	const char* p = env;
	while (true)
	{
		char* end = NULL;
		const unsigned long k = std::strtoul(p, &end, 10);
		unsigned long W = 0;
		if ((end!=p) && (*end==':'))
		{
			const char* w = end+1;
			W = std::strtoul(w, &end, 10);
			if (end==w)
				end = NULL;
		}
		if ((end==NULL) || (end==p) || (k<2) || (k>(UINT_MAX-1)/2) || (W>UINT_MAX) || ((*end!='\0') && (*end!=',')))
			throw std::invalid_argument(std::string(name) + " must be 0, or k or k:W (k at least 2), several separated by commas");
		tiles.push_back(std::make_pair((unsigned int)k, (unsigned int)std::max<unsigned long>(W, 2*k+1)));
		if (*end=='\0')
			return true;
		p = end+1;
	}
}

double Engine_NumberFromEnv(const char* name, double fallback, bool zero)
{
	const char* env = std::getenv(name);
	if ((env==NULL) || (env[0]=='\0'))
		return fallback;
	char* end = NULL;
	const double v = std::strtod(env, &end);
	if ((end==env) || (*end!='\0') || (v<0) || ((v==0) && !zero))
		throw std::invalid_argument(std::string(name) + " must be a number " + (zero ? "of 0 or more" : "greater than 0"));
	return v;
}
