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

//! The tile borders along an axis of \a N lines, tiles of \a width lines: 0, width, ..., N, a short last tile folded into the one before
static std::vector<unsigned int> TileBorders(unsigned int N, unsigned int width)
{
	const unsigned int tiles = Engine_BlockedPlan::TileCount(N, width);
	std::vector<unsigned int> border;
	for (unsigned int i=0; i<tiles; ++i)
		border.push_back(i*width);
	border.push_back(N);
	return border;
}

namespace
{
//! A piece of one axis of N lines in a block of k timesteps: per timestep t its voltage and current ranges
struct Piece
{
	enum Kind {CORE, WEDGE, DRIFT};
	Kind kind;
	unsigned int a, b;   //!< a core or a parallelogram: the tile [a, b); a wedge: the border a (b unused)

	//! the lines [lo, hi) of the voltages or the currents (one line less: the last line of an axis has none) at timestep t
	void Range(unsigned int N, unsigned int t, bool volt, unsigned int& lo, unsigned int& hi) const
	{
		const unsigned int c = volt ? 0 : 1;
		switch (kind)
		{
		case CORE:    // shrinks by a line per timestep on a side that borders another tile
			lo = (a>0) ? a+t : 0;
			hi = (b<N) ? b-t-c : N-c;
			break;
		case WEDGE:   // grows into the gap the cores on either side left
			lo = a-t-c;
			hi = a+t;
			break;
		case DRIFT:   // moves down by a line per timestep, the first from 0, the last to the end
			lo = (a>0) ? a-t-c : 0;
			hi = (b<N) ? b-t-c : N-c;
			break;
		}
	}
};

//! the pieces of \a kind along an axis with the tile borders \a border: a core or a parallelogram per tile, a wedge per inner border
std::vector<Piece> Pieces(const std::vector<unsigned int>& border, Piece::Kind kind)
{
	std::vector<Piece> p;
	if (kind==Piece::WEDGE)
		for (size_t i=1; i+1<border.size(); ++i)
			p.push_back(Piece{kind, border[i], 0});
	else
		for (size_t i=0; i+1<border.size(); ++i)
			p.push_back(Piece{kind, border[i], border[i+1]});
	return p;
}

//! the steps of the product of \a s along \a S and \a m along \a M through the k timesteps of a block, added to \a steps; returns their cells
double AddSteps(std::vector<Engine_Step>& steps, const unsigned int numLines[3], int S, int M, const Piece& s, const Piece& m, unsigned int k)
{
	double cells = 0;
	for (unsigned int t=0; t<k; ++t)
		for (int volt=1; volt>=0; --volt)
		{
			Engine_Step st;
			st.volt = volt;
			st.t = t;
			for (int a=0; a<3; ++a)
			{
				st.box.start[a] = 0;
				st.box.stop[a] = numLines[a] - (volt ? 0 : 1);
			}
			s.Range(numLines[S], t, volt, st.box.start[S], st.box.stop[S]);
			m.Range(numLines[M], t, volt, st.box.start[M], st.box.stop[M]);
			if (st.box.Empty())
				continue;
			steps.push_back(st);
			cells += (double)(st.box.stop[0]-st.box.start[0])*(st.box.stop[1]-st.box.start[1])*(st.box.stop[2]-st.box.start[2]);
		}
	return cells;
}
}

Engine_BlockedPlan Engine_BlockedPlan::Parallelograms(const unsigned int numLines[3], int S, int M,
                                                      unsigned int widthS, unsigned int widthM, unsigned int k)
{
	Engine_BlockedPlan plan;
	plan.k = k;
	const std::vector<unsigned int> bS = TileBorders(numLines[S], std::max(widthS, 2*k+1));
	const std::vector<unsigned int> bM = TileBorders(numLines[M], std::max(widthM, 2*k+1));
	plan.tiles[0] = bS.size()-1;
	plan.tiles[1] = bM.size()-1;
	std::vector<Engine_Step> steps;
	for (const Piece& s : Pieces(bS, Piece::DRIFT))
		for (const Piece& m : Pieces(bM, Piece::DRIFT))
			AddSteps(steps, numLines, S, M, s, m, k);
	plan.phases.push_back(std::vector< std::vector<Engine_Step> >(1, steps));
	return plan;
}

Engine_BlockedPlan Engine_BlockedPlan::Tiles(const unsigned int numLines[3], unsigned int numThreads, int S, int M,
                                             unsigned int widthS, unsigned int widthM, unsigned int k)
{
	Engine_BlockedPlan plan;
	plan.k = k;
	const std::vector<unsigned int> bS = TileBorders(numLines[S], std::max(widthS, 2*k+1));
	const std::vector<unsigned int> bM = TileBorders(numLines[M], std::max(widthM, 2*k+1));
	plan.tiles[0] = bS.size()-1;
	plan.tiles[1] = bM.size()-1;
	numThreads = std::max(1u, numThreads);
	for (int phase=0; phase<4; ++phase)
	{
		// cores x cores, wedges x cores, cores x wedges, wedges x wedges: per piece its steps and cells
		std::vector< std::vector<Engine_Step> > pieces;
		std::vector<double> cells;
		double total = 0;
		for (const Piece& s : Pieces(bS, (phase&1) ? Piece::WEDGE : Piece::CORE))
			for (const Piece& m : Pieces(bM, (phase&2) ? Piece::WEDGE : Piece::CORE))
			{
				std::vector<Engine_Step> steps;
				const double n = AddSteps(steps, numLines, S, M, s, m, k);
				if (steps.empty())
					continue;
				pieces.push_back(steps);
				cells.push_back(n);
				total += n;
			}
		if (pieces.empty())
			continue;
		// dealt in order (neighbouring pieces on one thread), each thread up to its share of the cells
		std::vector< std::vector<Engine_Step> > threads(numThreads);
		double done = 0;
		unsigned int th = 0;
		for (size_t p=0; p<pieces.size(); ++p)
		{
			// on to the next thread once this one has its share, judged by the middle of the piece
			while ((th+1<numThreads) && (done + cells[p]/2 > total*(th+1)/numThreads))
				++th;
			threads[th].insert(threads[th].end(), pieces[p].begin(), pieces[p].end());
			done += cells[p];
		}
		plan.phases.push_back(threads);
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
