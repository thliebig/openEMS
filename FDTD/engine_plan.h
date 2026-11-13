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

#include <cstddef>
#include <utility>
#include <vector>

//! A box of mesh lines, [start, stop) per axis
struct Engine_Box
{
	unsigned int start[3];
	unsigned int stop[3];

	bool Empty() const {return (start[0]>=stop[0]) || (start[1]>=stop[1]) || (start[2]>=stop[2]);}
	bool Contains(unsigned int x, unsigned int y, unsigned int z) const
	{
		return (x>=start[0]) && (x<stop[0]) && (y>=start[1]) && (y<stop[1]) && (z>=start[2]) && (z<stop[2]);
	}
};

//! Whether \a box (all cells if NULL) holds the cell (x, y, z): the filter of the extensions' box hooks
inline bool Engine_Box_Contains(const Engine_Box* box, unsigned int x, unsigned int y, unsigned int z)
{
	return (box==NULL) || box->Contains(x, y, z);
}

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
	/*!
	  \a first: per axis the first line updated (default 0), e.g. the outer grid of a
	  cylindrical multi-grid leaves the lines its sub-grid covers (Operator::CalcStartStopLines())
	  */
	static Engine_Plan Flat(const unsigned int numLines[3], unsigned int numThreads, int splitAxis, const unsigned int* first=NULL);
};

//! One half step of one thread in a blocked plan: the voltages or currents of a box at timestep t of the block
struct Engine_Step
{
	bool volt;
	unsigned int t;
	Engine_Box box;
};

//! Several timesteps per cache-sized tile (temporal blocking)
/*!
  A block advances the whole domain by k timesteps in phases; a barrier ends
  each phase, within a phase each thread walks its steps in order. Every cell
  is updated exactly once per half step, its neighbours at the right time
  levels (tests/test_engine_plan.cpp). The extensions run per step on its box
  (Engine_Extension::SupportsBoxes()).
  */
struct Engine_BlockedPlan
{
	unsigned int k = 0;   //!< the timesteps of a block
	unsigned int tiles = 0;   //!< the tiles along the split axis
	std::vector< std::vector< std::vector<Engine_Step> > > phases;   //!< [phase][thread] its steps

	unsigned int GetNumThreads() const {return phases.empty() ? 0 : phases[0].size();}

	//! Trapezoids of \a width lines along \a splitAxis, all threads on one tile, the tile split along \a threadAxis
	/*!
	  The cores (the tiles, shrinking by a line per timestep on each side that
	  borders another tile) first, then the wedges between them, which grow
	  into the gaps; each a phase per half step. Needs width >= 2k+1; a short
	  last tile is folded into the one before.
	  */
	static Engine_BlockedPlan Trapezoids(const unsigned int numLines[3], unsigned int numThreads, int splitAxis, int threadAxis,
	                                     unsigned int width, unsigned int k);
	//! Parallelograms of \a width lines along \a splitAxis, one tile after the other, all threads on one tile
	/*!
	  A tile moves down by a line per half step: V_t on [a-2t, b-2t), I_t on
	  [a-2t-1, b-2t-1). It reads only lines of the tile before it, which that
	  one has finished up to there, so no wedges are needed, but the tiles run
	  in order: a phase per half step of a tile, half the phases of
	  Trapezoids(). The first tile starts at line 0 at every step, the last
	  reaches the end of the axis. Needs width >= 2k+1; a short last tile is
	  folded into the one before.
	  */
	static Engine_BlockedPlan Parallelograms(const unsigned int numLines[3], unsigned int numThreads, int splitAxis, int threadAxis,
	                                         unsigned int width, unsigned int k);
};

//! The blocked plans an engine offers its thread pool (Engine::BuildBlockedPlans())
/*!
  Each candidate is a set of plans, one per k by decreasing k (the processing
  between two calls of the engine may leave fewer timesteps than the largest
  k). The candidates are ordered, e.g. by decreasing tile size. With a trial
  window the pool runs, at the start of the run, candidates for that many
  seconds each (after a warm-up of \a start): \a start, then on to the next
  as long as each is faster than the one before; if already the one after
  \a start is slower, from \a start towards the first. Then the flat sweep,
  and it keeps the fastest: the blocked and the flat sweep compute the same,
  so it may switch between two calls. Without a window it runs \a start.
  */
struct Engine_BlockedChoice
{
	std::vector< std::vector<Engine_BlockedPlan> > candidates;
	unsigned int start = 0;    //!< the candidate expected fastest
	double trialSeconds = 0;   //!< 0: no trials, \a start
};

//! The tiles of the temporal blocking the environment variable \a name sets, as (k, W): "0" none, "k:W" or "k:W,k:W,..." (k timesteps per block, at least 2; tiles of W lines, at least 2k+1, "k" alone the narrowest)
/*!
  For OPENEMS_CPU_TBLOCK and OPENEMS_GPU_TBLOCK. Returns false if the variable is
  unset or empty (the engine chooses), true with \a tiles empty for "0". Throws
  std::invalid_argument for anything else.
  */
bool Engine_TilesFromEnv(const char* name, std::vector< std::pair<unsigned int, unsigned int> >& tiles);
//! The number the environment variable \a name sets, \a fallback if it is unset or empty; throws std::invalid_argument if it is no number, negative, or 0 and not \a zero
double Engine_NumberFromEnv(const char* name, double fallback, bool zero);

#endif // ENGINE_PLAN_H
