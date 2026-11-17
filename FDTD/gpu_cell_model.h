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

#ifndef GPU_CELL_MODEL_H
#define GPU_CELL_MODEL_H

#include <algorithm>
#include <cstdint>
#include <vector>
#include "kernels/cell_update.h"
#include "engine_plan.h"

class CellModel;

//! The cell model of a grid expanded for a GPU: per node and field an index into the combinations of the three component sets
/*!
  Built from a CellModel with the fast axis z, split x, mid y (the ArrayNIJK layout of the GPU engine). Per field
  (0 voltages, 1 currents): the index of each node into \a combos, its state
  slot (UPML and dispersive nodes), the poles and profiles of all three
  components, and the mesh factors. The kernels: kernels/cell_update.h.
  */
struct GPU_CellModel
{
	struct Field
	{
		std::vector<uint32_t> index;      //!< per node its combination
		std::vector<CM_Combo> combos;
		std::vector<uint32_t> slot;       //!< per node its state slot, 0xffffffff none; empty if no node has state
		uint64_t numSlots = 0;
		unsigned int slotFloats = 0;      //!< floats per state slot
		std::vector<CM_Pole> poles;
		std::vector<float> profile[5];
		std::vector<float> gxy[3];    //!< gx*gy per line (x, y), see CM_Field
		std::vector<float> gz[3];
	};
	unsigned int numLines[3];
	Field field[2];

	//! Expand \a model (built along z, split x, mid y); throws if it was built along other axes
	explicit GPU_CellModel(const CellModel& model, const unsigned int numLines[3]);

	//! Whether a device stores the node index in 16 bit: both fields have at most 65536 combinations
	bool Index16() const {return (field[0].combos.size()<=65536) && (field[1].combos.size()<=65536);}
	//! The device memory a main update reads per node, on average: the fields, the node index, the state slot and the state (the tile size of the temporal blocking)
	double NodeBytes() const;
};

//! The nodes of \a box on a grid of \a nx x \a ny x \a nz lines as the main update kernels take them; false if there are none
inline bool GPU_UpdateBox(const Engine_Box& box, unsigned int nx, unsigned int ny, unsigned int nz, CM_Box& b)
{
	b.x0 = box.start[0];
	b.x1 = (std::min)(box.stop[0], nx);
	b.y0 = box.start[1];
	b.y1 = (std::min)(box.stop[1], ny);
	b.z0 = box.start[2];
	b.z1 = (std::min)(box.stop[2], nz);
	return (b.x0<b.x1) && (b.y0<b.y1) && (b.z0<b.z1);
}

#endif // GPU_CELL_MODEL_H
