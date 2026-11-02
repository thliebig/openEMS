/*
*	Copyright (C) 2026 Sean Mollet (sean@malmoset.com)
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

#ifndef MULTIGRID_INTERPOLATION_H
#define MULTIGRID_INTERPOLATION_H

#include <vector>

//! The interpolation of a cylindrical multi-grid level from its sub-grid (Operator_CylinderMultiGrid::GetInterpolation())
/*!
  The sub-grid has the same r and z lines as the base grid (up to its size) and
  every other alpha line. For each alpha line a of the base grid and direction n
  (0: r and z components, 1: alpha component), a field of the base grid is
  interpolated from the sub-grid alpha lines pos_2p[n][a] and pos_2pp[n][a] with
  the weights w_2p[n][a] and w_2pp[n][a] (_v: voltages, _i: currents).
  Independent of the engine: each coupling (Engine_SubGridCoupling,
  GPU_MultiGridLink) takes it into its own layout.
  */
struct MultiGridInterpolation
{
	unsigned int split_pos = 0;   //!< first base grid r line not covered by the sub-grid
	std::vector<unsigned int> pos_v_2p[2], pos_v_2pp[2], pos_i_2p[2], pos_i_2pp[2];
	std::vector<float> w_v_2p[2], w_v_2pp[2], w_i_2p[2], w_i_2pp[2];
};

#endif // MULTIGRID_INTERPOLATION_H
