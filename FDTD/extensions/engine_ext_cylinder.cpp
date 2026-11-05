/*
*	Copyright (C) 2010 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#include "FDTD/engine.h"
#include "engine_ext_cylinder.h"
#include "operator_ext_cylinder.h"
#include "FDTD/engine_sse.h"
#include "engine_extension_dispatcher.h"

#include <algorithm>
#include <stdexcept>

Engine_Ext_Cylinder::Engine_Ext_Cylinder(Operator_Ext_Cylinder* op_ext) : Engine_Extension(op_ext)
{
	cyl_Op = op_ext;

	CC_closedAlpha = op_ext->CC_closedAlpha;
	CC_R0_included = op_ext->CC_R0_included;

	for (int n=0; n<3; ++n)
		numLines[n] = op_ext->m_Op->GetNumberOfLines(n,true);

	//this cylindrical extension should be executed first?
	m_Priority = ENG_EXT_PRIO_CYLINDER;
}

// the flat hooks: the whole domain as the box
void Engine_Ext_Cylinder::DoPostVoltageUpdates()
{
	const Engine_Box all = {{0, 0, 0}, {numLines[0], numLines[1], numLines[2]}};
	ENG_DISPATCH_ARGS(BoxImpl, true, &all);
}

void Engine_Ext_Cylinder::DoPostCurrentUpdates()
{
	const Engine_Box all = {{0, 0, 0}, {numLines[0], numLines[1], numLines[2]}};
	ENG_DISPATCH_ARGS(BoxImpl, false, &all);
}

// The r and z lines of the box, which must span the whole alpha axis: the r=0 axis sums
// the currents around it, the closed mesh copies the first alpha plane from the last
// (voltages) and the other way round (currents).
template <typename EngType>
void Engine_Ext_Cylinder::BoxImpl(EngType* eng, bool volt, const Engine_Box* box)
{
	if (CC_closedAlpha==false) return;
	if ((box->start[1]!=0) || (box->stop[1]<numLines[1]))
		throw std::logic_error("Engine_Ext_Cylinder: a box does not span the whole alpha axis");
	const unsigned int last_A_Line = numLines[1]-2;
	unsigned int pos[3];
	if (volt)
	{
		const unsigned int z0 = box->start[2], z1 = std::min(box->stop[2], numLines[2]);
		if (CC_R0_included && (box->start[0]==0))
		{
			for (pos[2]=z0; pos[2]<z1; ++pos[2])
			{
				FDTD_FLOAT v = eng->EngType::GetVolt(2,0,0,pos[2])*cyl_Op->vv_R0[pos[2]];
				// the alpha lines last_A_Line and last_A_Line+1 repeat the lines 0 and 1: each
				// alpha current around the axis once, as in the capacitance of vi_R0
				for (pos[1]=0; pos[1]<last_A_Line; ++pos[1])
					v += cyl_Op->vi_R0[pos[2]] * eng->EngType::GetCurr(1,0,pos[1],pos[2]);
				eng->EngType::SetVolt(2,0,0,pos[2], v);
			}
			for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
				for (pos[2]=z0; pos[2]<z1; ++pos[2])
				{
					eng->EngType::SetVolt(1,0,pos[1],pos[2], 0); //no voltage in alpha-direction at r=0
					eng->EngType::SetVolt(2,0,pos[1],pos[2], eng->EngType::GetVolt(2,0,0,pos[2]) );
				}
		}
		for (pos[0]=box->start[0]; pos[0]<std::min(box->stop[0], numLines[0]); ++pos[0])
			for (pos[2]=z0; pos[2]<z1; ++pos[2])
			{
				eng->EngType::SetVolt(0,pos[0],0,pos[2], eng->EngType::GetVolt(0,pos[0],last_A_Line,pos[2]) );
				eng->EngType::SetVolt(2,pos[0],0,pos[2], eng->EngType::GetVolt(2,pos[0],last_A_Line,pos[2]) );
			}
	}
	else
	{
		for (pos[0]=box->start[0]; pos[0]<std::min(box->stop[0], numLines[0]-1); ++pos[0])
			for (pos[2]=box->start[2]; pos[2]<std::min(box->stop[2], numLines[2]-1); ++pos[2])
			{
				eng->EngType::SetCurr(0,pos[0],last_A_Line,pos[2], eng->EngType::GetCurr(0,pos[0],0,pos[2]) );
				eng->EngType::SetCurr(2,pos[0],last_A_Line,pos[2], eng->EngType::GetCurr(2,pos[0],0,pos[2]) );
			}
	}
}

void Engine_Ext_Cylinder::DoPostVoltageUpdatesBox(const Engine_Box& box, unsigned int ts)
{
	(void)ts;
	ENG_DISPATCH_ARGS(BoxImpl, true, &box);
}

void Engine_Ext_Cylinder::DoPostCurrentUpdatesBox(const Engine_Box& box, unsigned int ts)
{
	(void)ts;
	ENG_DISPATCH_ARGS(BoxImpl, false, &box);
}
