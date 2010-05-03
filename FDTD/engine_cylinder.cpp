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

#include "engine_cylinder.h"
#include "engine_extension.h"
#include "operator_extension.h"

Engine_Cylinder* Engine_Cylinder::New(const Operator_Cylinder* op)
{
	Engine_Cylinder* e = new Engine_Cylinder(op);
	e->Init();
	return e;
}

Engine_Cylinder::Engine_Cylinder(const Operator_Cylinder* op) : Engine(op)
{
	cyl_Op = op;
	if (cyl_Op->GetClosedAlpha())
	{
		++numLines[1]; //necessary for dobled voltage and current line in alpha-dir, operator will return one smaller for correct post-processing
	}
}

inline void Engine_Cylinder::CloseAlphaVoltages()
{
	unsigned int pos[3];
	// copy voltages from last alpha-plane to first
	unsigned int last_A_Line = numLines[1]-1;
	for (pos[0]=0;pos[0]<numLines[0];++pos[0])
	{
		for (pos[2]=0;pos[2]<numLines[2];++pos[2])
		{
			volt[0][pos[0]][0][pos[2]] = volt[0][pos[0]][last_A_Line][pos[2]];
			volt[1][pos[0]][0][pos[2]] = volt[1][pos[0]][last_A_Line][pos[2]];
			volt[2][pos[0]][0][pos[2]] = volt[2][pos[0]][last_A_Line][pos[2]];
		}
	}
}

void Engine_Cylinder::R0IncludeVoltages()
{
	unsigned int pos[3];
	pos[0] = 0;
	for (pos[2]=0;pos[2]<numLines[2];++pos[2])
	{
		volt[2][0][0][pos[2]] *= cyl_Op->vv_R0[pos[2]];
		for (pos[1]=0;pos[1]<numLines[1]-cyl_Op->GetClosedAlpha();++pos[1])
		{
				volt[2][0][0][pos[2]] += cyl_Op->vi_R0[pos[2]] *  curr[1][0][pos[1]][pos[2]];
		}
	}
	for (pos[1]=0;pos[1]<numLines[1];++pos[1])
	{
		for (pos[2]=0;pos[2]<numLines[2];++pos[2])
		{
			volt[1][0][pos[1]][pos[2]] = 0; //no voltage in alpha-direction at r=0
			volt[2][0][pos[1]][pos[2]] = volt[2][0][0][pos[2]];
		}
	}
}

inline void Engine_Cylinder::CloseAlphaCurrents()
{
	unsigned int pos[3];
	// copy currents from first alpha-plane to last
	for (pos[0]=0;pos[0]<numLines[0]-1;++pos[0])
	{
		unsigned int last_A_Line = numLines[1]-1;
		for (pos[2]=0;pos[2]<numLines[2]-1;++pos[2])
		{
			curr[0][pos[0]][last_A_Line][pos[2]] = curr[0][pos[0]][0][pos[2]];
			curr[1][pos[0]][last_A_Line][pos[2]] = curr[1][pos[0]][0][pos[2]];
			curr[2][pos[0]][last_A_Line][pos[2]] = curr[2][pos[0]][0][pos[2]];
		}
	}
}

bool Engine_Cylinder::IterateTS(unsigned int iterTS)
{
	if (cyl_Op->GetClosedAlpha()==false)
		return Engine::IterateTS(iterTS);

	for (unsigned int iter=0;iter<iterTS;++iter)
	{
		//voltage updates with extensions
		for (size_t n=0;n<m_Eng_exts.size();++n)
			m_Eng_exts.at(n)->DoPreVoltageUpdates();

		UpdateVoltages();

		for (size_t n=0;n<m_Eng_exts.size();++n)
			m_Eng_exts.at(n)->DoPostVoltageUpdates();
		for (size_t n=0;n<m_Eng_exts.size();++n)
			m_Eng_exts.at(n)->Apply2Voltages();

		if (cyl_Op->GetR0Included())
			R0IncludeVoltages();

		ApplyVoltageExcite();

		CloseAlphaVoltages();

		//current updates with extensions
		for (size_t n=0;n<m_Eng_exts.size();++n)
			m_Eng_exts.at(n)->DoPreCurrentUpdates();

		UpdateCurrents();

		for (size_t n=0;n<m_Eng_exts.size();++n)
			m_Eng_exts.at(n)->DoPostCurrentUpdates();
		for (size_t n=0;n<m_Eng_exts.size();++n)
			m_Eng_exts.at(n)->Apply2Current();

		ApplyCurrentExcite();

		CloseAlphaCurrents();

		++numTS;
	}

	return true;
}
