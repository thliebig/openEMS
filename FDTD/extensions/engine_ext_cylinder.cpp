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

Engine_Ext_Cylinder::Engine_Ext_Cylinder(Operator_Ext_Cylinder* op_ext) : Engine_Extension(op_ext)
{
	cyl_Op = op_ext;
	m_Eng_SSE = NULL;

	CC_closedAlpha = op_ext->CC_closedAlpha;
	CC_R0_included = op_ext->CC_R0_included;

	for (int n=0; n<3; ++n)
		numLines[n] = op_ext->m_Op->GetNumberOfLines(n,true);

	//this cylindrical extension should be executed first?
	m_Priority = ENG_EXT_PRIO_CYLINDER;
}

void Engine_Ext_Cylinder::SetEngine(Engine* eng)
{
	Engine_Extension::SetEngine(eng);
	m_Eng_SSE = dynamic_cast<Engine_sse*>(m_Eng);
}

void Engine_Ext_Cylinder::DoPostVoltageUpdates()
{
	if (CC_closedAlpha==false) return;

	if (CC_R0_included)
	{
		unsigned int pos[3];
		pos[0] = 0;
		FDTD_FLOAT volt=0;
		for (pos[2]=0; pos[2]<numLines[2]; ++pos[2])
		{
			volt =m_Eng_SSE->Engine_sse::GetVolt(2,0,0,pos[2])*cyl_Op->vv_R0[pos[2]];
			for (pos[1]=0; pos[1]<numLines[1]-1; ++pos[1])
				volt +=cyl_Op->vi_R0[pos[2]] *  m_Eng_SSE->Engine_sse::GetCurr(1,0,pos[1],pos[2]);
			m_Eng_SSE->Engine_sse::SetVolt(2,0,0,pos[2], volt);
		}

		for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
		{
			for (pos[2]=0; pos[2]<numLines[2]; ++pos[2])
			{
				m_Eng_SSE->Engine_sse::SetVolt(1,0,pos[1],pos[2], 0); //no voltage in alpha-direction at r=0
				m_Eng_SSE->Engine_sse::SetVolt(2,0,pos[1],pos[2], m_Eng_SSE->Engine_sse::GetVolt(2,0,0,pos[2]) );
			}
		}
	}

	//close alpha
	unsigned int pos[3];
	// copy tangential voltages from last alpha-plane to first
	unsigned int last_A_Line = numLines[1]-2;
	for (pos[0]=0; pos[0]<numLines[0]; ++pos[0])
	{
		for (pos[2]=0; pos[2]<numLines[2]; ++pos[2])
		{
			m_Eng_SSE->Engine_sse::SetVolt(0,pos[0],0,pos[2], m_Eng_SSE->Engine_sse::GetVolt(0,pos[0],last_A_Line,pos[2]) );
			m_Eng_SSE->Engine_sse::SetVolt(2,pos[0],0,pos[2], m_Eng_SSE->Engine_sse::GetVolt(2,pos[0],last_A_Line,pos[2]) );
		}
	}
}

void Engine_Ext_Cylinder::DoPostCurrentUpdates()
{
	if (CC_closedAlpha==false) return;

	//close alpha
	unsigned int pos[3];
	// copy tangential currents from first alpha-plane to last
	for (pos[0]=0; pos[0]<numLines[0]-1; ++pos[0])
	{
		unsigned int last_A_Line = numLines[1]-2;
		for (pos[2]=0; pos[2]<numLines[2]-1; ++pos[2])
		{
			m_Eng_SSE->Engine_sse::SetCurr(0,pos[0],last_A_Line,pos[2], m_Eng_SSE->Engine_sse::GetCurr(0,pos[0],0,pos[2]) );
			m_Eng_SSE->Engine_sse::SetCurr(2,pos[0],last_A_Line,pos[2], m_Eng_SSE->Engine_sse::GetCurr(2,pos[0],0,pos[2]) );
		}
	}
}
