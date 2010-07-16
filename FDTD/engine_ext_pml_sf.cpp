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

#include "engine_ext_pml_sf.h"
#include "operator_ext_pml_sf.h"
#include "engine.h"
#include "engine_sse.h"
#include "tools/array_ops.h"

/************************************************  Engine_Ext_PML_SF  **************************************************************************/
Engine_Ext_PML_SF::Engine_Ext_PML_SF(Operator_Ext_PML_SF* op_ext) : Engine_Extension(op_ext)
{
	m_Op_PML_SF = op_ext;

	volt[0] = Create_N_3DArray<FDTD_FLOAT>(m_Op_PML_SF->m_numLines);
	volt[1] = Create_N_3DArray<FDTD_FLOAT>(m_Op_PML_SF->m_numLines);
	curr[0] = Create_N_3DArray<FDTD_FLOAT>(m_Op_PML_SF->m_numLines);
	curr[1] = Create_N_3DArray<FDTD_FLOAT>(m_Op_PML_SF->m_numLines);
}

Engine_Ext_PML_SF::~Engine_Ext_PML_SF()
{
	Delete_N_3DArray<FDTD_FLOAT>(volt[0],m_Op_PML_SF->m_numLines);
	volt[0]=NULL;
	Delete_N_3DArray<FDTD_FLOAT>(volt[1],m_Op_PML_SF->m_numLines);
	volt[1]=NULL;
	Delete_N_3DArray<FDTD_FLOAT>(curr[0],m_Op_PML_SF->m_numLines);
	curr[0]=NULL;
	Delete_N_3DArray<FDTD_FLOAT>(curr[1],m_Op_PML_SF->m_numLines);
	curr[1]=NULL;
}

void Engine_Ext_PML_SF::UpdateVoltages(unsigned int startX, unsigned int numX)
{
	unsigned int pos[3];
	bool shift[3];

	pos[0] = startX;
	//voltage updates
	for (unsigned int posX=0;posX<numX;++posX)
	{
		shift[0]=pos[0];
		for (pos[1]=0;pos[1]<m_Op_PML_SF->m_numLines[1];++pos[1])
		{
			shift[1]=pos[1];
			for (pos[2]=0;pos[2]<m_Op_PML_SF->m_numLines[2];++pos[2])
			{
				shift[2]=pos[2];
				//do the updates here
				//for x
				volt[0][0][pos[0]][pos[1]][pos[2]] *= m_Op_PML_SF->vv[0][0][pos[0]][pos[1]][pos[2]];
				volt[0][0][pos[0]][pos[1]][pos[2]] += m_Op_PML_SF->vi[0][0][pos[0]][pos[1]][pos[2]] * ( curr[0][2][pos[0]][pos[1]][pos[2]]		    + curr[1][2][pos[0]][pos[1]][pos[2]]			- curr[0][2][pos[0]][pos[1]-shift[1]][pos[2]]	- curr[1][2][pos[0]][pos[1]-shift[1]][pos[2]]);

				volt[1][0][pos[0]][pos[1]][pos[2]] *= m_Op_PML_SF->vv[1][0][pos[0]][pos[1]][pos[2]];
				volt[1][0][pos[0]][pos[1]][pos[2]] += m_Op_PML_SF->vi[1][0][pos[0]][pos[1]][pos[2]] * ( curr[0][1][pos[0]][pos[1]][pos[2]-shift[2]] + curr[1][1][pos[0]][pos[1]][pos[2]-shift[2]]	- curr[0][1][pos[0]][pos[1]][pos[2]]			- curr[1][1][pos[0]][pos[1]][pos[2]]);

				//for y
				volt[0][1][pos[0]][pos[1]][pos[2]] *= m_Op_PML_SF->vv[0][1][pos[0]][pos[1]][pos[2]];
				volt[0][1][pos[0]][pos[1]][pos[2]] += m_Op_PML_SF->vi[0][1][pos[0]][pos[1]][pos[2]] * ( curr[0][0][pos[0]][pos[1]][pos[2]]			+ curr[1][0][pos[0]][pos[1]][pos[2]]			- curr[0][0][pos[0]][pos[1]][pos[2]-shift[2]]	- curr[1][0][pos[0]][pos[1]][pos[2]-shift[2]]);

				volt[1][1][pos[0]][pos[1]][pos[2]] *= m_Op_PML_SF->vv[1][1][pos[0]][pos[1]][pos[2]];
				volt[1][1][pos[0]][pos[1]][pos[2]] += m_Op_PML_SF->vi[1][1][pos[0]][pos[1]][pos[2]] * ( curr[0][2][pos[0]-shift[0]][pos[1]][pos[2]]	+ curr[1][2][pos[0]-shift[0]][pos[1]][pos[2]]	- curr[0][2][pos[0]][pos[1]][pos[2]]			- curr[1][2][pos[0]][pos[1]][pos[2]]);

				//for z
				volt[0][2][pos[0]][pos[1]][pos[2]] *= m_Op_PML_SF->vv[0][2][pos[0]][pos[1]][pos[2]];
				volt[0][2][pos[0]][pos[1]][pos[2]] += m_Op_PML_SF->vi[0][2][pos[0]][pos[1]][pos[2]] * ( curr[0][1][pos[0]][pos[1]][pos[2]]			+ curr[1][1][pos[0]][pos[1]][pos[2]]			- curr[0][1][pos[0]-shift[0]][pos[1]][pos[2]]	- curr[1][1][pos[0]-shift[0]][pos[1]][pos[2]]);

				volt[1][2][pos[0]][pos[1]][pos[2]] *= m_Op_PML_SF->vv[1][2][pos[0]][pos[1]][pos[2]];
				volt[1][2][pos[0]][pos[1]][pos[2]] += m_Op_PML_SF->vi[1][2][pos[0]][pos[1]][pos[2]] * ( curr[0][0][pos[0]][pos[1]-shift[1]][pos[2]]	+ curr[1][0][pos[0]][pos[1]-shift[1]][pos[2]]	- curr[0][0][pos[0]][pos[1]][pos[2]]			- curr[1][0][pos[0]][pos[1]][pos[2]]);
			}
		}
		++pos[0];
	}
}

void Engine_Ext_PML_SF::UpdateCurrents(unsigned int startX, unsigned int numX)
{
	unsigned int pos[3];
	pos[0] = startX;
	for (unsigned int posX=0;posX<numX;++posX)
	{
		for (pos[1]=0;pos[1]<m_Op_PML_SF->m_numLines[1]-1;++pos[1])
		{
			for (pos[2]=0;pos[2]<m_Op_PML_SF->m_numLines[2]-1;++pos[2])
			{
				//do the updates here
				//for x
				curr[0][0][pos[0]][pos[1]][pos[2]] *= m_Op_PML_SF->ii[0][0][pos[0]][pos[1]][pos[2]];
				curr[0][0][pos[0]][pos[1]][pos[2]] += m_Op_PML_SF->iv[0][0][pos[0]][pos[1]][pos[2]] * ( volt[0][2][pos[0]][pos[1]][pos[2]]		+ volt[1][2][pos[0]][pos[1]][pos[2]]	- volt[0][2][pos[0]][pos[1]+1][pos[2]]	- volt[1][2][pos[0]][pos[1]+1][pos[2]]);

				curr[1][0][pos[0]][pos[1]][pos[2]] *= m_Op_PML_SF->ii[1][0][pos[0]][pos[1]][pos[2]];
				curr[1][0][pos[0]][pos[1]][pos[2]] += m_Op_PML_SF->iv[1][0][pos[0]][pos[1]][pos[2]] * ( volt[0][1][pos[0]][pos[1]][pos[2]+1]	+ volt[1][1][pos[0]][pos[1]][pos[2]+1]	- volt[0][1][pos[0]][pos[1]][pos[2]]	- volt[1][1][pos[0]][pos[1]][pos[2]]);

//				cerr << volt[0][0][pos[0]][pos[1]][pos[2]] << " ";
//				cerr << volt[0][1][pos[0]][pos[1]][pos[2]] << " ";
//				cerr << volt[0][2][pos[0]][pos[1]][pos[2]] << endl;

				//for y
				curr[0][1][pos[0]][pos[1]][pos[2]] *= m_Op_PML_SF->ii[0][1][pos[0]][pos[1]][pos[2]];
				curr[0][1][pos[0]][pos[1]][pos[2]] += m_Op_PML_SF->iv[0][1][pos[0]][pos[1]][pos[2]] * ( volt[0][0][pos[0]][pos[1]][pos[2]]		+ volt[1][0][pos[0]][pos[1]][pos[2]]	- volt[0][0][pos[0]][pos[1]][pos[2]+1]	- volt[1][0][pos[0]][pos[1]][pos[2]+1]);

				curr[1][1][pos[0]][pos[1]][pos[2]] *= m_Op_PML_SF->ii[1][1][pos[0]][pos[1]][pos[2]];
				curr[1][1][pos[0]][pos[1]][pos[2]] += m_Op_PML_SF->iv[1][1][pos[0]][pos[1]][pos[2]] * ( volt[0][2][pos[0]+1][pos[1]][pos[2]]	+ volt[1][2][pos[0]+1][pos[1]][pos[2]]	- volt[0][2][pos[0]][pos[1]][pos[2]]	- volt[1][2][pos[0]][pos[1]][pos[2]]);

				//for z
				curr[0][2][pos[0]][pos[1]][pos[2]] *= m_Op_PML_SF->ii[0][2][pos[0]][pos[1]][pos[2]];
				curr[0][2][pos[0]][pos[1]][pos[2]] += m_Op_PML_SF->iv[0][2][pos[0]][pos[1]][pos[2]] * ( volt[0][1][pos[0]][pos[1]][pos[2]]		+ volt[1][1][pos[0]][pos[1]][pos[2]]	- volt[0][1][pos[0]+1][pos[1]][pos[2]]	- volt[1][1][pos[0]+1][pos[1]][pos[2]]);

				curr[1][2][pos[0]][pos[1]][pos[2]] *= m_Op_PML_SF->ii[1][2][pos[0]][pos[1]][pos[2]];
				curr[1][2][pos[0]][pos[1]][pos[2]] += m_Op_PML_SF->iv[1][2][pos[0]][pos[1]][pos[2]] * ( volt[0][0][pos[0]][pos[1]+1][pos[2]]	+ volt[1][0][pos[0]][pos[1]+1][pos[2]]	- volt[0][0][pos[0]][pos[1]][pos[2]]	- volt[1][0][pos[0]][pos[1]][pos[2]]);
			}
		}
		++pos[0];
	}
}

/************************************************  Engine_Ext_PML_SF_Plane  **************************************************************************/
Engine_Ext_PML_SF_Plane::Engine_Ext_PML_SF_Plane(Operator_Ext_PML_SF_Plane* op_ext) : Engine_Ext_PML_SF(op_ext)
{
	m_Op_PML_SF_PL = op_ext;
}

Engine_Ext_PML_SF_Plane::~Engine_Ext_PML_SF_Plane()
{
}

void Engine_Ext_PML_SF_Plane::Apply2Voltages()
{
	unsigned int pos[] = {0,0,0};
	unsigned int pml_pos[] = {0,0,0};
	unsigned int m_ny = m_Op_PML_SF_PL->m_ny;
	unsigned int m_nyP = m_Op_PML_SF_PL->m_nyP;
	unsigned int m_nyPP = m_Op_PML_SF_PL->m_nyPP;

	pos[m_ny] = 0;
	pml_pos[m_ny] = m_Op_PML_SF_PL->m_numLines[m_ny]-1;
	// copy currents data from main engine to pml engine (lowest main line to highest pml line)
	if (m_Op_PML_SF_PL->m_top==false)
	{
		for (pos[m_nyP]=0;pos[m_nyP]<m_Op_PML_SF_PL->m_numLines[m_nyP]-1;++pos[m_nyP])
		{
			pml_pos[m_nyP] = pos[m_nyP];
			for (pos[m_nyPP]=0;pos[m_nyPP]<m_Op_PML_SF_PL->m_numLines[m_nyPP]-1;++pos[m_nyPP])
			{
				pml_pos[m_nyPP] = pos[m_nyPP];
				curr[0][0][pml_pos[0]][pml_pos[1]][pml_pos[2]] = m_Eng->GetCurr(0,pos);
				curr[1][0][pml_pos[0]][pml_pos[1]][pml_pos[2]] = 0;
				curr[0][1][pml_pos[0]][pml_pos[1]][pml_pos[2]] = m_Eng->GetCurr(1,pos);
				curr[1][1][pml_pos[0]][pml_pos[1]][pml_pos[2]] = 0;
				curr[0][2][pml_pos[0]][pml_pos[1]][pml_pos[2]] = m_Eng->GetCurr(2,pos);
				curr[1][2][pml_pos[0]][pml_pos[1]][pml_pos[2]] = 0;
			}
		}
	}

	// do the voltage updates
	UpdateVoltages(0,m_Op_PML_SF->m_numLines[0]);

	if (m_Op_PML_SF_PL->m_top==false)
	{
		// copy voltage data from pml engine to main engine (highest pml line to lowest main line)
		pos[m_ny] = 0;
		pml_pos[m_ny] = m_Op_PML_SF_PL->m_numLines[m_ny]-1;
		for (pos[m_nyP]=0;pos[m_nyP]<m_Op_PML_SF_PL->m_numLines[m_nyP];++pos[m_nyP])
		{
			pml_pos[m_nyP] = pos[m_nyP];
			for (pos[m_nyPP]=0;pos[m_nyPP]<m_Op_PML_SF_PL->m_numLines[m_nyPP];++pos[m_nyPP])
			{
				pml_pos[m_nyPP] = pos[m_nyPP];
				m_Eng->GetVolt(0,pos) = volt[0][0][pml_pos[0]][pml_pos[1]][pml_pos[2]] + volt[1][0][pml_pos[0]][pml_pos[1]][pml_pos[2]];
				m_Eng->GetVolt(1,pos) = volt[0][1][pml_pos[0]][pml_pos[1]][pml_pos[2]] + volt[1][1][pml_pos[0]][pml_pos[1]][pml_pos[2]];
				m_Eng->GetVolt(2,pos) = volt[0][2][pml_pos[0]][pml_pos[1]][pml_pos[2]] + volt[1][2][pml_pos[0]][pml_pos[1]][pml_pos[2]];
			}
		}
	}
}

void Engine_Ext_PML_SF_Plane::Apply2Current()
{
	unsigned int pos[] = {0,0,0};
	unsigned int pml_pos[] = {0,0,0};
	unsigned int m_ny = m_Op_PML_SF_PL->m_ny;
	unsigned int m_nyP = m_Op_PML_SF_PL->m_nyP;
	unsigned int m_nyPP = m_Op_PML_SF_PL->m_nyPP;

	pos[m_ny] = m_Op_PML_SF_PL->m_LineNr;
	pml_pos[m_ny] = 0;
	// copy voltage data from main engine to pml engine (highest main line to lowest pml line)
	if (m_Op_PML_SF_PL->m_top)
	{
		for (pos[m_nyP]=0;pos[m_nyP]<m_Op_PML_SF_PL->m_numLines[m_nyP];++pos[m_nyP])
		{
			pml_pos[m_nyP] = pos[m_nyP];
			for (pos[m_nyPP]=0;pos[m_nyPP]<m_Op_PML_SF_PL->m_numLines[m_nyPP];++pos[m_nyPP])
			{
				pml_pos[m_nyPP] = pos[m_nyPP];
				volt[0][0][pml_pos[0]][pml_pos[1]][pml_pos[2]] = m_Eng->GetVolt(0,pos);
				volt[1][0][pml_pos[0]][pml_pos[1]][pml_pos[2]] = 0;
				volt[0][1][pml_pos[0]][pml_pos[1]][pml_pos[2]] = m_Eng->GetVolt(1,pos);
				volt[1][1][pml_pos[0]][pml_pos[1]][pml_pos[2]] = 0;
				volt[0][2][pml_pos[0]][pml_pos[1]][pml_pos[2]] = m_Eng->GetVolt(2,pos);
				volt[1][2][pml_pos[0]][pml_pos[1]][pml_pos[2]] = 0;
			}
		}
	}

	UpdateCurrents(0,m_Op_PML_SF->m_numLines[0]-1);

	pos[m_ny] = m_Op_PML_SF_PL->m_LineNr;
	pml_pos[m_ny] = 0;
	// copy currents data from pml engine to main engine (lowest pml line to highest main line)
	if (m_Op_PML_SF_PL->m_top)
	{
		for (pos[m_nyP]=0;pos[m_nyP]<m_Op_PML_SF_PL->m_numLines[m_nyP]-1;++pos[m_nyP])
		{
			pml_pos[m_nyP] = pos[m_nyP];
			for (pos[m_nyPP]=0;pos[m_nyPP]<m_Op_PML_SF_PL->m_numLines[m_nyPP]-1;++pos[m_nyPP])
			{
				pml_pos[m_nyPP] = pos[m_nyPP];

				m_Eng->GetCurr(0,pos) = curr[0][0][pml_pos[0]][pml_pos[1]][pml_pos[2]] + curr[1][0][pml_pos[0]][pml_pos[1]][pml_pos[2]];
				m_Eng->GetCurr(1,pos) = curr[0][1][pml_pos[0]][pml_pos[1]][pml_pos[2]] + curr[1][1][pml_pos[0]][pml_pos[1]][pml_pos[2]];
				m_Eng->GetCurr(2,pos) = curr[0][2][pml_pos[0]][pml_pos[1]][pml_pos[2]] + curr[1][2][pml_pos[0]][pml_pos[1]][pml_pos[2]];
			}
		}
	}
}
