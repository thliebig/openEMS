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

#include "operator_ext_mur_abc.h"
#include "engine_ext_mur_abc.h"

#include "tools/array_ops.h"

Operator_Ext_Mur_ABC::Operator_Ext_Mur_ABC(Operator* op) : Operator_Extension(op)
{
	m_Mur_Coeff = NULL;
	m_ny = -1;
	m_nyP = -1;
	m_nyPP = -1;
	m_LineNr = 0;
	m_LineNr_Shift = 0;
	m_Mur_Coeff = 0;
}

Operator_Ext_Mur_ABC::~Operator_Ext_Mur_ABC()
{
}

void Operator_Ext_Mur_ABC::SetDirection(int ny, bool top_ny)
{
	if ((ny<0) || (ny>2))
		return;
	m_ny = ny;
	m_nyP = (ny+1)%3;
	m_nyPP = (ny+2)%3;
	if (!top_ny)
	{
		m_LineNr = 0;
		m_LineNr_Shift = 1;
	}
	else
	{
		m_LineNr = m_Op->GetNumberOfLines(m_ny)-1;
		m_LineNr_Shift = m_Op->GetNumberOfLines(m_ny) - 2;
	}

	m_numLines[0] = m_Op->GetNumberOfLines(m_nyP);
	m_numLines[1] = m_Op->GetNumberOfLines(m_nyPP);
}

bool Operator_Ext_Mur_ABC::BuildExtension()
{
	if (m_ny<0)
	{
		cerr << "Operator_Ext_Mur_ABC::BuildExtension: Warning, Extension not initialized! Use SetDirection!! Abort build!!" << endl;
		return false;
	}
	double dT = m_Op->GetTimestep();
	unsigned int pos[] = {0,0,0};
	pos[m_ny] = m_LineNr;
	double delta = fabs(m_Op->GetMeshDelta(m_ny,pos));
	m_Mur_Coeff = (__C0__ * dT - delta) / (__C0__ * dT + delta);
//	cerr << "Operator_Ext_Mur_ABC::BuildExtension(): " << m_Mur_Coeff << endl;
	return true;
}

Engine_Extension* Operator_Ext_Mur_ABC::CreateEngineExtention()
{
	Engine_Ext_Mur_ABC* eng_ext = new Engine_Ext_Mur_ABC(this);
	return eng_ext;
}
