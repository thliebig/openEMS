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

#include "operator_ext_dispersive.h"

#include "tools/array_ops.h"

Operator_Ext_Dispersive::Operator_Ext_Dispersive(Operator* op)  : Operator_Extension(op)
{
	m_curr_ADE_On = false;
	m_volt_ADE_On = false;

	m_LM_pos[0]=NULL;
	m_LM_pos[1]=NULL;
	m_LM_pos[2]=NULL;
}

Operator_Ext_Dispersive::~Operator_Ext_Dispersive()
{
	delete[] m_LM_pos[0];
	delete[] m_LM_pos[1];
	delete[] m_LM_pos[2];

	m_LM_pos[0]=NULL;
	m_LM_pos[1]=NULL;
	m_LM_pos[2]=NULL;
}

void Operator_Ext_Dispersive::ShowStat(ostream &ostr)  const
{
	Operator_Extension::ShowStat(ostr);
	string On_Off[2] = {"Off", "On"};
	ostr << " Active cells\t\t: " << 	m_LM_Count << endl;
	ostr << " Voltage ADE is \t: " << On_Off[m_volt_ADE_On] << endl;
	ostr << " Current ADE is \t: " << On_Off[m_curr_ADE_On] << endl;
}
