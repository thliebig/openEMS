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

using namespace std;

Operator_Ext_Dispersive::Operator_Ext_Dispersive(Operator* op) : Operator_Extension(op)
{
	m_curr_ADE_On = NULL;
	m_volt_ADE_On = NULL;

	m_LM_pos=NULL;
	m_curr_ADE_On=NULL;
	m_volt_ADE_On=NULL;

	m_Order = 0;
}

Operator_Ext_Dispersive::Operator_Ext_Dispersive(Operator* op, Operator_Ext_Dispersive* op_ext) : Operator_Extension(op,op_ext)
{
	m_curr_ADE_On = NULL;
	m_volt_ADE_On = NULL;

	m_LM_pos=NULL;
	m_curr_ADE_On=NULL;
	m_volt_ADE_On=NULL;

	m_Order = 0;
}

Operator_Ext_Dispersive::~Operator_Ext_Dispersive()
{
	delete[] m_curr_ADE_On;
	delete[] m_volt_ADE_On;
	m_curr_ADE_On=NULL;
	m_volt_ADE_On=NULL;

	for (int n=0;n<m_Order;++n)
	{
		delete[] m_LM_pos[n][0];
		delete[] m_LM_pos[n][1];
		delete[] m_LM_pos[n][2];
	}
	delete[] m_LM_pos;
	m_LM_pos=NULL;
	m_Order=0;
	m_LM_Count.clear();
}

void Operator_Ext_Dispersive::ShowStat(ostream &ostr)  const
{
	Operator_Extension::ShowStat(ostr);
	string On_Off[2] = {"Off", "On"};
	ostr << " Max. Dispersion Order N = " << m_Order << endl;
	for (int i=0;i<m_Order;++i)
	{
		ostr << " N=" << i << ":\t Active cells\t\t: " << 	m_LM_Count.at(i) << endl;
		ostr << " N=" << i << ":\t Voltage ADE is \t: " << On_Off[m_volt_ADE_On[i]] << endl;
		ostr << " N=" << i << ":\t Current ADE is \t: " << On_Off[m_curr_ADE_On[i]] << endl;
	}
}
