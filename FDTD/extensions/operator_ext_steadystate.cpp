/*
*	Copyright (C) 2015 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#include "operator_ext_steadystate.h"
#include "engine_ext_steadystate.h"

Operator_Ext_SteadyState::Operator_Ext_SteadyState(Operator* op, double period): Operator_Extension(op)
{
	this->Reset();
	m_T_period = period;
}

Operator_Ext_SteadyState::Operator_Ext_SteadyState(Operator* op, Operator_Ext_SteadyState* op_ext): Operator_Extension(op, op_ext)
{
	this->Reset();
	m_T_period = op_ext->m_T_period;
}

Operator_Ext_SteadyState::~Operator_Ext_SteadyState()
{
}

Operator_Extension* Operator_Ext_SteadyState::Clone(Operator* op)
{
	//disable cloning, only the main operator needs to have a steady state detection
	UNUSED(op);
	return NULL;
}

bool Operator_Ext_SteadyState::BuildExtension()
{
	double dT = m_Op->GetTimestep();
	m_TS_period = round(m_T_period/dT);
	m_T_period  = m_TS_period*dT;
	return true;
}

void Operator_Ext_SteadyState::Reset()
{
	for (int n=0;n<3;++n)
	{
		m_E_probe_pos[n].clear();
		m_H_probe_pos[n].clear();
	}
	m_E_probe_dir.clear();
	m_H_probe_dir.clear();
	m_T_period  = 0;
	m_TS_period = 0;
}

bool Operator_Ext_SteadyState::Add_E_Probe(unsigned int pos[3], int dir)
{
	if ((dir<0) || (dir>2))
		return false;
	for (int n=0;n<3;++n)
		if (pos[n]>=m_Op->GetNumberOfLines(n))
			return false;
	for (int n=0;n<3;++n)
		m_E_probe_pos[n].push_back(pos[n]);
	m_E_probe_dir.push_back(dir);
	return true;
}

bool Operator_Ext_SteadyState::Add_H_Probe(unsigned int pos[3], int dir)
{
	if ((dir<0) || (dir>2))
		return false;
	for (int n=0;n<3;++n)
		if (pos[n]>=m_Op->GetNumberOfLines(n))
			return false;
	for (int n=0;n<3;++n)
		m_H_probe_pos[n].push_back(pos[n]);
	m_H_probe_dir.push_back(dir);
	return true;
}

Engine_Extension* Operator_Ext_SteadyState::CreateEngineExtention()
{
	m_Eng_Ext = new Engine_Ext_SteadyState(this);
	return m_Eng_Ext;
}

void Operator_Ext_SteadyState::ShowStat(ostream &ostr)  const
{
	Operator_Extension::ShowStat(ostr);
	cout << "Period time (s): "  << m_T_period << "\t Period TS: " << m_TS_period << endl;
	cout << "Number of E probes\t: " << m_E_probe_dir.size() << endl;
	cout << "Number of H probes\t: " << m_H_probe_dir.size() << endl;
}

