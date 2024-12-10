/*
*	Copyright (C) 2011 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#include "engine_ext_excitation.h"
#include "operator_ext_excitation.h"
#include "FDTD/engine_sse.h"

Engine_Ext_Excitation::Engine_Ext_Excitation(Operator_Ext_Excitation* op_ext) : Engine_Extension(op_ext)
{
	m_Op_Exc = op_ext;
	m_Priority = ENG_EXT_PRIO_EXCITATION;
}

Engine_Ext_Excitation::~Engine_Ext_Excitation()
{

}

template <typename EngType>
void Engine_Ext_Excitation::Apply2VoltagesImpl(EngType* eng)
{
	//soft voltage excitation here (E-field excite)
	int exc_pos;
	unsigned int ny;
	unsigned int pos[3];
	int numTS = m_Eng->GetNumberOfTimesteps();
	unsigned int length = m_Op_Exc->m_Exc->GetLength();
	FDTD_FLOAT* exc_volt =  m_Op_Exc->m_Exc->GetVoltageSignal();

	int p = numTS+1;
	if (m_Op_Exc->m_Exc->GetSignalPeriod()>0)
		p = int(m_Op_Exc->m_Exc->GetSignalPeriod()/m_Op_Exc->m_Exc->GetTimestep());

	for (unsigned int n=0; n<m_Op_Exc->Volt_Count; ++n)
	{
		exc_pos = numTS - (int)m_Op_Exc->Volt_delay[n];
		exc_pos *= (exc_pos>0);
		exc_pos %= p;
		exc_pos *= (exc_pos<(int)length);
		ny = m_Op_Exc->Volt_dir[n];
		pos[0]=m_Op_Exc->Volt_index[0][n];
		pos[1]=m_Op_Exc->Volt_index[1][n];
		pos[2]=m_Op_Exc->Volt_index[2][n];
		eng->EngType::SetVolt(ny,pos, eng->EngType::GetVolt(ny,pos) + m_Op_Exc->Volt_amp[n]*exc_volt[exc_pos]);
	}
}

void Engine_Ext_Excitation::Apply2Voltages()
{
	ENG_DISPATCH(Apply2VoltagesImpl);
}

template <typename EngType>
void Engine_Ext_Excitation::Apply2CurrentImpl(EngType* eng)
{
	//soft current excitation here (H-field excite)
	int exc_pos;
	unsigned int ny;
	unsigned int pos[3];
	int numTS = m_Eng->GetNumberOfTimesteps();
	unsigned int length = m_Op_Exc->m_Exc->GetLength();
	FDTD_FLOAT* exc_curr =  m_Op_Exc->m_Exc->GetCurrentSignal();

	int p = numTS+1;
	if (m_Op_Exc->m_Exc->GetSignalPeriod()>0)
		p = int(m_Op_Exc->m_Exc->GetSignalPeriod()/m_Op_Exc->m_Exc->GetTimestep());

	for (unsigned int n=0; n<m_Op_Exc->Curr_Count; ++n)
	{
		exc_pos = numTS - (int)m_Op_Exc->Curr_delay[n];
		exc_pos *= (exc_pos>0);
		exc_pos %= p;
		exc_pos *= (exc_pos<(int)length);
		ny = m_Op_Exc->Curr_dir[n];
		pos[0]=m_Op_Exc->Curr_index[0][n];
		pos[1]=m_Op_Exc->Curr_index[1][n];
		pos[2]=m_Op_Exc->Curr_index[2][n];
		eng->EngType::SetCurr(ny,pos, eng->EngType::GetCurr(ny,pos) + m_Op_Exc->Curr_amp[n]*exc_curr[exc_pos]);
	}
}

void Engine_Ext_Excitation::Apply2Current()
{
	ENG_DISPATCH(Apply2CurrentImpl);
}
