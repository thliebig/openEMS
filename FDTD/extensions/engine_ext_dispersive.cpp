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

#include "engine_ext_dispersive.h"
#include "operator_ext_dispersive.h"
#include "FDTD/engine_sse.h"

Engine_Ext_Dispersive::Engine_Ext_Dispersive(Operator_Ext_Dispersive* op_ext_disp) : Engine_Extension(op_ext_disp)
{
	m_Op_Ext_Disp = op_ext_disp;
	for (int n=0;n<3;++n)
	{
		if (m_Op_Ext_Disp->m_curr_ADE_On==true)
		{
			curr_ADE[n] = new FDTD_FLOAT[m_Op_Ext_Disp->m_LM_Count];
			for (unsigned int i=0;i<m_Op_Ext_Disp->m_LM_Count;++i)
				curr_ADE[n][i]=0.0;
		}
		else
			curr_ADE[n] = NULL;
		if (m_Op_Ext_Disp->m_volt_ADE_On==true)
		{
			volt_ADE[n] = new FDTD_FLOAT[m_Op_Ext_Disp->m_LM_Count];
			for (unsigned int i=0;i<m_Op_Ext_Disp->m_LM_Count;++i)
				volt_ADE[n][i]=0.0;
		}
		else
			volt_ADE[n] = NULL;
	}
}

Engine_Ext_Dispersive::~Engine_Ext_Dispersive()
{
	for (int n=0;n<3;++n)
	{
		delete[] curr_ADE[n];
		curr_ADE[n] = NULL;
		delete[] volt_ADE[n];
		volt_ADE[n] = NULL;
	}
}

void Engine_Ext_Dispersive::Apply2Voltages()
{
	if (m_Op_Ext_Disp->m_volt_ADE_On==false) return;

	unsigned int **pos = m_Op_Ext_Disp->m_LM_pos;

	//switch for different engine types to access faster inline engine functions
	switch (m_Eng->GetType())
	{
	case Engine::BASIC:
		{
			for (unsigned int i=0;i<m_Op_Ext_Disp->m_LM_Count;++i)
			{
				m_Eng->Engine::SetVolt(0,pos[0][i],pos[1][i],pos[2][i], m_Eng->Engine::GetVolt(0,pos[0][i],pos[1][i],pos[2][i]) - volt_ADE[0][i]);
				m_Eng->Engine::SetVolt(1,pos[0][i],pos[1][i],pos[2][i], m_Eng->Engine::GetVolt(1,pos[0][i],pos[1][i],pos[2][i]) - volt_ADE[1][i]);
				m_Eng->Engine::SetVolt(2,pos[0][i],pos[1][i],pos[2][i], m_Eng->Engine::GetVolt(2,pos[0][i],pos[1][i],pos[2][i]) - volt_ADE[2][i]);
			}
			break;
		}
	case Engine::SSE:
		{
			Engine_sse* eng_sse = (Engine_sse*)m_Eng;
			for (unsigned int i=0;i<m_Op_Ext_Disp->m_LM_Count;++i)
			{
				eng_sse->Engine_sse::SetVolt(0,pos[0][i],pos[1][i],pos[2][i], eng_sse->Engine_sse::GetVolt(0,pos[0][i],pos[1][i],pos[2][i]) - volt_ADE[0][i]);
				eng_sse->Engine_sse::SetVolt(1,pos[0][i],pos[1][i],pos[2][i], eng_sse->Engine_sse::GetVolt(1,pos[0][i],pos[1][i],pos[2][i]) - volt_ADE[1][i]);
				eng_sse->Engine_sse::SetVolt(2,pos[0][i],pos[1][i],pos[2][i], eng_sse->Engine_sse::GetVolt(2,pos[0][i],pos[1][i],pos[2][i]) - volt_ADE[2][i]);
			}
			break;
		}
	default:
		for (unsigned int i=0;i<m_Op_Ext_Disp->m_LM_Count;++i)
		{
			m_Eng->SetVolt(0,pos[0][i],pos[1][i],pos[2][i], m_Eng->GetVolt(0,pos[0][i],pos[1][i],pos[2][i]) - volt_ADE[0][i]);
			m_Eng->SetVolt(1,pos[0][i],pos[1][i],pos[2][i], m_Eng->GetVolt(1,pos[0][i],pos[1][i],pos[2][i]) - volt_ADE[1][i]);
			m_Eng->SetVolt(2,pos[0][i],pos[1][i],pos[2][i], m_Eng->GetVolt(2,pos[0][i],pos[1][i],pos[2][i]) - volt_ADE[2][i]);
		}
		break;
	}
}

void Engine_Ext_Dispersive::Apply2Current()
{
	if (m_Op_Ext_Disp->m_curr_ADE_On==false) return;

	unsigned int **pos = m_Op_Ext_Disp->m_LM_pos;

	//switch for different engine types to access faster inline engine functions
	switch (m_Eng->GetType())
	{
	case Engine::BASIC:
		{
			for (unsigned int i=0;i<m_Op_Ext_Disp->m_LM_Count;++i)
			{
				m_Eng->Engine::SetCurr(0,pos[0][i],pos[1][i],pos[2][i], m_Eng->Engine::GetCurr(0,pos[0][i],pos[1][i],pos[2][i]) - curr_ADE[0][i]);
				m_Eng->Engine::SetCurr(1,pos[0][i],pos[1][i],pos[2][i], m_Eng->Engine::GetCurr(1,pos[0][i],pos[1][i],pos[2][i]) - curr_ADE[1][i]);
				m_Eng->Engine::SetCurr(2,pos[0][i],pos[1][i],pos[2][i], m_Eng->Engine::GetCurr(2,pos[0][i],pos[1][i],pos[2][i]) - curr_ADE[2][i]);
			}
			break;
		}
	case Engine::SSE:
		{
			Engine_sse* eng_sse = (Engine_sse*)m_Eng;
			for (unsigned int i=0;i<m_Op_Ext_Disp->m_LM_Count;++i)
			{
				eng_sse->Engine_sse::SetCurr(0,pos[0][i],pos[1][i],pos[2][i], eng_sse->Engine_sse::GetCurr(0,pos[0][i],pos[1][i],pos[2][i]) - curr_ADE[0][i]);
				eng_sse->Engine_sse::SetCurr(1,pos[0][i],pos[1][i],pos[2][i], eng_sse->Engine_sse::GetCurr(1,pos[0][i],pos[1][i],pos[2][i]) - curr_ADE[1][i]);
				eng_sse->Engine_sse::SetCurr(2,pos[0][i],pos[1][i],pos[2][i], eng_sse->Engine_sse::GetCurr(2,pos[0][i],pos[1][i],pos[2][i]) - curr_ADE[2][i]);
			}
			break;
		}
	default:
		for (unsigned int i=0;i<m_Op_Ext_Disp->m_LM_Count;++i)
		{
			m_Eng->SetCurr(0,pos[0][i],pos[1][i],pos[2][i], m_Eng->GetCurr(0,pos[0][i],pos[1][i],pos[2][i]) - curr_ADE[0][i]);
			m_Eng->SetCurr(1,pos[0][i],pos[1][i],pos[2][i], m_Eng->GetCurr(1,pos[0][i],pos[1][i],pos[2][i]) - curr_ADE[1][i]);
			m_Eng->SetCurr(2,pos[0][i],pos[1][i],pos[2][i], m_Eng->GetCurr(2,pos[0][i],pos[1][i],pos[2][i]) - curr_ADE[2][i]);
		}
		break;
	}

}
