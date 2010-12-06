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

#include "engine_ext_lorentzmaterial.h"
#include "operator_ext_lorentzmaterial.h"
#include "FDTD/engine_sse.h"

Engine_Ext_LorentzMaterial::Engine_Ext_LorentzMaterial(Operator_Ext_LorentzMaterial* op_ext_lorentz) : Engine_Ext_Dispersive(op_ext_lorentz)
{
	m_Op_Ext_Lor = op_ext_lorentz;
}

Engine_Ext_LorentzMaterial::~Engine_Ext_LorentzMaterial()
{

}

void Engine_Ext_LorentzMaterial::DoPreVoltageUpdates()
{
	if (m_Op_Ext_Lor->m_volt_ADE_On==false) return;

	unsigned int **pos = m_Op_Ext_Lor->m_LM_pos;

	//switch for different engine types to access faster inline engine functions
	switch (m_Eng->GetType())
	{
	case Engine::BASIC:
		{
			for (unsigned int i=0; i<m_Op_Ext_Lor->m_LM_Count; ++i)
			{
				volt_ADE[0][i] *= m_Op_Ext_Lor->v_int_ADE[0][i];
				volt_ADE[0][i] += m_Op_Ext_Lor->v_ext_ADE[0][i] * m_Eng->Engine::GetVolt(0,pos[0][i],pos[1][i],pos[2][i]);

				volt_ADE[1][i] *= m_Op_Ext_Lor->v_int_ADE[1][i];
				volt_ADE[1][i] += m_Op_Ext_Lor->v_ext_ADE[1][i] * m_Eng->Engine::GetVolt(1,pos[0][i],pos[1][i],pos[2][i]);

				volt_ADE[2][i] *= m_Op_Ext_Lor->v_int_ADE[2][i];
				volt_ADE[2][i] += m_Op_Ext_Lor->v_ext_ADE[2][i] * m_Eng->Engine::GetVolt(2,pos[0][i],pos[1][i],pos[2][i]);
			}
			break;
		}
	case Engine::SSE:
		{
			Engine_sse* eng_sse = (Engine_sse*)m_Eng;
			for (unsigned int i=0; i<m_Op_Ext_Lor->m_LM_Count; ++i)
			{
				volt_ADE[0][i] *= m_Op_Ext_Lor->v_int_ADE[0][i];
				volt_ADE[0][i] += m_Op_Ext_Lor->v_ext_ADE[0][i] * eng_sse->Engine_sse::GetVolt(0,pos[0][i],pos[1][i],pos[2][i]);

				volt_ADE[1][i] *= m_Op_Ext_Lor->v_int_ADE[1][i];
				volt_ADE[1][i] += m_Op_Ext_Lor->v_ext_ADE[1][i] * eng_sse->Engine_sse::GetVolt(1,pos[0][i],pos[1][i],pos[2][i]);

				volt_ADE[2][i] *= m_Op_Ext_Lor->v_int_ADE[2][i];
				volt_ADE[2][i] += m_Op_Ext_Lor->v_ext_ADE[2][i] * eng_sse->Engine_sse::GetVolt(2,pos[0][i],pos[1][i],pos[2][i]);
			}
			break;
		}
	default:
		for (unsigned int i=0; i<m_Op_Ext_Lor->m_LM_Count; ++i)
		{
			volt_ADE[0][i] *= m_Op_Ext_Lor->v_int_ADE[0][i];
			volt_ADE[0][i] += m_Op_Ext_Lor->v_ext_ADE[0][i] * m_Eng->GetVolt(0,pos[0][i],pos[1][i],pos[2][i]);

			volt_ADE[1][i] *= m_Op_Ext_Lor->v_int_ADE[1][i];
			volt_ADE[1][i] += m_Op_Ext_Lor->v_ext_ADE[1][i] * m_Eng->GetVolt(1,pos[0][i],pos[1][i],pos[2][i]);

			volt_ADE[2][i] *= m_Op_Ext_Lor->v_int_ADE[2][i];
			volt_ADE[2][i] += m_Op_Ext_Lor->v_ext_ADE[2][i] * m_Eng->GetVolt(2,pos[0][i],pos[1][i],pos[2][i]);
		}
		break;
	}
}

void Engine_Ext_LorentzMaterial::DoPreCurrentUpdates()
{
	if (m_Op_Ext_Lor->m_curr_ADE_On==false) return;

	unsigned int **pos = m_Op_Ext_Lor->m_LM_pos;

	//switch for different engine types to access faster inline engine functions
	switch (m_Eng->GetType())
	{
	case Engine::BASIC:
		{
			for (unsigned int i=0; i<m_Op_Ext_Lor->m_LM_Count; ++i)
			{
				curr_ADE[0][i] *= m_Op_Ext_Lor->i_int_ADE[0][i];
				curr_ADE[0][i] += m_Op_Ext_Lor->i_ext_ADE[0][i] * m_Eng->Engine::GetCurr(0,pos[0][i],pos[1][i],pos[2][i]);

				curr_ADE[1][i] *= m_Op_Ext_Lor->i_int_ADE[1][i];
				curr_ADE[1][i] += m_Op_Ext_Lor->i_ext_ADE[1][i] * m_Eng->Engine::GetCurr(1,pos[0][i],pos[1][i],pos[2][i]);

				curr_ADE[2][i] *= m_Op_Ext_Lor->i_int_ADE[2][i];
				curr_ADE[2][i] += m_Op_Ext_Lor->i_ext_ADE[2][i] * m_Eng->Engine::GetCurr(2,pos[0][i],pos[1][i],pos[2][i]);
			}
			break;
		}
	case Engine::SSE:
		{
			Engine_sse* eng_sse = (Engine_sse*)m_Eng;
			for (unsigned int i=0; i<m_Op_Ext_Lor->m_LM_Count; ++i)
			{
				curr_ADE[0][i] *= m_Op_Ext_Lor->i_int_ADE[0][i];
				curr_ADE[0][i] += m_Op_Ext_Lor->i_ext_ADE[0][i] * eng_sse->Engine_sse::GetCurr(0,pos[0][i],pos[1][i],pos[2][i]);

				curr_ADE[1][i] *= m_Op_Ext_Lor->i_int_ADE[1][i];
				curr_ADE[1][i] += m_Op_Ext_Lor->i_ext_ADE[1][i] * eng_sse->Engine_sse::GetCurr(1,pos[0][i],pos[1][i],pos[2][i]);

				curr_ADE[2][i] *= m_Op_Ext_Lor->i_int_ADE[2][i];
				curr_ADE[2][i] += m_Op_Ext_Lor->i_ext_ADE[2][i] * eng_sse->Engine_sse::GetCurr(2,pos[0][i],pos[1][i],pos[2][i]);
			}
			break;
		}
	default:
		for (unsigned int i=0; i<m_Op_Ext_Lor->m_LM_Count; ++i)
		{
			curr_ADE[0][i] *= m_Op_Ext_Lor->i_int_ADE[0][i];
			curr_ADE[0][i] += m_Op_Ext_Lor->i_ext_ADE[0][i] * m_Eng->GetCurr(0,pos[0][i],pos[1][i],pos[2][i]);

			curr_ADE[1][i] *= m_Op_Ext_Lor->i_int_ADE[1][i];
			curr_ADE[1][i] += m_Op_Ext_Lor->i_ext_ADE[1][i] * m_Eng->GetCurr(1,pos[0][i],pos[1][i],pos[2][i]);

			curr_ADE[2][i] *= m_Op_Ext_Lor->i_int_ADE[2][i];
			curr_ADE[2][i] += m_Op_Ext_Lor->i_ext_ADE[2][i] * m_Eng->GetCurr(2,pos[0][i],pos[1][i],pos[2][i]);
		}
		break;
	}
}

