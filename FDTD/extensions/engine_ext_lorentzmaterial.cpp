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
	m_Order = m_Op_Ext_Lor->GetDispersionOrder();
	int order = m_Op_Ext_Lor->m_Order;

	curr_Lor_ADE = new FDTD_FLOAT**[order];
	volt_Lor_ADE = new FDTD_FLOAT**[order];
	for (int o=0;o<order;++o)
	{
		curr_Lor_ADE[o] = new FDTD_FLOAT*[3];
		volt_Lor_ADE[o] = new FDTD_FLOAT*[3];
		for (int n=0; n<3; ++n)
		{
			if (m_Op_Ext_Lor->m_curr_Lor_ADE_On[o]==true)
			{
				curr_Lor_ADE[o][n] = new FDTD_FLOAT[m_Op_Ext_Lor->m_LM_Count[o]];
				for (unsigned int i=0; i<m_Op_Ext_Lor->m_LM_Count[o]; ++i)
					curr_Lor_ADE[o][n][i]=0.0;
			}
			else
				curr_Lor_ADE[o][n] = NULL;

			if (m_Op_Ext_Lor->m_volt_Lor_ADE_On[o]==true)
			{
				volt_Lor_ADE[o][n] = new FDTD_FLOAT[m_Op_Ext_Lor->m_LM_Count[o]];
				for (unsigned int i=0; i<m_Op_Ext_Lor->m_LM_Count[o]; ++i)
					volt_Lor_ADE[o][n][i]=0.0;
			}
			else
				volt_Lor_ADE[o][n] = NULL;
		}
	}
}

Engine_Ext_LorentzMaterial::~Engine_Ext_LorentzMaterial()
{
	if (curr_Lor_ADE==NULL && volt_Lor_ADE==NULL)
		return;

	for (int o=0;o<m_Op_Ext_Lor->m_Order;++o)
	{
		for (int n=0; n<3; ++n)
		{
			delete[] curr_Lor_ADE[o][n];
			delete[] volt_Lor_ADE[o][n];
		}
		delete[] curr_Lor_ADE[o];
		delete[] volt_Lor_ADE[o];
	}
	delete[] curr_Lor_ADE;
	curr_Lor_ADE=NULL;

	delete[] volt_Lor_ADE;
	volt_Lor_ADE=NULL;
}

void Engine_Ext_LorentzMaterial::DoPreVoltageUpdates()
{
	for (int o=0;o<m_Order;++o)
	{
		if (m_Op_Ext_Lor->m_volt_ADE_On[o]==false) continue;

		unsigned int **pos = m_Op_Ext_Lor->m_LM_pos[o];

		if (m_Op_Ext_Lor->m_volt_Lor_ADE_On[o])
		{
			//switch for different engine types to access faster inline engine functions
			switch (m_Eng->GetType())
			{
			case Engine::BASIC:
			{
				for (unsigned int i=0; i<m_Op_Ext_Lor->m_LM_Count.at(o); ++i)
				{
					volt_Lor_ADE[o][0][i]+=m_Op_Ext_Lor->v_Lor_ADE[o][0][i]*volt_ADE[o][0][i];
					volt_ADE[o][0][i] *= m_Op_Ext_Lor->v_int_ADE[o][0][i];
					volt_ADE[o][0][i] += m_Op_Ext_Lor->v_ext_ADE[o][0][i] * (m_Eng->Engine::GetVolt(0,pos[0][i],pos[1][i],pos[2][i])-volt_Lor_ADE[o][0][i]);

					volt_Lor_ADE[o][1][i]+=m_Op_Ext_Lor->v_Lor_ADE[o][1][i]*volt_ADE[o][1][i];
					volt_ADE[o][1][i] *= m_Op_Ext_Lor->v_int_ADE[o][1][i];
					volt_ADE[o][1][i] += m_Op_Ext_Lor->v_ext_ADE[o][1][i] * (m_Eng->Engine::GetVolt(1,pos[0][i],pos[1][i],pos[2][i])-volt_Lor_ADE[o][2][i]);

					volt_Lor_ADE[o][2][i]+=m_Op_Ext_Lor->v_Lor_ADE[o][2][i]*volt_ADE[o][2][i];
					volt_ADE[o][2][i] *= m_Op_Ext_Lor->v_int_ADE[o][2][i];
					volt_ADE[o][2][i] += m_Op_Ext_Lor->v_ext_ADE[o][2][i] * (m_Eng->Engine::GetVolt(2,pos[0][i],pos[1][i],pos[2][i])-volt_Lor_ADE[o][2][i]);
				}
				break;
			}
			case Engine::SSE:
			{
				Engine_sse* eng_sse = (Engine_sse*)m_Eng;
				for (unsigned int i=0; i<m_Op_Ext_Lor->m_LM_Count.at(o); ++i)
				{
					volt_Lor_ADE[o][0][i]+=m_Op_Ext_Lor->v_Lor_ADE[o][0][i]*volt_ADE[o][0][i];
					volt_ADE[o][0][i] *= m_Op_Ext_Lor->v_int_ADE[o][0][i];
					volt_ADE[o][0][i] += m_Op_Ext_Lor->v_ext_ADE[o][0][i] * (eng_sse->Engine_sse::GetVolt(0,pos[0][i],pos[1][i],pos[2][i])-volt_Lor_ADE[o][0][i]);

					volt_Lor_ADE[o][1][i]+=m_Op_Ext_Lor->v_Lor_ADE[o][1][i]*volt_ADE[o][1][i];
					volt_ADE[o][1][i] *= m_Op_Ext_Lor->v_int_ADE[o][1][i];
					volt_ADE[o][1][i] += m_Op_Ext_Lor->v_ext_ADE[o][1][i] * (eng_sse->Engine_sse::GetVolt(1,pos[0][i],pos[1][i],pos[2][i])-volt_Lor_ADE[o][1][i]);

					volt_Lor_ADE[o][2][i]+=m_Op_Ext_Lor->v_Lor_ADE[o][2][i]*volt_ADE[o][2][i];
					volt_ADE[o][2][i] *= m_Op_Ext_Lor->v_int_ADE[o][2][i];
					volt_ADE[o][2][i] += m_Op_Ext_Lor->v_ext_ADE[o][2][i] * (eng_sse->Engine_sse::GetVolt(2,pos[0][i],pos[1][i],pos[2][i])-volt_Lor_ADE[o][2][i]);
				}
				break;
			}
			default:
				for (unsigned int i=0; i<m_Op_Ext_Lor->m_LM_Count.at(o); ++i)
				{
					volt_Lor_ADE[o][0][i]+=m_Op_Ext_Lor->v_Lor_ADE[o][0][i]*volt_ADE[o][0][i];
					volt_ADE[o][0][i] *= m_Op_Ext_Lor->v_int_ADE[o][0][i];
					volt_ADE[o][0][i] += m_Op_Ext_Lor->v_ext_ADE[o][0][i] * (m_Eng->GetVolt(0,pos[0][i],pos[1][i],pos[2][i])-volt_Lor_ADE[o][0][i]);

					volt_Lor_ADE[o][1][i]+=m_Op_Ext_Lor->v_Lor_ADE[o][1][i]*volt_ADE[o][1][i];
					volt_ADE[o][1][i] *= m_Op_Ext_Lor->v_int_ADE[o][1][i];
					volt_ADE[o][1][i] += m_Op_Ext_Lor->v_ext_ADE[o][1][i] * (m_Eng->GetVolt(1,pos[0][i],pos[1][i],pos[2][i])-volt_Lor_ADE[o][1][i]);

					volt_Lor_ADE[o][2][i]+=m_Op_Ext_Lor->v_Lor_ADE[o][2][i]*volt_ADE[o][2][i];
					volt_ADE[o][2][i] *= m_Op_Ext_Lor->v_int_ADE[o][2][i];
					volt_ADE[o][2][i] += m_Op_Ext_Lor->v_ext_ADE[o][2][i] * (m_Eng->GetVolt(2,pos[0][i],pos[1][i],pos[2][i])-volt_Lor_ADE[o][2][i]);
				}
				break;
			}
		}
		else
		{
			//switch for different engine types to access faster inline engine functions
			switch (m_Eng->GetType())
			{
			case Engine::BASIC:
			{
				for (unsigned int i=0; i<m_Op_Ext_Lor->m_LM_Count.at(o); ++i)
				{
					volt_ADE[o][0][i] *= m_Op_Ext_Lor->v_int_ADE[o][0][i];
					volt_ADE[o][0][i] += m_Op_Ext_Lor->v_ext_ADE[o][0][i] * m_Eng->Engine::GetVolt(0,pos[0][i],pos[1][i],pos[2][i]);

					volt_ADE[o][1][i] *= m_Op_Ext_Lor->v_int_ADE[o][1][i];
					volt_ADE[o][1][i] += m_Op_Ext_Lor->v_ext_ADE[o][1][i] * m_Eng->Engine::GetVolt(1,pos[0][i],pos[1][i],pos[2][i]);

					volt_ADE[o][2][i] *= m_Op_Ext_Lor->v_int_ADE[o][2][i];
					volt_ADE[o][2][i] += m_Op_Ext_Lor->v_ext_ADE[o][2][i] * m_Eng->Engine::GetVolt(2,pos[0][i],pos[1][i],pos[2][i]);
				}
				break;
			}
			case Engine::SSE:
			{
				Engine_sse* eng_sse = (Engine_sse*)m_Eng;
				for (unsigned int i=0; i<m_Op_Ext_Lor->m_LM_Count.at(o); ++i)
				{
					volt_ADE[o][0][i] *= m_Op_Ext_Lor->v_int_ADE[o][0][i];
					volt_ADE[o][0][i] += m_Op_Ext_Lor->v_ext_ADE[o][0][i] * eng_sse->Engine_sse::GetVolt(0,pos[0][i],pos[1][i],pos[2][i]);

					volt_ADE[o][1][i] *= m_Op_Ext_Lor->v_int_ADE[o][1][i];
					volt_ADE[o][1][i] += m_Op_Ext_Lor->v_ext_ADE[o][1][i] * eng_sse->Engine_sse::GetVolt(1,pos[0][i],pos[1][i],pos[2][i]);

					volt_ADE[o][2][i] *= m_Op_Ext_Lor->v_int_ADE[o][2][i];
					volt_ADE[o][2][i] += m_Op_Ext_Lor->v_ext_ADE[o][2][i] * eng_sse->Engine_sse::GetVolt(2,pos[0][i],pos[1][i],pos[2][i]);
				}
				break;
			}
			default:
				for (unsigned int i=0; i<m_Op_Ext_Lor->m_LM_Count.at(o); ++i)
				{
					volt_ADE[o][0][i] *= m_Op_Ext_Lor->v_int_ADE[o][0][i];
					volt_ADE[o][0][i] += m_Op_Ext_Lor->v_ext_ADE[o][0][i] * m_Eng->GetVolt(0,pos[0][i],pos[1][i],pos[2][i]);

					volt_ADE[o][1][i] *= m_Op_Ext_Lor->v_int_ADE[o][1][i];
					volt_ADE[o][1][i] += m_Op_Ext_Lor->v_ext_ADE[o][1][i] * m_Eng->GetVolt(1,pos[0][i],pos[1][i],pos[2][i]);

					volt_ADE[o][2][i] *= m_Op_Ext_Lor->v_int_ADE[o][2][i];
					volt_ADE[o][2][i] += m_Op_Ext_Lor->v_ext_ADE[o][2][i] * m_Eng->GetVolt(2,pos[0][i],pos[1][i],pos[2][i]);
				}
				break;
			}
		}
	}
}

void Engine_Ext_LorentzMaterial::DoPreCurrentUpdates()
{
	for (int o=0;o<m_Order;++o)
	{
		if (m_Op_Ext_Lor->m_curr_ADE_On[o]==false) continue;

		unsigned int **pos = m_Op_Ext_Lor->m_LM_pos[o];

		if (m_Op_Ext_Lor->m_curr_Lor_ADE_On[o])
		{
			//switch for different engine types to access faster inline engine functions
			switch (m_Eng->GetType())
			{
			case Engine::BASIC:
			{
				for (unsigned int i=0; i<m_Op_Ext_Lor->m_LM_Count.at(o); ++i)
				{
					curr_Lor_ADE[o][0][i]+=m_Op_Ext_Lor->i_Lor_ADE[o][0][i]*curr_ADE[o][0][i];
					curr_ADE[o][0][i] *= m_Op_Ext_Lor->i_int_ADE[o][0][i];
					curr_ADE[o][0][i] += m_Op_Ext_Lor->i_ext_ADE[o][0][i] * (m_Eng->Engine::GetCurr(0,pos[0][i],pos[1][i],pos[2][i])-curr_Lor_ADE[o][0][i]);

					curr_Lor_ADE[o][1][i]+=m_Op_Ext_Lor->i_Lor_ADE[o][1][i]*curr_ADE[o][1][i];
					curr_ADE[o][1][i] *= m_Op_Ext_Lor->i_int_ADE[o][1][i];
					curr_ADE[o][1][i] += m_Op_Ext_Lor->i_ext_ADE[o][1][i] * (m_Eng->Engine::GetCurr(1,pos[0][i],pos[1][i],pos[2][i])-curr_Lor_ADE[o][1][i]);

					curr_Lor_ADE[o][2][i]+=m_Op_Ext_Lor->i_Lor_ADE[o][2][i]*curr_ADE[o][2][i];
					curr_ADE[o][2][i] *= m_Op_Ext_Lor->i_int_ADE[o][2][i];
					curr_ADE[o][2][i] += m_Op_Ext_Lor->i_ext_ADE[o][2][i] * (m_Eng->Engine::GetCurr(2,pos[0][i],pos[1][i],pos[2][i])-curr_Lor_ADE[o][2][i]);
				}
				break;
			}
			case Engine::SSE:
			{
				Engine_sse* eng_sse = (Engine_sse*)m_Eng;
				for (unsigned int i=0; i<m_Op_Ext_Lor->m_LM_Count.at(o); ++i)
				{
					curr_Lor_ADE[o][0][i]+=m_Op_Ext_Lor->i_Lor_ADE[o][0][i]*curr_ADE[o][0][i];
					curr_ADE[o][0][i] *= m_Op_Ext_Lor->i_int_ADE[o][0][i];
					curr_ADE[o][0][i] += m_Op_Ext_Lor->i_ext_ADE[o][0][i] * (eng_sse->Engine_sse::GetCurr(0,pos[0][i],pos[1][i],pos[2][i])-curr_Lor_ADE[o][0][i]);

					curr_Lor_ADE[o][1][i]+=m_Op_Ext_Lor->i_Lor_ADE[o][1][i]*curr_ADE[o][1][i];
					curr_ADE[o][1][i] *= m_Op_Ext_Lor->i_int_ADE[o][1][i];
					curr_ADE[o][1][i] += m_Op_Ext_Lor->i_ext_ADE[o][1][i] * (eng_sse->Engine_sse::GetCurr(1,pos[0][i],pos[1][i],pos[2][i])-curr_Lor_ADE[o][1][i]);

					curr_Lor_ADE[o][2][i]+=m_Op_Ext_Lor->i_Lor_ADE[o][2][i]*curr_ADE[o][2][i];
					curr_ADE[o][2][i] *= m_Op_Ext_Lor->i_int_ADE[o][2][i];
					curr_ADE[o][2][i] += m_Op_Ext_Lor->i_ext_ADE[o][2][i] * (eng_sse->Engine_sse::GetCurr(2,pos[0][i],pos[1][i],pos[2][i])-curr_Lor_ADE[o][2][i]);
				}
				break;
			}
			default:
				for (unsigned int i=0; i<m_Op_Ext_Lor->m_LM_Count.at(o); ++i)
				{
					curr_Lor_ADE[o][0][i]+=m_Op_Ext_Lor->i_Lor_ADE[o][0][i]*curr_ADE[o][0][i];
					curr_ADE[o][0][i] *= m_Op_Ext_Lor->i_int_ADE[o][0][i];
					curr_ADE[o][0][i] += m_Op_Ext_Lor->i_ext_ADE[o][0][i] * (m_Eng->GetCurr(0,pos[0][i],pos[1][i],pos[2][i])-curr_Lor_ADE[o][0][i]);

					curr_Lor_ADE[o][1][i]+=m_Op_Ext_Lor->i_Lor_ADE[o][1][i]*curr_ADE[o][1][i];
					curr_ADE[o][1][i] *= m_Op_Ext_Lor->i_int_ADE[o][1][i];
					curr_ADE[o][1][i] += m_Op_Ext_Lor->i_ext_ADE[o][1][i] * (m_Eng->GetCurr(1,pos[0][i],pos[1][i],pos[2][i])-curr_Lor_ADE[o][1][i]);

					curr_Lor_ADE[o][2][i]+=m_Op_Ext_Lor->i_Lor_ADE[o][2][i]*curr_ADE[o][2][i];
					curr_ADE[o][2][i] *= m_Op_Ext_Lor->i_int_ADE[o][2][i];
					curr_ADE[o][2][i] += m_Op_Ext_Lor->i_ext_ADE[o][2][i] * (m_Eng->GetCurr(2,pos[0][i],pos[1][i],pos[2][i])-curr_Lor_ADE[o][2][i]);
				}
				break;
			}
		}
		else
		{
			//switch for different engine types to access faster inline engine functions
			switch (m_Eng->GetType())
			{
			case Engine::BASIC:
			{
				for (unsigned int i=0; i<m_Op_Ext_Lor->m_LM_Count.at(o); ++i)
				{
					curr_ADE[o][0][i] *= m_Op_Ext_Lor->i_int_ADE[o][0][i];
					curr_ADE[o][0][i] += m_Op_Ext_Lor->i_ext_ADE[o][0][i] * m_Eng->Engine::GetCurr(0,pos[0][i],pos[1][i],pos[2][i]);

					curr_ADE[o][1][i] *= m_Op_Ext_Lor->i_int_ADE[o][1][i];
					curr_ADE[o][1][i] += m_Op_Ext_Lor->i_ext_ADE[o][1][i] * m_Eng->Engine::GetCurr(1,pos[0][i],pos[1][i],pos[2][i]);

					curr_ADE[o][2][i] *= m_Op_Ext_Lor->i_int_ADE[o][2][i];
					curr_ADE[o][2][i] += m_Op_Ext_Lor->i_ext_ADE[o][2][i] * m_Eng->Engine::GetCurr(2,pos[0][i],pos[1][i],pos[2][i]);
				}
				break;
			}
			case Engine::SSE:
			{
				Engine_sse* eng_sse = (Engine_sse*)m_Eng;
				for (unsigned int i=0; i<m_Op_Ext_Lor->m_LM_Count.at(o); ++i)
				{
					curr_ADE[o][0][i] *= m_Op_Ext_Lor->i_int_ADE[o][0][i];
					curr_ADE[o][0][i] += m_Op_Ext_Lor->i_ext_ADE[o][0][i] * eng_sse->Engine_sse::GetCurr(0,pos[0][i],pos[1][i],pos[2][i]);

					curr_ADE[o][1][i] *= m_Op_Ext_Lor->i_int_ADE[o][1][i];
					curr_ADE[o][1][i] += m_Op_Ext_Lor->i_ext_ADE[o][1][i] * eng_sse->Engine_sse::GetCurr(1,pos[0][i],pos[1][i],pos[2][i]);

					curr_ADE[o][2][i] *= m_Op_Ext_Lor->i_int_ADE[o][2][i];
					curr_ADE[o][2][i] += m_Op_Ext_Lor->i_ext_ADE[o][2][i] * eng_sse->Engine_sse::GetCurr(2,pos[0][i],pos[1][i],pos[2][i]);
				}
				break;
			}
			default:
				for (unsigned int i=0; i<m_Op_Ext_Lor->m_LM_Count.at(o); ++i)
				{
					curr_ADE[o][0][i] *= m_Op_Ext_Lor->i_int_ADE[o][0][i];
					curr_ADE[o][0][i] += m_Op_Ext_Lor->i_ext_ADE[o][0][i] * m_Eng->GetCurr(0,pos[0][i],pos[1][i],pos[2][i]);

					curr_ADE[o][1][i] *= m_Op_Ext_Lor->i_int_ADE[o][1][i];
					curr_ADE[o][1][i] += m_Op_Ext_Lor->i_ext_ADE[o][1][i] * m_Eng->GetCurr(1,pos[0][i],pos[1][i],pos[2][i]);

					curr_ADE[o][2][i] *= m_Op_Ext_Lor->i_int_ADE[o][2][i];
					curr_ADE[o][2][i] += m_Op_Ext_Lor->i_ext_ADE[o][2][i] * m_Eng->GetCurr(2,pos[0][i],pos[1][i],pos[2][i]);
				}
				break;
			}
		}
	}
}

