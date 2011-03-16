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

#include "engine_ext_upml.h"
#include "operator_ext_upml.h"
#include "FDTD/engine.h"
#include "FDTD/engine_sse.h"
#include "tools/array_ops.h"
#include "tools/useful.h"

Engine_Ext_UPML::Engine_Ext_UPML(Operator_Ext_UPML* op_ext) : Engine_Extension(op_ext)
{
	m_Op_UPML = op_ext;

	//this ABC extension should be executed first!
	m_Priority = ENG_EXT_PRIO_UPML;

	volt_flux = Create_N_3DArray<FDTD_FLOAT>(m_Op_UPML->m_numLines);
	curr_flux = Create_N_3DArray<FDTD_FLOAT>(m_Op_UPML->m_numLines);

	SetNumberOfThreads(1);
}

Engine_Ext_UPML::~Engine_Ext_UPML()
{
	Delete_N_3DArray<FDTD_FLOAT>(volt_flux,m_Op_UPML->m_numLines);
	volt_flux=NULL;
	Delete_N_3DArray<FDTD_FLOAT>(curr_flux,m_Op_UPML->m_numLines);
	curr_flux=NULL;
}

void Engine_Ext_UPML::SetNumberOfThreads(int nrThread)
{
	Engine_Extension::SetNumberOfThreads(nrThread);

	m_numX = AssignJobs2Threads(m_Op_UPML->m_numLines[0],m_NrThreads,false);
	m_start.resize(m_NrThreads,0);
	m_start.at(0)=0;
	for (size_t n=1; n<m_numX.size(); ++n)
		m_start.at(n) = m_start.at(n-1) + m_numX.at(n-1);
}


void Engine_Ext_UPML::DoPreVoltageUpdates(int threadID)
{
	if (m_Eng==NULL)
		return;

	if (threadID>=m_NrThreads)
		return;

	unsigned int pos[3];
	unsigned int loc_pos[3];
	FDTD_FLOAT f_help;
	switch (m_Eng->GetType())
	{
	case Engine::BASIC:
		{
			for (unsigned int lineX=0; lineX<m_numX.at(threadID); ++lineX)
			{
				loc_pos[0]=lineX+m_start.at(threadID);
				pos[0] = loc_pos[0] + m_Op_UPML->m_StartPos[0];
				for (loc_pos[1]=0; loc_pos[1]<m_Op_UPML->m_numLines[1]; ++loc_pos[1])
				{
					pos[1] = loc_pos[1] + m_Op_UPML->m_StartPos[1];
					for (loc_pos[2]=0; loc_pos[2]<m_Op_UPML->m_numLines[2]; ++loc_pos[2])
					{
						pos[2] = loc_pos[2] + m_Op_UPML->m_StartPos[2];

						f_help = m_Op_UPML->vv[0][loc_pos[0]][loc_pos[1]][loc_pos[2]]   * m_Eng->Engine::GetVolt(0,pos)
						         - m_Op_UPML->vvfo[0][loc_pos[0]][loc_pos[1]][loc_pos[2]] * volt_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]];
						m_Eng->Engine::SetVolt(0,pos, volt_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]]);
						volt_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]] = f_help;

						f_help = m_Op_UPML->vv[1][loc_pos[0]][loc_pos[1]][loc_pos[2]]   * m_Eng->Engine::GetVolt(1,pos)
						         - m_Op_UPML->vvfo[1][loc_pos[0]][loc_pos[1]][loc_pos[2]] * volt_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]];
						m_Eng->Engine::SetVolt(1,pos, volt_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]]);
						volt_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]] = f_help;

						f_help = m_Op_UPML->vv[2][loc_pos[0]][loc_pos[1]][loc_pos[2]]   * m_Eng->Engine::GetVolt(2,pos)
						         - m_Op_UPML->vvfo[2][loc_pos[0]][loc_pos[1]][loc_pos[2]] * volt_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]];
						m_Eng->Engine::SetVolt(2,pos, volt_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]]);
						volt_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]] = f_help;
					}
				}
			}
			break;
		}
	case Engine::SSE:
		{
			Engine_sse* eng_sse = (Engine_sse*) m_Eng;
			for (unsigned int lineX=0; lineX<m_numX.at(threadID); ++lineX)
			{
				loc_pos[0]=lineX+m_start.at(threadID);
				pos[0] = loc_pos[0] + m_Op_UPML->m_StartPos[0];
				for (loc_pos[1]=0; loc_pos[1]<m_Op_UPML->m_numLines[1]; ++loc_pos[1])
				{
					pos[1] = loc_pos[1] + m_Op_UPML->m_StartPos[1];
					for (loc_pos[2]=0; loc_pos[2]<m_Op_UPML->m_numLines[2]; ++loc_pos[2])
					{
						pos[2] = loc_pos[2] + m_Op_UPML->m_StartPos[2];

						f_help = m_Op_UPML->vv[0][loc_pos[0]][loc_pos[1]][loc_pos[2]]   * eng_sse->Engine_sse::GetVolt(0,pos)
						         - m_Op_UPML->vvfo[0][loc_pos[0]][loc_pos[1]][loc_pos[2]] * volt_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]];
						eng_sse->Engine_sse::SetVolt(0,pos, volt_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]]);
						volt_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]] = f_help;

						f_help = m_Op_UPML->vv[1][loc_pos[0]][loc_pos[1]][loc_pos[2]]   * eng_sse->Engine_sse::GetVolt(1,pos)
						         - m_Op_UPML->vvfo[1][loc_pos[0]][loc_pos[1]][loc_pos[2]] * volt_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]];
						eng_sse->Engine_sse::SetVolt(1,pos, volt_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]]);
						volt_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]] = f_help;

						f_help = m_Op_UPML->vv[2][loc_pos[0]][loc_pos[1]][loc_pos[2]]   * eng_sse->Engine_sse::GetVolt(2,pos)
						         - m_Op_UPML->vvfo[2][loc_pos[0]][loc_pos[1]][loc_pos[2]] * volt_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]];
						eng_sse->Engine_sse::SetVolt(2,pos, volt_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]]);
						volt_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]] = f_help;
					}
				}
			}
			break;
		}
	default:
		{
			for (unsigned int lineX=0; lineX<m_numX.at(threadID); ++lineX)
			{
				loc_pos[0]=lineX+m_start.at(threadID);
				pos[0] = loc_pos[0] + m_Op_UPML->m_StartPos[0];
				for (loc_pos[1]=0; loc_pos[1]<m_Op_UPML->m_numLines[1]; ++loc_pos[1])
				{
					pos[1] = loc_pos[1] + m_Op_UPML->m_StartPos[1];
					for (loc_pos[2]=0; loc_pos[2]<m_Op_UPML->m_numLines[2]; ++loc_pos[2])
					{
						pos[2] = loc_pos[2] + m_Op_UPML->m_StartPos[2];

						f_help = m_Op_UPML->vv[0][loc_pos[0]][loc_pos[1]][loc_pos[2]]   * m_Eng->GetVolt(0,pos)
						         - m_Op_UPML->vvfo[0][loc_pos[0]][loc_pos[1]][loc_pos[2]] * volt_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]];
						m_Eng->SetVolt(0,pos, volt_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]]);
						volt_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]] = f_help;

						f_help = m_Op_UPML->vv[1][loc_pos[0]][loc_pos[1]][loc_pos[2]]   * m_Eng->GetVolt(1,pos)
						         - m_Op_UPML->vvfo[1][loc_pos[0]][loc_pos[1]][loc_pos[2]] * volt_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]];
						m_Eng->SetVolt(1,pos, volt_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]]);
						volt_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]] = f_help;

						f_help = m_Op_UPML->vv[2][loc_pos[0]][loc_pos[1]][loc_pos[2]]   * m_Eng->GetVolt(2,pos)
						         - m_Op_UPML->vvfo[2][loc_pos[0]][loc_pos[1]][loc_pos[2]] * volt_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]];
						m_Eng->SetVolt(2,pos, volt_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]]);
						volt_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]] = f_help;
					}
				}
			}
			break;
		}
	}

}

void Engine_Ext_UPML::DoPostVoltageUpdates(int threadID)
{
	if (m_Eng==NULL)
		return;
	if (threadID>=m_NrThreads)
		return;

	unsigned int pos[3];
	unsigned int loc_pos[3];
	FDTD_FLOAT f_help;

	switch (m_Eng->GetType())
	{
	case Engine::BASIC:
		{
			for (unsigned int lineX=0; lineX<m_numX.at(threadID); ++lineX)
			{
				loc_pos[0]=lineX+m_start.at(threadID);
				pos[0] = loc_pos[0] + m_Op_UPML->m_StartPos[0];
				for (loc_pos[1]=0; loc_pos[1]<m_Op_UPML->m_numLines[1]; ++loc_pos[1])
				{
					pos[1] = loc_pos[1] + m_Op_UPML->m_StartPos[1];
					for (loc_pos[2]=0; loc_pos[2]<m_Op_UPML->m_numLines[2]; ++loc_pos[2])
					{
						pos[2] = loc_pos[2] + m_Op_UPML->m_StartPos[2];

						f_help = volt_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]];
						volt_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]] = m_Eng->Engine::GetVolt(0,pos);
						m_Eng->Engine::SetVolt(0,pos, f_help + m_Op_UPML->vvfn[0][loc_pos[0]][loc_pos[1]][loc_pos[2]] * volt_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]]);

						f_help = volt_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]];
						volt_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]] = m_Eng->Engine::GetVolt(1,pos);
						m_Eng->Engine::SetVolt(1,pos, f_help + m_Op_UPML->vvfn[1][loc_pos[0]][loc_pos[1]][loc_pos[2]] * volt_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]]);

						f_help = volt_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]];
						volt_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]] = m_Eng->Engine::GetVolt(2,pos);
						m_Eng->Engine::SetVolt(2,pos, f_help + m_Op_UPML->vvfn[2][loc_pos[0]][loc_pos[1]][loc_pos[2]] * volt_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]]);
					}
				}
			}
			break;
		}
	case Engine::SSE:
		{
			Engine_sse* eng_sse = (Engine_sse*) m_Eng;
			for (unsigned int lineX=0; lineX<m_numX.at(threadID); ++lineX)
			{
				loc_pos[0]=lineX+m_start.at(threadID);
				pos[0] = loc_pos[0] + m_Op_UPML->m_StartPos[0];
				for (loc_pos[1]=0; loc_pos[1]<m_Op_UPML->m_numLines[1]; ++loc_pos[1])
				{
					pos[1] = loc_pos[1] + m_Op_UPML->m_StartPos[1];
					for (loc_pos[2]=0; loc_pos[2]<m_Op_UPML->m_numLines[2]; ++loc_pos[2])
					{
						pos[2] = loc_pos[2] + m_Op_UPML->m_StartPos[2];

						f_help = volt_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]];
						volt_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]] = eng_sse->Engine_sse::GetVolt(0,pos);
						eng_sse->Engine_sse::SetVolt(0,pos, f_help + m_Op_UPML->vvfn[0][loc_pos[0]][loc_pos[1]][loc_pos[2]] * volt_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]]);

						f_help = volt_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]];
						volt_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]] = eng_sse->Engine_sse::GetVolt(1,pos);
						eng_sse->Engine_sse::SetVolt(1,pos, f_help + m_Op_UPML->vvfn[1][loc_pos[0]][loc_pos[1]][loc_pos[2]] * volt_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]]);

						f_help = volt_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]];
						volt_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]] = eng_sse->Engine_sse::GetVolt(2,pos);
						eng_sse->Engine_sse::SetVolt(2,pos, f_help + m_Op_UPML->vvfn[2][loc_pos[0]][loc_pos[1]][loc_pos[2]] * volt_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]]);
					}
				}
			}
			break;
		}
	default:
		{
			for (unsigned int lineX=0; lineX<m_numX.at(threadID); ++lineX)
			{
				loc_pos[0]=lineX+m_start.at(threadID);
				pos[0] = loc_pos[0] + m_Op_UPML->m_StartPos[0];
				for (loc_pos[1]=0; loc_pos[1]<m_Op_UPML->m_numLines[1]; ++loc_pos[1])
				{
					pos[1] = loc_pos[1] + m_Op_UPML->m_StartPos[1];
					for (loc_pos[2]=0; loc_pos[2]<m_Op_UPML->m_numLines[2]; ++loc_pos[2])
					{
						pos[2] = loc_pos[2] + m_Op_UPML->m_StartPos[2];

						f_help = volt_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]];
						volt_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]] = m_Eng->GetVolt(0,pos);
						m_Eng->SetVolt(0,pos, f_help + m_Op_UPML->vvfn[0][loc_pos[0]][loc_pos[1]][loc_pos[2]] * volt_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]]);

						f_help = volt_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]];
						volt_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]] = m_Eng->GetVolt(1,pos);
						m_Eng->SetVolt(1,pos, f_help + m_Op_UPML->vvfn[1][loc_pos[0]][loc_pos[1]][loc_pos[2]] * volt_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]]);

						f_help = volt_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]];
						volt_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]] = m_Eng->GetVolt(2,pos);
						m_Eng->SetVolt(2,pos, f_help + m_Op_UPML->vvfn[2][loc_pos[0]][loc_pos[1]][loc_pos[2]] * volt_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]]);
					}
				}
			}
			break;
		}
	}

}

void Engine_Ext_UPML::DoPreCurrentUpdates(int threadID)
{
	if (m_Eng==NULL)
		return;
	if (threadID>=m_NrThreads)
		return;

	unsigned int pos[3];
	unsigned int loc_pos[3];
	FDTD_FLOAT f_help;

	switch (m_Eng->GetType())
	{
	case Engine::BASIC:
		{
			for (unsigned int lineX=0; lineX<m_numX.at(threadID); ++lineX)
			{
				loc_pos[0]=lineX+m_start.at(threadID);
				pos[0] = loc_pos[0] + m_Op_UPML->m_StartPos[0];
				for (loc_pos[1]=0; loc_pos[1]<m_Op_UPML->m_numLines[1]; ++loc_pos[1])
				{
					pos[1] = loc_pos[1] + m_Op_UPML->m_StartPos[1];
					for (loc_pos[2]=0; loc_pos[2]<m_Op_UPML->m_numLines[2]; ++loc_pos[2])
					{
						pos[2] = loc_pos[2] + m_Op_UPML->m_StartPos[2];

						f_help = m_Op_UPML->ii[0][loc_pos[0]][loc_pos[1]][loc_pos[2]]   * m_Eng->Engine::GetCurr(0,pos)
						         - m_Op_UPML->iifo[0][loc_pos[0]][loc_pos[1]][loc_pos[2]] * curr_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]];
						m_Eng->Engine::SetCurr(0,pos, curr_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]]);
						curr_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]] = f_help;

						f_help = m_Op_UPML->ii[1][loc_pos[0]][loc_pos[1]][loc_pos[2]]   * m_Eng->Engine::GetCurr(1,pos)
						         - m_Op_UPML->iifo[1][loc_pos[0]][loc_pos[1]][loc_pos[2]] * curr_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]];
						m_Eng->Engine::SetCurr(1,pos, curr_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]]);
						curr_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]] = f_help;

						f_help = m_Op_UPML->ii[2][loc_pos[0]][loc_pos[1]][loc_pos[2]]   * m_Eng->Engine::GetCurr(2,pos)
						         - m_Op_UPML->iifo[2][loc_pos[0]][loc_pos[1]][loc_pos[2]] * curr_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]];
						m_Eng->Engine::SetCurr(2,pos, curr_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]]);
						curr_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]] = f_help;
					}
				}
			}
			break;
		}
	case Engine::SSE:
		{
			Engine_sse* eng_sse = (Engine_sse*) m_Eng;
			for (unsigned int lineX=0; lineX<m_numX.at(threadID); ++lineX)
			{
				loc_pos[0]=lineX+m_start.at(threadID);
				pos[0] = loc_pos[0] + m_Op_UPML->m_StartPos[0];
				for (loc_pos[1]=0; loc_pos[1]<m_Op_UPML->m_numLines[1]; ++loc_pos[1])
				{
					pos[1] = loc_pos[1] + m_Op_UPML->m_StartPos[1];
					for (loc_pos[2]=0; loc_pos[2]<m_Op_UPML->m_numLines[2]; ++loc_pos[2])
					{
						pos[2] = loc_pos[2] + m_Op_UPML->m_StartPos[2];

						f_help = m_Op_UPML->ii[0][loc_pos[0]][loc_pos[1]][loc_pos[2]]   * eng_sse->Engine_sse::GetCurr(0,pos)
						         - m_Op_UPML->iifo[0][loc_pos[0]][loc_pos[1]][loc_pos[2]] * curr_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]];
						eng_sse->Engine_sse::SetCurr(0,pos, curr_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]]);
						curr_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]] = f_help;

						f_help = m_Op_UPML->ii[1][loc_pos[0]][loc_pos[1]][loc_pos[2]]   * eng_sse->Engine_sse::GetCurr(1,pos)
						         - m_Op_UPML->iifo[1][loc_pos[0]][loc_pos[1]][loc_pos[2]] * curr_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]];
						eng_sse->Engine_sse::SetCurr(1,pos, curr_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]]);
						curr_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]] = f_help;

						f_help = m_Op_UPML->ii[2][loc_pos[0]][loc_pos[1]][loc_pos[2]]   * eng_sse->Engine_sse::GetCurr(2,pos)
						         - m_Op_UPML->iifo[2][loc_pos[0]][loc_pos[1]][loc_pos[2]] * curr_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]];
						eng_sse->Engine_sse::SetCurr(2,pos, curr_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]]);
						curr_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]] = f_help;

					}
				}
			}
			break;
		}
	default:
		{
			for (unsigned int lineX=0; lineX<m_numX.at(threadID); ++lineX)
			{
				loc_pos[0]=lineX+m_start.at(threadID);
				pos[0] = loc_pos[0] + m_Op_UPML->m_StartPos[0];
				for (loc_pos[1]=0; loc_pos[1]<m_Op_UPML->m_numLines[1]; ++loc_pos[1])
				{
					pos[1] = loc_pos[1] + m_Op_UPML->m_StartPos[1];
					for (loc_pos[2]=0; loc_pos[2]<m_Op_UPML->m_numLines[2]; ++loc_pos[2])
					{
						pos[2] = loc_pos[2] + m_Op_UPML->m_StartPos[2];

						f_help = m_Op_UPML->ii[0][loc_pos[0]][loc_pos[1]][loc_pos[2]]   * m_Eng->GetCurr(0,pos)
						         - m_Op_UPML->iifo[0][loc_pos[0]][loc_pos[1]][loc_pos[2]] * curr_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]];
						m_Eng->SetCurr(0,pos, curr_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]]);
						curr_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]] = f_help;

						f_help = m_Op_UPML->ii[1][loc_pos[0]][loc_pos[1]][loc_pos[2]]   * m_Eng->GetCurr(1,pos)
						         - m_Op_UPML->iifo[1][loc_pos[0]][loc_pos[1]][loc_pos[2]] * curr_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]];
						m_Eng->SetCurr(1,pos, curr_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]]);
						curr_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]] = f_help;

						f_help = m_Op_UPML->ii[2][loc_pos[0]][loc_pos[1]][loc_pos[2]]   * m_Eng->GetCurr(2,pos)
						         - m_Op_UPML->iifo[2][loc_pos[0]][loc_pos[1]][loc_pos[2]] * curr_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]];
						m_Eng->SetCurr(2,pos, curr_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]]);
						curr_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]] = f_help;
					}
				}
			}
			break;
		}
	}
}

void Engine_Ext_UPML::DoPostCurrentUpdates(int threadID)
{
	if (m_Eng==NULL)
		return;
	if (threadID>=m_NrThreads)
		return;

	unsigned int pos[3];
	unsigned int loc_pos[3];
	FDTD_FLOAT f_help;

	switch (m_Eng->GetType())
	{
	case Engine::BASIC:
		{
			for (unsigned int lineX=0; lineX<m_numX.at(threadID); ++lineX)
			{
				loc_pos[0]=lineX+m_start.at(threadID);
				pos[0] = loc_pos[0] + m_Op_UPML->m_StartPos[0];
				for (loc_pos[1]=0; loc_pos[1]<m_Op_UPML->m_numLines[1]; ++loc_pos[1])
				{
					pos[1] = loc_pos[1] + m_Op_UPML->m_StartPos[1];
					for (loc_pos[2]=0; loc_pos[2]<m_Op_UPML->m_numLines[2]; ++loc_pos[2])
					{
						pos[2] = loc_pos[2] + m_Op_UPML->m_StartPos[2];

						f_help = curr_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]];
						curr_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]] = m_Eng->Engine::GetCurr(0,pos);
						m_Eng->Engine::SetCurr(0,pos, f_help + m_Op_UPML->iifn[0][loc_pos[0]][loc_pos[1]][loc_pos[2]] * curr_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]]);

						f_help = curr_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]];
						curr_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]] = m_Eng->Engine::GetCurr(1,pos);
						m_Eng->Engine::SetCurr(1,pos, f_help + m_Op_UPML->iifn[1][loc_pos[0]][loc_pos[1]][loc_pos[2]] * curr_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]]);

						f_help = curr_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]];
						curr_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]] = m_Eng->Engine::GetCurr(2,pos);
						m_Eng->Engine::SetCurr(2,pos, f_help + m_Op_UPML->iifn[2][loc_pos[0]][loc_pos[1]][loc_pos[2]] * curr_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]]);
					}
				}
			}
			break;
		}
	case Engine::SSE:
		{
			Engine_sse* eng_sse = (Engine_sse*) m_Eng;
			for (unsigned int lineX=0; lineX<m_numX.at(threadID); ++lineX)
			{
				loc_pos[0]=lineX+m_start.at(threadID);
				pos[0] = loc_pos[0] + m_Op_UPML->m_StartPos[0];
				for (loc_pos[1]=0; loc_pos[1]<m_Op_UPML->m_numLines[1]; ++loc_pos[1])
				{
					pos[1] = loc_pos[1] + m_Op_UPML->m_StartPos[1];
					for (loc_pos[2]=0; loc_pos[2]<m_Op_UPML->m_numLines[2]; ++loc_pos[2])
					{
						pos[2] = loc_pos[2] + m_Op_UPML->m_StartPos[2];

						f_help = curr_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]];
						curr_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]] = eng_sse->Engine_sse::GetCurr(0,pos);
						eng_sse->Engine_sse::SetCurr(0,pos, f_help + m_Op_UPML->iifn[0][loc_pos[0]][loc_pos[1]][loc_pos[2]] * curr_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]]);

						f_help = curr_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]];
						curr_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]] = eng_sse->Engine_sse::GetCurr(1,pos);
						eng_sse->Engine_sse::SetCurr(1,pos, f_help + m_Op_UPML->iifn[1][loc_pos[0]][loc_pos[1]][loc_pos[2]] * curr_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]]);

						f_help = curr_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]];
						curr_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]] = eng_sse->Engine_sse::GetCurr(2,pos);
						eng_sse->Engine_sse::SetCurr(2,pos, f_help + m_Op_UPML->iifn[2][loc_pos[0]][loc_pos[1]][loc_pos[2]] * curr_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]]);
					}
				}
			}
			break;
		}
	default:
		{
			for (unsigned int lineX=0; lineX<m_numX.at(threadID); ++lineX)
			{
				loc_pos[0]=lineX+m_start.at(threadID);
				pos[0] = loc_pos[0] + m_Op_UPML->m_StartPos[0];
				for (loc_pos[1]=0; loc_pos[1]<m_Op_UPML->m_numLines[1]; ++loc_pos[1])
				{
					pos[1] = loc_pos[1] + m_Op_UPML->m_StartPos[1];
					for (loc_pos[2]=0; loc_pos[2]<m_Op_UPML->m_numLines[2]; ++loc_pos[2])
					{
						pos[2] = loc_pos[2] + m_Op_UPML->m_StartPos[2];

						f_help = curr_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]];
						curr_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]] = m_Eng->GetCurr(0,pos);
						m_Eng->SetCurr(0,pos, f_help + m_Op_UPML->iifn[0][loc_pos[0]][loc_pos[1]][loc_pos[2]] * curr_flux[0][loc_pos[0]][loc_pos[1]][loc_pos[2]]);

						f_help = curr_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]];
						curr_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]] = m_Eng->GetCurr(1,pos);
						m_Eng->SetCurr(1,pos, f_help + m_Op_UPML->iifn[1][loc_pos[0]][loc_pos[1]][loc_pos[2]] * curr_flux[1][loc_pos[0]][loc_pos[1]][loc_pos[2]]);

						f_help = curr_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]];
						curr_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]] = m_Eng->GetCurr(2,pos);
						m_Eng->SetCurr(2,pos, f_help + m_Op_UPML->iifn[2][loc_pos[0]][loc_pos[1]][loc_pos[2]] * curr_flux[2][loc_pos[0]][loc_pos[1]][loc_pos[2]]);
					}
				}
			}
			break;
		}
	}
}
