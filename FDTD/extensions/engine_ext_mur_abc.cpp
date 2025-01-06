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

#include "engine_ext_mur_abc.h"
#include "operator_ext_mur_abc.h"
#include "FDTD/engine.h"
#include "FDTD/engine_sse.h"
#include "tools/array_ops.h"
#include "tools/useful.h"
#include "operator_ext_excitation.h"

Engine_Ext_Mur_ABC::Engine_Ext_Mur_ABC(Operator_Ext_Mur_ABC* op_ext) :
	Engine_Extension(op_ext),
	m_Mur_Coeff_nyP (op_ext->m_Mur_Coeff_nyP),
	m_Mur_Coeff_nyPP(op_ext->m_Mur_Coeff_nyPP),
	m_volt_nyP ("volt_nyP",  op_ext->m_numLines),
	m_volt_nyPP("volt_nyPP", op_ext->m_numLines)
{
	m_Op_mur = op_ext;
	m_numLines[0] = m_Op_mur->m_numLines[0];
	m_numLines[1] = m_Op_mur->m_numLines[1];
	m_ny = m_Op_mur->m_ny;
	m_nyP = m_Op_mur->m_nyP;
	m_nyPP = m_Op_mur->m_nyPP;
	m_LineNr = m_Op_mur->m_LineNr;
	m_LineNr_Shift = m_Op_mur->m_LineNr_Shift;

	//find if some excitation is on this mur-abc and find the max length of this excite, so that the abc can start after the excitation is done...
	int maxDelay=-1;
	Operator_Ext_Excitation* Exc_ext = m_Op_mur->m_Op->GetExcitationExtension();
	for (unsigned int n=0; n<Exc_ext->GetVoltCount(); ++n)
	{
		if ( ((Exc_ext->Volt_dir[n]==m_nyP) || (Exc_ext->Volt_dir[n]==m_nyPP)) && (Exc_ext->Volt_index[m_ny][n]==m_LineNr) )
		{
			if ((int)Exc_ext->Volt_delay[n]>maxDelay)
				maxDelay = (int)Exc_ext->Volt_delay[n];
		}
	}
	m_start_TS = 0;
	if (maxDelay>=0)
	{
		m_start_TS = maxDelay + m_Op_mur->m_Op->GetExcitationSignal()->GetLength() + 10; //give it some extra timesteps, for the excitation to travel at least one cell away
		cerr << "Engine_Ext_Mur_ABC::Engine_Ext_Mur_ABC: Warning: Excitation inside the Mur-ABC #" <<  m_ny << "-" << (int)(m_LineNr>0) << " found!!!!  Mur-ABC will be switched on after excitation is done at " << m_start_TS << " timesteps!!! " << endl;
	}

	SetNumberOfThreads(1);
}

Engine_Ext_Mur_ABC::~Engine_Ext_Mur_ABC()
{
}


void Engine_Ext_Mur_ABC::SetNumberOfThreads(int nrThread)
{
	Engine_Extension::SetNumberOfThreads(nrThread);

	m_numX = AssignJobs2Threads(m_numLines[0],m_NrThreads,false);
	m_start.resize(m_NrThreads,0);
	m_start.at(0)=0;
	for (size_t n=1; n<m_numX.size(); ++n)
		m_start.at(n) = m_start.at(n-1) + m_numX.at(n-1);
}


template <typename EngType>
void Engine_Ext_Mur_ABC::DoPreVoltageUpdatesImpl(EngType* eng, int threadID)
{
	if (IsActive()==false) return;
	if (m_Eng==NULL) return;
	if (threadID>=m_NrThreads)
		return;

	// See detailed comments in operator_ext_mur_abc.h, not repeated here.
	unsigned int pos[] = {0,0,0};
	unsigned int pos_shift[] = {0,0,0};
	pos[m_ny] = m_LineNr;
	pos_shift[m_ny] = m_LineNr_Shift;

	for (unsigned int i = m_start.at(threadID); i < m_start.at(threadID) + m_numX.at(threadID); i++)
	{
		pos[m_nyP] = i;
		pos_shift[m_nyP] = i;

		for (unsigned int j = 0; j < m_numLines[1]; j++)
		{
			pos[m_nyPP] = j;
			pos_shift[m_nyPP] = j;

			m_volt_nyP[i][j] = eng->EngType::GetVolt(m_nyP, pos_shift) -
				m_Op_mur->m_Mur_Coeff_nyP[i][j] * eng->EngType::GetVolt(m_nyP, pos);

			m_volt_nyPP[i][j] = eng->EngType::GetVolt(m_nyPP, pos_shift) -
				m_Op_mur->m_Mur_Coeff_nyPP[i][j] * eng->EngType::GetVolt(m_nyPP, pos);
		}
	}
}

void Engine_Ext_Mur_ABC::DoPreVoltageUpdates(int threadID)
{
	ENG_DISPATCH_ARGS(DoPreVoltageUpdatesImpl, threadID);
}

template <typename EngType>
void Engine_Ext_Mur_ABC::DoPostVoltageUpdatesImpl(EngType* eng, int threadID)
{
	if (IsActive()==false) return;
	if (m_Eng==NULL) return;
	if (threadID>=m_NrThreads)
		return;
	unsigned int pos_shift[] = {0,0,0};
	pos_shift[m_ny] = m_LineNr_Shift;

	for (unsigned int i = m_start.at(threadID); i < m_start.at(threadID) + m_numX.at(threadID); i++)
	{
		pos_shift[m_nyP] = i;

		for (unsigned int j = 0; j < m_numLines[1]; j++)
		{
			pos_shift[m_nyPP] = j;

			m_volt_nyP[i][j] +=
				m_Op_mur->m_Mur_Coeff_nyP[i][j] * eng->EngType::GetVolt(m_nyP, pos_shift);
			m_volt_nyPP[i][j] +=
				m_Op_mur->m_Mur_Coeff_nyPP[i][j] * eng->EngType::GetVolt(m_nyPP, pos_shift);
		}
	}
}

void Engine_Ext_Mur_ABC::DoPostVoltageUpdates(int threadID)
{
	ENG_DISPATCH_ARGS(DoPostVoltageUpdatesImpl, threadID);
}

template <typename EngType>
void Engine_Ext_Mur_ABC::Apply2VoltagesImpl(EngType* eng, int threadID)
{
	if (IsActive()==false) return;
	if (threadID>=m_NrThreads)
		return;
	if (m_Eng==NULL) return;
	unsigned int pos[] = {0,0,0};
	pos[m_ny] = m_LineNr;

	for (unsigned int i = m_start.at(threadID); i < m_start.at(threadID) + m_numX.at(threadID); i++)
	{
		pos[m_nyP] = i;

		for (unsigned int j = 0; j < m_numLines[1]; j++)
		{
			pos[m_nyPP] = j;

			eng->EngType::SetVolt(m_nyP, pos, m_volt_nyP[i][j]);
			eng->EngType::SetVolt(m_nyPP, pos, m_volt_nyPP[i][j]);
		}
	}
}

void Engine_Ext_Mur_ABC::Apply2Voltages(int threadID)
{
	ENG_DISPATCH_ARGS(Apply2VoltagesImpl, threadID);
}
