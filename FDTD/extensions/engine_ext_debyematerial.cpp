/*
*	Copyright (C) 2026 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#include "engine_ext_debyematerial.h"
#include "operator_ext_debyematerial.h"
#include "FDTD/engine_sse.h"

Engine_Ext_DebyeMaterial::Engine_Ext_DebyeMaterial(Operator_Ext_DebyeMaterial* op_ext_debye) : Engine_Ext_Dispersive(op_ext_debye)
{
	m_Op_Ext_Deb = op_ext_debye;
	// the base allocated volt_ADE[0], which carries the joint correction
	unsigned int count = m_Op_Ext_Deb->m_LM_Count.at(0);

	volt_pole_ADE = new FDTD_FLOAT**[m_Op_Ext_Deb->m_PoleCount];
	for (int o=0; o<m_Op_Ext_Deb->m_PoleCount; ++o)
	{
		volt_pole_ADE[o] = new FDTD_FLOAT*[3];
		for (int n=0; n<3; ++n)
		{
			volt_pole_ADE[o][n] = new FDTD_FLOAT[count];
			for (unsigned int i=0; i<count; ++i)
				volt_pole_ADE[o][n][i] = 0.0;
		}
	}

	volt_pre_ADE = new FDTD_FLOAT*[3];
	for (int n=0; n<3; ++n)
	{
		volt_pre_ADE[n] = new FDTD_FLOAT[count];
		for (unsigned int i=0; i<count; ++i)
			volt_pre_ADE[n][i] = 0.0;
	}
}

Engine_Ext_DebyeMaterial::~Engine_Ext_DebyeMaterial()
{
	if (volt_pole_ADE!=NULL)
	{
		for (int o=0; o<m_Op_Ext_Deb->m_PoleCount; ++o)
		{
			for (int n=0; n<3; ++n)
				delete[] volt_pole_ADE[o][n];
			delete[] volt_pole_ADE[o];
		}
		delete[] volt_pole_ADE;
		volt_pole_ADE = NULL;
	}

	if (volt_pre_ADE!=NULL)
		for (int n=0; n<3; ++n)
			delete[] volt_pre_ADE[n];
	delete[] volt_pre_ADE;
	volt_pre_ADE = NULL;
}

template <typename EngType>
void Engine_Ext_DebyeMaterial::DoPreVoltageUpdatesImpl(EngType* eng)
{
	unsigned int **pos = m_Op_Ext_Deb->m_LM_pos[0];
	int poles = m_Op_Ext_Deb->m_PoleCount;

	for (unsigned int i=0; i<m_Op_Ext_Deb->m_LM_Count.at(0); ++i)
	{
		for (int n=0; n<3; ++n)
		{
			FDTD_FLOAT volt = eng->EngType::GetVolt(n,pos[0][i],pos[1][i],pos[2][i]);
			FDTD_FLOAT pre = 0.0;
			for (int o=0; o<poles; ++o)
			{
				// as much of the pole's new state as V^n already fixes
				FDTD_FLOAT W = m_Op_Ext_Deb->v_relax_ADE[o][n][i]*volt_pole_ADE[o][n][i]
				             + m_Op_Ext_Deb->v_drive_ADE[o][n][i]*volt;
				pre += W - volt_pole_ADE[o][n][i];
				volt_pole_ADE[o][n][i] = W;
			}
			volt_pre_ADE[n][i] = pre;
		}
	}
}

void Engine_Ext_DebyeMaterial::DoPreVoltageUpdates()
{
	ENG_DISPATCH(DoPreVoltageUpdatesImpl);
}

template <typename EngType>
void Engine_Ext_DebyeMaterial::DoPostVoltageUpdatesImpl(EngType* eng)
{
	unsigned int **pos = m_Op_Ext_Deb->m_LM_pos[0];
	int poles = m_Op_Ext_Deb->m_PoleCount;

	for (unsigned int i=0; i<m_Op_Ext_Deb->m_LM_Count.at(0); ++i)
	{
		for (int n=0; n<3; ++n)
		{
			FDTD_FLOAT v_raw = eng->EngType::GetVolt(n,pos[0][i],pos[1][i],pos[2][i]);
			// solved for the V^(n+1) the poles themselves depend on
			FDTD_FLOAT volt = (v_raw - volt_pre_ADE[n][i]) * m_Op_Ext_Deb->v_solve_ADE[n][i];
			volt_ADE[0][n][i] = v_raw - volt;   // what Apply2Voltages() subtracts
			// V^(n+1) is now known, so the pole states can be finished
			for (int o=0; o<poles; ++o)
				volt_pole_ADE[o][n][i] += m_Op_Ext_Deb->v_drive_ADE[o][n][i]*volt;
		}
	}
}

void Engine_Ext_DebyeMaterial::DoPostVoltageUpdates()
{
	ENG_DISPATCH(DoPostVoltageUpdatesImpl);
}
