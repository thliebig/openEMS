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
#include "engine.h"
#include "tools/array_ops.h"

Engine_Ext_Mur_ABC::Engine_Ext_Mur_ABC(Operator_Ext_Mur_ABC* op_ext) : Engine_Extension(op_ext)
{
	m_Op_mur = op_ext;
	m_numLines[0] = m_Op_mur->m_numLines[0];
	m_numLines[1] = m_Op_mur->m_numLines[1];
	m_ny = m_Op_mur->m_ny;
	m_nyP = m_Op_mur->m_nyP;
	m_nyPP = m_Op_mur->m_nyPP;
	m_LineNr = m_Op_mur->m_LineNr;
	m_LineNr_Shift = m_Op_mur->m_LineNr_Shift;

	m_Mur_Coeff_nyP = m_Op_mur->m_Mur_Coeff_nyP;
	m_Mur_Coeff_nyPP = m_Op_mur->m_Mur_Coeff_nyPP;

	m_volt_nyP = Create2DArray(m_numLines);
	m_volt_nyPP = Create2DArray(m_numLines);

	//find if some excitation is on this mur-abc and find the max length of this excite, so that the abc can start after the excitation is done...
	int maxDelay=-1;
	for (unsigned int n=0;n<m_Op_mur->m_Op->Exc->E_Count;++n)
	{
		if ( ((m_Op_mur->m_Op->Exc->E_dir[n]==m_nyP) || (m_Op_mur->m_Op->Exc->E_dir[n]==m_nyPP)) && (m_Op_mur->m_Op->Exc->E_index[m_ny][n]==m_LineNr) )
		{
			if ((int)m_Op_mur->m_Op->Exc->E_delay[n]>maxDelay)
				maxDelay = (int)m_Op_mur->m_Op->Exc->E_delay[n];
		}
	}
	m_start_TS = 0;
	if (maxDelay>=0)
	{
		m_start_TS = maxDelay + m_Op_mur->m_Op->Exc->Length + 10; //give it some extra timesteps, for the excitation to travel at least one cell away
		cerr << "Engine_Ext_Mur_ABC::Engine_Ext_Mur_ABC: Warning: Excitation inside the Mur-ABC #" <<  m_ny << "-" << (int)(m_LineNr>0) << " found!!!!  Mur-ABC will be switched on after excitation is done at " << m_start_TS << " timesteps!!! " << endl;
	}
}

Engine_Ext_Mur_ABC::~Engine_Ext_Mur_ABC()
{
	Delete2DArray(m_volt_nyP,m_numLines);
	m_volt_nyP = NULL;
	Delete2DArray(m_volt_nyPP,m_numLines);
	m_volt_nyPP = NULL;
}

void Engine_Ext_Mur_ABC::DoPreVoltageUpdates()
{
	if (IsActive()==false) return;
	if (m_Eng==NULL) return;
	unsigned int pos[] = {0,0,0};
	unsigned int pos_shift[] = {0,0,0};
	pos[m_ny] = m_LineNr;
	pos_shift[m_ny] = m_LineNr_Shift;

	for (pos[m_nyP]=0;pos[m_nyP]<m_numLines[0];++pos[m_nyP])
	{
		pos_shift[m_nyP] = pos[m_nyP];
		for (pos[m_nyPP]=0;pos[m_nyPP]<m_numLines[1];++pos[m_nyPP])
		{
			pos_shift[m_nyPP] = pos[m_nyPP];
			m_volt_nyP[pos[m_nyP]][pos[m_nyPP]] = m_Eng->GetVolt(m_nyP,pos_shift) - m_Op_mur->m_Mur_Coeff_nyP[pos[m_nyP]][pos[m_nyPP]] * m_Eng->GetVolt(m_nyP,pos);
			m_volt_nyPP[pos[m_nyP]][pos[m_nyPP]] = m_Eng->GetVolt(m_nyPP,pos_shift) - m_Op_mur->m_Mur_Coeff_nyPP[pos[m_nyP]][pos[m_nyPP]] * m_Eng->GetVolt(m_nyPP,pos);
		}
	}
}

void Engine_Ext_Mur_ABC::DoPostVoltageUpdates()
{
	if (IsActive()==false) return;
	if (m_Eng==NULL) return;
	unsigned int pos[] = {0,0,0};
	unsigned int pos_shift[] = {0,0,0};
	pos[m_ny] = m_LineNr;
	pos_shift[m_ny] = m_LineNr_Shift;

	for (pos[m_nyP]=0;pos[m_nyP]<m_numLines[0];++pos[m_nyP])
	{
		pos_shift[m_nyP] = pos[m_nyP];
		for (pos[m_nyPP]=0;pos[m_nyPP]<m_numLines[1];++pos[m_nyPP])
		{
			pos_shift[m_nyPP] = pos[m_nyPP];
			m_volt_nyP[pos[m_nyP]][pos[m_nyPP]] += m_Op_mur->m_Mur_Coeff_nyP[pos[m_nyP]][pos[m_nyPP]] * m_Eng->GetVolt(m_nyP,pos_shift);
			m_volt_nyPP[pos[m_nyP]][pos[m_nyPP]] += m_Op_mur->m_Mur_Coeff_nyPP[pos[m_nyP]][pos[m_nyPP]] * m_Eng->GetVolt(m_nyPP,pos_shift);
		}
	}
}

void Engine_Ext_Mur_ABC::Apply2Voltages()
{
	if (IsActive()==false) return;
	if (m_Eng==NULL) return;
	unsigned int pos[] = {0,0,0};
	pos[m_ny] = m_LineNr;

	for (pos[m_nyP]=0;pos[m_nyP]<m_numLines[0];++pos[m_nyP])
	{
		for (pos[m_nyPP]=0;pos[m_nyPP]<m_numLines[1];++pos[m_nyPP])
		{
			m_Eng->GetVolt(m_nyP,pos) = m_volt_nyP[pos[m_nyP]][pos[m_nyPP]];
			m_Eng->GetVolt(m_nyPP,pos) = m_volt_nyPP[pos[m_nyP]][pos[m_nyPP]];
		}
	}
}
