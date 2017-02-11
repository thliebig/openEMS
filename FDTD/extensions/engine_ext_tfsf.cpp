/*
*	Copyright (C) 2012 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#include "engine_ext_tfsf.h"
#include "operator_ext_tfsf.h"
#include "FDTD/engine_sse.h"

Engine_Ext_TFSF::Engine_Ext_TFSF(Operator_Ext_TFSF* op_ext) : Engine_Extension(op_ext)
{
	m_Op_TFSF = op_ext;
	m_Priority = ENG_EXT_PRIO_TFSF;

	m_DelayLookup = new unsigned int[m_Op_TFSF->m_maxDelay+1];
}

Engine_Ext_TFSF::~Engine_Ext_TFSF()
{
	delete[] m_DelayLookup;
	m_DelayLookup = NULL;
}

void Engine_Ext_TFSF::DoPostVoltageUpdates()
{
	unsigned int numTS = m_Eng->GetNumberOfTimesteps();
	unsigned int length = m_Op_TFSF->m_Exc->GetLength();

	int p = int(m_Op_TFSF->m_Exc->GetSignalPeriod()/m_Op_TFSF->m_Exc->GetTimestep());

	for (unsigned int n=0;n<=m_Op_TFSF->m_maxDelay;++n)
	{
		if ( numTS < n )
			m_DelayLookup[n]=0;
		else if ((numTS-n >= length) && (p==0))
			m_DelayLookup[n]=0;
		else
			m_DelayLookup[n] = numTS - n;
		if (p>0)
			m_DelayLookup[n] = (m_DelayLookup[n] % p);
	}

	//get the current signal since an H-field is added ...
	FDTD_FLOAT* signal =  m_Op_TFSF->m_Exc->GetCurrentSignal();

	int nP,nPP;
	unsigned int ui_pos;
	unsigned int pos[3];
	for (int n=0;n<3;++n)
	{
		nP = (n+1)%3;
		nPP = (n+2)%3;

		// lower plane
		pos[nP] = m_Op_TFSF->m_Start[nP];
		ui_pos = 0;
		if (m_Op_TFSF->m_ActiveDir[n][0])
		{

			for (unsigned int i=0;i<m_Op_TFSF->m_numLines[nP];++i)
			{
				pos[nPP] = m_Op_TFSF->m_Start[nPP];
				for (unsigned int j=0;j<m_Op_TFSF->m_numLines[nPP];++j)
				{
					// current updates
					pos[n] = m_Op_TFSF->m_Start[n];

					m_Eng->SetVolt(nP,pos, m_Eng->GetVolt(nP,pos)
								   + (1.0-m_Op_TFSF->m_VoltDelayDelta[n][0][0][ui_pos])*m_Op_TFSF->m_VoltAmp[n][0][0][ui_pos]*signal[m_DelayLookup[  m_Op_TFSF->m_VoltDelay[n][0][0][ui_pos]]]
								   +      m_Op_TFSF->m_VoltDelayDelta[n][0][0][ui_pos] *m_Op_TFSF->m_VoltAmp[n][0][0][ui_pos]*signal[m_DelayLookup[1+m_Op_TFSF->m_VoltDelay[n][0][0][ui_pos]]] );

					m_Eng->SetVolt(nPP,pos, m_Eng->GetVolt(nPP,pos)
								   + (1.0-m_Op_TFSF->m_VoltDelayDelta[n][0][1][ui_pos])*m_Op_TFSF->m_VoltAmp[n][0][1][ui_pos]*signal[m_DelayLookup[  m_Op_TFSF->m_VoltDelay[n][0][1][ui_pos]]]
								   +      m_Op_TFSF->m_VoltDelayDelta[n][0][1][ui_pos] *m_Op_TFSF->m_VoltAmp[n][0][1][ui_pos]*signal[m_DelayLookup[1+m_Op_TFSF->m_VoltDelay[n][0][1][ui_pos]]] );

					++pos[nPP];
					++ui_pos;
				}
				++pos[nP];
			}
		}

		// upper plane
		pos[nP] = m_Op_TFSF->m_Start[nP];
		ui_pos = 0;
		if (m_Op_TFSF->m_ActiveDir[n][1])
		{

			for (unsigned int i=0;i<m_Op_TFSF->m_numLines[nP];++i)
			{
				pos[nPP] = m_Op_TFSF->m_Start[nPP];
				for (unsigned int j=0;j<m_Op_TFSF->m_numLines[nPP];++j)
				{
					// current updates
					pos[n] = m_Op_TFSF->m_Stop[n];

					m_Eng->SetVolt(nP,pos, m_Eng->GetVolt(nP,pos)
								   + (1.0-m_Op_TFSF->m_VoltDelayDelta[n][1][0][ui_pos])*m_Op_TFSF->m_VoltAmp[n][1][0][ui_pos]*signal[m_DelayLookup[  m_Op_TFSF->m_VoltDelay[n][1][0][ui_pos]]]
								   +      m_Op_TFSF->m_VoltDelayDelta[n][1][0][ui_pos] *m_Op_TFSF->m_VoltAmp[n][1][0][ui_pos]*signal[m_DelayLookup[1+m_Op_TFSF->m_VoltDelay[n][1][0][ui_pos]]] );

					m_Eng->SetVolt(nPP,pos, m_Eng->GetVolt(nPP,pos)
								   + (1.0-m_Op_TFSF->m_VoltDelayDelta[n][1][1][ui_pos])*m_Op_TFSF->m_VoltAmp[n][1][1][ui_pos]*signal[m_DelayLookup[  m_Op_TFSF->m_VoltDelay[n][1][1][ui_pos]]]
								   +      m_Op_TFSF->m_VoltDelayDelta[n][1][1][ui_pos] *m_Op_TFSF->m_VoltAmp[n][1][1][ui_pos]*signal[m_DelayLookup[1+m_Op_TFSF->m_VoltDelay[n][1][1][ui_pos]]] );

					++pos[nPP];
					++ui_pos;
				}
				++pos[nP];
			}
		}
	}
}

void Engine_Ext_TFSF::DoPostCurrentUpdates()
{
	unsigned int numTS = m_Eng->GetNumberOfTimesteps();
	unsigned int length = m_Op_TFSF->m_Exc->GetLength();

	int p = int(m_Op_TFSF->m_Exc->GetSignalPeriod()/m_Op_TFSF->m_Exc->GetTimestep());

	for (unsigned int n=0;n<m_Op_TFSF->m_maxDelay;++n)
	{
		if ( numTS < n )
			m_DelayLookup[n]=0;
		else if ((numTS-n >= length) && (p==0))
			m_DelayLookup[n]=0;
		else
			m_DelayLookup[n] = numTS - n;
		if (p>0)
			m_DelayLookup[n] = (m_DelayLookup[n] % p);
	}

	//get the current signal since an E-field is added ...
	FDTD_FLOAT* signal =  m_Op_TFSF->m_Exc->GetVoltageSignal();

	int nP,nPP;
	unsigned int ui_pos;
	unsigned int pos[3];
	for (int n=0;n<3;++n)
	{
		if (!m_Op_TFSF->m_ActiveDir[n][0] && !m_Op_TFSF->m_ActiveDir[n][1])
			continue;

		nP = (n+1)%3;
		nPP = (n+2)%3;

		// lower plane
		pos[nP] = m_Op_TFSF->m_Start[nP];
		ui_pos = 0;
		if (m_Op_TFSF->m_ActiveDir[n][0])
		{
			for (unsigned int i=0;i<m_Op_TFSF->m_numLines[nP];++i)
			{
				pos[nPP] = m_Op_TFSF->m_Start[nPP];
				for (unsigned int j=0;j<m_Op_TFSF->m_numLines[nPP];++j)
				{
					// current updates
					pos[n] = m_Op_TFSF->m_Start[n]-1;

					m_Eng->SetCurr(nP,pos, m_Eng->GetCurr(nP,pos)
								   + (1.0-m_Op_TFSF->m_CurrDelayDelta[n][0][0][ui_pos])*m_Op_TFSF->m_CurrAmp[n][0][0][ui_pos]*signal[m_DelayLookup[  m_Op_TFSF->m_CurrDelay[n][0][0][ui_pos]]]
								   +      m_Op_TFSF->m_CurrDelayDelta[n][0][0][ui_pos] *m_Op_TFSF->m_CurrAmp[n][0][0][ui_pos]*signal[m_DelayLookup[1+m_Op_TFSF->m_CurrDelay[n][0][0][ui_pos]]] );

					m_Eng->SetCurr(nPP,pos, m_Eng->GetCurr(nPP,pos)
								   + (1.0-m_Op_TFSF->m_CurrDelayDelta[n][0][1][ui_pos])*m_Op_TFSF->m_CurrAmp[n][0][1][ui_pos]*signal[m_DelayLookup[  m_Op_TFSF->m_CurrDelay[n][0][1][ui_pos]]]
								   +      m_Op_TFSF->m_CurrDelayDelta[n][0][1][ui_pos] *m_Op_TFSF->m_CurrAmp[n][0][1][ui_pos]*signal[m_DelayLookup[1+m_Op_TFSF->m_CurrDelay[n][0][1][ui_pos]]] );

					++pos[nPP];
					++ui_pos;
				}
				++pos[nP];
			}
		}

		// upper plane
		pos[nP] = m_Op_TFSF->m_Start[nP];
		ui_pos = 0;
		if (m_Op_TFSF->m_ActiveDir[n][1])
		{
			for (unsigned int i=0;i<m_Op_TFSF->m_numLines[nP];++i)
			{
				pos[nPP] = m_Op_TFSF->m_Start[nPP];
				for (unsigned int j=0;j<m_Op_TFSF->m_numLines[nPP];++j)
				{
					// current updates
					pos[n] = m_Op_TFSF->m_Stop[n];

					m_Eng->SetCurr(nP,pos, m_Eng->GetCurr(nP,pos)
								   + (1.0-m_Op_TFSF->m_CurrDelayDelta[n][1][0][ui_pos])*m_Op_TFSF->m_CurrAmp[n][1][0][ui_pos]*signal[m_DelayLookup[  m_Op_TFSF->m_CurrDelay[n][1][0][ui_pos]]]
								   +      m_Op_TFSF->m_CurrDelayDelta[n][1][0][ui_pos] *m_Op_TFSF->m_CurrAmp[n][1][0][ui_pos]*signal[m_DelayLookup[1+m_Op_TFSF->m_CurrDelay[n][1][0][ui_pos]]] );

					m_Eng->SetCurr(nPP,pos, m_Eng->GetCurr(nPP,pos)
								   + (1.0-m_Op_TFSF->m_CurrDelayDelta[n][1][1][ui_pos])*m_Op_TFSF->m_CurrAmp[n][1][1][ui_pos]*signal[m_DelayLookup[  m_Op_TFSF->m_CurrDelay[n][1][1][ui_pos]]]
								   +      m_Op_TFSF->m_CurrDelayDelta[n][1][1][ui_pos] *m_Op_TFSF->m_CurrAmp[n][1][1][ui_pos]*signal[m_DelayLookup[1+m_Op_TFSF->m_CurrDelay[n][1][1][ui_pos]]] );

					++pos[nPP];
					++ui_pos;
				}
				++pos[nP];
			}
		}
	}
}
