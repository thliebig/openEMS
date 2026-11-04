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

#include <algorithm>
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
								   + (1.0f-m_Op_TFSF->m_VoltDelayDelta[n][0][0][ui_pos])*m_Op_TFSF->m_VoltAmp[n][0][0][ui_pos]*signal[m_DelayLookup[  m_Op_TFSF->m_VoltDelay[n][0][0][ui_pos]]]
								   +      m_Op_TFSF->m_VoltDelayDelta[n][0][0][ui_pos] *m_Op_TFSF->m_VoltAmp[n][0][0][ui_pos]*signal[m_DelayLookup[1+m_Op_TFSF->m_VoltDelay[n][0][0][ui_pos]]] );

					m_Eng->SetVolt(nPP,pos, m_Eng->GetVolt(nPP,pos)
								   + (1.0f-m_Op_TFSF->m_VoltDelayDelta[n][0][1][ui_pos])*m_Op_TFSF->m_VoltAmp[n][0][1][ui_pos]*signal[m_DelayLookup[  m_Op_TFSF->m_VoltDelay[n][0][1][ui_pos]]]
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
								   + (1.0f-m_Op_TFSF->m_VoltDelayDelta[n][1][0][ui_pos])*m_Op_TFSF->m_VoltAmp[n][1][0][ui_pos]*signal[m_DelayLookup[  m_Op_TFSF->m_VoltDelay[n][1][0][ui_pos]]]
								   +      m_Op_TFSF->m_VoltDelayDelta[n][1][0][ui_pos] *m_Op_TFSF->m_VoltAmp[n][1][0][ui_pos]*signal[m_DelayLookup[1+m_Op_TFSF->m_VoltDelay[n][1][0][ui_pos]]] );

					m_Eng->SetVolt(nPP,pos, m_Eng->GetVolt(nPP,pos)
								   + (1.0f-m_Op_TFSF->m_VoltDelayDelta[n][1][1][ui_pos])*m_Op_TFSF->m_VoltAmp[n][1][1][ui_pos]*signal[m_DelayLookup[  m_Op_TFSF->m_VoltDelay[n][1][1][ui_pos]]]
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
								   + (1.0f-m_Op_TFSF->m_CurrDelayDelta[n][0][0][ui_pos])*m_Op_TFSF->m_CurrAmp[n][0][0][ui_pos]*signal[m_DelayLookup[  m_Op_TFSF->m_CurrDelay[n][0][0][ui_pos]]]
								   +      m_Op_TFSF->m_CurrDelayDelta[n][0][0][ui_pos] *m_Op_TFSF->m_CurrAmp[n][0][0][ui_pos]*signal[m_DelayLookup[1+m_Op_TFSF->m_CurrDelay[n][0][0][ui_pos]]] );

					m_Eng->SetCurr(nPP,pos, m_Eng->GetCurr(nPP,pos)
								   + (1.0f-m_Op_TFSF->m_CurrDelayDelta[n][0][1][ui_pos])*m_Op_TFSF->m_CurrAmp[n][0][1][ui_pos]*signal[m_DelayLookup[  m_Op_TFSF->m_CurrDelay[n][0][1][ui_pos]]]
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
								   + (1.0f-m_Op_TFSF->m_CurrDelayDelta[n][1][0][ui_pos])*m_Op_TFSF->m_CurrAmp[n][1][0][ui_pos]*signal[m_DelayLookup[  m_Op_TFSF->m_CurrDelay[n][1][0][ui_pos]]]
								   +      m_Op_TFSF->m_CurrDelayDelta[n][1][0][ui_pos] *m_Op_TFSF->m_CurrAmp[n][1][0][ui_pos]*signal[m_DelayLookup[1+m_Op_TFSF->m_CurrDelay[n][1][0][ui_pos]]] );

					m_Eng->SetCurr(nPP,pos, m_Eng->GetCurr(nPP,pos)
								   + (1.0f-m_Op_TFSF->m_CurrDelayDelta[n][1][1][ui_pos])*m_Op_TFSF->m_CurrAmp[n][1][1][ui_pos]*signal[m_DelayLookup[  m_Op_TFSF->m_CurrDelay[n][1][1][ui_pos]]]
								   +      m_Op_TFSF->m_CurrDelayDelta[n][1][1][ui_pos] *m_Op_TFSF->m_CurrAmp[n][1][1][ui_pos]*signal[m_DelayLookup[1+m_Op_TFSF->m_CurrDelay[n][1][1][ui_pos]]] );

					++pos[nPP];
					++ui_pos;
				}
				++pos[nP];
			}
		}
	}
}

// the position in the excitation signal of a delay of d timesteps at timestep ts, as in m_DelayLookup
static inline unsigned int Delay(unsigned int ts, unsigned int d, unsigned int length, int p)
{
	unsigned int v;
	if (ts<d)
		v = 0;
	else if ((ts-d >= length) && (p==0))
		v = 0;
	else
		v = ts-d;
	if (p>0)
		v = v % p;
	return v;
}

// The planes of the flat hooks restricted to the box, in the same order, so a cell on
// an edge of the TF/SF box gets its additions in the same order too.
void Engine_Ext_TFSF::BoxImpl(bool volt, const Engine_Box& box, unsigned int ts)
{
	const Operator_Ext_TFSF* op = m_Op_TFSF;
	Engine* eng = m_Eng;
	const unsigned int length = op->m_Exc->GetLength();
	const int p = int(op->m_Exc->GetSignalPeriod()/op->m_Exc->GetTimestep());
	// a voltage gets the H-field signal and the other way round
	const FDTD_FLOAT* signal = volt ? op->m_Exc->GetCurrentSignal() : op->m_Exc->GetVoltageSignal();
	unsigned int pos[3];
	for (int n=0;n<3;++n)
	{
		const int nP = (n+1)%3, nPP = (n+2)%3;
		for (int side=0; side<2; ++side)
		{
			if (!op->m_ActiveDir[n][side])
				continue;
			pos[n] = (side==1) ? op->m_Stop[n] : (volt ? op->m_Start[n] : op->m_Start[n]-1);
			if ((pos[n]<box.start[n]) || (pos[n]>=box.stop[n]))
				continue;
			const unsigned int i0 = std::max(op->m_Start[nP], box.start[nP]);
			const unsigned int i1 = std::min(op->m_Start[nP]+op->m_numLines[nP], box.stop[nP]);
			const unsigned int j0 = std::max(op->m_Start[nPP], box.start[nPP]);
			const unsigned int j1 = std::min(op->m_Start[nPP]+op->m_numLines[nPP], box.stop[nPP]);
			for (pos[nP]=i0; pos[nP]<i1; ++pos[nP])
				for (pos[nPP]=j0; pos[nPP]<j1; ++pos[nPP])
				{
					const unsigned int ui_pos = (pos[nP]-op->m_Start[nP])*op->m_numLines[nPP] + (pos[nPP]-op->m_Start[nPP]);
					for (int c=0; c<2; ++c)
					{
						const int ny = (c==0) ? nP : nPP;
						const unsigned int d = volt ? op->m_VoltDelay[n][side][c][ui_pos] : op->m_CurrDelay[n][side][c][ui_pos];
						const FDTD_FLOAT dd = volt ? op->m_VoltDelayDelta[n][side][c][ui_pos] : op->m_CurrDelayDelta[n][side][c][ui_pos];
						const FDTD_FLOAT amp = volt ? op->m_VoltAmp[n][side][c][ui_pos] : op->m_CurrAmp[n][side][c][ui_pos];
						// the evaluation order of the flat hooks: (field + first) + second
						const FDTD_FLOAT s0 = signal[Delay(ts, d, length, p)], s1 = signal[Delay(ts, d+1, length, p)];
						if (volt)
							eng->SetVolt(ny, pos, eng->GetVolt(ny, pos) + (1.0f-dd)*amp*s0 + dd*amp*s1);
						else
							eng->SetCurr(ny, pos, eng->GetCurr(ny, pos) + (1.0f-dd)*amp*s0 + dd*amp*s1);
					}
				}
		}
	}
}

void Engine_Ext_TFSF::DoPostVoltageUpdatesBox(const Engine_Box& box, unsigned int ts)
{
	BoxImpl(true, box, ts);
}

void Engine_Ext_TFSF::DoPostCurrentUpdatesBox(const Engine_Box& box, unsigned int ts)
{
	BoxImpl(false, box, ts);
}
