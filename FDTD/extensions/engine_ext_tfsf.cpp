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
#include <climits>
#include "operator_ext_tfsf.h"
#include "FDTD/engine_sse.h"

Engine_Ext_TFSF::Engine_Ext_TFSF(Operator_Ext_TFSF* op_ext) : Engine_Extension(op_ext)
{
	m_Op_TFSF = op_ext;
	m_Priority = ENG_EXT_PRIO_TFSF;
}

Engine_Ext_TFSF::~Engine_Ext_TFSF()
{
}

// the flat hooks: the whole domain as the box
void Engine_Ext_TFSF::DoPostVoltageUpdates()
{
	const Engine_Box all = {{0, 0, 0}, {UINT_MAX, UINT_MAX, UINT_MAX}};   // the planes clip it
	PlanesImpl(true, all, m_Eng->GetNumberOfTimesteps());
}

void Engine_Ext_TFSF::DoPostCurrentUpdates()
{
	const Engine_Box all = {{0, 0, 0}, {UINT_MAX, UINT_MAX, UINT_MAX}};   // the planes clip it
	PlanesImpl(false, all, m_Eng->GetNumberOfTimesteps());
}

// the position in the excitation signal of a delay of d timesteps at timestep ts
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

// The planes of the TF/SF box in the box \a box: the lower and the upper plane normal to each
// axis, per cell both components, so a cell on an edge of the TF/SF box gets its additions in
// the same order in the flat and in the box hooks.
void Engine_Ext_TFSF::PlanesImpl(bool volt, const Engine_Box& box, unsigned int ts)
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
	PlanesImpl(true, box, ts);
}

void Engine_Ext_TFSF::DoPostCurrentUpdatesBox(const Engine_Box& box, unsigned int ts)
{
	PlanesImpl(false, box, ts);
}
