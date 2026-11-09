/*
*	Copyright (C) 2026 Ivan Iudice
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

#include "engine_ext_mpi.h"
#include "operator_ext_mpi.h"
#include "engine_extension_dispatcher.h"
#include "tools/mpi_domain.h"

#include <chrono>

Engine_Ext_MPI::Engine_Ext_MPI(Operator_Ext_MPI* op_ext) : Engine_Extension(op_ext)
{
	m_Op_Ext_MPI = op_ext;
	SetPriority(ENG_EXT_PRIO_MPI);
	for (int n=0; n<3; ++n)
	{
		size_t size = 2 * (size_t)op_ext->m_numLines[(n+1)%3] * op_ext->m_numLines[(n+2)%3];
		m_SendBuf[n].resize(size);
		m_RecvBuf[n].resize(size);
	}
}

Engine_Ext_MPI::~Engine_Ext_MPI()
{
}

template <typename EngType>
void Engine_Ext_MPI::PackPlane(EngType* eng, bool voltages, int n, unsigned int pos, float* buf) const
{
	const unsigned int* numLines = m_Op_Ext_MPI->m_numLines;
	int nP  = (n+1)%3;
	int nPP = (n+2)%3;
	unsigned int p[3];
	p[n] = pos;
	size_t i = 0;
	for (p[nP]=0; p[nP]<numLines[nP]; ++p[nP])
		for (p[nPP]=0; p[nPP]<numLines[nPP]; ++p[nPP])
		{
			if (voltages)
			{
				buf[i++] = eng->EngType::GetVolt(nP,  p[0], p[1], p[2]);
				buf[i++] = eng->EngType::GetVolt(nPP, p[0], p[1], p[2]);
			}
			else
			{
				buf[i++] = eng->EngType::GetCurr(nP,  p[0], p[1], p[2]);
				buf[i++] = eng->EngType::GetCurr(nPP, p[0], p[1], p[2]);
			}
		}
}

template <typename EngType>
void Engine_Ext_MPI::UnpackPlane(EngType* eng, bool voltages, int n, unsigned int pos, const float* buf) const
{
	const unsigned int* numLines = m_Op_Ext_MPI->m_numLines;
	int nP  = (n+1)%3;
	int nPP = (n+2)%3;
	unsigned int p[3];
	p[n] = pos;
	size_t i = 0;
	for (p[nP]=0; p[nP]<numLines[nP]; ++p[nP])
		for (p[nPP]=0; p[nPP]<numLines[nPP]; ++p[nPP])
		{
			if (voltages)
			{
				eng->EngType::SetVolt(nP,  p[0], p[1], p[2], buf[i++]);
				eng->EngType::SetVolt(nPP, p[0], p[1], p[2], buf[i++]);
			}
			else
			{
				eng->EngType::SetCurr(nP,  p[0], p[1], p[2], buf[i++]);
				eng->EngType::SetCurr(nPP, p[0], p[1], p[2], buf[i++]);
			}
		}
}

template <typename EngType>
void Engine_Ext_MPI::ExchangeImpl(EngType* eng, bool voltages)
{
	typedef std::chrono::steady_clock clock;
	clock::duration pack(0), comm(0);
	const unsigned int* numLines = m_Op_Ext_MPI->m_numLines;
	for (int n=0; n<3; ++n)
	{
		// voltages: the lower rank owns the tangential voltages on the split line, they go up;
		// currents: the upper rank owns the tangential currents above it, they go down.
		// The split line is the last but one line of the lower block, the second of the upper.
		int to   = voltages ? m_Op_Ext_MPI->m_NeighborUp[n]   : m_Op_Ext_MPI->m_NeighborDown[n];
		int from = voltages ? m_Op_Ext_MPI->m_NeighborDown[n] : m_Op_Ext_MPI->m_NeighborUp[n];
		if ((to<0) && (from<0))
			continue;
		unsigned int sendPos = voltages ? numLines[n]-2 : 1;
		unsigned int recvPos = voltages ? 1 : numLines[n]-2;
		clock::time_point t0 = clock::now();
		if (to>=0)
			PackPlane(eng, voltages, n, sendPos, m_SendBuf[n].data());
		clock::time_point t1 = clock::now();
		g_mpi.SendRecv(m_SendBuf[n].data(), to, m_RecvBuf[n].data(), from, (int)m_SendBuf[n].size(), (voltages ? 100 : 200) + n);
		clock::time_point t2 = clock::now();
		if (from>=0)
			UnpackPlane(eng, voltages, n, recvPos, m_RecvBuf[n].data());
		pack += (t1-t0) + (clock::now()-t2);
		comm += t2-t1;
	}
	g_mpi.AddHaloTime(std::chrono::duration<double>(pack).count(), std::chrono::duration<double>(comm).count());
}

void Engine_Ext_MPI::Apply2Voltages()
{
	ENG_DISPATCH_ARGS(ExchangeImpl, true);
}

void Engine_Ext_MPI::Apply2Current()
{
	ENG_DISPATCH_ARGS(ExchangeImpl, false);
}
