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

#include "engine_interface_fdtd.h"

Engine_Interface_FDTD::Engine_Interface_FDTD() : Engine_Interface_Base()
{
	m_Op = NULL;
	m_Eng = NULL;
}

Engine_Interface_FDTD::Engine_Interface_FDTD(Operator* op, Engine* eng) : Engine_Interface_Base()
{
	SetOperator(op);
	SetEngine(eng);
}

Engine_Interface_FDTD::~Engine_Interface_FDTD()
{
}

double* Engine_Interface_FDTD::GetEField(const unsigned int* pos, double* out) const
{
	unsigned int iPos[] = {pos[0],pos[1],pos[2]};
	int nP,nPP;
	double delta;
	switch (m_InterpolType)
	{
	default:
	case NO_INTERPOLATION:
		for (int n=0;n<3;++n)
		{
			delta =  m_Op->GetEdgeLength(n,pos,false);
			if (delta)
				out[n] = m_Eng->GetVolt(n,pos) / delta;
			else
				out[n] = 0.0;
		}
		break;
	case NODE_INTERPOLATE:
		for (int n=0;n<3;++n)
		{
			delta = m_Op->GetMeshDelta(n,iPos);
			out[n] = m_Eng->GetVolt(n,iPos);
			if (delta==0)
			{
				out[n]=0;
				continue;
			}
			if (pos[n]==0)
			{
				out[n] /= (delta * 2.0); //make it consistant with upper PEC boundary
				continue;
			}
			--iPos[n];
			double deltaDown = m_Op->GetMeshDelta(n,iPos);
			double deltaRel = delta / (delta+deltaDown);
			out[n] = out[n]*(1.0-deltaRel)/delta + (double)m_Eng->GetVolt(n,iPos)/deltaDown*deltaRel;
			++iPos[n];
		}
		break;
	case CELL_INTERPOLATE:
		for (int n=0;n<3;++n)
		{
			nP = (n+1)%3;
			nPP = (n+2)%3;
			if ((pos[0]==m_Op->GetNumberOfLines(0)-1) || (pos[1]==m_Op->GetNumberOfLines(1)-1) || (pos[2]==m_Op->GetNumberOfLines(2)-1))
			{
				out[n] = 0; //electric field outside the field domain is always zero
				continue;
			}
			delta = m_Op->GetMeshDelta(n,iPos);
			if (delta)
				out[n]=m_Eng->GetVolt(n,iPos)/delta;
			++iPos[nP];
			delta = m_Op->GetMeshDelta(n,iPos);
			if (delta)
				out[n]+=m_Eng->GetVolt(n,iPos)/delta;
			++iPos[nPP];
			delta = m_Op->GetMeshDelta(n,iPos);
			if (delta)
				out[n]+=m_Eng->GetVolt(n,iPos)/delta;
			--iPos[nP];
			delta = m_Op->GetMeshDelta(n,iPos);
			if (delta)
				out[n]+=m_Eng->GetVolt(n,iPos)/delta;
			--iPos[nPP];
			out[n]/=4;
		}
		break;
	}
	return out;
}

double* Engine_Interface_FDTD::GetHField(const unsigned int* pos, double* out) const
{
	unsigned int iPos[] = {pos[0],pos[1],pos[2]};
	int nP,nPP;
	double delta;
	switch (m_InterpolType)
	{
	default:
	case NO_INTERPOLATION:
		out[0] = m_Eng->GetCurr(0,pos) / m_Op->GetEdgeLength(0,pos,true);
		out[1] = m_Eng->GetCurr(1,pos) / m_Op->GetEdgeLength(1,pos,true);
		out[2] = m_Eng->GetCurr(2,pos) / m_Op->GetEdgeLength(2,pos,true);
		break;
	case NODE_INTERPOLATE:
		for (int n=0;n<3;++n)
		{
			nP = (n+1)%3;
			nPP = (n+2)%3;
			if ((pos[0]==m_Op->GetNumberOfLines(0)-1) || (pos[1]==m_Op->GetNumberOfLines(1)-1) || (pos[2]==m_Op->GetNumberOfLines(2)-1) || (pos[nP]==0) || (pos[nPP]==0))
			{
				out[n] = 0;
				continue;
			}
			out[n]=m_Eng->GetCurr(n,iPos)/m_Op->GetMeshDelta(n,iPos,true);
			--iPos[nP];
			out[n]+=m_Eng->GetCurr(n,iPos)/m_Op->GetMeshDelta(n,iPos,true);
			--iPos[nPP];
			out[n]+=m_Eng->GetCurr(n,iPos)/m_Op->GetMeshDelta(n,iPos,true);
			++iPos[nP];
			out[n]+=m_Eng->GetCurr(n,iPos)/m_Op->GetMeshDelta(n,iPos,true);
			++iPos[nPP];
			out[n]/=4;
		}
		break;
	case CELL_INTERPOLATE:
		for (int n=0;n<3;++n)
		{
			delta = m_Op->GetMeshDelta(n,iPos,true);
			out[n] = m_Eng->GetCurr(n,iPos);
			if ((pos[n]>=m_Op->GetNumberOfLines(n)-1))
			{
				out[n] = 0; //magnetic field on the outer boundaries is always zero
				continue;
			}
			++iPos[n];
			double deltaUp = m_Op->GetMeshDelta(n,iPos,true);
			double deltaRel = delta / (delta+deltaUp);
			out[n] = out[n]*(1.0-deltaRel)/delta + (double)m_Eng->GetCurr(n,iPos)/deltaUp*deltaRel;
			--iPos[n];
		}
		break;
	}

	return out;
}
