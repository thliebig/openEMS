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
#include <stdexcept>

using std::cerr;
using std::endl;

Engine_Interface_FDTD::Engine_Interface_FDTD(Operator* op) : Engine_Interface_Base(op)
{
	if (op==NULL)
		throw std::runtime_error("Engine_Interface_FDTD::Engine_Interface_FDTD: Error: Operator is not set!");
	m_Op = op;
	m_Eng = m_Op->GetEngine();
	if (m_Eng==NULL)
		throw std::runtime_error("Engine_Interface_FDTD::Engine_Interface_FDTD: Error: Engine is not set!");
	UpdateClosedMesh();
}

void Engine_Interface_FDTD::SetFDTDOperator(Operator* op)
{
	SetOperator(op);
	m_Op=op;
	UpdateClosedMesh();
}

void Engine_Interface_FDTD::UpdateClosedMesh()
{
	m_ClosedMesh = false;
	for (int n=0; n<3; ++n)
	{
		m_Closed[n] = m_Op->GetClosedMeshLines(n, m_ClosedFirst[n], m_ClosedLast[n]);
		m_ClosedMesh |= m_Closed[n];
	}
}

unsigned int Engine_Interface_FDTD::GetSeamLine(int n, unsigned int line) const
{
	if (!m_Closed[n])
		return line;
	// node interpolation reads the lower neighbours, cell interpolation the upper ones
	if ((m_InterpolType==NODE_INTERPOLATE) && (line==m_ClosedFirst[n]))
		return m_ClosedLast[n];
	if ((m_InterpolType==CELL_INTERPOLATE) && (line==m_ClosedLast[n]))
		return m_ClosedFirst[n];
	return line;
}

Engine_Interface_FDTD::~Engine_Interface_FDTD()
{
}

double* Engine_Interface_FDTD::GetEField(const unsigned int* pos, double* out) const
{
	return GetRawInterpolatedField(pos, out, 0);
}

double* Engine_Interface_FDTD::GetJField(const unsigned int* pos, double* out) const
{
	return GetRawInterpolatedField(pos, out, 1);
}

double* Engine_Interface_FDTD::GetDField(const unsigned int* pos, double* out) const
{
	return GetRawInterpolatedField(pos, out, 3);
}

double* Engine_Interface_FDTD::GetRotHField(const unsigned int* pos, double* out) const
{
	return GetRawInterpolatedField(pos, out, 2);
}

double* Engine_Interface_FDTD::GetRawInterpolatedField(const unsigned int* pos, double* out, int type) const
{
	InterpolateField(pos, out, [this, type](unsigned int n, const unsigned int* p) {return GetRawField(n, p, type);});
	return out;
}

double* Engine_Interface_FDTD::GetHField(const unsigned int* pos, double* out) const
{
	return GetRawInterpolatedDualField(pos, out, 0);
}

double* Engine_Interface_FDTD::GetBField(const unsigned int* pos, double* out) const
{
	return GetRawInterpolatedDualField(pos, out, 1);
}

double Engine_Interface_FDTD::GetRawDualField(unsigned int n, const unsigned int* pos, int type) const
{
	double value = m_Eng->GetCurr(n,pos[0],pos[1],pos[2]);
	double delta = m_Op->GetEdgeLength(n,pos,true);
	if ((type==0) && (delta))
		return value/delta;
	if ((type==1) && (m_Op->m_mueR_ptr) && (delta))
	{
		ArrayLib::ArrayNIJK<float>& m_mueR = *m_Op->m_mueR_ptr;
		return value * m_mueR(n, pos[0], pos[1], pos[2]) / delta;
	}
	return 0.0;
}

double* Engine_Interface_FDTD::GetRawInterpolatedDualField(const unsigned int* pos, double* out, int type) const
{
	InterpolateDualField(pos, out, [this, type](unsigned int n, const unsigned int* p) {return GetRawDualField(n, p, type);});
	return out;
}

double Engine_Interface_FDTD::CalcVoltageIntegral(const unsigned int* start, const unsigned int* stop) const
{
	if (((start[0]!=stop[0]) + (start[1]!=stop[1]) + (start[2]!=stop[2]))!=1)
	{
		cerr << "Engine_Interface_FDTD::CalcVoltageIntegral: Error, only a 1D/line integration is allowed" << endl;
		return 0;
	}
	//cerr << "CalcVoltageIntegral" << start[0] << ", " << start[1] << ", " << start[2] << " -> " << stop[0] << ", " << stop[1] << ", " << stop[2] << ", " << endl;
	double result=0;
	for (int n=0; n<3; ++n)
	{
		if (start[n]<stop[n])
		{
			unsigned int pos[3]={start[0],start[1],start[2]};
			for (; pos[n]<stop[n]; ++pos[n])
				result += m_Eng->GetVolt(n,pos[0],pos[1],pos[2]);

		}
		else
		{
			unsigned int pos[3]={stop[0],stop[1],stop[2]};
			for (; pos[n]<start[n]; ++pos[n])
				result -= m_Eng->GetVolt(n,pos[0],pos[1],pos[2]);
		}
	}
	return result;
}

double Engine_Interface_FDTD::GetRawField(unsigned int n, const unsigned int* pos, int type) const
{
	double value = m_Eng->GetVolt(n,pos[0],pos[1],pos[2]);
	double delta = m_Op->GetEdgeLength(n,pos);
	if ((type==0) && (delta))
		return value/delta;
	if ((type==1) && (m_Op->m_kappa_ptr) && (delta)) {
		ArrayLib::ArrayNIJK<float>& kappa = *m_Op->m_kappa_ptr;
		return value * kappa(n, pos[0], pos[1], pos[2]) / delta;
	}
	if ((type==3) && (m_Op->m_epsR_ptr) && (delta)) {
		ArrayLib::ArrayNIJK<float>& epsR = *m_Op->m_epsR_ptr;
		return value * epsR(n, pos[0], pos[1], pos[2]) / delta;
	}
	if (type==2) //calc rot(H)
	{
		int nP = (n+1)%3;
		int nPP = (n+2)%3;
		unsigned int locPos[] = {pos[0],pos[1],pos[2]};
		double area = m_Op->GetEdgeArea(n,pos);
		value  = m_Eng->GetCurr(nPP,pos);
		value -= m_Eng->GetCurr(nP,pos);
		if (pos[nPP]>0)
		{
			--locPos[nPP];
			value += m_Eng->GetCurr(nP,locPos);
			++locPos[nPP];
		}
		if (pos[nP]>0)
		{
			--locPos[nP];
			value -= m_Eng->GetCurr(nPP,locPos);
		}
		return value/area;
	}

	return 0.0;
}

double Engine_Interface_FDTD::CalcFastEnergy() const
{
	double E_energy=0.0;
	double H_energy=0.0;

	unsigned int pos[3];
	if (m_Eng->GetType()==Engine::BASIC)
	{
		for (pos[0]=0; pos[0]<m_Op->GetNumberOfLines(0)-1; ++pos[0])
		{
			for (pos[1]=0; pos[1]<m_Op->GetNumberOfLines(1)-1; ++pos[1])
			{
				for (pos[2]=0; pos[2]<m_Op->GetNumberOfLines(2)-1; ++pos[2])
				{
					E_energy+=m_Eng->Engine::GetVolt(0,pos[0],pos[1],pos[2]) * m_Eng->Engine::GetVolt(0,pos[0],pos[1],pos[2]);
					E_energy+=m_Eng->Engine::GetVolt(1,pos[0],pos[1],pos[2]) * m_Eng->Engine::GetVolt(1,pos[0],pos[1],pos[2]);
					E_energy+=m_Eng->Engine::GetVolt(2,pos[0],pos[1],pos[2]) * m_Eng->Engine::GetVolt(2,pos[0],pos[1],pos[2]);

					H_energy+=m_Eng->Engine::GetCurr(0,pos[0],pos[1],pos[2]) * m_Eng->Engine::GetCurr(0,pos[0],pos[1],pos[2]);
					H_energy+=m_Eng->Engine::GetCurr(1,pos[0],pos[1],pos[2]) * m_Eng->Engine::GetCurr(1,pos[0],pos[1],pos[2]);
					H_energy+=m_Eng->Engine::GetCurr(2,pos[0],pos[1],pos[2]) * m_Eng->Engine::GetCurr(2,pos[0],pos[1],pos[2]);
				}
			}
		}
	}
	else
	{
		for (pos[0]=0; pos[0]<m_Op->GetNumberOfLines(0)-1; ++pos[0])
		{
			for (pos[1]=0; pos[1]<m_Op->GetNumberOfLines(1)-1; ++pos[1])
			{
				for (pos[2]=0; pos[2]<m_Op->GetNumberOfLines(2)-1; ++pos[2])
				{
					E_energy+=m_Eng->GetVolt(0,pos[0],pos[1],pos[2]) * m_Eng->GetVolt(0,pos[0],pos[1],pos[2]);
					E_energy+=m_Eng->GetVolt(1,pos[0],pos[1],pos[2]) * m_Eng->GetVolt(1,pos[0],pos[1],pos[2]);
					E_energy+=m_Eng->GetVolt(2,pos[0],pos[1],pos[2]) * m_Eng->GetVolt(2,pos[0],pos[1],pos[2]);

					H_energy+=m_Eng->GetCurr(0,pos[0],pos[1],pos[2]) * m_Eng->GetCurr(0,pos[0],pos[1],pos[2]);
					H_energy+=m_Eng->GetCurr(1,pos[0],pos[1],pos[2]) * m_Eng->GetCurr(1,pos[0],pos[1],pos[2]);
					H_energy+=m_Eng->GetCurr(2,pos[0],pos[1],pos[2]) * m_Eng->GetCurr(2,pos[0],pos[1],pos[2]);
				}
			}
		}
	}
	return EPS0*E_energy + MUE0*H_energy;
}
