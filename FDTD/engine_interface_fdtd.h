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

#ifndef ENGINE_INTERFACE_FDTD_H
#define ENGINE_INTERFACE_FDTD_H

#include <cmath>

#include "Common/engine_interface_base.h"
#include "operator.h"
#include "engine.h"

class Engine_Interface_FDTD : public Engine_Interface_Base
{
public:
	Engine_Interface_FDTD(Operator* op);
	virtual ~Engine_Interface_FDTD();

	//! Set the FDTD operator
	virtual void SetFDTDOperator(Operator* op);
	//! Set the FDTD engine
	virtual void SetFDTDEngine(Engine* eng) {m_Eng=eng;}

	//! Get the FDTD engine in case direct access is needed. Direct access is not recommended!
	const Engine* GetFDTDEngine() const {return m_Eng;}
	//! Get the FDTD operator in case direct access is needed. Direct access is not recommended!
	const Operator* GetFDTDOperator() const {return m_Op;}

	virtual double* GetEField(const unsigned int* pos, double* out) const;
	virtual double* GetHField(const unsigned int* pos, double* out) const;
	virtual double* GetJField(const unsigned int* pos, double* out) const;
	virtual double* GetRotHField(const unsigned int* pos, double* out) const;
	virtual double* GetDField(const unsigned int* pos, double* out) const;
	virtual double* GetBField(const unsigned int* pos, double* out) const;

	virtual double CalcVoltageIntegral(const unsigned int* start, const unsigned int* stop) const;

	virtual double GetTime(bool dualTime=false) const {return ((double)m_Eng->GetNumberOfTimesteps() + (double)dualTime*0.5)*m_Op->GetTimestep();};
	virtual unsigned int GetNumberOfTimesteps() const {return m_Eng->GetNumberOfTimesteps();}

	virtual double CalcFastEnergy() const;

	//! The interpolated E field at \a pos from raw(n, pos), the voltage of edge n at pos divided by its length
	/*!
	  GetEField() with raw = GetRawField(); the dumps of the GPU engine read the raw values
	  from field snapshots and frequency domain sums instead, with the same interpolation.
	  */
	template <typename T, typename Raw>
	void InterpolateField(const unsigned int* pos, T* out, Raw raw) const
	{
		unsigned int seamPos[3];
		if (m_ClosedMesh)
		{
			for (int n=0; n<3; ++n)
				seamPos[n] = GetSeamLine(n, pos[n]);
			pos = seamPos;
		}
		unsigned int iPos[] = {pos[0],pos[1],pos[2]};
		int nP,nPP;
		double delta;
		switch (m_InterpolType)
		{
		default:
		case NO_INTERPOLATION:
			for (int n=0; n<3; ++n)
				out[n] = raw(n, pos);
			break;
		case NODE_INTERPOLATE:
			for (int n=0; n<3; ++n)
			{
				if (pos[n]==m_Op->GetNumberOfLines(n, true)-1)  // use only the "lower value" at the upper bound
				{
					--iPos[n];
					out[n] = raw(n, iPos);
					++iPos[n];
					continue;
				}
				delta = m_Op->GetEdgeLength(n,iPos);
				out[n] = raw(n, iPos);
				if (delta==0)
				{
					out[n]=0;
					continue;
				}
				if (pos[n]==0) // use only the "upper value" at the lower bound
					continue;
				--iPos[n];
				double deltaDown = m_Op->GetEdgeLength(n,iPos);
				double deltaRel = delta / (delta+deltaDown);
				out[n] = out[n]*(1.0-deltaRel) + raw(n, iPos)*deltaRel;
				++iPos[n];
			}
			break;
		case CELL_INTERPOLATE:
			for (int n=0; n<3; ++n)
			{
				nP = (n+1)%3;
				nPP = (n+2)%3;
				if ((pos[0]==m_Op->GetNumberOfLines(0,true)-1) || (pos[1]==m_Op->GetNumberOfLines(1,true)-1) || (pos[2]==m_Op->GetNumberOfLines(2,true)-1))
				{
					out[n] = 0; //electric field outside the field domain is always zero
					continue;
				}
				out[n]=raw(n, iPos);
				++iPos[nP];
				out[n]+=raw(n, iPos);
				++iPos[nPP];
				out[n]+=raw(n, iPos);
				--iPos[nP];
				out[n]+=raw(n, iPos);
				--iPos[nPP];
				out[n]/=4;
			}
			break;
		}
	}
	//! The interpolated H field at \a pos from raw(n, pos), the current of edge n at pos divided by its dual length, see InterpolateField()
	template <typename T, typename Raw>
	void InterpolateDualField(const unsigned int* pos, T* out, Raw raw) const
	{
		unsigned int seamPos[3];
		if (m_ClosedMesh)
		{
			for (int n=0; n<3; ++n)
				seamPos[n] = GetSeamLine(n, pos[n]);
			pos = seamPos;
		}
		unsigned int iPos[] = {pos[0],pos[1],pos[2]};
		int nP,nPP;
		double delta;
		switch (m_InterpolType)
		{
		default:
		case NO_INTERPOLATION:
			out[0] = raw(0, pos);
			out[1] = raw(1, pos);
			out[2] = raw(2, pos);
			break;
		case NODE_INTERPOLATE:
			for (int n=0; n<3; ++n)
			{
				nP = (n+1)%3;
				nPP = (n+2)%3;
				if ((pos[0]==m_Op->GetNumberOfLines(0,true)-1) || (pos[1]==m_Op->GetNumberOfLines(1,true)-1) || (pos[2]==m_Op->GetNumberOfLines(2,true)-1) || (pos[nP]==0) || (pos[nPP]==0))
				{
					out[n] = 0;
					continue;
				}
				out[n] = raw(n, iPos);
				--iPos[nP];
				out[n]+= raw(n, iPos);
				--iPos[nPP];
				out[n]+= raw(n, iPos);
				++iPos[nP];
				out[n]+= raw(n, iPos);
				++iPos[nPP];
				out[n]/=4;
			}
			break;
		case CELL_INTERPOLATE:
			for (int n=0; n<3; ++n)
			{
				delta = m_Op->GetEdgeLength(n,iPos,true);
				out[n] = raw(n, iPos);
				if ((pos[n]>=m_Op->GetNumberOfLines(n,true)-1))
				{
					out[n] = 0; //magnetic field on the outer boundaries is always zero
					continue;
				}
				++iPos[n];
				double deltaUp = m_Op->GetEdgeLength(n,iPos,true);
				double deltaRel = delta / (delta+deltaUp);
				out[n] = out[n]*(1.0-deltaRel) + raw(n, iPos)*deltaRel;
				--iPos[n];
			}
			break;
		}
	}

protected:
	Operator* m_Op;
	Engine* m_Eng;

	//! The closed (periodic) directions of the mesh, see Operator_Base::GetClosedMeshLines()
	bool m_ClosedMesh;
	bool m_Closed[3];
	unsigned int m_ClosedFirst[3], m_ClosedLast[3];
	void UpdateClosedMesh();
	//! The line of direction \a n an interpolated field at \a line is evaluated at: on the seam of a closed mesh its other side, where the interpolation has the neighbours
	unsigned int GetSeamLine(int n, unsigned int line) const;

	//! Internal method to get an interpolated field of a given type. (0: E, 1: J, 2: rotH, 3: D)
	virtual double* GetRawInterpolatedField(const unsigned int* pos, double* out, int type) const;
	//! Internal method to get a raw field of a given type. (0: E, 1: J, 2: rotH, 3: D)
	virtual double GetRawField(unsigned int n, const unsigned int* pos, int type) const;

	//! Internal method to get an interpolated dual field of a given type. (0: H, 1: B)
	virtual double* GetRawInterpolatedDualField(const unsigned int* pos, double* out, int type) const;
	//! Internal method to get a raw dual field of a given type. (0: H, 1: B)
	virtual double GetRawDualField(unsigned int n, const unsigned int* pos, int type) const;
};

#endif // ENGINE_INTERFACE_FDTD_H
