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
	virtual void SetFDTDOperator(Operator* op) {SetOperator(op); m_Op=op;}
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

protected:
	Operator* m_Op;
	Engine* m_Eng;

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
