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

#ifndef ENGINE_SSE_H
#define ENGINE_SSE_H

#include "operator.h"
#include "engine.h"

class Engine_sse : public Engine
{
public:
	static Engine_sse* New(const Operator* op);
	virtual ~Engine_sse();

	virtual void Init();
	virtual void Reset();

	//!Iterate a number of timesteps
	virtual bool IterateTS(unsigned int iterTS);

	virtual unsigned int GetNumberOfTimesteps() {return numTS;};

//	virtual f4vector**** GetVoltages() {return volt;};
	virtual FDTD_FLOAT**** GetCurrents() {return curr;};

protected:
	Engine_sse(const Operator* op);
	const Operator* Op;

	virtual void UpdateVoltages();
	virtual void ApplyVoltageExcite();
	virtual void UpdateCurrents();
	virtual void ApplyCurrentExcite();

	unsigned int numLines[3];

	f4vector**** volt_;
	FDTD_FLOAT**** curr;
	unsigned int numTS;
};

#endif // ENGINE_SSE_H
