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

#include "engine.h"
#include "operator_sse.h"

class Engine_sse : public Engine
{
public:
	static Engine_sse* New(const Operator_sse* op);
	virtual ~Engine_sse();

	virtual void Init();
	virtual void Reset();

	//!Iterate a number of timesteps
	virtual bool IterateTS(unsigned int iterTS);

	virtual unsigned int GetNumberOfTimesteps() {return numTS;};

	inline virtual FDTD_FLOAT GetVolt( unsigned int n, unsigned int x, unsigned int y, unsigned int z ) const { return volt_[n][x][y][z/4].f[z%4]; }
	inline virtual FDTD_FLOAT GetCurr( unsigned int n, unsigned int x, unsigned int y, unsigned int z ) const { return curr_[n][x][y][z/4].f[z%4]; }

protected:
	Engine_sse(const Operator_sse* op);
	const Operator_sse* Op;

	virtual void UpdateVoltages();
	virtual void ApplyVoltageExcite();
	virtual void UpdateCurrents();
	virtual void ApplyCurrentExcite();

	unsigned int numLines[3];

	f4vector**** volt_;
	f4vector**** curr_;
	unsigned int numTS;
};

#endif // ENGINE_SSE_H
