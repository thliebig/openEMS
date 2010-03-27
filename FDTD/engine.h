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

#ifndef ENGINE_H
#define ENGINE_H

#include "operator.h"

class Engine
{
	friend class Processing;
	friend class ProcessVoltage;
	friend class ProcessCurrent;
	friend class ProcessFields;
	friend class ProcessFieldsTD;
public:
	static Engine* createEngine(Operator* op);
	virtual ~Engine();

	virtual void Init();
	virtual void Reset();

	//!Iterate a number of timesteps
	virtual bool IterateTS(unsigned int iterTS);

	virtual unsigned int GetNumberOfTimesteps() {return numTS;};

protected:
	Engine(Operator* op);
	Operator* Op;

	FDTD_FLOAT**** volt;
	FDTD_FLOAT**** curr;
	unsigned int numTS;
};

#endif // ENGINE_H
