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

	virtual unsigned int GetNumberOfTimesteps() {return numTS;};

	//this access functions muss be overloaded by any new engine using a different storage model
	inline virtual FDTD_FLOAT GetVolt( unsigned int n, unsigned int x, unsigned int y, unsigned int z )	const { return f4_volt[n][x][y][z%numVectors].f[z/numVectors]; }
	inline virtual FDTD_FLOAT GetVolt( unsigned int n, const unsigned int pos[3] )						const { return f4_volt[n][pos[0]][pos[1]][pos[2]%numVectors].f[pos[2]/numVectors]; }
	inline virtual FDTD_FLOAT GetCurr( unsigned int n, unsigned int x, unsigned int y, unsigned int z )	const { return f4_curr[n][x][y][z%numVectors].f[z/numVectors]; }
	inline virtual FDTD_FLOAT GetCurr( unsigned int n, const unsigned int pos[3] )						const { return f4_curr[n][pos[0]][pos[1]][pos[2]%numVectors].f[pos[2]/numVectors]; }

	inline virtual void SetVolt( unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value)	{ f4_volt[n][x][y][z%numVectors].f[z/numVectors]=value; }
	inline virtual void SetVolt( unsigned int n, const unsigned int pos[3], FDTD_FLOAT value )						{ f4_volt[n][pos[0]][pos[1]][pos[2]%numVectors].f[pos[2]/numVectors]=value; }
	inline virtual void SetCurr( unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value)	{ f4_curr[n][x][y][z%numVectors].f[z/numVectors]=value; }
	inline virtual void SetCurr( unsigned int n, const unsigned int pos[3], FDTD_FLOAT value )						{ f4_curr[n][pos[0]][pos[1]][pos[2]%numVectors].f[pos[2]/numVectors]=value; }

protected:
	Engine_sse(const Operator_sse* op);
	const Operator_sse* Op;

	virtual void UpdateVoltages(unsigned int startX, unsigned int numX);
	virtual void UpdateCurrents(unsigned int startX, unsigned int numX);

	unsigned int numVectors;

public: //public access to the sse arrays for efficient extensions access... use careful...
	f4vector**** f4_volt;
	f4vector**** f4_curr;
};

#endif // ENGINE_SSE_H
