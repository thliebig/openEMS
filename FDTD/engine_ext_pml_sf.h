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

#ifndef ENGINE_EXT_PML_SF_H
#define ENGINE_EXT_PML_SF_H

#include "engine_extension.h"
#include "engine.h"
#include "operator.h"

class Operator_Ext_PML_SF;
class Operator_Ext_PML_SF_Plane;

//! Split field pml engine base class
class Engine_Ext_PML_SF : public Engine_Extension
{
public:
	virtual ~Engine_Ext_PML_SF();

	inline virtual FDTD_FLOAT GetVolt( unsigned int n, unsigned int x, unsigned int y, unsigned int z ) { return volt[0][n][x][y][z]+volt[1][n][x][y][z]; }
	inline virtual FDTD_FLOAT GetCurr( unsigned int n, unsigned int x, unsigned int y, unsigned int z ) { return curr[0][n][x][y][z]+curr[1][n][x][y][z]; }

protected:
	Engine_Ext_PML_SF(Operator_Ext_PML_SF* op_ext);

	void UpdateVoltages(unsigned int startX, unsigned int numX);
	void UpdateCurrents(unsigned int startX, unsigned int numX);

	Operator_Ext_PML_SF* m_Op_PML_SF;

	//split field voltages and currents
	FDTD_FLOAT**** volt[2];
	FDTD_FLOAT**** curr[2];
};

//! Split field pml plane engine class
class Engine_Ext_PML_SF_Plane : public Engine_Ext_PML_SF
{
public:
	Engine_Ext_PML_SF_Plane(Operator_Ext_PML_SF_Plane* op_ext);
	virtual ~Engine_Ext_PML_SF_Plane();

	virtual void Apply2Voltages();
	virtual void Apply2Current();

protected:
	Operator_Ext_PML_SF_Plane* m_Op_PML_SF_PL;
};

#endif // ENGINE_EXT_PML_SF_H
