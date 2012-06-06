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

#ifndef ENGINE_CYLINDER_H
#define ENGINE_CYLINDER_H

#include "FDTD/engine.h"
#include "engine_extension.h"
#include "FDTD/operator_cylinder.h"

class Operator_Ext_Cylinder;
class Engine_sse;

class Engine_Ext_Cylinder : public Engine_Extension
{
public:
	Engine_Ext_Cylinder(Operator_Ext_Cylinder* op_ext);

	virtual void DoPostVoltageUpdates();

	virtual void DoPostCurrentUpdates();

	virtual void SetEngine(Engine* eng);

protected:
	Operator_Ext_Cylinder* cyl_Op;
	Engine_sse* m_Eng_SSE;

	unsigned int numLines[3];

	bool CC_closedAlpha;
	bool CC_R0_included;
};

#endif // ENGINE_CYLINDER_H
