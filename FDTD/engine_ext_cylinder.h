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

#include "engine.h"
#include "engine_extension.h"
#include "operator_cylinder.h"

class Operator_Ext_Cylinder;

class Engine_Ext_Cylinder : public Engine_Extension
{
public:
	Engine_Ext_Cylinder(Operator_Ext_Cylinder* op_ext);

	virtual void Apply2Voltages();

	virtual void Apply2Current();

protected:
	Operator_Ext_Cylinder* cyl_Op;

	unsigned int numLines[3];

	bool CC_closedAlpha;
	bool CC_R0_included;
};

#endif // ENGINE_CYLINDER_H
