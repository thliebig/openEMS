/*
*	Copyright (C) 2012 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#include "operator_cylinder.h"
#include "engine_multithread.h"

class Engine_Cylinder : public Engine_Multithread
{
public:
	static Engine_Cylinder* New(const Operator_Cylinder* op, unsigned int numThreads = 0);
	virtual ~Engine_Cylinder();

protected:
	Engine_Cylinder(const Operator_Cylinder* op);

	const Operator_Cylinder* m_Op_Cyl;
};

#endif // ENGINE_CYLINDER_H
