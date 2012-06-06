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

#include "engine_cylinder.h"

Engine_Cylinder::Engine_Cylinder(const Operator_Cylinder* op) : Engine_Multithread(op)
{
	m_Op_Cyl = op;
}

Engine_Cylinder::~Engine_Cylinder()
{

}

Engine_Cylinder* Engine_Cylinder::New(const Operator_Cylinder* op, unsigned int numThreads)
{
	cout << "Create FDTD engine (cylindrical mesh using sse compression + multithreading)" << endl;
	Engine_Cylinder* e = new Engine_Cylinder(op);
	e->setNumThreads(numThreads);
	e->Init();
	return e;
}

