/*
*	Copyright (C) 2011 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#include "engine_interface_sse_fdtd.h"
#include "operator_cylinder.h"

#ifndef ENGINE_INTERFACE_CYLINDRICAL_FDTD_H
#define ENGINE_INTERFACE_CYLINDRICAL_FDTD_H

class Engine_Interface_Cylindrical_FDTD : public Engine_Interface_SSE_FDTD
{
public:
	Engine_Interface_Cylindrical_FDTD(Operator_sse* op);
	virtual ~Engine_Interface_Cylindrical_FDTD();

	virtual double* GetHField(const unsigned int* pos, double* out) const;

protected:
	Operator_Cylinder* m_Op_Cyl;

	//! Internal method to get an interpolated field of a given type. (0: E, 1: J, 2: rotH)
	virtual double* GetRawInterpolatedField(const unsigned int* pos, double* out, int type) const;
};

#endif // ENGINE_INTERFACE_CYLINDRICAL_FDTD_H
