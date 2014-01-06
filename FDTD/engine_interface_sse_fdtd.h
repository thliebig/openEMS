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

#ifndef ENGINE_INTERFACE_SSE_FDTD_H
#define ENGINE_INTERFACE_SSE_FDTD_H

#include "engine_interface_fdtd.h"
#include "operator_sse.h"
#include "engine_sse.h"

class Engine_Interface_SSE_FDTD : public Engine_Interface_FDTD
{
public:
	Engine_Interface_SSE_FDTD(Operator_sse* op);
	virtual ~Engine_Interface_SSE_FDTD();

	virtual double CalcFastEnergy() const;

protected:
	Operator_sse* m_Op_SSE;
	Engine_sse* m_Eng_SSE;
};

#endif // ENGINE_INTERFACE_SSE_FDTD_H
