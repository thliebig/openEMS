/*
*	Copyright (C) 2010 Sebastian Held (Sebastian.Held@gmx.de)
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

#ifndef OPERATOR_SSE_H
#define OPERATOR_SSE_H

#include "operator.h"
#include "tools/array_ops.h"

class Operator_sse : public Operator
{
public:
	//! Create a new operator
	static Operator_sse* New();
	virtual ~Operator_sse();

	virtual int CalcECOperator();

	virtual void ApplyMagneticBC(bool* dirs);

protected:
	//! use New() for creating a new Operator
	Operator_sse();

	virtual void Init();
	virtual void Reset();
	virtual void InitOperator();

	// engine/post-proc needs access
public:
	f4vector**** vv_; //calc new voltage from old voltage
	f4vector**** vi_; //calc new voltage from old current
	f4vector**** iv_; //calc new current from old current
	f4vector**** ii_; //calc new current from old voltage
};

#endif // OPERATOR_SSE_H
