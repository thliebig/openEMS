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

#ifndef OPERATOR_CYLINDER_H
#define OPERATOR_CYLINDER_H

#include "operator.h"

class Operator_Cylinder : public Operator
{
public:
	static Operator_Cylinder* New();
	virtual ~Operator_Cylinder();
	
	virtual bool SetGeometryCSX(ContinuousStructure* geo);

	virtual void ApplyElectricBC(bool* dirs);
	virtual void ApplyMagneticBC(bool* dirs);

	virtual void Reset();

protected:
	Operator_Cylinder();
	virtual void Init();

	virtual inline void Calc_ECOperatorPos(int n, unsigned int* pos);

	bool CC_closedAlpha;
	bool CC_R0_included;

	virtual bool Calc_ECPos(int n, unsigned int* pos, double* inEC);
	virtual bool Calc_EffMatPos(int n, unsigned int* pos, double* inMat);
};

#endif // OPERATOR_CYLINDER_H
