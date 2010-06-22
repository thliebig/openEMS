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

//! define the base class for the cylindrical coordinate FDTD
#define __OP_CYLINDER_BASE_CLASS__ Operator_Multithread

#include "operator_multithread.h"

//! This class creates an operator for a cylindrical FDTD.
/*!
This class creates an operator for a cylindrical FDTD. No special engine is necessary,
all special cases e.g. a closed alpha mesh or an included r=0 case is treated by an operator/engine extension \sa operator_ext_cylinder.
*/
class Operator_Cylinder : public __OP_CYLINDER_BASE_CLASS__
{
public:
	static Operator_Cylinder* New(unsigned int numThreads = 0);
	virtual ~Operator_Cylinder();

	virtual bool SetGeometryCSX(ContinuousStructure* geo);

	virtual void ApplyElectricBC(bool* dirs);
	virtual void ApplyMagneticBC(bool* dirs);

	virtual unsigned int GetNumberOfLines(int ny) const;

	//! Get the name for the given direction: 0 -> rho, 1 -> alpha, 2 -> z
	virtual string GetDirName(int ny) const;

	//! Get the mesh delta times the grid delta for a 3D position, including radius corrected alpha-mesh width
	virtual double GetMeshDelta(int n, const int* pos, bool dualMesh=false) const;

	bool GetClosedAlpha() const {return CC_closedAlpha;}
	bool GetR0Included() const {return CC_R0_included;}

	virtual void AddExtension(Operator_Extension* op_ext);

	virtual bool Calc_ECPos(int n, unsigned int* pos, double* inEC);
	virtual bool Calc_EffMatPos(int n, unsigned int* pos, double* inMat);

protected:
	Operator_Cylinder();
	virtual void Init();
	virtual void InitOperator();
	virtual void Reset();

	//Calc timestep only internal use
	virtual double CalcTimestep();

	bool CC_closedAlpha;
	bool CC_R0_included;
};

#endif // OPERATOR_CYLINDER_H
