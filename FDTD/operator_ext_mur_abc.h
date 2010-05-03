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

#ifndef OPERATOR_EXT_MUR_ABC_H
#define OPERATOR_EXT_MUR_ABC_H

#include "operator.h"
#include "operator_extension.h"

class Operator_Ext_Mur_ABC : public Operator_Extension
{
	friend class Engine_Ext_Mur_ABC;
public:
	Operator_Ext_Mur_ABC(Operator* op);
	~Operator_Ext_Mur_ABC();

	//! Define the direction of this ABC: ny=0,1,2 -> x,y,z and if at bottom_ny -> e.g. x=0 or x=end
	void SetDirection(int ny, bool top_ny);

	virtual bool BuildExtension();

	virtual Engine_Extension* CreateEngineExtention();

	virtual bool IsCylinderCoordsSave() {return false;} //not yet save to use with cylinder-coords

	virtual string GetExtensionName() {return string("Mur ABC extension");}

protected:
	int m_ny;
	int m_nyP,m_nyPP;
	unsigned int m_LineNr;
	int m_LineNr_Shift;

	unsigned int m_numLines[2];

	FDTD_FLOAT** m_Mur_Coeff_nyP;
	FDTD_FLOAT** m_Mur_Coeff_nyPP;
};

#endif // OPERATOR_EXT_MUR_ABC_H
