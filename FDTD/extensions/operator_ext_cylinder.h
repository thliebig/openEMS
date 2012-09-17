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

#ifndef OPERATOR_EXT_CYLINDER_H
#define OPERATOR_EXT_CYLINDER_H

#include "operator_extension.h"
#include "FDTD/operator.h"

class Operator_Cylinder;

class Operator_Ext_Cylinder : public Operator_Extension
{
	friend class Engine_Ext_Cylinder;
	friend class Operator_Ext_LorentzMaterial;
public:
	Operator_Ext_Cylinder(Operator_Cylinder* op);
	~Operator_Ext_Cylinder();

	virtual bool BuildExtension();

	virtual Engine_Extension* CreateEngineExtention();

	virtual bool IsCylinderCoordsSave(bool closedAlpha, bool R0_included) const {UNUSED(closedAlpha); UNUSED(R0_included); return true;}
	virtual bool IsCylindricalMultiGridSave(bool child) const {UNUSED(child); return true;}

	// FIXME, this extension is not save or unknown to be save to use with MPI
	virtual bool IsMPISave() const {return false;}

	virtual std::string GetExtensionName() const {return std::string("Extension for the Cylinder-Coords Operator");}

	virtual void ShowStat(ostream &ostr) const;

protected:
	Operator_Cylinder* m_Op_Cyl;

	bool CC_closedAlpha;
	bool CC_R0_included;

	//special EC operator for R0
	FDTD_FLOAT* vv_R0; //calc new voltage from old voltage
	FDTD_FLOAT* vi_R0; //calc new voltage from old current

};

#endif // OPERATOR_EXT_CYLINDER_H
