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

#ifndef OPERATOR_EXT_LORENTZMATERIAL_H
#define OPERATOR_EXT_LORENTZMATERIAL_H

#include "operator.h"
#include "operator_ext_dispersive.h"

class Operator_Ext_LorentzMaterial : public Operator_Ext_Dispersive
{
	friend class Engine_Ext_LorentzMaterial;
public:
	Operator_Ext_LorentzMaterial(Operator* op);
	virtual ~Operator_Ext_LorentzMaterial();

	virtual bool BuildExtension();

	virtual Engine_Extension* CreateEngineExtention();

	virtual bool IsCylinderCoordsSave() {return false;}

	virtual std::string GetExtensionName() {return std::string("Lorentz Dispersive Material Extension");}

protected:

	FDTD_FLOAT *v_int_ADE[3];
	FDTD_FLOAT *v_ext_ADE[3];
	FDTD_FLOAT *i_int_ADE[3];
	FDTD_FLOAT *i_ext_ADE[3];
};

#endif // OPERATOR_EXT_LORENTZMATERIAL_H
