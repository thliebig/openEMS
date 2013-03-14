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

#include "FDTD/operator.h"
#include "operator_ext_dispersive.h"

class Operator_Ext_LorentzMaterial : public Operator_Ext_Dispersive
{
	friend class Engine_Ext_LorentzMaterial;
public:
	Operator_Ext_LorentzMaterial(Operator* op);
	virtual ~Operator_Ext_LorentzMaterial();

	virtual Operator_Extension* Clone(Operator* op);

	virtual bool BuildExtension();

	virtual Engine_Extension* CreateEngineExtention();

	virtual bool IsCylinderCoordsSave(bool closedAlpha, bool R0_included) const {UNUSED(closedAlpha); UNUSED(R0_included); return true;}
	virtual bool IsCylindricalMultiGridSave(bool child) const {UNUSED(child); return true;}
	virtual bool IsMPISave() const {return true;}

	virtual string GetExtensionName() const {return string("Drude/Lorentz Dispersive Material Extension");}

	virtual void ShowStat(ostream &ostr) const;

protected:
	//! Copy constructor
	Operator_Ext_LorentzMaterial(Operator* op, Operator_Ext_LorentzMaterial* op_ext);

	//ADE update coefficients, array setup: coeff[N_order][direction][mesh_pos_index]
	FDTD_FLOAT ***v_int_ADE;
	FDTD_FLOAT ***v_ext_ADE;
	FDTD_FLOAT ***i_int_ADE;
	FDTD_FLOAT ***i_ext_ADE;

	bool *m_curr_Lor_ADE_On;
	bool *m_volt_Lor_ADE_On;

	FDTD_FLOAT ***v_Lor_ADE;
	FDTD_FLOAT ***i_Lor_ADE;
};

#endif // OPERATOR_EXT_LORENTZMATERIAL_H
