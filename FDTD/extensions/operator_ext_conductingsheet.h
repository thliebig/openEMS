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

#ifndef OPERATOR_EXT_CONDUCTINGSHEET_H
#define OPERATOR_EXT_CONDUCTINGSHEET_H

#include "FDTD/operator.h"
#include "operator_ext_lorentzmaterial.h"

/*!
  FDTD extension for a conducting sheet model as described in:
  Lauer, A.; Wolff, I.; , "A conducting sheet model for efficient wide band FDTD analysis of planar waveguides and circuits," Microwave Symposium Digest, 1999 IEEE MTT-S International , vol.4, no., pp.1589-1592 vol.4, 1999
  doi: 10.1109/MWSYM.1999.780262
  URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=780262&isnumber=16934
  */
class Operator_Ext_ConductingSheet : public Operator_Ext_LorentzMaterial
{
public:
	Operator_Ext_ConductingSheet(Operator* op, double f_max);

	virtual Operator_Extension* Clone(Operator* op);

	virtual bool BuildExtension();

	virtual bool IsCylinderCoordsSave(bool closedAlpha, bool R0_included) const {UNUSED(closedAlpha); UNUSED(R0_included); return true;}
	virtual bool IsCylindricalMultiGridSave(bool child) const {UNUSED(child); return true;}
	virtual bool IsMPISave() const {return true;}

	virtual string GetExtensionName() const {return string("Conducting Sheet Extension");}

protected:
	//! Copy constructor
	Operator_Ext_ConductingSheet(Operator* op, Operator_Ext_ConductingSheet* op_ext);
	double m_f_max;
};

#endif // OPERATOR_EXT_CONDUCTINGSHEET_H
