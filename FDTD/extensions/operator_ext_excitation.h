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

#ifndef OPERATOR_EXT_EXCITATION_H
#define OPERATOR_EXT_EXCITATION_H

#include "operator_extension.h"
#include "FDTD/operator.h"

class Excitation;

class Operator_Ext_Excitation : public Operator_Extension
{
	friend class Engine_Ext_Excitation;
public:
	Operator_Ext_Excitation(Operator* op, Excitation* exc);
	~Operator_Ext_Excitation();

	virtual Operator_Extension* Clone(Operator* op) {UNUSED(op);return NULL;}

	virtual bool BuildExtension();

	virtual Engine_Extension* CreateEngineExtention();

	virtual bool IsCylinderCoordsSave() const {return true;}
	virtual bool IsCylindricalMultiGridSave(bool child) const {UNUSED(child); return true;}

	virtual string GetExtensionName() const {return string("Excitation Extension");}

	virtual void ShowStat(ostream &ostr) const;

protected:
	Operator_Ext_Excitation(Operator* op, Operator_Ext_Excitation* op_ext);

	Excitation* m_Exc;
};

#endif // OPERATOR_EXT_EXCITATION_H
