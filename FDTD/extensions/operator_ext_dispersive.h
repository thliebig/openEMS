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

#ifndef OPERATOR_EXT_DISPERSIVE_H
#define OPERATOR_EXT_DISPERSIVE_H

//#include "operator.h"
#include "operator_extension.h"
#include "vector"

//! Abstract base class for all dispersive material models, based on an ADE (additional differential equation)
class Operator_Ext_Dispersive : public Operator_Extension
{
	friend class Engine_Ext_Dispersive;
public:
	virtual ~Operator_Ext_Dispersive();

	virtual int GetDispersionOrder() {return m_Order;}

	virtual std::string GetExtensionName() const {return std::string("Dispersive Material Abstract Base class");}

	virtual void ShowStat(std::ostream &ostr) const;

protected:
	Operator_Ext_Dispersive(Operator* op);
	//! Copy constructor
	Operator_Ext_Dispersive(Operator* op, Operator_Ext_Dispersive* op_ext);

	//! Dispersive order
	int m_Order;

	//! Dispersive material count
	std::vector<unsigned int> m_LM_Count;
	//! Index with dispersive material
	// Array setup: m_LM_pos[N_order][direction][mesh_pos]
	unsigned int ***m_LM_pos;

	bool *m_curr_ADE_On;
	bool *m_volt_ADE_On;
};

#endif // OPERATOR_EXT_DISPERSIVE_H
