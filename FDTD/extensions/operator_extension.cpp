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

#include "operator_extension.h"
#include "FDTD/operator.h"

Operator_Extension::Operator_Extension(Operator* op)
{
	m_Op = op;
}

Operator_Extension::~Operator_Extension()
{
}

Operator_Extension::Operator_Extension(Operator* op, Operator_Extension* op_ext)
{
	UNUSED(op_ext);
	m_Op = op;
}

void Operator_Extension::ShowStat(ostream &ostr) const
{
	ostr << "--- " << GetExtensionName() << " ---" << endl;
}

Operator_Extension* Operator_Extension::Clone(Operator* op)
{
	if (dynamic_cast<Operator_Extension*>(this)==NULL)
		return NULL;
	return new Operator_Extension(op, this);
}
