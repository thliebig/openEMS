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

#include "engine_interface_base.h"
#include "string"

Engine_Interface_Base::Engine_Interface_Base(Operator_Base* base_op)
{
	m_Op_Base = base_op;
	m_InterpolType = NO_INTERPOLATION;
}

std::string Engine_Interface_Base::GetInterpolationNameByType(InterpolationType mode)
{
	switch (mode)
	{
	case NO_INTERPOLATION:
		return std::string("None");
	case NODE_INTERPOLATE:
		return std::string("Node");
	case CELL_INTERPOLATE:
		return std::string("Cell");
	}
	return std::string();
}
