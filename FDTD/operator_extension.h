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

#ifndef OPERATOR_EXTENSION_H
#define OPERATOR_EXTENSION_H

class Operator;
class Engine_Extension;

//! Abstract base-class for all operator extensions
class Operator_Extension
{
public:
	virtual bool BuildExtension() {return true;}

	virtual Engine_Extension* CreateEngineExtention() {return 0;}

protected:
	Operator_Extension(Operator* op);
	Operator* m_Op;
};

#endif // OPERATOR_EXTENSION_H
