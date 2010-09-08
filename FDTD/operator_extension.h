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

#include <string>

#include <iostream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;


class Operator;
class Engine_Extension;

//! Abstract base-class for all operator extensions
class Operator_Extension
{
	friend class Engine_Extension;
public:
	//! Create a clone of this extension, will return NULL if this is impossible
	/*!
		Create a clone of this extension, will return NULL if this is impossible (e.g. derived extension has no clone method and copy-constructor)...
		BuildExtension has to be called separatly!
	*/
	virtual Operator_Extension* Clone(Operator* op);

	virtual bool BuildExtension() {return true;}

	virtual Engine_Extension* CreateEngineExtention() {return 0;}

	virtual bool IsCylinderCoordsSave() const {return false;}

	virtual std::string GetExtensionName() const {return std::string("Abstract Operator Extension Base Class");}

	virtual void ShowStat(ostream &ostr) const;

protected:
	Operator_Extension(Operator* op);
	//! Copy constructor, returns NULL if extension cannot be copied...
	Operator_Extension(Operator* op, Operator_Extension* op_ext);
	Operator* m_Op;
};

#endif // OPERATOR_EXTENSION_H
