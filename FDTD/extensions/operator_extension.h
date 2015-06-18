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

#include "tools/global.h"

class Operator;
class Operator_Cylinder;
class Engine_Extension;

//! Abstract base-class for all operator extensions
class Operator_Extension
{
	friend class Engine_Extension;
public:
	virtual ~Operator_Extension();

	//! Create a clone of this extension, will return NULL if this is impossible
	/*!
		Create a clone of this extension, will return NULL if this is impossible (e.g. derived extension has no clone method and copy-constructor)...
		BuildExtension has to be called separatly!
	*/
	virtual Operator_Extension* Clone(Operator* op) {UNUSED(op); return NULL;}

	virtual bool BuildExtension() {return true;}

	virtual Engine_Extension* CreateEngineExtention() {return 0;}
	virtual Engine_Extension* GetEngineExtention() {return m_Eng_Ext;}

	//! The cylindrical operator will check whether the extension is save to use. Default is false. Derive this method to override.
	virtual bool IsCylinderCoordsSave(bool closedAlpha, bool R0_included) const {UNUSED(closedAlpha); UNUSED(R0_included); return false;}

	//! The cylindrical multi grid operator will check whether the extension is save to use. Default is false. Derive this method to override.
	virtual bool IsCylindricalMultiGridSave(bool child) const {UNUSED(child); return false;}

	//! The MPI operator (if enabled) will check whether the extension is compatible with MPI. Default is false. Derive this method to override.
	virtual bool IsMPISave() const {return false;}

	virtual std::string GetExtensionName() const {return std::string("Abstract Operator Extension Base Class");}

	virtual void ShowStat(std::ostream &ostr) const;

	virtual bool IsActive() const {return m_Active;}
	virtual void SetActive(bool active=true) {m_Active=active;}

	virtual void Init() {}
	virtual void Reset() {}

protected:
	Operator_Extension(Operator* op);
	//! Copy constructor
	Operator_Extension(Operator* op, Operator_Extension* op_ext);

	bool m_Active;

	//FDTD Operator
	Operator* m_Op;
	Engine_Extension* m_Eng_Ext;

	//Cylindrical FDTD Operator (not NULL if a cylindrical FDTD is used)
	Operator_Cylinder* m_Op_Cyl;
	bool m_CC_R0_included;
};

#endif // OPERATOR_EXTENSION_H
