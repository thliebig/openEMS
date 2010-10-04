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

#ifndef ENGINE_EXT_UPML_H
#define ENGINE_EXT_UPML_H

#include "engine_extension.h"
#include "engine.h"
#include "operator.h"

class Operator_Ext_UPML;

class Engine_Ext_UPML : public Engine_Extension
{
public:
	Engine_Ext_UPML(Operator_Ext_UPML* op_ext);
	virtual ~Engine_Ext_UPML();

	virtual void DoPreVoltageUpdates();
	virtual void DoPostVoltageUpdates();

	virtual void DoPreCurrentUpdates();
	virtual void DoPostCurrentUpdates();

protected:
	Operator_Ext_UPML* m_Op_UPML;

	FDTD_FLOAT**** volt;
	FDTD_FLOAT**** curr;
	FDTD_FLOAT**** volt_flux;
	FDTD_FLOAT**** curr_flux;
};

#endif // ENGINE_EXT_UPML_H
