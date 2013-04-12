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

#ifndef ENGINE_EXT_TFSF_H
#define ENGINE_EXT_TFSF_H

#include "engine_extension.h"

class Operator_Ext_TFSF;

class Engine_Ext_TFSF : public Engine_Extension
{
public:
	Engine_Ext_TFSF(Operator_Ext_TFSF* op_ext);
	virtual ~Engine_Ext_TFSF();

	virtual void DoPostVoltageUpdates();
	virtual void DoPostCurrentUpdates();

protected:
	Operator_Ext_TFSF* m_Op_TFSF;

	unsigned int* m_DelayLookup;
};

#endif // ENGINE_EXT_TFSF_H
