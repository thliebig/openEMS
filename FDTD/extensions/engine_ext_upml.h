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
#include "FDTD/engine.h"
#include "FDTD/operator.h"

class Operator_Ext_UPML;

class Engine_Ext_UPML : public Engine_Extension
{
public:
	Engine_Ext_UPML(Operator_Ext_UPML* op_ext);
	virtual ~Engine_Ext_UPML();

	virtual void SetNumberOfThreads(int nrThread);

	virtual void DoPreVoltageUpdates() {Engine_Ext_UPML::DoPreVoltageUpdates(0);};
	virtual void DoPreVoltageUpdates(int threadID);
	virtual void DoPostVoltageUpdates() {Engine_Ext_UPML::DoPostVoltageUpdates(0);};
	virtual void DoPostVoltageUpdates(int threadID);

	virtual void DoPreCurrentUpdates() {Engine_Ext_UPML::DoPreCurrentUpdates(0);};
	virtual void DoPreCurrentUpdates(int threadID);
	virtual void DoPostCurrentUpdates() {Engine_Ext_UPML::DoPostCurrentUpdates(0);};
	virtual void DoPostCurrentUpdates(int threadID);

protected:
	Operator_Ext_UPML* m_Op_UPML;

	vector<unsigned int> m_start;
	vector<unsigned int> m_numX;

	FDTD_FLOAT**** volt_flux;
	FDTD_FLOAT**** curr_flux;
};

#endif // ENGINE_EXT_UPML_H
