/*
*	Copyright (C) 2015 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#ifndef ENGINE_EXT_STEADYSTATE_H
#define ENGINE_EXT_STEADYSTATE_H

#include "engine_extension.h"
#include "FDTD/engine.h"
#include "FDTD/operator.h"

class Operator_Ext_SteadyState;
class Engine_Interface_FDTD;

class Engine_Ext_SteadyState : public Engine_Extension
{
public:
	Engine_Ext_SteadyState(Operator_Ext_SteadyState* op_ext);
	virtual ~Engine_Ext_SteadyState();

	virtual void Apply2Voltages();
	virtual void Apply2Current();

	void SetEngineInterface(Engine_Interface_FDTD* eng_if) {m_Eng_Interface=eng_if;}
	double GetLastDiff() {return m_last_max_diff;}

protected:
	Operator_Ext_SteadyState* m_Op_SS;
	double m_last_max_diff;
	vector<double*> m_E_records;
	vector<double*> m_H_records;

	double last_total_energy;
	Engine_Interface_FDTD* m_Eng_Interface;
};


#endif // ENGINE_EXT_STEADYSTATE_H
