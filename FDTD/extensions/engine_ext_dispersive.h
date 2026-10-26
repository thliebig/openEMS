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

#ifndef ENGINE_EXT_DISPERSIVE_H
#define ENGINE_EXT_DISPERSIVE_H

#include "engine_extension.h"
#include "FDTD/engine.h"
#include "FDTD/operator.h"
#include "engine_extension_dispatcher.h"

class Operator_Ext_Dispersive;

class Engine_Ext_Dispersive : public Engine_Extension
{
public:
	Engine_Ext_Dispersive(Operator_Ext_Dispersive* op_ext_disp);
	virtual ~Engine_Ext_Dispersive();

	virtual void SetNumberOfThreads(int nrThread);

	virtual void Apply2Voltages() {Engine_Ext_Dispersive::Apply2Voltages(0);}
	virtual void Apply2Voltages(int threadID);
	virtual void Apply2Current() {Engine_Ext_Dispersive::Apply2Current(0);}
	virtual void Apply2Current(int threadID);

protected:
	template <typename EngType>
	void Apply2VoltagesImpl(EngType* eng, int threadID);

	template <typename EngType>
	void Apply2CurrentImpl(EngType* eng, int threadID);

	Operator_Ext_Dispersive* m_Op_Ext_Disp;

	//! Thread t updates the entries [m_start[o][t], m_start[o][t+1]) of order o.
	//! All orders are split at the same mesh positions, so that each field value
	//! is only written by one thread.
	std::vector<std::vector<size_t>> m_start;

	//! Dispersive order
	int m_Order;

	//! ADE currents
	// Array setup: curr_ADE[N_order][direction][mesh_pos]
	FDTD_FLOAT ***curr_ADE;

	//! ADE voltages
	// Array setup: volt_ADE[N_order][direction][mesh_pos]
	FDTD_FLOAT ***volt_ADE;
};

#endif // ENGINE_EXT_DISPERSIVE_H
