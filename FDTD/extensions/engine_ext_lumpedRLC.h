/*
*	Copyright (C) 2023 Gadi Lahav (gadi@rfwithcare.com)
*	Copyright (C) 2026 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#ifndef ENGINE_EXT_LUMPEDRLC_H
#define ENGINE_EXT_LUMPEDRLC_H

#include "engine_extension.h"
#include "engine_ext_cellindex.h"
#include "FDTD/engine.h"
#include "FDTD/operator.h"
#include "engine_extension_dispatcher.h"

class Operator_Ext_LumpedRLC;

class Engine_Ext_LumpedRLC : public Engine_Extension
{
	friend class Operator_Ext_LumpedRLC;
	friend class Operator;
	friend class ContinuousStructure;

public:

	Engine_Ext_LumpedRLC(Operator_Ext_LumpedRLC *op_ext_RLC);
	virtual ~Engine_Ext_LumpedRLC();

	virtual void DoPreVoltageUpdates();
	virtual void Apply2Voltages();

	virtual bool SupportsBoxes() const {return true;}
	virtual void PrepareBoxes(const unsigned int numLines[3]);
	virtual void DoPreVoltageUpdatesBox(const Engine_Box& box, unsigned int ts);
	virtual void Apply2VoltagesBox(const Engine_Box& box, unsigned int ts);

protected:
	//! the elements in \a box (all with NULL) at timestep \a ts
	void DoPreVoltageUpdatesImpl(unsigned int ts, const Engine_Box* box);
	template <typename EngType>
	void Apply2VoltagesImpl(EngType* eng, unsigned int ts, const Engine_Box* box);
	Engine_Ext_CellIndex m_Index;   //!< the elements by box

	Operator_Ext_LumpedRLC* m_Op_Ext_RLC;

	// Auxilliary containers

	// Array setup: volt_C_ADE[mesh_pos]
	FDTD_FLOAT *v_Il;		// Container for current on inductor- Parallel RLC

	// Containers for nodal vd and J of the last three timesteps, the one of timestep ts at [ts%3]
	FDTD_FLOAT **v_Vdn;
	FDTD_FLOAT **v_Jn;

};

#endif // ENGINE_EXT_LUMPEDRLC_H
