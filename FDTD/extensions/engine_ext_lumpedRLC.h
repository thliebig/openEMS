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

	//! Metal diamond offload descriptors. The operator owns the coefficients and
	//! positions; these accessors give Engine_Metal a compact, read-only view
	//! without exposing the internal arrays.
	bool MetalRLCOffloadSupported() const;
	unsigned int MetalRLCCount() const;
	int MetalRLCDir(unsigned int index) const;
	unsigned int MetalRLCPos(unsigned int axis, unsigned int index) const;
	FDTD_FLOAT MetalRLCDJdV(unsigned int index) const;
	FDTD_FLOAT MetalRLCAV(unsigned int index) const;
	FDTD_FLOAT MetalRLCAQ(unsigned int index) const;
	FDTD_FLOAT MetalRLCAJ(unsigned int index) const;
	FDTD_FLOAT MetalRLCVcd(unsigned int index) const;
	FDTD_FLOAT MetalRLCVvd(unsigned int index) const;
	//! Product i2v*ilv, the explicit parallel-inductor current coefficient.
	FDTD_FLOAT MetalRLCIlCoeff(unsigned int index) const;
	FDTD_FLOAT MetalRLCDtHalf() const;

protected:
	template <typename EngType>
	void Apply2VoltagesImpl(EngType* eng);

	Operator_Ext_LumpedRLC* m_Op_Ext_RLC;

	// Auxilliary containers

	// Array setup: volt_C_ADE[mesh_pos]
	FDTD_FLOAT *v_Il;		// Container for current on inductor- Parallel RLC
	FDTD_FLOAT *v_q;		// Container for the series charge q

	FDTD_FLOAT **v_Vdn;		// Container for nodal vd at [n],[n-1],[n-2]
	FDTD_FLOAT **v_Jn;		// Container for nodal J at [n],[n-1],[n-2]

};

#endif // ENGINE_EXT_LUMPEDRLC_H
