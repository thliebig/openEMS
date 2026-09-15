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

#ifndef ENGINE_EXT_LORENTZMATERIAL_H
#define ENGINE_EXT_LORENTZMATERIAL_H

#include "engine_ext_dispersive.h"

class Operator_Ext_LorentzMaterial;

class Engine_Ext_LorentzMaterial : public Engine_Ext_Dispersive
{
public:
	Engine_Ext_LorentzMaterial(Operator_Ext_LorentzMaterial* op_ext_lorentz);
	virtual ~Engine_Ext_LorentzMaterial();

	virtual void DoPreVoltageUpdates();

	virtual void DoPreCurrentUpdates();

	//! The Metal engine can run the plain volt-ADE scheme (e.g. the conducting-sheet
	//! model) on the GPU. Models that need extra Lorentz flux states or ADE
	//! currents keep the CPU path.
	bool MetalADEOffloadSupported() const;

	int MetalADEOrder() const { return m_Order; }
	bool MetalADEVoltOn(int order) const;
	unsigned int MetalADECount(int order) const;
	const unsigned int* MetalADEPos(int order, int dir) const;
	const FDTD_FLOAT* MetalADEVoltInt(int order, int dir) const;
	const FDTD_FLOAT* MetalADEVoltExt(int order, int dir) const;

protected:
	template <typename EngType>
	void DoPreVoltageUpdatesImpl(EngType* eng);

	template <typename EngType>
	void DoPreCurrentUpdatesImpl(EngType* eng);

	Operator_Ext_LorentzMaterial* m_Op_Ext_Lor;

	//! ADE Lorentz voltages
	// Array setup: volt_Lor_ADE[N_order][direction][mesh_pos]
	FDTD_FLOAT ***volt_Lor_ADE;

	//! ADE Lorentz currents
	// Array setup: curr_Lor_ADE[N_order][direction][mesh_pos]
	FDTD_FLOAT ***curr_Lor_ADE;

};

#endif // ENGINE_EXT_LORENTZMATERIAL_H
