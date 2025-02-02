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

#ifndef ENGINE_EXT_ABSORBING_BC_H
#define ENGINE_EXT_ABSORBING_BC_H

#include "engine_extension.h"
#include "FDTD/engine.h"
#include "FDTD/operator.h"
#include "engine_extension_dispatcher.h"

class Operator_Ext_Absorbing_BC;

class Engine_Ext_Absorbing_BC : public Engine_Extension
{
public:
	Engine_Ext_Absorbing_BC(Operator_Ext_Absorbing_BC* op_ext);
	virtual ~Engine_Ext_Absorbing_BC();

	virtual void SetNumberOfThreads(int nrThread);

	virtual void DoPreVoltageUpdates() {Engine_Ext_Absorbing_BC::DoPreVoltageUpdates(0);}
	virtual void DoPreVoltageUpdates(int threadID);
	virtual void DoPostVoltageUpdates() {Engine_Ext_Absorbing_BC::DoPostVoltageUpdates(0);}
	virtual void DoPostVoltageUpdates(int threadID);
	virtual void Apply2Voltages() {Engine_Ext_Absorbing_BC::Apply2Voltages(0);}
	virtual void Apply2Voltages(int threadID);

	virtual void DoPreCurrentUpdates() {Engine_Ext_Absorbing_BC::DoPreCurrentUpdates(0);}
	virtual void DoPreCurrentUpdates(int threadID);
	virtual void DoPostCurrentUpdates() {Engine_Ext_Absorbing_BC::DoPostCurrentUpdates(0);}
	virtual void DoPostCurrentUpdates(int threadID);
	virtual void Apply2Current() {Engine_Ext_Absorbing_BC::Apply2Current(0);}
	virtual void Apply2Current(int threadID);

protected:
	Operator_Ext_Absorbing_BC* m_Op_ABC;

	template <typename EngType>
	void DoPreVoltageUpdatesImpl(EngType* eng, int threadID);

	template <typename EngType>
	void DoPostVoltageUpdatesImpl(EngType* eng, int threadID);

	template <typename EngType>
	void Apply2VoltagesImpl(EngType* eng, int threadID);

	template <typename EngType>
	void DoPreCurrentUpdatesImpl(EngType* eng, int threadID);

	template <typename EngType>
	void DoPostCurrentUpdatesImpl(EngType* eng, int threadID);

	template <typename EngType>
	void Apply2CurrentImpl(EngType* eng, int threadID);

	inline bool IsActive() {if (m_Eng->GetNumberOfTimesteps() < m_start_TS) return false; return true;}
	unsigned int m_start_TS;

	unsigned int 	m_numPrims;
	unsigned int	*m_numCells;
	unsigned int	**m_posStart;
	unsigned int	**m_posStop;
	unsigned int	*m_pos_ny0_I;
	unsigned int	*m_pos_ny0_shift_V;
	unsigned int	*m_pos_ny0_shift_I;
	unsigned int	**m_dir;

	vector<unsigned int>	v_primsPerThread;
	vector<unsigned int>	v_threadStartPrim;

	// Coefficients predefined in the operator
	FDTD_FLOAT		***m_K1;
	FDTD_FLOAT		***m_K2;

	// Containers for E-fields
	FDTD_FLOAT		**m_V_ny1; // n+1 direction
	FDTD_FLOAT		**m_V_ny2; // n+2 direction

	// Containers for H-fields calculated with the ABC
	FDTD_FLOAT		**m_I_ny1; // n+1 direction
	FDTD_FLOAT		**m_I_ny2; // n+2 direction

	// Containers for H-fields calculated from the curl equations
	FDTD_FLOAT		**m_Ic_ny1; // n+1 direction
	FDTD_FLOAT		**m_Ic_ny2; // n+2 direction

};

#endif // ENGINE_EXT_ABSORBING_BC_H
