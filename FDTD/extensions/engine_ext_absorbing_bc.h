/*
*	Copyright (C) 2023-2025 Gadi Lahav (gadi@rfwithcare.com), 2010 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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
#include "tools/arraylib/array_ij.h"

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

	unsigned int	m_numLines[2];
	unsigned int	m_posStart[3];
	unsigned int	m_posStop[3];
	unsigned int	m_pos_ny0_I;
	unsigned int	m_pos_ny0_shift_V;
	unsigned int	m_pos_ny0_shift_I;

	int				m_ny,
					m_nyP,
					m_nyPP;

	int				m_ABCtype;					// Declared as integer here, because of the forward declaration formatting.

	std::vector<unsigned int>	m_threadStartLine;
	std::vector<unsigned int>	m_linesPerThread;

	ArrayLib::ArrayIJ<FDTD_FLOAT>&	m_K1_nyP;	// Copy of first set of coefficients, direction n + 1
	ArrayLib::ArrayIJ<FDTD_FLOAT>&	m_K1_nyPP;	// Copy of second set of coefficients, direction n + 2
	ArrayLib::ArrayIJ<FDTD_FLOAT>&	m_K2_nyP;	// Copy of first set of coefficients, direction n + 1
	ArrayLib::ArrayIJ<FDTD_FLOAT>&	m_K2_nyPP;	// Copy of second set of coefficients, direction n + 2

	ArrayLib::ArrayIJ<FDTD_FLOAT>	m_V_nyP;	// Storage for voltage, direction n + 1
	ArrayLib::ArrayIJ<FDTD_FLOAT>	m_V_nyPP;	// Storage for voltage, direction n + 2
	ArrayLib::ArrayIJ<FDTD_FLOAT>	m_I_nyP;	// Storage for currents, direction n + 1
	ArrayLib::ArrayIJ<FDTD_FLOAT>	m_I_nyPP;	// Storage for currents, direction n + 2


};

#endif // ENGINE_EXT_ABSORBING_BC_H
