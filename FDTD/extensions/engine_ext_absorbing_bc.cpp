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

#include "engine_ext_absorbing_bc.h"
#include "operator_ext_absorbing_bc.h"
#include "FDTD/engine.h"
#include "FDTD/engine_sse.h"
#include "tools/array_ops.h"
#include "tools/useful.h"
#include "operator_ext_excitation.h"

#include "CSPropAbsorbingBC.h"

#define SAFE_DELETE(v) 	\
	if (v != NULL) 		\
	{ 					\
		delete[] v;		\
		v = NULL;		\
	}					\

Engine_Ext_Absorbing_BC::Engine_Ext_Absorbing_BC(Operator_Ext_Absorbing_BC* op_ext) : Engine_Extension(op_ext)
{

	m_Op_ABC = op_ext;

	m_K1 = m_Op_ABC->m_K1;
	m_K2 = m_Op_ABC->m_K2;

	m_numCells = m_Op_ABC->m_numCells;
	m_numPrims = m_Op_ABC->m_numPrims;

	m_start_TS = 0;


	// Initialize shifted positions
	m_pos_ny0_shift_V = new uint[m_Op_ABC->m_numPrims];

	// For magnetic fields, two containers are necessary
	m_pos_ny0_I = new uint[m_Op_ABC->m_numPrims];
	m_pos_ny0_shift_I = new uint[m_Op_ABC->m_numPrims];

	m_posStart = m_Op_ABC->m_sheetX0;
	m_posStop = m_Op_ABC->m_sheetX1;

	m_dir = m_Op_ABC->m_dir;

	// Initialize shifted location
	int normDir;
	uint cDir,cStartPos;
	for (uint primIdx = 0 ; primIdx < m_Op_ABC->m_numPrims ; primIdx++)
	{
		// Change according to normal direction
		normDir = m_Op_ABC->m_normDir[primIdx];
		cDir = abs(normDir) - 1;
		normDir = normDir >= 0 ? 1 : -1;

		cStartPos = m_posStart[primIdx][cDir];

		// This shift is also valid for magnetic field. The half cell
		// shift is always to the right (positive direction)
		m_pos_ny0_shift_V[primIdx] = cStartPos + normDir;

		// For SABC, also initialize magnetic field containers
		if ((CSPropAbsorbingBC::ABCtype)(m_Op_ABC->m_boundaryTypes[primIdx]) == CSPropAbsorbingBC::MUR_1ST_SA)
		{
			normDir = normDir >= 0 ? 0 : -1;
			if (m_pos_ny0_I)
				m_pos_ny0_I[primIdx] = cStartPos + normDir;

			// Convert normal direction indicator to fit with dual grid
			normDir = normDir >= 0 ? 1 : -2;

			m_pos_ny0_shift_I[primIdx] = cStartPos + normDir;
		}
		else
		{
			m_pos_ny0_I = 0;
			m_pos_ny0_shift_I = 0;
		}

	}

	m_V_ny1 = new FDTD_FLOAT*[m_Op_ABC->m_numPrims];
	m_V_ny2 = new FDTD_FLOAT*[m_Op_ABC->m_numPrims];

	m_I_ny1 = new FDTD_FLOAT*[m_Op_ABC->m_numPrims];
	m_I_ny2 = new FDTD_FLOAT*[m_Op_ABC->m_numPrims];

	m_Ic_ny1 = new FDTD_FLOAT*[m_Op_ABC->m_numPrims];
	m_Ic_ny2 = new FDTD_FLOAT*[m_Op_ABC->m_numPrims];

	for (uint primIdx = 0 ; primIdx < m_Op_ABC->m_numPrims ; primIdx++)
	{
		m_V_ny1[primIdx] = new FDTD_FLOAT[m_numCells[primIdx]];
		m_V_ny2[primIdx] = new FDTD_FLOAT[m_numCells[primIdx]];

		// Allocate if this is a SABC (maybe in the future something else
		if ((CSPropAbsorbingBC::ABCtype)(m_Op_ABC->m_boundaryTypes[primIdx]) == CSPropAbsorbingBC::MUR_1ST_SA)
		{
			m_I_ny1[primIdx] = new FDTD_FLOAT[m_numCells[primIdx]];
			m_I_ny2[primIdx] = new FDTD_FLOAT[m_numCells[primIdx]];

			m_Ic_ny1[primIdx] = new FDTD_FLOAT[m_numCells[primIdx]];
			m_Ic_ny2[primIdx] = new FDTD_FLOAT[m_numCells[primIdx]];
		}
		else
		{
			m_I_ny1[primIdx] = NULL;
			m_I_ny2[primIdx] = NULL;
			m_I_ny1[primIdx] = NULL;
			m_I_ny2[primIdx] = NULL;
		}

	}


//	//find if some excitation is on this mur-abc and find the max length of this excite, so that the abc can start after the excitation is done...
//	int maxDelay=-1;
//
//	Operator_Ext_Excitation* Exc_ext = m_Op_ABC->m_Op->GetExcitationExtension();
//
//	for (unsigned int n=0; n<Exc_ext->GetVoltCount(); ++n)
//	{
//		if ( ((Exc_ext->Volt_dir[n] == m_nyP) || (Exc_ext->Volt_dir[n]==m_nyPP)) && (Exc_ext->Volt_index[m_ny][n]==m_LineNr) )
//		{
//			if ((int)Exc_ext->Volt_delay[n]>maxDelay)
//				maxDelay = (int)Exc_ext->Volt_delay[n];
//		}
//	}
//	m_start_TS = 0;
//	if (maxDelay >= 0)
//	{
//		m_start_TS = maxDelay + m_Op_mur->m_Op->GetExcitationSignal()->GetLength() + 10; //give it some extra timesteps, for the excitation to travel at least one cell away
//		cerr << "Engine_Ext_Mur_ABC::Engine_Ext_Mur_ABC: Warning: Excitation inside the Mur-ABC #" <<  m_ny << "-" << (int)(m_LineNr>0) << " found!!!!  Mur-ABC will be switched on after excitation is done at " << m_start_TS << " timesteps!!! " << endl;
//	}

	// Da fuq?
	SetNumberOfThreads(m_Op_ABC->m_Nthreads);
}

Engine_Ext_Absorbing_BC::~Engine_Ext_Absorbing_BC()
{
	for (uint primIdx = 0 ; primIdx < m_Op_ABC->m_numPrims ; primIdx++)
	{
		delete[] m_V_ny1[primIdx];
		delete[] m_V_ny2[primIdx];

		SAFE_DELETE(m_I_ny1[primIdx]);
		SAFE_DELETE(m_I_ny2[primIdx]);

		SAFE_DELETE(m_Ic_ny1[primIdx]);
		SAFE_DELETE(m_Ic_ny2[primIdx]);

	}

	delete[] m_V_ny1;
	delete[] m_V_ny2;

	delete[] m_I_ny1;
	delete[] m_I_ny2;

	delete[] m_Ic_ny1;
	delete[] m_Ic_ny2;

	delete[] m_pos_ny0_shift_V;
	delete[] m_pos_ny0_shift_I;
	delete[] m_pos_ny0_I;
}

void Engine_Ext_Absorbing_BC::SetNumberOfThreads(int nrThread)
{
	Engine_Extension::SetNumberOfThreads(nrThread);

	// This command assigns the number of jobs (primitives) handled by each thread
	v_primsPerThread = AssignJobs2Threads(m_numPrims,m_NrThreads,false);

	// Basically cumsum. Starting point of each thread.
	v_threadStartPrim.resize(m_NrThreads,0);
	v_threadStartPrim.at(0) = 0;
	for (size_t threadIdx = 1; threadIdx < v_threadStartPrim.size(); threadIdx++)
		v_threadStartPrim.at(threadIdx) = v_threadStartPrim.at(threadIdx - 1) + v_primsPerThread.at(threadIdx - 1);
}

template <typename EngType>
void Engine_Ext_Absorbing_BC::DoPreVoltageUpdatesImpl(EngType* eng, int threadID)
{
	if (IsActive()==false) return;

	if (m_Eng==NULL) return;

	if (threadID >= m_NrThreads)
		return;

	uint pos[] = {0,0,0};
	uint pos0[] = {0,0,0};
	uint pos1[] = {0,0,0};
	uint pos_shift[] = {0,0,0};

	uint dir[] = {0,0,0};
	uint primIdx = 0;

	uint cellCtr;

	for (uint primCtr = 0 ; primCtr < v_primsPerThread.at(threadID) ; primCtr++)
	{
		// current primitive index for this thread
		primIdx = primCtr + v_threadStartPrim.at(threadID);
		for (uint dimIdx = 0 ; dimIdx < 3 ; dimIdx++)
		{
			// Position #1
			pos[dimIdx]			= m_posStart[primIdx][dimIdx];
			pos0[dimIdx]		= m_posStart[primIdx][dimIdx];
			pos1[dimIdx]		= m_posStop[primIdx][dimIdx];
			dir[dimIdx]			= m_dir[primIdx][dimIdx];
		}

		// Store shifted location in this container
		pos_shift[dir[0]] = m_pos_ny0_shift_V[primIdx];

		cellCtr = 0;
		for (pos[dir[1]] = pos0[dir[1]] ; pos[dir[1]] < pos1[dir[1]] ; pos[dir[1]]++)
		{
			pos_shift[dir[1]] = pos[dir[1]];
			for (pos[dir[2]] = pos0[dir[2]] ; pos[dir[2]] < pos1[dir[2]] ; pos[dir[2]]++)
			{
				pos_shift[dir[2]] = pos[dir[2]];
				m_V_ny1[primIdx][cellCtr] = eng->EngType::GetVolt(dir[1],pos_shift) - m_K1[primIdx][dir[0]][cellCtr] * eng->EngType::GetVolt(dir[1],pos);
				m_V_ny2[primIdx][cellCtr] = eng->EngType::GetVolt(dir[2],pos_shift) - m_K1[primIdx][dir[0]][cellCtr] * eng->EngType::GetVolt(dir[2],pos);

				cellCtr++;
			}
		}

	}
}

void Engine_Ext_Absorbing_BC::DoPreVoltageUpdates(int threadID)
{
	ENG_DISPATCH_ARGS(DoPreVoltageUpdatesImpl, threadID);
}

template <typename EngType>
void Engine_Ext_Absorbing_BC::DoPostVoltageUpdatesImpl(EngType* eng, int threadID)
{
	if (IsActive()==false) return;

	if (m_Eng==NULL) return;

	if (threadID >= m_NrThreads)
		return;

	uint pos[] = {0,0,0};
	uint pos0[] = {0,0,0};
	uint pos1[] = {0,0,0};
	uint pos_shift[] = {0,0,0};

	uint dir[] = {0,0,0};
	uint primIdx = 0;

	uint cellCtr;

	// Distribute jobs between threads
	for (uint primCtr = 0 ; primCtr < v_primsPerThread.at(threadID) ; primCtr++)
	{
		// current primitive index for this thread
		primIdx = primCtr + v_threadStartPrim.at(threadID);
		for (uint dimIdx = 0; dimIdx < 3 ; dimIdx++)
		{
			// Position #1
			pos[dimIdx]			= m_posStart[primIdx][dimIdx];
			pos0[dimIdx]		= m_posStart[primIdx][dimIdx];
			pos1[dimIdx]		= m_posStop[primIdx][dimIdx];
			dir[dimIdx]			= m_dir[primIdx][dimIdx];
		}

		// Store shifted location in this container
		pos_shift[dir[0]] = m_pos_ny0_shift_V[primIdx];

		cellCtr = 0;
		for (pos[dir[1]] = pos0[dir[1]] ; pos[dir[1]] < pos1[dir[1]] ; pos[dir[1]]++)
		{
			pos_shift[dir[1]] = pos[dir[1]];
			for (pos[dir[2]] = pos0[dir[2]] ; pos[dir[2]] < pos1[dir[2]] ; pos[dir[2]]++)
			{
				pos_shift[dir[2]] = pos[dir[2]];
				m_V_ny1[primIdx][cellCtr] += m_K1[primIdx][dir[0]][cellCtr] * eng->EngType::GetVolt(dir[1],pos_shift);
				m_V_ny2[primIdx][cellCtr] += m_K1[primIdx][dir[0]][cellCtr] * eng->EngType::GetVolt(dir[2],pos_shift);

				cellCtr++;
			}
		}

	}
}

void Engine_Ext_Absorbing_BC::DoPostVoltageUpdates(int threadID)
{
	ENG_DISPATCH_ARGS(DoPostVoltageUpdatesImpl, threadID);
}

template <typename EngType>
void Engine_Ext_Absorbing_BC::Apply2VoltagesImpl(EngType* eng, int threadID)
{
	// if (IsActive()==false) return;

	if (m_Eng==NULL) return;

	if (threadID >= m_NrThreads)
		return;

	uint pos[] = {0,0,0};
	uint pos0[] = {0,0,0};
	uint pos1[] = {0,0,0};

	uint dir[] = {0,0,0};
	uint primIdx = 0;

	uint cellCtr;

	for (uint primCtr = 0 ; primCtr < v_primsPerThread.at(threadID) ; primCtr++)
	{
		// current primitive index for this thread
		primIdx = primCtr + v_threadStartPrim.at(threadID);
		for (uint dimIdx = 0 ; dimIdx < 3 ; dimIdx++)
		{
			// Position #1
			pos[dimIdx]			= m_posStart[primIdx][dimIdx];
			pos0[dimIdx]		= m_posStart[primIdx][dimIdx];
			pos1[dimIdx]		= m_posStop[primIdx][dimIdx];
			dir[dimIdx]			= m_dir[primIdx][dimIdx];
		}

		// Store shifted location in this container
		cellCtr = 0;
		for (pos[dir[1]] = pos0[dir[1]] ; pos[dir[1]] < pos1[dir[1]] ; pos[dir[1]]++)
		{
			for (pos[dir[2]] = pos0[dir[2]] ; pos[dir[2]] < pos1[dir[2]] ; pos[dir[2]]++)
			{
				eng->EngType::SetVolt(dir[1],pos, m_V_ny1[primIdx][cellCtr]);
				eng->EngType::SetVolt(dir[2],pos, m_V_ny2[primIdx][cellCtr]);

				cellCtr++;
			}
		}

	}
}

void Engine_Ext_Absorbing_BC::Apply2Voltages(int threadID)
{
	ENG_DISPATCH_ARGS(Apply2VoltagesImpl, threadID);
}

template <typename EngType>
void Engine_Ext_Absorbing_BC::DoPreCurrentUpdatesImpl(EngType* eng, int threadID)
{
	// if (IsActive()==false) return;

	if (m_Eng==NULL) return;

	if (threadID >= m_NrThreads)
		return;

	uint pos[] = {0,0,0};
	uint pos0[] = {0,0,0};
	uint pos1[] = {0,0,0};
	uint pos_shift[] = {0,0,0};

	uint dir[] = {0,0,0};
	uint primIdx = 0;

	uint cellCtr;

	for (uint primCtr = 0 ; primCtr < v_primsPerThread.at(threadID) ; primCtr++)
	{

		// current primitive index for this thread
		primIdx = primCtr + v_threadStartPrim.at(threadID);

		// If this isn't the appropriate boundary type, move on to the next primitive
		if ((CSPropAbsorbingBC::ABCtype)(m_Op_ABC->m_boundaryTypes[primIdx]) != CSPropAbsorbingBC::MUR_1ST_SA)
			continue;

		for (uint dimIdx = 0 ; dimIdx < 3 ; dimIdx++)
		{
			// Position #1
			pos[dimIdx]			= m_posStart[primIdx][dimIdx];
			pos0[dimIdx]		= m_posStart[primIdx][dimIdx];
			pos1[dimIdx]		= m_posStop[primIdx][dimIdx] - 1;	// Deduct 1 because of the dual grid
			dir[dimIdx]			= m_dir[primIdx][dimIdx];
		}

		// Convert start position to H field start position
		pos[dir[0]] = m_pos_ny0_I[primIdx];
		// Store shifted location in this container
		pos_shift[dir[0]] = m_pos_ny0_shift_I[primIdx];
		//switch for different engine types to access faster inline engine functions

		cellCtr = 0;
		for (pos[dir[1]] = pos0[dir[1]] ; pos[dir[1]] < pos1[dir[1]] ; pos[dir[1]]++)
		{
			pos_shift[dir[1]] = pos[dir[1]];
			for (pos[dir[2]] = pos0[dir[2]] ; pos[dir[2]] < pos1[dir[2]] ; pos[dir[2]]++)
			{
				pos_shift[dir[2]] = pos[dir[2]];
				m_I_ny1[primIdx][cellCtr] = eng->EngType::GetCurr(dir[1],pos_shift) - m_K1[primIdx][dir[0]][cellCtr] * eng->EngType::GetCurr(dir[1],pos);
				m_I_ny2[primIdx][cellCtr] = eng->EngType::GetCurr(dir[2],pos_shift) - m_K1[primIdx][dir[0]][cellCtr] * eng->EngType::GetCurr(dir[2],pos);

				cellCtr++;
			}
		}

	}
}

void Engine_Ext_Absorbing_BC::DoPreCurrentUpdates(int threadID)
{
	ENG_DISPATCH_ARGS(DoPreCurrentUpdatesImpl, threadID);
}

template <typename EngType>
void Engine_Ext_Absorbing_BC::DoPostCurrentUpdatesImpl(EngType* eng, int threadID)
{
	// if (IsActive()==false) return;

	if (m_Eng==NULL) return;

	if (threadID >= m_NrThreads)
		return;

	uint pos[] = {0,0,0};
	uint pos0[] = {0,0,0};
	uint pos1[] = {0,0,0};
	uint pos_shift[] = {0,0,0};

	uint dir[] = {0,0,0};
	uint primIdx = 0;

	uint cellCtr;

	for (uint primCtr = 0 ; primCtr < v_primsPerThread.at(threadID) ; primCtr++)
	{
		// current primitive index for this thread
		primIdx = primCtr + v_threadStartPrim.at(threadID);

		// If this isn't the appropriate boundary type, move on to the next primitive
		if ((CSPropAbsorbingBC::ABCtype)(m_Op_ABC->m_boundaryTypes[primIdx]) != CSPropAbsorbingBC::MUR_1ST_SA)
			continue;

		for (uint dimIdx = 0 ; dimIdx < 3 ; dimIdx++)
		{
			// Position #1
			pos[dimIdx]			= m_posStart[primIdx][dimIdx];
			pos0[dimIdx]		= m_posStart[primIdx][dimIdx];
			pos1[dimIdx]		= m_posStop[primIdx][dimIdx] - 1; // Deduct 1 because of the H-field dual grid
			dir[dimIdx]			= m_dir[primIdx][dimIdx];
		}

		// Convert start position to H field start position
		pos[dir[0]] = m_pos_ny0_I[primIdx];
		// Store shifted location in this container
		pos_shift[dir[0]] = m_pos_ny0_shift_I[primIdx];


		cellCtr = 0;
		for (pos[dir[1]] = pos0[dir[1]] ; pos[dir[1]] < pos1[dir[1]] ; pos[dir[1]]++)
		{
			pos_shift[dir[1]] = pos[dir[1]];
			for (pos[dir[2]] = pos0[dir[2]] ; pos[dir[2]] < pos1[dir[2]] ; pos[dir[2]]++)
			{
				pos_shift[dir[2]] = pos[dir[2]];
				m_I_ny1[primIdx][cellCtr] += m_K1[primIdx][dir[0]][cellCtr] * eng->EngType::GetCurr(dir[1],pos_shift);
				m_I_ny2[primIdx][cellCtr] += m_K1[primIdx][dir[0]][cellCtr] * eng->EngType::GetCurr(dir[2],pos_shift);

				cellCtr++;
			}
		}

	}
}

void Engine_Ext_Absorbing_BC::DoPostCurrentUpdates(int threadID)
{
	ENG_DISPATCH_ARGS(DoPostCurrentUpdatesImpl, threadID);
}

template <typename EngType>
void Engine_Ext_Absorbing_BC::Apply2CurrentImpl(EngType* eng, int threadID)
{
	// if (IsActive()==false) return;

	if (m_Eng==NULL) return;

	if (threadID >= m_NrThreads)
		return;

	uint pos[] = {0,0,0};
	uint pos0[] = {0,0,0};
	uint pos1[] = {0,0,0};

	uint dir[] = {0,0,0};
	uint primIdx = 0;

	uint cellCtr;

	for (uint primCtr = 0 ; primCtr < v_primsPerThread.at(threadID) ; primCtr++)
	{
		// current primitive index for this thread
		primIdx = primCtr + v_threadStartPrim.at(threadID);

		// If this isn't the appropriate boundary type, move on to the next primitive
		if ((CSPropAbsorbingBC::ABCtype)(m_Op_ABC->m_boundaryTypes[primIdx]) != CSPropAbsorbingBC::MUR_1ST_SA)
			continue;

		for (uint dimIdx = 0 ; dimIdx < 3 ; dimIdx++)
		{
			// Position #1
			pos[dimIdx]			= m_posStart[primIdx][dimIdx];
			pos0[dimIdx]		= m_posStart[primIdx][dimIdx];
			pos1[dimIdx]		= m_posStop[primIdx][dimIdx] - 1;  // Deduct 1 because of the H-field dual grid
			dir[dimIdx]			= m_dir[primIdx][dimIdx];
		}

		// Convert start position to H field start position
		pos[dir[0]] = m_pos_ny0_I[primIdx];

		// First, load the field values calculated by the curl equations
		cellCtr = 0;
		for (pos[dir[1]] = pos0[dir[1]] ; pos[dir[1]] < pos1[dir[1]] ; pos[dir[1]]++)
			for (pos[dir[2]] = pos0[dir[2]] ; pos[dir[2]] < pos1[dir[2]] ; pos[dir[2]]++)
			{
				m_Ic_ny1[primIdx][cellCtr] = eng->EngType::GetCurr(dir[1],pos);
				m_Ic_ny2[primIdx][cellCtr] = eng->EngType::GetCurr(dir[2],pos);
				cellCtr++;
			}

	}
}

void Engine_Ext_Absorbing_BC::Apply2Current(int threadID)
{
	ENG_DISPATCH_ARGS(Apply2CurrentImpl, threadID);
}
