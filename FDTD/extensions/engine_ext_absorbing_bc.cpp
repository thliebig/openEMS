/*
*	Copyright (C) 2023-2025 Gadi Lahav (gadi@rfwithcare.com), Thorsten Liebig (Thorsten.Liebig@gmx.de)
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


/* This version of absorbing boundary conditions is based on this article:
 *
 * Betz, Vaughn Timothy, and R. Mittra. "Absorbing boundary conditions for the finite-difference time-domain analysis of guided-wave structures." Coordinated Science Laboratory Report no. UILU-ENG-93-2243 (1993).
 *
 * After some trial and error, it was discovered that the simplest and most efficient implementations
 * are:
 * 1. Mur first order boundary conditions
 * 2. First order Mur with "super-absorption".
 * Later I discovered that the latter is equivalent to the so-called "Surface impedance boundary
 * conditions" (SIBC).
 */


#include "engine_ext_absorbing_bc.h"
#include "operator_ext_absorbing_bc.h"

#include <algorithm>
#include <stdexcept>
#include "FDTD/engine.h"
#include "FDTD/engine_sse.h"
#include "tools/array_ops.h"
#include "tools/useful.h"
#include "operator_ext_excitation.h"

Engine_Ext_Absorbing_BC::Engine_Ext_Absorbing_BC(Operator_Ext_Absorbing_BC* op_ext) :
	Engine_Extension(op_ext),
	m_K1_nyP(op_ext->m_K1_nyP),
	m_K1_nyPP(op_ext->m_K1_nyPP),
	m_K2_nyP(op_ext->m_K2_nyP),
	m_K2_nyPP(op_ext->m_K2_nyPP)
{

	m_Op_ABC = op_ext;
	m_ABCtype = int(m_Op_ABC->m_ABCtype);

	for (unsigned int dimIdx = 0 ; dimIdx < 3 ; dimIdx++)
	{
		m_posStart[dimIdx] = m_Op_ABC->m_sheetX0[dimIdx];
		m_posStop[dimIdx] = m_Op_ABC->m_sheetX1[dimIdx];
	}

	m_ny	= m_Op_ABC->m_ny;
	m_nyP 	= m_Op_ABC->m_nyP;
	m_nyPP 	= m_Op_ABC->m_nyPP;

	m_numLines[0] = m_Op_ABC->m_numLines[0];
	m_numLines[1] = m_Op_ABC->m_numLines[1];

	bool normalSignPositive = m_Op_ABC->m_normalSignPositive;

	m_start_TS = 0;

	// Initialize shifted position for V
	m_pos_ny0_shift_V = m_posStart[m_ny] + (normalSignPositive  ? 1 : -1);

	// Initialize shifted position for I. Different for super-absorption
	m_pos_ny0_I = m_posStart[m_ny] + (normalSignPositive ? 0 : -1);
	m_pos_ny0_shift_I = m_posStart[m_ny] + (normalSignPositive ? 1 : -2);

	m_V_nyP.Init("volt_nyP",m_numLines);
	m_V_nyPP.Init("volt_nyPP",m_numLines);
	m_I_nyP.Init("curr_nyP",m_numLines);
	m_I_nyPP.Init("curr_nyPP",m_numLines);

	// One thread per boundary
	SetNumberOfThreads(1);
}

Engine_Ext_Absorbing_BC::~Engine_Ext_Absorbing_BC()
{
}

void Engine_Ext_Absorbing_BC::SetNumberOfThreads(int nrThread)
{
	Engine_Extension::SetNumberOfThreads(nrThread);

	// This command assigns the number of jobs (primitives) handled by each thread
	m_linesPerThread = AssignJobs2Threads(m_numLines[0],m_NrThreads,false);

	// Basically cumulative sum. Starting point of each thread.
	m_threadStartLine.resize(m_NrThreads,0);
	m_threadStartLine.at(0) = 0;
	for (size_t n = 1; n < m_linesPerThread.size(); ++n)
		m_threadStartLine.at(n) = m_threadStartLine.at(n - 1) + m_linesPerThread.at(n - 1);
}

// The first order Mur boundary condition is based on the following step:
// The field at time step n, location i, is u
// E(i,n + 1) = E(i + s,n) + K1*[E(i,n + 1) - E(i,n)]
// where s is the +-1 shift, depending on the direction, and
// K1 = (vp*Dt - Dx)/(vp*Dt + Dx)
//
// Super-absorption (the currents, MUR_1ST_SA only):
//
// 1. Re-iterate the Mur B.C.
// Hsa(i,n + 1) = H(i + s,n) + K1*[H(i,n + 1) - H(i,n)]
//
// 2. Update the H(i,n + 1) as such
// H(i,n + 1) = (K2*Hsa(i,n + 1) + Hc(i,n + 1))/(K2 + 1)
// Where Hsa(i,n + 1) is the field calculated by the boundary condition, and
// Hc(i,n + 1) is the field calculated by the FDTD step.
// and K2 = vp*Dt/Dx
//
// The cells [i0, i1) x [j0, j1) of the sheet (its own indices): hook 0 pre, 1 post, 2 apply, of the
// voltages or currents. The flat hooks run it on the lines of their thread, the box hooks on the part
// in the box.
template <typename EngType>
void Engine_Ext_Absorbing_BC::SheetImpl(EngType* eng, bool volt, int hook, unsigned int i0, unsigned int i1, unsigned int j0, unsigned int j1)
{
	unsigned int pos[] = {0,0,0};
	unsigned int pos_shift[] = {0,0,0};
	pos[m_ny] = volt ? m_posStart[m_ny] : m_pos_ny0_I;
	pos_shift[m_ny] = volt ? m_pos_ny0_shift_V : m_pos_ny0_shift_I;
	for (unsigned int i=i0; i<i1; ++i)
	{
		pos_shift[m_nyP] = pos[m_nyP] = m_posStart[m_nyP] + i;
		for (unsigned int j=j0; j<j1; ++j)
		{
			pos_shift[m_nyPP] = pos[m_nyPP] = m_posStart[m_nyPP] + j;
			if (volt)
			{
				if (hook==0)
				{
					// E(i + s,n) - K1*E(i,n)
					m_V_nyP (i,j) = eng->EngType::GetVolt(m_nyP ,pos_shift) - m_Op_ABC->m_K1_nyP (i,j) * eng->EngType::GetVolt(m_nyP ,pos);
					m_V_nyPP(i,j) = eng->EngType::GetVolt(m_nyPP,pos_shift) - m_Op_ABC->m_K1_nyPP(i,j) * eng->EngType::GetVolt(m_nyPP,pos);
				}
				else if (hook==1)
				{
					// E(i + s,n) - K1*E(i,n) + K1*E(i + s,n) = E(i + s,n) + [K1*E(i + s,n) - K1*E(i,n)]
					m_V_nyP (i,j) += m_K1_nyP (i,j) * eng->EngType::GetVolt(m_nyP ,pos_shift);
					m_V_nyPP(i,j) += m_K1_nyPP(i,j) * eng->EngType::GetVolt(m_nyPP,pos_shift);
				}
				else
				{
					// E(i,n + 1) = E(i + s,n) + [K1*E(i + s,n) - K1*E(i,n)]
					eng->EngType::SetVolt(m_nyP ,pos, m_V_nyP (i,j));
					eng->EngType::SetVolt(m_nyPP,pos, m_V_nyPP(i,j));
				}
			}
			else
			{
				if (hook==0)
				{
					// H(i + s,n) - K1*H(i,n)
					m_I_nyP (i,j) = eng->EngType::GetCurr(m_nyP ,pos_shift) - m_K1_nyP (i,j)*eng->EngType::GetCurr(m_nyP ,pos);
					m_I_nyPP(i,j) = eng->EngType::GetCurr(m_nyPP,pos_shift) - m_K1_nyPP(i,j)*eng->EngType::GetCurr(m_nyPP,pos);
				}
				else if (hook==1)
				{
					// H(i + s,n) + K1*[H(i,n + 1) - H(i,n)]
					m_I_nyP (i,j) += m_K1_nyP (i,j)*eng->EngType::GetCurr(m_nyP ,pos_shift);
					m_I_nyPP(i,j) += m_K1_nyPP(i,j)*eng->EngType::GetCurr(m_nyPP,pos_shift);
				}
				else
				{
					// H(i + s,n) = (Hsa*K2 + Hc)/(1 + K2)
					eng->EngType::SetCurr(m_nyP ,pos, (m_I_nyP (i,j)*m_K2_nyP (i,j) + eng->EngType::GetCurr(m_nyP ,pos))/(m_K2_nyP (i,j) + 1.0f));
					eng->EngType::SetCurr(m_nyPP,pos, (m_I_nyPP(i,j)*m_K2_nyPP(i,j) + eng->EngType::GetCurr(m_nyPP,pos))/(m_K2_nyPP(i,j) + 1.0f));
				}
			}
		}
	}
}

template <typename EngType>
void Engine_Ext_Absorbing_BC::FlatImpl(EngType* eng, bool volt, int hook, int threadID)
{
	if (IsActive()==false) return;
	if (m_Eng==NULL) return;
	if (threadID >= m_NrThreads)
		return;
	if (!volt && ((Operator_Ext_Absorbing_BC::ABCtype)(m_ABCtype) != Operator_Ext_Absorbing_BC::MUR_1ST_SA))
		return;
	const unsigned int i0 = m_threadStartLine.at(threadID);
	unsigned int i1 = i0 + m_linesPerThread.at(threadID);
	if (volt)
		SheetImpl(eng, true, hook, i0, i1, 0, m_numLines[1]);
	else   // the currents: a line less (dual grid)
		SheetImpl(eng, false, hook, i0, std::min(i1, m_numLines[0] - 1), 0, m_numLines[1] - 1);
}

void Engine_Ext_Absorbing_BC::DoPreVoltageUpdates(int threadID)
{
	ENG_DISPATCH_ARGS(FlatImpl, true, 0, threadID);
}

void Engine_Ext_Absorbing_BC::DoPostVoltageUpdates(int threadID)
{
	ENG_DISPATCH_ARGS(FlatImpl, true, 1, threadID);
}

void Engine_Ext_Absorbing_BC::Apply2Voltages(int threadID)
{
	ENG_DISPATCH_ARGS(FlatImpl, true, 2, threadID);
}

void Engine_Ext_Absorbing_BC::DoPreCurrentUpdates(int threadID)
{
	ENG_DISPATCH_ARGS(FlatImpl, false, 0, threadID);
}

void Engine_Ext_Absorbing_BC::DoPostCurrentUpdates(int threadID)
{
	ENG_DISPATCH_ARGS(FlatImpl, false, 1, threadID);
}

void Engine_Ext_Absorbing_BC::Apply2Current(int threadID)
{
	ENG_DISPATCH_ARGS(FlatImpl, false, 2, threadID);
}

// The box hooks: the box must hold the line next to the sheet as well (the post hook
// reads it after the update of the box).
template <typename EngType>
void Engine_Ext_Absorbing_BC::BoxImpl(EngType* eng, bool volt, int hook, const Engine_Box* box, unsigned int ts)
{
	if (ts<m_start_TS)
		return;
	if (!volt && ((Operator_Ext_Absorbing_BC::ABCtype)(m_ABCtype) != Operator_Ext_Absorbing_BC::MUR_1ST_SA))
		return;
	const unsigned int line  = volt ? m_posStart[m_ny] : m_pos_ny0_I;
	const unsigned int shift = volt ? m_pos_ny0_shift_V : m_pos_ny0_shift_I;
	if ((line<box->start[m_ny]) || (line>=box->stop[m_ny]))
		return;
	if ((shift<box->start[m_ny]) || (shift>=box->stop[m_ny]))
		throw std::logic_error("Engine_Ext_Absorbing_BC: a box holds the sheet but not the line next to it");
	// the sheet's own indices (i, j) of the box; the currents have a line less (dual grid)
	const unsigned int n0 = volt ? m_numLines[0] : m_numLines[0]-1;
	const unsigned int n1 = volt ? m_numLines[1] : m_numLines[1]-1;
	auto clip = [](unsigned int lo, unsigned int hi, unsigned int start, unsigned int n, unsigned int& a, unsigned int& b) {
		a = (lo>start) ? lo-start : 0;
		b = (hi>start) ? std::min(hi-start, n) : 0;
	};
	unsigned int i0, i1, j0, j1;
	clip(box->start[m_nyP],  box->stop[m_nyP],  m_posStart[m_nyP],  n0, i0, i1);
	clip(box->start[m_nyPP], box->stop[m_nyPP], m_posStart[m_nyPP], n1, j0, j1);
	SheetImpl(eng, volt, hook, i0, i1, j0, j1);
}

void Engine_Ext_Absorbing_BC::DoPreVoltageUpdatesBox(const Engine_Box& box, unsigned int ts)
{
	ENG_DISPATCH_ARGS(BoxImpl, true, 0, &box, ts);
}

void Engine_Ext_Absorbing_BC::DoPostVoltageUpdatesBox(const Engine_Box& box, unsigned int ts)
{
	ENG_DISPATCH_ARGS(BoxImpl, true, 1, &box, ts);
}

void Engine_Ext_Absorbing_BC::Apply2VoltagesBox(const Engine_Box& box, unsigned int ts)
{
	ENG_DISPATCH_ARGS(BoxImpl, true, 2, &box, ts);
}

void Engine_Ext_Absorbing_BC::DoPreCurrentUpdatesBox(const Engine_Box& box, unsigned int ts)
{
	ENG_DISPATCH_ARGS(BoxImpl, false, 0, &box, ts);
}

void Engine_Ext_Absorbing_BC::DoPostCurrentUpdatesBox(const Engine_Box& box, unsigned int ts)
{
	ENG_DISPATCH_ARGS(BoxImpl, false, 1, &box, ts);
}

void Engine_Ext_Absorbing_BC::Apply2CurrentBox(const Engine_Box& box, unsigned int ts)
{
	ENG_DISPATCH_ARGS(BoxImpl, false, 2, &box, ts);
}
