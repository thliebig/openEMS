/*
*	Copyright (C) 2010-2026 Thorsten Liebig (Thorsten.Liebig@gmx.de)
*	Copyright (C) 2026 Sean Mollet (sean@malmoset.com)
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

#include "engine_subgrid.h"
#include "engine_sse.h"
#include "engine_gpu.h"
#include "operator_cylindermultigrid.h"

#include <iostream>
#include <stdexcept>

// ----------------------------------------------------------------------------

Engine_SubGridCoupling* Engine_SubGridCoupling::New(Engine* base, Engine* sub, const Operator_CylinderMultiGrid* op)
{
	if (base->GetType()!=sub->GetType())
		throw std::logic_error("Engine_SubGridCoupling::New: the grid and its sub-grid run different engine types");
	switch (base->GetType())
	{
	case Engine::SSE:
		return new SubGrid_Coupling_SSE(static_cast<Engine_sse*>(base), static_cast<Engine_sse*>(sub), op);
	case Engine::BASIC:
		return new SubGrid_Coupling_Basic(base, sub, op);
	case Engine::GPU:
		return new SubGrid_Coupling_GPU(static_cast<Engine_GPU*>(base), static_cast<Engine_GPU*>(sub), op);
	default:
		throw std::logic_error("Engine_SubGridCoupling::New: no coupling for this engine type");
	}
}

// ----------------------------------------------------------------------------

SubGrid_Coupling_Basic::SubGrid_Coupling_Basic(Engine* base, Engine* sub, const Operator_CylinderMultiGrid* op)
{
	m_Base = base;
	m_Sub = sub;
	m_Interpol = op->GetInterpolation();
	for (int n=0; n<3; ++n)
		m_numLines[n] = op->GetNumberOfLines(n, true);
}

void SubGrid_Coupling_Basic::SyncVoltages()
{
	SyncVoltages(*m_Base->volt_ptr, *m_Sub->volt_ptr, m_Interpol, m_numLines);
}

void SubGrid_Coupling_Basic::SyncCurrents()
{
	if (m_Interpol.split_pos>=2)
		Interpolate(*m_Base->curr_ptr, *m_Sub->curr_ptr, m_Interpol, true, m_Interpol.split_pos-2, m_numLines);
}

void SubGrid_Coupling_Basic::Interpolate()
{
	for (unsigned int r=0; r+1<m_Interpol.split_pos; ++r)
		Interpolate(*m_Base->volt_ptr, *m_Sub->volt_ptr, m_Interpol, false, r, m_numLines);
	for (unsigned int r=0; r+2<m_Interpol.split_pos; ++r)
		Interpolate(*m_Base->curr_ptr, *m_Sub->curr_ptr, m_Interpol, true, r, m_numLines);
}

// see SubGrid_Coupling_SSE::SyncVoltages()
void SubGrid_Coupling_Basic::SyncVoltages(const ArrayLib::ArrayNIJK<FDTD_FLOAT>& base, ArrayLib::ArrayNIJK<FDTD_FLOAT>& sub, const MultiGridInterpolation& ip, const unsigned int numLines[3])
{
	if (ip.split_pos<1)
		return;
	const unsigned int r = ip.split_pos-1;
	for (unsigned int a=0; a<numLines[1]-1; a+=2)
	{
		const unsigned int h = a/2;
		for (unsigned int z=0; z<numLines[2]; ++z)
		{
			sub(0, r, h, z) = 0;
			sub(2, r, h, z) = base(2, r, a, z);
			sub(1, r, h, z) = base(1, r, a, z) + base(1, r, a+1, z);
		}
	}
}

// see SubGrid_Coupling_SSE::InterpolVoltChild2Base() and InterpolCurrChild2Base()
void SubGrid_Coupling_Basic::Interpolate(ArrayLib::ArrayNIJK<FDTD_FLOAT>& base, const ArrayLib::ArrayNIJK<FDTD_FLOAT>& sub, const MultiGridInterpolation& ip, bool currents, unsigned int r, const unsigned int numLines[3])
{
	const std::vector<unsigned int>* pos_2p  = currents ? ip.pos_i_2p  : ip.pos_v_2p;
	const std::vector<unsigned int>* pos_2pp = currents ? ip.pos_i_2pp : ip.pos_v_2pp;
	const std::vector<float>* w_2p  = currents ? ip.w_i_2p  : ip.w_v_2p;
	const std::vector<float>* w_2pp = currents ? ip.w_i_2pp : ip.w_v_2pp;
	for (unsigned int a=0; a<numLines[1]; ++a)
	{
		for (unsigned int z=0; z<numLines[2]; ++z)
		{
			//r and z direction
			base(0, r, a, z) = w_2p[0][a]*sub(0, r, pos_2p[0][a], z) + w_2pp[0][a]*sub(0, r, pos_2pp[0][a], z);
			base(2, r, a, z) = w_2p[0][a]*sub(2, r, pos_2p[0][a], z) + w_2pp[0][a]*sub(2, r, pos_2pp[0][a], z);
			//alpha direction
			base(1, r, a, z) = w_2p[1][a]*sub(1, r, pos_2p[1][a], z) + w_2pp[1][a]*sub(1, r, pos_2pp[1][a], z);
		}
	}
}

// ----------------------------------------------------------------------------

SubGrid_Coupling_GPU::SubGrid_Coupling_GPU(Engine_GPU* base, Engine_GPU* sub, const Operator_CylinderMultiGrid* op)
	: SubGrid_Coupling_Basic(base, sub, op)
{
	m_Link = NULL;
	// all levels on the device with a device coupling, or all in the host fallback
	if (!base->FieldsOnHost() && !sub->FieldsOnHost())
		m_Link = base->GetBackend()->CreateMultiGridLink(sub->GetBackend(), m_Interpol);
	if (m_Link==NULL)
	{
		if (!base->FieldsOnHost() || !sub->FieldsOnHost())
			std::cerr << "SubGrid_Coupling_GPU: Warning: multi-grid level " << op->GetMultiGridLevel()
			          << " and its sub-grid run in the host fallback" << std::endl;
		base->ForceHostFallback();
		sub->ForceHostFallback();
	}
}

SubGrid_Coupling_GPU::~SubGrid_Coupling_GPU()
{
	delete m_Link;
}

void SubGrid_Coupling_GPU::DropLink()
{
	delete m_Link;
	m_Link = NULL;
}

void SubGrid_Coupling_GPU::SyncVoltages()
{
	if (m_Link)
		m_Link->SyncVoltages();
	else
		SubGrid_Coupling_Basic::SyncVoltages();
}

void SubGrid_Coupling_GPU::SyncCurrents()
{
	if (m_Link)
		m_Link->SyncCurrents();
	else
		SubGrid_Coupling_Basic::SyncCurrents();
}

void SubGrid_Coupling_GPU::Interpolate()
{
	if (m_Link)
		m_Link->InterpolateToBase();
	else
		SubGrid_Coupling_Basic::Interpolate();
}

// ----------------------------------------------------------------------------

SubGrid_Coupling_SSE::SubGrid_Coupling_SSE(Engine_sse* base, Engine_sse* sub, const Operator_CylinderMultiGrid* op)
{
	m_Base = base;
	m_Sub = sub;
	m_Interpol = op->GetInterpolation();
	for (int n=0; n<3; ++n)
		m_numLines[n] = op->GetNumberOfLines(n, true);
	m_numVectors = base->GetNumberOfVectors();

	// the weights in all four lanes, for the f4vector fields
	const std::vector<float>* w[4] = {m_Interpol.w_v_2p, m_Interpol.w_v_2pp, m_Interpol.w_i_2p, m_Interpol.w_i_2pp};
	f4weights* f4w[4] = {m_w_v_2p, m_w_v_2pp, m_w_i_2p, m_w_i_2pp};
	for (int k=0; k<4; ++k)
		for (int n=0; n<2; ++n)
		{
			f4w[k][n].resize(w[k][n].size());
			for (size_t a=0; a<w[k][n].size(); ++a)
				for (int l=0; l<4; ++l)
					f4w[k][n][a].f[l] = w[k][n][a];
		}
}

void SubGrid_Coupling_SSE::SyncVoltages()
{

	//interpolate voltages from base engine to child engine...
	unsigned int pos[3];
	pos[0] = m_Interpol.split_pos-1;
	unsigned int pos1_half = 0;
	f4vector v_null;
	v_null.f[0] = 0;
	v_null.f[1] = 0;
	v_null.f[2] = 0;
	v_null.f[3] = 0;
	for (pos[1]=0; pos[1]<m_numLines[1]-1; pos[1]+=2)
	{
		pos1_half = pos[1]/2;
		for (pos[2]=0; pos[2]<m_numVectors; ++pos[2])
		{
			ArrayLib::ArrayNIJK<f4vector>& mgEng_f4_volt = *m_Base->f4_volt_ptr;
			ArrayLib::ArrayNIJK<f4vector>& innerEng_f4_volt = *m_Sub->f4_volt_ptr;

			//r - direczion
			innerEng_f4_volt(0, pos[0], pos1_half, pos[2]).v = v_null.v;

			//z - direction
			innerEng_f4_volt(2, pos[0], pos1_half, pos[2]).v = mgEng_f4_volt(2, pos[0], pos[1], pos[2]).v;

			//alpha - direction
			innerEng_f4_volt(1, pos[0], pos1_half, pos[2]).v  = mgEng_f4_volt(1, pos[0], pos[1]  , pos[2]).v;
			innerEng_f4_volt(1, pos[0], pos1_half, pos[2]).v += mgEng_f4_volt(1, pos[0], pos[1]+1, pos[2]).v;
		}
	}
}

void SubGrid_Coupling_SSE::SyncCurrents()
{
	InterpolCurrChild2Base(m_Interpol.split_pos-2);
}

void SubGrid_Coupling_SSE::Interpolate()
{
	//interpolate child data to base mesh...
	for (unsigned int n=0; n<m_Interpol.split_pos-1; ++n)
		InterpolVoltChild2Base(n);
	for (unsigned int n=0; n<m_Interpol.split_pos-2; ++n)
		InterpolCurrChild2Base(n);
}

void SubGrid_Coupling_SSE::InterpolVoltChild2Base(unsigned int rPos)
{
	//interpolate voltages from child engine to the base engine...
	ArrayLib::ArrayNIJK<f4vector>& f4_volt = *m_Base->f4_volt_ptr;
	ArrayLib::ArrayNIJK<f4vector>& innerEngine_f4_volt = *m_Sub->f4_volt_ptr;

	unsigned int pos[3];
	pos[0] = rPos;
	for (pos[1]=0; pos[1]<m_numLines[1]; ++pos[1])
	{
		for (pos[2]=0; pos[2]<m_numVectors; ++pos[2])
		{
			//r - direction
			f4_volt(0, pos[0], pos[1], pos[2]).v  = m_w_v_2p[0][pos[1]].v * innerEngine_f4_volt(0, pos[0], m_Interpol.pos_v_2p[0][pos[1]], pos[2]).v
													+ m_w_v_2pp[0][pos[1]].v * innerEngine_f4_volt(0, pos[0], m_Interpol.pos_v_2pp[0][pos[1]], pos[2]).v;

			//z - direction
			f4_volt(2, pos[0], pos[1], pos[2]).v  = m_w_v_2p[0][pos[1]].v * innerEngine_f4_volt(2, pos[0], m_Interpol.pos_v_2p[0][pos[1]], pos[2]).v
													+ m_w_v_2pp[0][pos[1]].v * innerEngine_f4_volt(2, pos[0], m_Interpol.pos_v_2pp[0][pos[1]], pos[2]).v;

			//alpha - direction
			f4_volt(1, pos[0], pos[1], pos[2]).v  = m_w_v_2p[1][pos[1]].v * innerEngine_f4_volt(1, pos[0], m_Interpol.pos_v_2p[1][pos[1]], pos[2]).v
													+ m_w_v_2pp[1][pos[1]].v * innerEngine_f4_volt(1, pos[0], m_Interpol.pos_v_2pp[1][pos[1]], pos[2]).v;
		}
	}
}

void SubGrid_Coupling_SSE::InterpolCurrChild2Base(unsigned int rPos)
{
	//interpolate voltages from child engine to the base engine...
	ArrayLib::ArrayNIJK<f4vector>& f4_curr = *m_Base->f4_curr_ptr;
	ArrayLib::ArrayNIJK<f4vector>& innerEngine_f4_curr = *m_Sub->f4_curr_ptr;

	unsigned int pos[3];
	pos[0] = rPos;
	for (pos[1]=0; pos[1]<m_numLines[1]; ++pos[1])
	{
		for (pos[2]=0; pos[2]<m_numVectors; ++pos[2])
		{
			//r - direction
			f4_curr(0, pos[0], pos[1], pos[2]).v  = m_w_i_2p[0][pos[1]].v * innerEngine_f4_curr(0, pos[0], m_Interpol.pos_i_2p[0][pos[1]], pos[2]).v
													+ m_w_i_2pp[0][pos[1]].v * innerEngine_f4_curr(0, pos[0], m_Interpol.pos_i_2pp[0][pos[1]], pos[2]).v;

			//z - direction
			f4_curr(2, pos[0], pos[1], pos[2]).v  = m_w_i_2p[0][pos[1]].v * innerEngine_f4_curr(2, pos[0], m_Interpol.pos_i_2p[0][pos[1]], pos[2]).v
													+ m_w_i_2pp[0][pos[1]].v * innerEngine_f4_curr(2, pos[0], m_Interpol.pos_i_2pp[0][pos[1]], pos[2]).v;

			//alpha - direction
			f4_curr(1, pos[0], pos[1], pos[2]).v  = m_w_i_2p[1][pos[1]].v * innerEngine_f4_curr(1, pos[0], m_Interpol.pos_i_2p[1][pos[1]], pos[2]).v
													+ m_w_i_2pp[1][pos[1]].v * innerEngine_f4_curr(1, pos[0], m_Interpol.pos_i_2pp[1][pos[1]], pos[2]).v;
		}
	}
}
