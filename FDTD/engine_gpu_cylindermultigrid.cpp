/*
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

#include <stdexcept>

#include "engine_gpu_cylindermultigrid.h"
#include "operator_cylindermultigrid.h"

using std::cout;
using std::cerr;
using std::endl;

Engine_GPU_CylinderMultiGrid* Engine_GPU_CylinderMultiGrid::New(const Operator_CylinderMultiGrid* op, const std::string& backend, GPU_Backend* parent)
{
	cout << "Create FDTD engine (GPU, cylindrical multi-grid level " << op->GetMultiGridLevel() << ")" << endl;
	Engine_GPU_CylinderMultiGrid* e = new Engine_GPU_CylinderMultiGrid(op, backend, parent);
	e->Init();
	return e;
}

Engine_GPU_CylinderMultiGrid::Engine_GPU_CylinderMultiGrid(const Operator_CylinderMultiGrid* op, const std::string& backend, GPU_Backend* parent)
	: Engine_GPU(op, backend, parent)
{
	Op_CMG = op;
	m_InnerEngine = NULL;
	m_Link = NULL;
}

Engine_GPU_CylinderMultiGrid::~Engine_GPU_CylinderMultiGrid()
{
	Reset();
}

void Engine_GPU_CylinderMultiGrid::Init()
{
	// this grid: backend, extensions, device or host fallback
	Engine_GPU::Init();

	m_Interpol.split_pos = Op_CMG->m_Split_Pos;
	for (int n=0; n<2; ++n)
		for (unsigned int a=0; a<numLines[1]; ++a)
		{
			m_Interpol.pos_v_2p[n].push_back(Op_CMG->m_interpol_pos_v_2p[n][a]);
			m_Interpol.pos_v_2pp[n].push_back(Op_CMG->m_interpol_pos_v_2pp[n][a]);
			m_Interpol.pos_i_2p[n].push_back(Op_CMG->m_interpol_pos_i_2p[n][a]);
			m_Interpol.pos_i_2pp[n].push_back(Op_CMG->m_interpol_pos_i_2pp[n][a]);
			// all four lanes of the SIMD coefficients are the same
			m_Interpol.w_v_2p[n].push_back(Op_CMG->f4_interpol_v_2p[n][a].f[0]);
			m_Interpol.w_v_2pp[n].push_back(Op_CMG->f4_interpol_v_2pp[n][a].f[0]);
			m_Interpol.w_i_2p[n].push_back(Op_CMG->f4_interpol_i_2p[n][a].f[0]);
			m_Interpol.w_i_2pp[n].push_back(Op_CMG->f4_interpol_i_2pp[n][a].f[0]);
		}

	// the sub-grid engine, created by its operator on a backend sharing our work stream
	Operator_Cylinder* inner_op = Op_CMG->GetInnerOperator();
	inner_op->SetGPUBackend(m_BackendName, m_Backend);
	m_InnerEngine = dynamic_cast<Engine_GPU*>(inner_op->CreateEngine());
	if (m_InnerEngine==NULL)
		throw std::runtime_error("Engine_GPU_CylinderMultiGrid::Init: the sub-grid engine is not a GPU engine");

	// all levels on the device with a device coupling, or all in the host fallback
	if (!m_FieldsOnHost && !m_InnerEngine->FieldsOnHost())
		m_Link = m_Backend->CreateMultiGridLink(m_InnerEngine->GetBackend(), m_Interpol);
	if (m_Link==NULL)
	{
		if (!m_FieldsOnHost || !m_InnerEngine->FieldsOnHost())
		{
			cout << "Engine_GPU: multi-grid level " << Op_CMG->GetMultiGridLevel() << " runs in the host fallback" << endl;
			cerr << "Engine_GPU_CylinderMultiGrid::Init: Warning: multi-grid level " << Op_CMG->GetMultiGridLevel()
			     << " and its sub-grid run in the host fallback" << endl;
		}
		ForceHostFallback();
	}
}

void Engine_GPU_CylinderMultiGrid::Reset()
{
	// the coupling refers to the backends of both grids
	delete m_Link;
	m_Link = NULL;
	delete m_InnerEngine;
	m_InnerEngine = NULL;
	Engine_GPU::Reset();
}

void Engine_GPU_CylinderMultiGrid::ForceHostFallback()
{
	delete m_Link;
	m_Link = NULL;
	Engine_GPU::ForceHostFallback();
	if (m_InnerEngine)
		m_InnerEngine->ForceHostFallback();
}

void Engine_GPU_CylinderMultiGrid::VoltageHalfStep()
{
	Engine_GPU::VoltageHalfStep();
	m_InnerEngine->VoltageHalfStep();
	if (m_Link)
		m_Link->SyncVoltages();
	else
		HostSyncVoltages();
}

void Engine_GPU_CylinderMultiGrid::CurrentHalfStep()
{
	Engine_GPU::CurrentHalfStep();
	m_InnerEngine->CurrentHalfStep();
	if (m_Link)
		m_Link->SyncCurrents();
	else if (m_Interpol.split_pos>=2)
		HostInterpolate(*curr_ptr, m_InnerEngine->HostCurrents(), true, m_Interpol.split_pos-2);
}

void Engine_GPU_CylinderMultiGrid::NextTimestep()
{
	Engine_GPU::NextTimestep();
	m_InnerEngine->NextTimestep();
}

void Engine_GPU_CylinderMultiGrid::FinishBatch()
{
	// the sub-grid first interpolates from its own sub-grids (see Engine_CylinderMultiGrid::IterateTS())
	m_InnerEngine->FinishBatch();
	if (m_Link)
		m_Link->InterpolateToBase();
	else
	{
		for (unsigned int r=0; r+1<m_Interpol.split_pos; ++r)
			HostInterpolate(*volt_ptr, m_InnerEngine->HostVoltages(), false, r);
		for (unsigned int r=0; r+2<m_Interpol.split_pos; ++r)
			HostInterpolate(*curr_ptr, m_InnerEngine->HostCurrents(), true, r);
	}
	Engine_GPU::FinishBatch();
}

// see Engine_Ext_CylinderMultiGrid::SyncVoltages()
void Engine_GPU_CylinderMultiGrid::HostSyncVoltages()
{
	if (m_Interpol.split_pos<1)
		return;
	ArrayLib::ArrayNIJK<FDTD_FLOAT>& base = *volt_ptr;
	ArrayLib::ArrayNIJK<FDTD_FLOAT>& sub = m_InnerEngine->HostVoltages();
	const unsigned int r = m_Interpol.split_pos-1;
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

// see Engine_CylinderMultiGrid::InterpolVoltChild2Base() and InterpolCurrChild2Base()
void Engine_GPU_CylinderMultiGrid::HostInterpolate(ArrayLib::ArrayNIJK<FDTD_FLOAT>& base, const ArrayLib::ArrayNIJK<FDTD_FLOAT>& sub, bool currents, unsigned int r)
{
	const std::vector<unsigned int>* pos_2p  = currents ? m_Interpol.pos_i_2p  : m_Interpol.pos_v_2p;
	const std::vector<unsigned int>* pos_2pp = currents ? m_Interpol.pos_i_2pp : m_Interpol.pos_v_2pp;
	const std::vector<float>* w_2p  = currents ? m_Interpol.w_i_2p  : m_Interpol.w_v_2p;
	const std::vector<float>* w_2pp = currents ? m_Interpol.w_i_2pp : m_Interpol.w_v_2pp;
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
