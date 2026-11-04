/*
*	Additional
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

#include "engine_ext_lumpedRLC.h"
#include "operator_ext_lumpedRLC.h"

#include "FDTD/engine_sse.h"

Engine_Ext_LumpedRLC::Engine_Ext_LumpedRLC(Operator_Ext_LumpedRLC* op_ext_RLC) : Engine_Extension(op_ext_RLC)
{
	// Local pointer of the operator.
	m_Op_Ext_RLC = op_ext_RLC;

	v_Vdn		= new FDTD_FLOAT*[3];
	v_Jn		= new FDTD_FLOAT*[3];
	v_Il		= NULL;

	// No additional allocations are required if there are no actual lumped elements.
	if (!(m_Op_Ext_RLC->RLC_count))
		return;

	// Initialize ADE containers for currents and voltages
	v_Il		= new FDTD_FLOAT[m_Op_Ext_RLC->RLC_count];

	for (size_t posIdx = 0 ; posIdx < m_Op_Ext_RLC->RLC_count ; ++posIdx)
		v_Il[posIdx] 	= 0.0;

	for (unsigned int k = 0 ; k < 3 ; k++)
	{
		v_Vdn[k] = new FDTD_FLOAT[m_Op_Ext_RLC->RLC_count];
		v_Jn[k] = new FDTD_FLOAT[m_Op_Ext_RLC->RLC_count];

		for (size_t posIdx = 0 ; posIdx < m_Op_Ext_RLC->RLC_count ; ++posIdx)
		{
			v_Jn[k][posIdx] = 0.0;
			v_Vdn[k][posIdx] = 0.0;
		}
	}

}

Engine_Ext_LumpedRLC::~Engine_Ext_LumpedRLC()
{
	// Only delete if values were allocated in the first place
	if (m_Op_Ext_RLC->RLC_count)
	{
		delete[] v_Il;

		for (unsigned int k = 0 ; k < 3 ; k++)
		{
			delete[] v_Vdn[k];
			delete[] v_Jn[k];
		}
	}

	delete[] v_Vdn;
	delete[] v_Jn;

	v_Il	= NULL;

	v_Vdn	= NULL;
	v_Jn	= NULL;

	m_Op_Ext_RLC = NULL;


}

bool Engine_Ext_LumpedRLC::InBox(size_t pIdx, const Engine_Box* box) const
{
	if (box==NULL)
		return true;
	for (int a=0; a<3; ++a)
	{
		const unsigned int p = m_Op_Ext_RLC->v_RLC_pos[a][pIdx];
		if ((p<box->start[a]) || (p>=box->stop[a]))
			return false;
	}
	return true;
}

// The rings hold timesteps ts, ts-1, ts-2 at [ts%3], [(ts+2)%3], [(ts+1)%3]: what the
// rotation of three pointers per timestep did, but independent of the calls, so the
// elements can run per box of a blocked plan.
void Engine_Ext_LumpedRLC::DoPreVoltageUpdatesImpl(unsigned int ts, const Engine_Box* box)
{
	if (!m_Op_Ext_RLC->RLC_count)
		return;
	const FDTD_FLOAT* Vd1 = v_Vdn[(ts+2)%3];

	// In pre-process, only update the parallel inductor current:
	for (size_t pIdx = 0 ; pIdx < m_Op_Ext_RLC->RLC_count ; pIdx++)
		if (InBox(pIdx, box))
			v_Il[pIdx] += (m_Op_Ext_RLC->v_RLC_i2v[pIdx])*(m_Op_Ext_RLC->v_RLC_ilv[pIdx])*Vd1[pIdx];
}

void Engine_Ext_LumpedRLC::DoPreVoltageUpdates()
{
	DoPreVoltageUpdatesImpl(m_Eng->GetNumberOfTimesteps(), NULL);
}

void Engine_Ext_LumpedRLC::DoPreVoltageUpdatesBox(const Engine_Box& box, unsigned int ts)
{
	DoPreVoltageUpdatesImpl(ts, &box);
}

template <typename EngType>
void Engine_Ext_LumpedRLC::Apply2VoltagesImpl(EngType* eng, unsigned int ts, const Engine_Box* box)
{
	if (!m_Op_Ext_RLC->RLC_count)
		return;

	unsigned int **pos = m_Op_Ext_RLC->v_RLC_pos;
	int *dir = m_Op_Ext_RLC->v_RLC_dir;
	FDTD_FLOAT* Vd0 = v_Vdn[ts%3];
	const FDTD_FLOAT* Vd2 = v_Vdn[(ts+1)%3];
	FDTD_FLOAT* J0 = v_Jn[ts%3];
	const FDTD_FLOAT* J1 = v_Jn[(ts+2)%3];
	const FDTD_FLOAT* J2 = v_Jn[(ts+1)%3];

	// Read engine calculated node voltage
	for (size_t pIdx = 0 ; pIdx < m_Op_Ext_RLC->RLC_count ; pIdx++)
		if (InBox(pIdx, box))
			Vd0[pIdx] = eng->EngType::GetVolt(dir[pIdx],pos[0][pIdx],pos[1][pIdx],pos[2][pIdx]);

	// Post process: Calculate node voltage with respect to the lumped RLC auxilliary quantity, J
	for (size_t pIdx = 0 ; pIdx < m_Op_Ext_RLC->RLC_count ; pIdx++)
	{
		if (!InBox(pIdx, box))
			continue;
		// Calculate updated node voltage, with series and parallel additions
		Vd0[pIdx] = 	(m_Op_Ext_RLC->v_RLC_vvd[pIdx])*(
							Vd0[pIdx] - v_Il[pIdx]						// Addition for Parallel inductor
							+
							(m_Op_Ext_RLC->v_RLC_vv2[pIdx])*Vd2[pIdx]	// Vd[n-2] addition
							+
							(m_Op_Ext_RLC->v_RLC_vj1[pIdx])*J1[pIdx]	// J[n-1] addition
							+
							(m_Op_Ext_RLC->v_RLC_vj2[pIdx])*J2[pIdx]);	// J[n-2] addition

		// Update J[0]
		J0[pIdx] =	(m_Op_Ext_RLC->v_RLC_ib0[pIdx])*(Vd0[pIdx] - Vd2[pIdx])
						-
						((m_Op_Ext_RLC->v_RLC_b1[pIdx])*(m_Op_Ext_RLC->v_RLC_ib0[pIdx]))*J1[pIdx]
						-
						((m_Op_Ext_RLC->v_RLC_b2[pIdx])*(m_Op_Ext_RLC->v_RLC_ib0[pIdx]))*J2[pIdx];
	}

	// Update node voltage
	for (size_t pIdx = 0 ; pIdx < m_Op_Ext_RLC->RLC_count ; pIdx++)
		if (InBox(pIdx, box))
			eng->EngType::SetVolt(dir[pIdx],pos[0][pIdx],pos[1][pIdx],pos[2][pIdx],Vd0[pIdx]);
}

void Engine_Ext_LumpedRLC::Apply2Voltages()
{
	ENG_DISPATCH_ARGS(Apply2VoltagesImpl, m_Eng->GetNumberOfTimesteps(), (const Engine_Box*)NULL);
}

void Engine_Ext_LumpedRLC::Apply2VoltagesBox(const Engine_Box& box, unsigned int ts)
{
	ENG_DISPATCH_ARGS(Apply2VoltagesImpl, ts, &box);
}
