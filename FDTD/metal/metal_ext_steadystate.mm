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

#include <cmath>

#include "metal_internal.h"
#include "FDTD/engine.h"
#include "FDTD/engine_interface_fdtd.h"
#include "FDTD/extensions/operator_ext_steadystate.h"
#include "FDTD/extensions/engine_ext_steadystate.h"

// Steady-state detection, see Engine_Ext_SteadyState. The probe voltages are
// recorded on the device every timestep; once per period the device is
// synchronized and the host compares the energy (read from the shared host
// mirror) and the recorded periods, and stores the result in the engine
// extension, where openEMS reads it (GetLastDiff()).
static const char* STEADYSTATE_SOURCE = R"MSL(
struct SSDParam { uint count; uint length; uint rel_pos; };

kernel void ssd_record(const device float* volt   [[buffer(0)]],
                       device float* records      [[buffer(1)]],
                       const device uint* index   [[buffer(2)]],
                       constant SSDParam& P       [[buffer(3)]],
                       uint n [[thread_position_in_grid]])
{
	if (n<P.count)
		records[n*P.length + P.rel_pos] = volt[index[n]];
}
)MSL";

class Metal_Ext_SteadyState : public GPU_Extension
{
public:
	Metal_Ext_SteadyState(GPU_Backend_Metal::Impl* impl, Operator_Ext_SteadyState* op_ext, Engine_Ext_SteadyState* eng_ext, Engine* eng);

	virtual void Apply2Voltages();

protected:
	void CheckPeriod(unsigned int rel_pos);

	GPU_Backend_Metal::Impl* d;
	Engine_Ext_SteadyState* m_Eng_SSD;
	Engine* m_Eng;
	unsigned int m_Period;   //!< timesteps per period
	unsigned int m_Count;    //!< number of probes
	id<MTLBuffer> m_Index;   //!< flat NIJK index of the probed voltages
	id<MTLBuffer> m_Records; //!< [probe][2*period]
};

Metal_Ext_SteadyState::Metal_Ext_SteadyState(GPU_Backend_Metal::Impl* impl, Operator_Ext_SteadyState* op_ext, Engine_Ext_SteadyState* eng_ext, Engine* eng)
{
	d = impl;
	m_Eng_SSD = eng_ext;
	m_Eng = eng;
	m_Period = op_ext->m_TS_period;
	m_Count = op_ext->m_E_probe_dir.size();

	std::vector<uint32_t> index(std::max(m_Count, 1u), 0);
	for (unsigned int n=0; n<m_Count; ++n)
		index[n] = ((op_ext->m_E_probe_dir.at(n)*d->dim.nx + op_ext->m_E_probe_pos[0].at(n))*d->dim.ny
		            + op_ext->m_E_probe_pos[1].at(n))*d->dim.nz + op_ext->m_E_probe_pos[2].at(n);
	m_Index = d->NewBuffer(index.size()*sizeof(uint32_t), index.data());
	m_Records = d->NewBuffer(std::max(m_Count, 1u)*2*(size_t)m_Period*sizeof(float));
	d->Pipeline(STEADYSTATE_SOURCE, "ssd_record");
}

void Metal_Ext_SteadyState::Apply2Voltages()
{
	const unsigned int p = m_Period;
	const unsigned int TS = m_Eng->GetNumberOfTimesteps();
	const unsigned int rel_pos = TS%(2*p);

	if (m_Count>0)
	{
		struct {uint32_t count, length, rel_pos;} param = {m_Count, 2*p, rel_pos};
		id<MTLComputePipelineState> pso = d->Pipeline(STEADYSTATE_SOURCE, "ssd_record");
		id<MTLComputeCommandEncoder> enc = d->Encoder();
		[enc setComputePipelineState:pso];
		[enc setBuffer:d->volt offset:0 atIndex:0];
		[enc setBuffer:m_Records offset:0 atIndex:1];
		[enc setBuffer:m_Index offset:0 atIndex:2];
		[enc setBytes:&param length:sizeof(param) atIndex:3];
		d->Dispatch(pso, m_Count);
	}

	if ((TS%p==0) && (TS>=2*p))
		CheckPeriod(rel_pos);
}

// the check of Engine_Ext_SteadyState, on the records in shared memory
void Metal_Ext_SteadyState::CheckPeriod(unsigned int rel_pos)
{
	const unsigned int p = m_Period;

	// finish all work so far, the host mirror (shared memory) and the records are up to date
	d->Flush();

	const float* records = static_cast<const float*>([m_Records contents]);
	m_Eng_SSD->CheckPeriod(rel_pos, [records, p](size_t n, unsigned int t) {return (double)records[n*2*(size_t)p + t];});
}

GPU_Extension* Metal_CreateExt_SteadyState(GPU_Backend_Metal::Impl* d, Engine_Extension* eng_ext, Engine* eng)
{
	Engine_Ext_SteadyState* ssd_ext = dynamic_cast<Engine_Ext_SteadyState*>(eng_ext);
	if (!ssd_ext)
		return NULL;
	Operator_Ext_SteadyState* op_ext = dynamic_cast<Operator_Ext_SteadyState*>(eng_ext->GetOperatorExtension());
	if (!op_ext)
		return NULL;
	// the energy is computed by the host on the host mirror, which is the (shared) device memory
	return new Metal_Ext_SteadyState(d, op_ext, ssd_ext, eng);
}
