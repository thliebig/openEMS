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

#include "hip/hip_runtime.h"
#include <cmath>

#include "hip_internal.h"
#include "FDTD/engine_gpu.h"
#include "FDTD/engine_interface_fdtd.h"
#include "FDTD/extensions/operator_ext_steadystate.h"
#include "FDTD/extensions/engine_ext_steadystate.h"

// Steady-state detection, see Engine_Ext_SteadyState. The probe voltages are
// recorded on the device every timestep; once per period the fields and the
// records are copied to the host, which compares the energy and the recorded
// periods, and stores the result in the engine extension, where openEMS reads
// it (GetLastDiff()).
struct SSDParam { unsigned int count; unsigned int length; unsigned int rel_pos; };

__global__ void ssd_record(const float* volt, float* records, const unsigned int* index, SSDParam P)
{
	const unsigned int n = blockIdx.x*blockDim.x + threadIdx.x;
	if (n<P.count)
		records[n*P.length + P.rel_pos] = volt[index[n]];
}

class HIP_Ext_SteadyState : public GPU_Extension
{
public:
	HIP_Ext_SteadyState(GPU_Backend_HIP::Impl* impl, Operator_Ext_SteadyState* op_ext, Engine_Ext_SteadyState* eng_ext, Engine_GPU* eng);

	virtual void Apply2Voltages();

protected:
	void CheckPeriod(unsigned int rel_pos);

	GPU_Backend_HIP::Impl* d;
	Engine_Ext_SteadyState* m_Eng_SSD;
	Engine_GPU* m_Eng;
	unsigned int m_Period;       //!< timesteps per period
	unsigned int m_Count;        //!< number of probes
	unsigned int* m_Index;       //!< flat NIJK index of the probed voltages
	float* m_Records;            //!< [probe][2*period]: the two periods CheckPeriod() compares
	std::vector<float> m_HostRecords;
};

HIP_Ext_SteadyState::HIP_Ext_SteadyState(GPU_Backend_HIP::Impl* impl, Operator_Ext_SteadyState* op_ext, Engine_Ext_SteadyState* eng_ext, Engine_GPU* eng)
{
	d = impl;
	m_Eng_SSD = eng_ext;
	m_Eng = eng;
	m_Period = op_ext->m_TS_period;
	m_Count = op_ext->m_E_probe_dir.size();

	std::vector<unsigned int> index(std::max(m_Count, 1u), 0);
	for (unsigned int n=0; n<m_Count; ++n)
		index[n] = ((op_ext->m_E_probe_dir.at(n)*d->dim.nx + op_ext->m_E_probe_pos[0].at(n))*d->dim.ny
		            + op_ext->m_E_probe_pos[1].at(n))*d->dim.nz + op_ext->m_E_probe_pos[2].at(n);
	m_Index = d->Alloc<unsigned int>(index.size(), index.data());
	m_HostRecords.resize(std::max(m_Count, 1u)*2*(size_t)m_Period);
	m_Records = d->Alloc<float>(m_HostRecords.size());
}

void HIP_Ext_SteadyState::Apply2Voltages()
{
	const unsigned int p = m_Period;
	const unsigned int TS = m_Eng->GetNumberOfTimesteps();
	const unsigned int rel_pos = TS%(2*p);

	if (m_Count>0)
	{
		SSDParam param = {m_Count, 2*p, rel_pos};
		HIP_Launch(d, "ssd_record", ssd_record, m_Count, 1, 1, (const float*)d->volt, m_Records, (const unsigned int*)m_Index, param);
	}

	if ((TS%p==0) && (TS>=2*p))
		CheckPeriod(rel_pos);
}

// same as the period check in Engine_Ext_SteadyState::Apply2Voltages()
void HIP_Ext_SteadyState::CheckPeriod(unsigned int rel_pos)
{
	const unsigned int p = m_Period;

	// the energy is computed by the host on the host mirror
	m_Eng->VoltagesToHost();
	m_Eng->CurrentsToHost();
	HIP_Check(hipMemcpyAsync(m_HostRecords.data(), m_Records, m_HostRecords.size()*sizeof(float), hipMemcpyDeviceToHost, d->Stream()), "copy the steady-state records");
	d->Flush();

	bool no_valid = true;
	double& last_max_diff = m_Eng_SSD->m_last_max_diff;
	last_max_diff = 0;
	double curr_total_energy = m_Eng_SSD->m_Eng_Interface->CalcFastEnergy();
	if (m_Eng_SSD->last_total_energy>0)
	{
		last_max_diff = std::fabs(curr_total_energy-m_Eng_SSD->last_total_energy)/m_Eng_SSD->last_total_energy;
		no_valid = false;
	}
	m_Eng_SSD->last_total_energy = curr_total_energy;

	unsigned int old_pos = 0;
	unsigned int new_pos = p;
	if (rel_pos<=p)
	{
		new_pos = 0;
		old_pos = p;
	}

	std::vector<double> curr_pow(m_Count, 0), diff_pow(m_Count, 0);
	double max_pow = 0;
	for (unsigned int n=0; n<m_Count; ++n)
	{
		const float* buf = m_HostRecords.data() + n*2*(size_t)p;
		for (unsigned int nt=0; nt<p; ++nt)
		{
			const double b_new = buf[nt+new_pos];
			const double b_old = buf[nt+old_pos];
			curr_pow[n] += b_new*b_new;
			diff_pow[n] += (b_old-b_new)*(b_old-b_new);
		}
		max_pow = std::max(max_pow, curr_pow[n]);
	}
	for (unsigned int n=0; n<m_Count; ++n)
	{
		if (curr_pow[n]>max_pow*1e-2)
		{
			last_max_diff = std::max(last_max_diff, diff_pow[n]/curr_pow[n]);
			no_valid = false;
		}
	}
	if ((no_valid) || (last_max_diff>1))
		last_max_diff = 1;
}

GPU_Extension* HIP_CreateExt_SteadyState(GPU_Backend_HIP::Impl* d, Engine_Extension* eng_ext, Engine* eng)
{
	Engine_Ext_SteadyState* ssd_ext = dynamic_cast<Engine_Ext_SteadyState*>(eng_ext);
	if (!ssd_ext)
		return NULL;
	Operator_Ext_SteadyState* op_ext = dynamic_cast<Operator_Ext_SteadyState*>(eng_ext->GetOperatorExtension());
	if (!op_ext)
		return NULL;
	Engine_GPU* gpu_eng = dynamic_cast<Engine_GPU*>(eng);
	if (!gpu_eng)
		return NULL;
	return new HIP_Ext_SteadyState(d, op_ext, ssd_ext, gpu_eng);
}
