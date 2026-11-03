/*
*	Copyright (C) 2026 Sean Mollet (sean@malmoset.com)
*	Copyright (C) 2026 Alexandru Iercosan-Lucaci (alex.lucaci@alacrity.ro)
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

#include "sycl_internal.h"
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

static inline void ssd_record(const sycl::nd_item<3>& it, const float* volt, float* records, const size_t* index, SSDParam P)
{
	const unsigned int n = it.get_group(2)*it.get_local_range(2) + it.get_local_id(2);
	if (n<P.count)
		records[n*P.length + P.rel_pos] = volt[index[n]];
}

class SYCL_Ext_SteadyState : public GPU_Extension
{
public:
	SYCL_Ext_SteadyState(GPU_Backend_SYCL::Impl* impl, Operator_Ext_SteadyState* op_ext, Engine_Ext_SteadyState* eng_ext, Engine_GPU* eng);

	virtual void Apply2Voltages();

protected:
	void CheckPeriod(unsigned int rel_pos);

	GPU_Backend_SYCL::Impl* d;
	Engine_Ext_SteadyState* m_Eng_SSD;
	Engine_GPU* m_Eng;
	unsigned int m_Period;       //!< timesteps per period
	unsigned int m_Count;        //!< number of probes
	size_t* m_Index;             //!< flat NIJK index of the probed voltages
	float* m_Records;            //!< [probe][2*period]: the two periods CheckPeriod() compares
	std::vector<float> m_HostRecords;
};

SYCL_Ext_SteadyState::SYCL_Ext_SteadyState(GPU_Backend_SYCL::Impl* impl, Operator_Ext_SteadyState* op_ext, Engine_Ext_SteadyState* eng_ext, Engine_GPU* eng)
{
	d = impl;
	m_Eng_SSD = eng_ext;
	m_Eng = eng;
	m_Period = op_ext->m_TS_period;
	m_Count = op_ext->m_E_probe_dir.size();

	std::vector<size_t> index(std::max(m_Count, 1u), 0);
	for (unsigned int n=0; n<m_Count; ++n)
		index[n] = (((size_t)op_ext->m_E_probe_dir.at(n)*d->dim.nx + op_ext->m_E_probe_pos[0].at(n))*d->dim.ny
		            + op_ext->m_E_probe_pos[1].at(n))*d->dim.nz + op_ext->m_E_probe_pos[2].at(n);
	m_Index = d->Alloc<size_t>(index.size(), index.data());
	m_HostRecords.resize(std::max(m_Count, 1u)*2*(size_t)m_Period);
	m_Records = d->Alloc<float>(m_HostRecords.size());
}

void SYCL_Ext_SteadyState::Apply2Voltages()
{
	const unsigned int p = m_Period;
	const unsigned int TS = m_Eng->GetNumberOfTimesteps();
	const unsigned int rel_pos = TS%(2*p);

	if (m_Count>0)
	{
		SSDParam param = {m_Count, 2*p, rel_pos};
		const float* volt = d->volt;
		float* records = m_Records;
		const size_t* index = m_Index;
		SYCL_Launch(d, "ssd_record", m_Count, 1, 1, [=](sycl::nd_item<3> it)
		{
			ssd_record(it, volt, records, index, param);
		});
	}

	if ((TS%p==0) && (TS>=2*p))
		CheckPeriod(rel_pos);
}

// the check of Engine_Ext_SteadyState, on the records downloaded from the device
void SYCL_Ext_SteadyState::CheckPeriod(unsigned int rel_pos)
{
	const unsigned int p = m_Period;

	// the energy comes from the GPU engine interface, computed on the device
	SYCL_Check("records download", [&] {d->Queue().memcpy(m_HostRecords.data(), m_Records, m_HostRecords.size()*sizeof(float));});
	d->Flush();

	const float* records = m_HostRecords.data();
	m_Eng_SSD->CheckPeriod(rel_pos, [records, p](size_t n, unsigned int t) {return (double)records[n*2*(size_t)p + t];});
}

GPU_Extension* SYCL_CreateExt_SteadyState(GPU_Backend_SYCL::Impl* d, Engine_Extension* eng_ext, Engine* eng)
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
	return new SYCL_Ext_SteadyState(d, op_ext, ssd_ext, gpu_eng);
}
