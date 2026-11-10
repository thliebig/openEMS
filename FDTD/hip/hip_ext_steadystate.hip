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

__global__ void ssd_record(const float* volt, float* records, const size_t* index, SSDParam P)
{
	const unsigned int n = blockIdx.x*blockDim.x + threadIdx.x;
	if (n<P.count)
		records[n*P.length + P.rel_pos] = volt[index[n]];
}

class HIP_Ext_SteadyState : public GPU_Extension
{
public:
	HIP_Ext_SteadyState(const std::vector<GPU_Backend_HIP::Impl*>& parts, Operator_Ext_SteadyState* op_ext, Engine_Ext_SteadyState* eng_ext, Engine_GPU* eng);

	virtual void Apply2Voltages();

protected:
	void CheckPeriod(unsigned int rel_pos);

	//! the probes of one part (see GPU_Backend_HIP::Impl::Owns()), recorded on its device
	struct Probes
	{
		GPU_Backend_HIP::Impl* d;
		unsigned int count;          //!< number of probes
		size_t* index;               //!< flat NIJK index of the probed voltages
		float* records;              //!< [probe][2*period]: the two periods CheckPeriod() compares
		std::vector<float> host;
	};

	Engine_Ext_SteadyState* m_Eng_SSD;
	Engine_GPU* m_Eng;
	unsigned int m_Period;       //!< timesteps per period
	std::vector<Probes> m_Parts;
	//! probe n of the operator extension is probe m_Where[n].second of part m_Where[n].first
	std::vector<std::pair<unsigned int, unsigned int>> m_Where;
};

HIP_Ext_SteadyState::HIP_Ext_SteadyState(const std::vector<GPU_Backend_HIP::Impl*>& parts, Operator_Ext_SteadyState* op_ext, Engine_Ext_SteadyState* eng_ext, Engine_GPU* eng)
{
	m_Eng_SSD = eng_ext;
	m_Eng = eng;
	m_Period = op_ext->m_TS_period;
	const unsigned int count = op_ext->m_E_probe_dir.size();
	m_Where.resize(count);

	// The check compares the probes of all parts at once: each part records the probes it owns,
	// CheckPeriod() downloads the records of all
	for (size_t q=0; q<parts.size(); ++q)
	{
		GPU_Backend_HIP::Impl* d = parts[q];
		Probes P;
		P.d = d;
		std::vector<size_t> index;
		for (unsigned int n=0; n<count; ++n)
		{
			if (!d->Owns(op_ext->m_E_probe_pos[0].at(n)))
				continue;
			m_Where[n] = std::make_pair((unsigned int)q, (unsigned int)index.size());
			index.push_back((((size_t)op_ext->m_E_probe_dir.at(n)*d->dim.nx + d->LocalX(op_ext->m_E_probe_pos[0].at(n)))*d->dim.ny
			                 + op_ext->m_E_probe_pos[1].at(n))*d->dim.nz + op_ext->m_E_probe_pos[2].at(n));
		}
		P.count = index.size();
		index.resize(std::max(P.count, 1u), 0);
		P.index = d->Alloc<size_t>(index.size(), index.data());
		P.host.resize(std::max(P.count, 1u)*2*(size_t)m_Period);
		P.records = d->Alloc<float>(P.host.size());
		m_Parts.push_back(P);
	}
}

void HIP_Ext_SteadyState::Apply2Voltages()
{
	const unsigned int p = m_Period;
	const unsigned int TS = m_Eng->GetNumberOfTimesteps();
	const unsigned int rel_pos = TS%(2*p);

	for (Probes& P : m_Parts)
		if (P.count>0)
		{
			SSDParam param = {P.count, 2*p, rel_pos};
			HIP_Launch(P.d, "ssd_record", ssd_record, P.count, 1, 1, (const float*)P.d->volt, P.records, (const size_t*)P.index, param);
		}

	if ((TS%p==0) && (TS>=2*p))
		CheckPeriod(rel_pos);
}

// the check of Engine_Ext_SteadyState, on the records downloaded from the device
void HIP_Ext_SteadyState::CheckPeriod(unsigned int rel_pos)
{
	const unsigned int p = m_Period;

	// the energy comes from the GPU engine interface, computed on the device
	for (Probes& P : m_Parts)
	{
		P.d->ctx->MakeCurrent();
		HIP_Check(hipMemcpyAsync(P.host.data(), P.records, P.host.size()*sizeof(float), hipMemcpyDeviceToHost, P.d->Stream()), "copy the steady-state records");
	}
	for (Probes& P : m_Parts)
		P.d->Flush();

	const std::vector<Probes>& parts = m_Parts;
	const std::vector<std::pair<unsigned int, unsigned int>>& where = m_Where;
	m_Eng_SSD->CheckPeriod(rel_pos, [&parts, &where, p](size_t n, unsigned int t)
	                       {return (double)parts[where[n].first].host[where[n].second*2*(size_t)p + t];});
}

GPU_Extension* HIP_CreateExt_SteadyStateParts(const std::vector<GPU_Backend_HIP::Impl*>& parts, Engine_Extension* eng_ext, Engine* eng)
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
	return new HIP_Ext_SteadyState(parts, op_ext, ssd_ext, gpu_eng);
}

GPU_Extension* HIP_CreateExt_SteadyState(GPU_Backend_HIP::Impl* d, Engine_Extension* eng_ext, Engine* eng)
{
	return HIP_CreateExt_SteadyStateParts(std::vector<GPU_Backend_HIP::Impl*>(1, d), eng_ext, eng);
}
