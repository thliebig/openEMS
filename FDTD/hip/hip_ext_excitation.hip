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
#include <set>

#include "hip_internal.h"
#include "FDTD/engine.h"
#include "FDTD/excitation.h"
#include "FDTD/extensions/operator_ext_excitation.h"
#include "FDTD/extensions/engine_ext_excitation.h"

// Soft voltage/current excitation, see Engine_Ext_Excitation.
// One thread per excitation point, or a single thread if points share an edge.
struct ExcitationParam { int numTS; int period; int length; unsigned int count; };

__device__ __forceinline__ void excite(float* field, const unsigned int* index, const float* amp,
                                       const unsigned int* delay, const float* signal, const ExcitationParam& P, unsigned int n)
{
	int exc_pos = P.numTS - (int)delay[n];
	exc_pos *= (exc_pos>0);
	exc_pos %= P.period;
	exc_pos *= (exc_pos<P.length);
	field[index[n]] = field[index[n]] + amp[n]*signal[exc_pos];
}

__global__ void excitation_apply(float* field, const unsigned int* index, const float* amp,
                                 const unsigned int* delay, const float* signal, ExcitationParam P)
{
	const unsigned int n = blockIdx.x*blockDim.x + threadIdx.x;
	if (n<P.count)
		excite(field, index, amp, delay, signal, P, n);
}

// same, in the order of the host, for excitation points sharing an edge
__global__ void excitation_apply_serial(float* field, const unsigned int* index, const float* amp,
                                        const unsigned int* delay, const float* signal, ExcitationParam P)
{
	if (blockIdx.x*blockDim.x + threadIdx.x != 0)
		return;
	for (unsigned int n=0; n<P.count; ++n)
		excite(field, index, amp, delay, signal, P, n);
}

class HIP_Ext_Excitation : public GPU_Extension
{
public:
	HIP_Ext_Excitation(GPU_Backend_HIP::Impl* impl, Operator_Ext_Excitation* op_ext, Engine* eng);

	virtual void Apply2Voltages() {Apply(m_Volt, d->volt);}
	virtual void Apply2Current() {Apply(m_Curr, d->curr);}

protected:
	struct Sources
	{
		unsigned int count;
		bool serial;            //!< points share an edge, apply them one after the other
		unsigned int* index;    //!< flat NIJK index of the excited edge
		float* amp;
		unsigned int* delay;
		float* signal;
	};

	void Setup(Sources& s, unsigned int count, unsigned int* const index[3], const unsigned short* dir,
	           const FDTD_FLOAT* amp, const unsigned int* delay, const FDTD_FLOAT* signal);
	void Apply(Sources& s, float* field);

	GPU_Backend_HIP::Impl* d;
	Operator_Ext_Excitation* m_Op_Exc;
	Engine* m_Eng;

	Sources m_Volt;
	Sources m_Curr;
};

HIP_Ext_Excitation::HIP_Ext_Excitation(GPU_Backend_HIP::Impl* impl, Operator_Ext_Excitation* op_ext, Engine* eng)
{
	d = impl;
	m_Op_Exc = op_ext;
	m_Eng = eng;
	Setup(m_Volt, op_ext->Volt_Count, op_ext->Volt_index, op_ext->Volt_dir, op_ext->Volt_amp, op_ext->Volt_delay, op_ext->m_Exc->GetVoltageSignal());
	Setup(m_Curr, op_ext->Curr_Count, op_ext->Curr_index, op_ext->Curr_dir, op_ext->Curr_amp, op_ext->Curr_delay, op_ext->m_Exc->GetCurrentSignal());

	// fused step (see update_fused): the currents next to the excited voltages are recomputed;
	// excited currents on the last mesh lines (never updated) would differ between its buffers
	for (unsigned int n=0; n<op_ext->Volt_Count; ++n)
		d->volt_modified.push_back(((op_ext->Volt_dir[n]*d->dim.nx + op_ext->Volt_index[0][n])*d->dim.ny + op_ext->Volt_index[1][n])*d->dim.nz + op_ext->Volt_index[2][n]);
	for (unsigned int n=0; n<op_ext->Curr_Count; ++n)
		if ((op_ext->Curr_index[0][n]+1>=d->dim.nx) || (op_ext->Curr_index[1][n]+1>=d->dim.ny) || (op_ext->Curr_index[2][n]+1>=d->dim.nz))
			d->fused_blockers++;
}

void HIP_Ext_Excitation::Setup(Sources& s, unsigned int count, unsigned int* const index[3], const unsigned short* dir,
                                const FDTD_FLOAT* amp, const unsigned int* delay, const FDTD_FLOAT* signal)
{
	s.count = count;
	s.serial = false;
	if (count==0)
		return;

	std::vector<unsigned int> flat(count);
	std::set<unsigned int> edges;
	for (unsigned int n=0; n<count; ++n)
	{
		flat[n] = ((dir[n]*d->dim.nx + index[0][n])*d->dim.ny + index[1][n])*d->dim.nz + index[2][n];
		if (!edges.insert(flat[n]).second)
			s.serial = true;
	}
	s.index = d->Alloc<unsigned int>(count, flat.data());
	s.amp = d->Alloc<float>(count, amp);
	s.delay = d->Alloc<unsigned int>(count, delay);
	s.signal = d->Alloc<float>(m_Op_Exc->m_Exc->GetLength(), signal);
}

void HIP_Ext_Excitation::Apply(Sources& s, float* field)
{
	if (s.count==0)
		return;

	ExcitationParam param;
	param.numTS = m_Eng->GetNumberOfTimesteps();
	param.period = param.numTS+1;
	if (m_Op_Exc->m_Exc->GetSignalPeriod()>0)
		param.period = int(m_Op_Exc->m_Exc->GetSignalPeriod()/m_Op_Exc->m_Exc->GetTimestep());
	param.length = m_Op_Exc->m_Exc->GetLength();
	param.count = s.count;

	if (s.serial)
		HIP_Launch(d, "excitation_apply_serial", excitation_apply_serial, 1, 1, 1,
		            field, (const unsigned int*)s.index, (const float*)s.amp, (const unsigned int*)s.delay, (const float*)s.signal, param);
	else
		HIP_Launch(d, "excitation_apply", excitation_apply, s.count, 1, 1,
		            field, (const unsigned int*)s.index, (const float*)s.amp, (const unsigned int*)s.delay, (const float*)s.signal, param);
}

GPU_Extension* HIP_CreateExt_Excitation(GPU_Backend_HIP::Impl* d, Engine_Extension* eng_ext, Engine* eng)
{
	if (!dynamic_cast<Engine_Ext_Excitation*>(eng_ext))
		return NULL;
	Operator_Ext_Excitation* op_ext = dynamic_cast<Operator_Ext_Excitation*>(eng_ext->GetOperatorExtension());
	if (!op_ext)
		return NULL;
	return new HIP_Ext_Excitation(d, op_ext, eng);
}
