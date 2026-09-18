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

#include <set>

#include "metal_internal.h"
#include "FDTD/engine.h"
#include "FDTD/excitation.h"
#include "FDTD/extensions/operator_ext_excitation.h"
#include "FDTD/extensions/engine_ext_excitation.h"

// Soft voltage/current excitation, see Engine_Ext_Excitation.
// One thread per excitation point, or a single thread if points share an edge.
static const char* EXCITATION_SOURCE = R"MSL(
struct ExcitationParam { int numTS; int period; int length; uint count; };

inline void excite(device float* field, const device uint* index, const device float* amp,
                   const device uint* delay, const device float* signal, constant ExcitationParam& P, uint n)
{
	int exc_pos = P.numTS - (int)delay[n];
	exc_pos *= (exc_pos>0);
	exc_pos %= P.period;
	exc_pos *= (exc_pos<P.length);
	field[index[n]] = field[index[n]] + amp[n]*signal[exc_pos];
}

kernel void excitation_apply(device float* field          [[buffer(0)]],
                             const device uint* index     [[buffer(1)]],
                             const device float* amp      [[buffer(2)]],
                             const device uint* delay     [[buffer(3)]],
                             const device float* signal   [[buffer(4)]],
                             constant ExcitationParam& P  [[buffer(5)]],
                             uint n [[thread_position_in_grid]])
{
	if (n<P.count)
		excite(field, index, amp, delay, signal, P, n);
}

// same, in the order of the host, for excitation points sharing an edge
kernel void excitation_apply_serial(device float* field          [[buffer(0)]],
                                    const device uint* index     [[buffer(1)]],
                                    const device float* amp      [[buffer(2)]],
                                    const device uint* delay     [[buffer(3)]],
                                    const device float* signal   [[buffer(4)]],
                                    constant ExcitationParam& P  [[buffer(5)]])
{
	for (uint n=0; n<P.count; ++n)
		excite(field, index, amp, delay, signal, P, n);
}
)MSL";

class Metal_Ext_Excitation : public GPU_Extension
{
public:
	Metal_Ext_Excitation(GPU_Backend_Metal::Impl* impl, Operator_Ext_Excitation* op_ext, Engine* eng);

	virtual void Apply2Voltages() {Apply(m_Volt, d->volt);}
	virtual void Apply2Current() {Apply(m_Curr, d->curr);}

protected:
	struct Sources
	{
		unsigned int count;
		bool serial;          //!< points share an edge, apply them one after the other
		id<MTLBuffer> index;  //!< flat NIJK index of the excited edge
		id<MTLBuffer> amp;
		id<MTLBuffer> delay;
		id<MTLBuffer> signal;
	};

	void Setup(Sources& s, unsigned int count, unsigned int* const index[3], const unsigned short* dir,
	           const FDTD_FLOAT* amp, const unsigned int* delay, const FDTD_FLOAT* signal);
	void Apply(Sources& s, id<MTLBuffer> field);

	GPU_Backend_Metal::Impl* d;
	Operator_Ext_Excitation* m_Op_Exc;
	Engine* m_Eng;

	Sources m_Volt;
	Sources m_Curr;
};

Metal_Ext_Excitation::Metal_Ext_Excitation(GPU_Backend_Metal::Impl* impl, Operator_Ext_Excitation* op_ext, Engine* eng)
{
	d = impl;
	m_Op_Exc = op_ext;
	m_Eng = eng;
	Setup(m_Volt, op_ext->Volt_Count, op_ext->Volt_index, op_ext->Volt_dir, op_ext->Volt_amp, op_ext->Volt_delay, op_ext->m_Exc->GetVoltageSignal());
	Setup(m_Curr, op_ext->Curr_Count, op_ext->Curr_index, op_ext->Curr_dir, op_ext->Curr_amp, op_ext->Curr_delay, op_ext->m_Exc->GetCurrentSignal());
	d->Pipeline(EXCITATION_SOURCE, "excitation_apply");
	d->Pipeline(EXCITATION_SOURCE, "excitation_apply_serial");
}

void Metal_Ext_Excitation::Setup(Sources& s, unsigned int count, unsigned int* const index[3], const unsigned short* dir,
                                 const FDTD_FLOAT* amp, const unsigned int* delay, const FDTD_FLOAT* signal)
{
	s.count = count;
	s.serial = false;
	if (count==0)
		return;

	std::vector<uint32_t> flat(count);
	std::set<uint32_t> edges;
	for (unsigned int n=0; n<count; ++n)
	{
		flat[n] = ((dir[n]*d->dim.nx + index[0][n])*d->dim.ny + index[1][n])*d->dim.nz + index[2][n];
		if (!edges.insert(flat[n]).second)
			s.serial = true;
	}
	s.index = d->NewBuffer(count*sizeof(uint32_t), flat.data());
	s.amp = d->NewBuffer(count*sizeof(float), amp);
	s.delay = d->NewBuffer(count*sizeof(uint32_t), delay);
	s.signal = d->NewBuffer(m_Op_Exc->m_Exc->GetLength()*sizeof(float), signal);
}

void Metal_Ext_Excitation::Apply(Sources& s, id<MTLBuffer> field)
{
	if (s.count==0)
		return;

	struct {int numTS; int period; int length; uint32_t count;} param;
	param.numTS = m_Eng->GetNumberOfTimesteps();
	param.period = param.numTS+1;
	if (m_Op_Exc->m_Exc->GetSignalPeriod()>0)
		param.period = int(m_Op_Exc->m_Exc->GetSignalPeriod()/m_Op_Exc->m_Exc->GetTimestep());
	param.length = m_Op_Exc->m_Exc->GetLength();
	param.count = s.count;

	id<MTLComputePipelineState> pso = d->Pipeline(EXCITATION_SOURCE, s.serial ? "excitation_apply_serial" : "excitation_apply");
	id<MTLComputeCommandEncoder> enc = d->Encoder();
	[enc setComputePipelineState:pso];
	[enc setBuffer:field offset:0 atIndex:0];
	[enc setBuffer:s.index offset:0 atIndex:1];
	[enc setBuffer:s.amp offset:0 atIndex:2];
	[enc setBuffer:s.delay offset:0 atIndex:3];
	[enc setBuffer:s.signal offset:0 atIndex:4];
	[enc setBytes:&param length:sizeof(param) atIndex:5];
	d->Dispatch(pso, s.serial ? 1 : s.count);
}

GPU_Extension* Metal_CreateExt_Excitation(GPU_Backend_Metal::Impl* d, Engine_Extension* eng_ext, Engine* eng)
{
	if (!dynamic_cast<Engine_Ext_Excitation*>(eng_ext))
		return NULL;
	Operator_Ext_Excitation* op_ext = dynamic_cast<Operator_Ext_Excitation*>(eng_ext->GetOperatorExtension());
	if (!op_ext)
		return NULL;
	return new Metal_Ext_Excitation(d, op_ext, eng);
}
