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

#include "metal_internal.h"
#include "FDTD/engine.h"
#include "FDTD/excitation.h"
#include "FDTD/extensions/operator_ext_tfsf.h"
#include "FDTD/extensions/engine_ext_tfsf.h"

// Total-field/scattered-field plane wave, see Engine_Ext_TFSF.
// One dispatch per active box plane (planes share edges, so they run one after the
// other in the host order), one thread per (i,j) of the plane, gid = (j, i).
// The host evaluates part of the sum in double precision, this runs in float.
static const char* TFSF_SOURCE = R"MSL(
struct TFSFParam
{
	uint numTS; uint length; int period;   // signal lookup
	uint n, nP, nPP;                       // plane normal and the two field components
	uint line, startP, startPP;            // plane position
	uint ni, nj;                           // plane size (along nP, nPP)
};

// signal index for a delay of k timesteps, same as Engine_Ext_TFSF::m_DelayLookup
inline uint tfsf_lookup(uint k, constant TFSFParam& P)
{
	uint v;
	if (P.numTS < k)
		v = 0;
	else if ((P.numTS-k >= P.length) && (P.period==0))
		v = 0;
	else
		v = P.numTS - k;
	if (P.period>0)
		v = v % (uint)P.period;
	return v;
}

kernel void tfsf_apply(device float* field         [[buffer(0)]],
                       const device uint* delay    [[buffer(1)]],
                       const device float* delta   [[buffer(2)]],
                       const device float* amp     [[buffer(3)]],
                       const device float* signal  [[buffer(4)]],
                       constant GridDim& N         [[buffer(5)]],
                       constant TFSFParam& P       [[buffer(6)]],
                       uint2 gid [[thread_position_in_grid]])
{
	const uint i = gid.y, j = gid.x;
	uint pos[3];
	pos[P.n] = P.line;
	pos[P.nP] = P.startP + i;
	pos[P.nPP] = P.startPP + j;
	const uint count = P.ni*P.nj;
	for (uint c=0; c<2; ++c)
	{
		const uint k = c*count + i*P.nj + j;
		const uint g = nijk(N, (c==0) ? P.nP : P.nPP, pos[0], pos[1], pos[2]);
		field[g] = field[g]
		         + (1.0f-delta[k])*amp[k]*signal[tfsf_lookup(delay[k], P)]
		         +       delta[k] *amp[k]*signal[tfsf_lookup(1+delay[k], P)];
	}
}
)MSL";

class Metal_Ext_TFSF : public GPU_Extension
{
public:
	Metal_Ext_TFSF(GPU_Backend_Metal::Impl* impl, Operator_Ext_TFSF* op_ext, Engine* eng);

	// the voltages get the current signal and vice versa (see Engine_Ext_TFSF)
	virtual void DoPostVoltageUpdates() {Apply(m_VoltPlanes, d->volt, m_CurrSignal);}
	virtual void DoPostCurrentUpdates() {Apply(m_CurrPlanes, d->curr, m_VoltSignal);}

protected:
	struct Plane
	{
		struct {uint32_t numTS, length; int32_t period; uint32_t n, nP, nPP, line, startP, startPP, ni, nj;} param;
		id<MTLBuffer> delay, delta, amp;   //!< [component][i*nj+j]
	};

	void SetupPlanes(std::vector<Plane>& planes, bool currents, Operator_Ext_TFSF* op_ext);
	void Apply(std::vector<Plane>& planes, id<MTLBuffer> field, id<MTLBuffer> signal);

	GPU_Backend_Metal::Impl* d;
	Operator_Ext_TFSF* m_Op_TFSF;
	Engine* m_Eng;

	std::vector<Plane> m_VoltPlanes, m_CurrPlanes;
	id<MTLBuffer> m_VoltSignal, m_CurrSignal;
};

Metal_Ext_TFSF::Metal_Ext_TFSF(GPU_Backend_Metal::Impl* impl, Operator_Ext_TFSF* op_ext, Engine* eng)
{
	d = impl;
	m_Op_TFSF = op_ext;
	m_Eng = eng;
	const size_t sig_bytes = op_ext->m_Exc->GetLength()*sizeof(float);
	m_VoltSignal = d->NewBuffer(sig_bytes, op_ext->m_Exc->GetVoltageSignal());
	m_CurrSignal = d->NewBuffer(sig_bytes, op_ext->m_Exc->GetCurrentSignal());
	SetupPlanes(m_VoltPlanes, false, op_ext);
	SetupPlanes(m_CurrPlanes, true, op_ext);
	d->Pipeline(TFSF_SOURCE, "tfsf_apply");
}

void Metal_Ext_TFSF::SetupPlanes(std::vector<Plane>& planes, bool currents, Operator_Ext_TFSF* op_ext)
{
	for (int n=0; n<3; ++n)
	{
		const int nP = (n+1)%3;
		const int nPP = (n+2)%3;
		for (int l=0; l<2; ++l)
		{
			if (!op_ext->m_ActiveDir[n][l])
				continue;
			Plane pl;
			pl.param.n = n;
			pl.param.nP = nP;
			pl.param.nPP = nPP;
			// voltages on the box faces, currents half a cell outside the lower face
			if (l==0)
				pl.param.line = currents ? op_ext->m_Start[n]-1 : op_ext->m_Start[n];
			else
				pl.param.line = op_ext->m_Stop[n];
			pl.param.startP = op_ext->m_Start[nP];
			pl.param.startPP = op_ext->m_Start[nPP];
			pl.param.ni = op_ext->m_numLines[nP];
			pl.param.nj = op_ext->m_numLines[nPP];

			const size_t count = (size_t)pl.param.ni*pl.param.nj;
			std::vector<uint32_t> delay(2*count);
			std::vector<float> delta(2*count), amp(2*count);
			for (int c=0; c<2; ++c)
			{
				const unsigned int* src_delay = currents ? op_ext->m_CurrDelay[n][l][c] : op_ext->m_VoltDelay[n][l][c];
				const FDTD_FLOAT* src_delta = currents ? op_ext->m_CurrDelayDelta[n][l][c] : op_ext->m_VoltDelayDelta[n][l][c];
				const FDTD_FLOAT* src_amp = currents ? op_ext->m_CurrAmp[n][l][c] : op_ext->m_VoltAmp[n][l][c];
				for (size_t k=0; k<count; ++k)
				{
					delay[c*count+k] = src_delay[k];
					delta[c*count+k] = src_delta[k];
					amp[c*count+k] = src_amp[k];
				}
			}
			pl.delay = d->NewBuffer(delay.size()*sizeof(uint32_t), delay.data());
			pl.delta = d->NewBuffer(delta.size()*sizeof(float), delta.data());
			pl.amp = d->NewBuffer(amp.size()*sizeof(float), amp.data());
			planes.push_back(pl);
		}
	}
}

void Metal_Ext_TFSF::Apply(std::vector<Plane>& planes, id<MTLBuffer> field, id<MTLBuffer> signal)
{
	const unsigned int numTS = m_Eng->GetNumberOfTimesteps();
	const int period = int(m_Op_TFSF->m_Exc->GetSignalPeriod()/m_Op_TFSF->m_Exc->GetTimestep());
	id<MTLComputePipelineState> pso = d->Pipeline(TFSF_SOURCE, "tfsf_apply");
	for (size_t p=0; p<planes.size(); ++p)
	{
		Plane& pl = planes[p];
		pl.param.numTS = numTS;
		pl.param.length = m_Op_TFSF->m_Exc->GetLength();
		pl.param.period = period;
		id<MTLComputeCommandEncoder> enc = d->Encoder();
		[enc setComputePipelineState:pso];
		[enc setBuffer:field offset:0 atIndex:0];
		[enc setBuffer:pl.delay offset:0 atIndex:1];
		[enc setBuffer:pl.delta offset:0 atIndex:2];
		[enc setBuffer:pl.amp offset:0 atIndex:3];
		[enc setBuffer:signal offset:0 atIndex:4];
		d->SetGridDim(5);
		[enc setBytes:&pl.param length:sizeof(pl.param) atIndex:6];
		d->Dispatch(pso, pl.param.nj, pl.param.ni);
	}
}

GPU_Extension* Metal_CreateExt_TFSF(GPU_Backend_Metal::Impl* d, Engine_Extension* eng_ext, Engine* eng)
{
	if (!dynamic_cast<Engine_Ext_TFSF*>(eng_ext))
		return NULL;
	Operator_Ext_TFSF* op_ext = dynamic_cast<Operator_Ext_TFSF*>(eng_ext->GetOperatorExtension());
	if (!op_ext)
		return NULL;
	return new Metal_Ext_TFSF(d, op_ext, eng);
}
