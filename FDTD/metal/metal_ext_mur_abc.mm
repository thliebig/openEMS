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
#include "FDTD/extensions/operator_ext_mur_abc.h"
#include "FDTD/extensions/engine_ext_mur_abc.h"

// First order Mur ABC on one boundary plane, see Engine_Ext_Mur_ABC.
// One thread per (i,j) of the plane, i along m_nyP, j along m_nyPP; gid = (j, i).
static const char* MUR_SOURCE = R"MSL(
struct MurParam { int ny, nyP, nyPP; uint line, line_shift; uint ni, nj; };

inline uint mur_index(constant GridDim& N, constant MurParam& P, uint n, uint line, uint i, uint j)
{
	uint pos[3];
	pos[P.ny] = line;
	pos[P.nyP] = i;
	pos[P.nyPP] = j;
	return nijk(N, n, pos[0], pos[1], pos[2]);
}

kernel void mur_pre(const device float* volt   [[buffer(0)]],
                    device float* v_nyP        [[buffer(1)]],
                    device float* v_nyPP       [[buffer(2)]],
                    const device float* c_nyP  [[buffer(3)]],
                    const device float* c_nyPP [[buffer(4)]],
                    constant GridDim& N        [[buffer(5)]],
                    constant MurParam& P       [[buffer(6)]],
                    uint2 gid [[thread_position_in_grid]])
{
	const uint i = gid.y, j = gid.x;
	const uint ij = i*P.nj + j;
	v_nyP[ij]  = volt[mur_index(N, P, P.nyP,  P.line_shift, i, j)] - c_nyP[ij]  * volt[mur_index(N, P, P.nyP,  P.line, i, j)];
	v_nyPP[ij] = volt[mur_index(N, P, P.nyPP, P.line_shift, i, j)] - c_nyPP[ij] * volt[mur_index(N, P, P.nyPP, P.line, i, j)];
}

kernel void mur_post(const device float* volt   [[buffer(0)]],
                     device float* v_nyP        [[buffer(1)]],
                     device float* v_nyPP       [[buffer(2)]],
                     const device float* c_nyP  [[buffer(3)]],
                     const device float* c_nyPP [[buffer(4)]],
                     constant GridDim& N        [[buffer(5)]],
                     constant MurParam& P       [[buffer(6)]],
                     uint2 gid [[thread_position_in_grid]])
{
	const uint i = gid.y, j = gid.x;
	const uint ij = i*P.nj + j;
	v_nyP[ij]  = v_nyP[ij]  + c_nyP[ij]  * volt[mur_index(N, P, P.nyP,  P.line_shift, i, j)];
	v_nyPP[ij] = v_nyPP[ij] + c_nyPP[ij] * volt[mur_index(N, P, P.nyPP, P.line_shift, i, j)];
}

kernel void mur_apply(device float* volt          [[buffer(0)]],
                      const device float* v_nyP   [[buffer(1)]],
                      const device float* v_nyPP  [[buffer(2)]],
                      constant GridDim& N         [[buffer(5)]],
                      constant MurParam& P        [[buffer(6)]],
                      uint2 gid [[thread_position_in_grid]])
{
	const uint i = gid.y, j = gid.x;
	const uint ij = i*P.nj + j;
	volt[mur_index(N, P, P.nyP,  P.line, i, j)] = v_nyP[ij];
	volt[mur_index(N, P, P.nyPP, P.line, i, j)] = v_nyPP[ij];
}
)MSL";

class Metal_Ext_Mur_ABC : public GPU_Extension
{
public:
	Metal_Ext_Mur_ABC(GPU_Backend_Metal::Impl* impl, Operator_Ext_Mur_ABC* op_ext, Engine_Ext_Mur_ABC* eng_ext, Engine* eng);

	virtual void DoPreVoltageUpdates()  {Run("mur_pre");}
	virtual void DoPostVoltageUpdates() {Run("mur_post");}
	virtual void Apply2Voltages()       {Run("mur_apply");}

protected:
	void Run(const char* kernel);

	GPU_Backend_Metal::Impl* d;
	Engine* m_Eng;
	unsigned int m_StartTS;   //!< the ABC is off until an excitation on its plane is done
	struct {int32_t ny, nyP, nyPP; uint32_t line, line_shift; uint32_t ni, nj;} m_Param;

	id<MTLBuffer> m_Volt_nyP, m_Volt_nyPP;
	id<MTLBuffer> m_Coeff_nyP, m_Coeff_nyPP;
};

Metal_Ext_Mur_ABC::Metal_Ext_Mur_ABC(GPU_Backend_Metal::Impl* impl, Operator_Ext_Mur_ABC* op_ext, Engine_Ext_Mur_ABC* eng_ext, Engine* eng)
{
	d = impl;
	m_Eng = eng;
	m_StartTS = eng_ext->GetStartTimestep();
	m_Param.ny = op_ext->m_ny;
	m_Param.nyP = op_ext->m_nyP;
	m_Param.nyPP = op_ext->m_nyPP;
	m_Param.line = op_ext->m_LineNr;
	m_Param.line_shift = op_ext->m_LineNr_Shift;
	m_Param.ni = op_ext->m_numLines[0];
	m_Param.nj = op_ext->m_numLines[1];

	const size_t bytes = (size_t)m_Param.ni*m_Param.nj*sizeof(FDTD_FLOAT);
	m_Volt_nyP  = d->NewBuffer(bytes);
	m_Volt_nyPP = d->NewBuffer(bytes);
	m_Coeff_nyP  = d->NewBuffer(bytes, op_ext->m_Mur_Coeff_nyP.data());
	m_Coeff_nyPP = d->NewBuffer(bytes, op_ext->m_Mur_Coeff_nyPP.data());

	d->Pipeline(MUR_SOURCE, "mur_pre");
	d->Pipeline(MUR_SOURCE, "mur_post");
	d->Pipeline(MUR_SOURCE, "mur_apply");
}

void Metal_Ext_Mur_ABC::Run(const char* kernel)
{
	if (m_Eng->GetNumberOfTimesteps()<m_StartTS)
		return;
	id<MTLComputePipelineState> pso = d->Pipeline(MUR_SOURCE, kernel);
	id<MTLComputeCommandEncoder> enc = d->Encoder();
	[enc setComputePipelineState:pso];
	[enc setBuffer:d->volt offset:0 atIndex:0];
	[enc setBuffer:m_Volt_nyP offset:0 atIndex:1];
	[enc setBuffer:m_Volt_nyPP offset:0 atIndex:2];
	[enc setBuffer:m_Coeff_nyP offset:0 atIndex:3];
	[enc setBuffer:m_Coeff_nyPP offset:0 atIndex:4];
	d->SetGridDim(5);
	[enc setBytes:&m_Param length:sizeof(m_Param) atIndex:6];
	d->Dispatch(pso, m_Param.nj, m_Param.ni);
}

GPU_Extension* Metal_CreateExt_Mur_ABC(GPU_Backend_Metal::Impl* d, Engine_Extension* eng_ext, Engine* eng)
{
	Engine_Ext_Mur_ABC* mur_ext = dynamic_cast<Engine_Ext_Mur_ABC*>(eng_ext);
	if (!mur_ext)
		return NULL;
	Operator_Ext_Mur_ABC* op_ext = dynamic_cast<Operator_Ext_Mur_ABC*>(eng_ext->GetOperatorExtension());
	if (!op_ext)
		return NULL;
	return new Metal_Ext_Mur_ABC(d, op_ext, mur_ext, eng);
}
