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
#include "FDTD/extensions/operator_ext_cylinder.h"
#include "FDTD/extensions/engine_ext_cylinder.h"

// Closed alpha direction and r=0 axis of cylindrical meshes, see Engine_Ext_Cylinder.
// The mesh directions are (r, alpha, z) = (x, y, z) of the grid.
static const char* CYLINDER_SOURCE = R"MSL(
// Ez on the axis: one thread per z line, summing the alpha currents in the host order
kernel void cyl_r0_ez(device float* volt          [[buffer(0)]],
                      const device float* curr    [[buffer(1)]],
                      const device float* vv_R0   [[buffer(2)]],
                      const device float* vi_R0   [[buffer(3)]],
                      constant GridDim& N         [[buffer(4)]],
                      uint z [[thread_position_in_grid]])
{
	if (z>=N.nz)
		return;
	float v = volt[nijk(N, 2, 0, 0, z)]*vv_R0[z];
	for (uint a=0; a<N.ny-1; ++a)
		v = v + vi_R0[z]*curr[nijk(N, 1, 0, a, z)];
	volt[nijk(N, 2, 0, 0, z)] = v;
}

// on the axis: no alpha voltage, the same Ez for all alpha; gid = (z, alpha)
kernel void cyl_r0_spread(device float* volt   [[buffer(0)]],
                          constant GridDim& N  [[buffer(4)]],
                          uint2 gid [[thread_position_in_grid]])
{
	const uint z = gid.x, a = gid.y;
	volt[nijk(N, 1, 0, a, z)] = 0;
	if (a>0)
		volt[nijk(N, 2, 0, a, z)] = volt[nijk(N, 2, 0, 0, z)];
}

// copy the tangential voltages from the last alpha plane to the first; gid = (z, r)
kernel void cyl_close_volt(device float* volt   [[buffer(0)]],
                           constant GridDim& N  [[buffer(4)]],
                           uint2 gid [[thread_position_in_grid]])
{
	const uint z = gid.x, r = gid.y, last = N.ny-2;
	volt[nijk(N, 0, r, 0, z)] = volt[nijk(N, 0, r, last, z)];
	volt[nijk(N, 2, r, 0, z)] = volt[nijk(N, 2, r, last, z)];
}

// copy the tangential currents from the first alpha plane to the last; gid = (z, r)
kernel void cyl_close_curr(device float* curr   [[buffer(0)]],
                           constant GridDim& N  [[buffer(4)]],
                           uint2 gid [[thread_position_in_grid]])
{
	const uint z = gid.x, r = gid.y, last = N.ny-2;
	curr[nijk(N, 0, r, last, z)] = curr[nijk(N, 0, r, 0, z)];
	curr[nijk(N, 2, r, last, z)] = curr[nijk(N, 2, r, 0, z)];
}
)MSL";

class Metal_Ext_Cylinder : public GPU_Extension
{
public:
	Metal_Ext_Cylinder(GPU_Backend_Metal::Impl* impl, Operator_Ext_Cylinder* op_ext);

	virtual void DoPostVoltageUpdates();
	virtual void DoPostCurrentUpdates();

protected:
	id<MTLComputeCommandEncoder> Begin(const char* kernel, id<MTLBuffer> field);

	GPU_Backend_Metal::Impl* d;
	bool m_ClosedAlpha;
	bool m_R0Included;
	id<MTLBuffer> m_vv_R0, m_vi_R0;
};

Metal_Ext_Cylinder::Metal_Ext_Cylinder(GPU_Backend_Metal::Impl* impl, Operator_Ext_Cylinder* op_ext)
{
	d = impl;
	m_ClosedAlpha = op_ext->CC_closedAlpha;
	m_R0Included = op_ext->CC_R0_included;
	if (m_ClosedAlpha && m_R0Included)
	{
		const size_t bytes = d->dim.nz*sizeof(FDTD_FLOAT);
		m_vv_R0 = d->NewBuffer(bytes, op_ext->vv_R0);
		m_vi_R0 = d->NewBuffer(bytes, op_ext->vi_R0);
	}
	d->Pipeline(CYLINDER_SOURCE, "cyl_r0_ez");
	d->Pipeline(CYLINDER_SOURCE, "cyl_r0_spread");
	d->Pipeline(CYLINDER_SOURCE, "cyl_close_volt");
	d->Pipeline(CYLINDER_SOURCE, "cyl_close_curr");
}

id<MTLComputeCommandEncoder> Metal_Ext_Cylinder::Begin(const char* kernel, id<MTLBuffer> field)
{
	id<MTLComputeCommandEncoder> enc = d->Encoder();
	[enc setComputePipelineState:d->Pipeline(CYLINDER_SOURCE, kernel)];
	[enc setBuffer:field offset:0 atIndex:0];
	d->SetGridDim(4);
	return enc;
}

void Metal_Ext_Cylinder::DoPostVoltageUpdates()
{
	if (!m_ClosedAlpha)
		return;
	if (m_R0Included)
	{
		id<MTLComputeCommandEncoder> enc = Begin("cyl_r0_ez", d->volt);
		[enc setBuffer:d->curr offset:0 atIndex:1];
		[enc setBuffer:m_vv_R0 offset:0 atIndex:2];
		[enc setBuffer:m_vi_R0 offset:0 atIndex:3];
		d->Dispatch(d->Pipeline(CYLINDER_SOURCE, "cyl_r0_ez"), d->dim.nz);

		Begin("cyl_r0_spread", d->volt);
		d->Dispatch(d->Pipeline(CYLINDER_SOURCE, "cyl_r0_spread"), d->dim.nz, d->dim.ny);
	}
	Begin("cyl_close_volt", d->volt);
	d->Dispatch(d->Pipeline(CYLINDER_SOURCE, "cyl_close_volt"), d->dim.nz, d->dim.nx);
}

void Metal_Ext_Cylinder::DoPostCurrentUpdates()
{
	if (!m_ClosedAlpha)
		return;
	Begin("cyl_close_curr", d->curr);
	d->Dispatch(d->Pipeline(CYLINDER_SOURCE, "cyl_close_curr"), d->dim.nz-1, d->dim.nx-1);
}

GPU_Extension* Metal_CreateExt_Cylinder(GPU_Backend_Metal::Impl* d, Engine_Extension* eng_ext, Engine* eng)
{
	UNUSED(eng);
	if (!dynamic_cast<Engine_Ext_Cylinder*>(eng_ext))
		return NULL;
	Operator_Ext_Cylinder* op_ext = dynamic_cast<Operator_Ext_Cylinder*>(eng_ext->GetOperatorExtension());
	if (!op_ext)
		return NULL;
	return new Metal_Ext_Cylinder(d, op_ext);
}
