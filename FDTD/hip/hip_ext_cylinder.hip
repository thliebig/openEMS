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
#include "hip_internal.h"
#include "FDTD/engine.h"
#include "FDTD/extensions/operator_ext_cylinder.h"
#include "FDTD/extensions/engine_extension.h"
// not engine_ext_cylinder.h: it pulls in Boost.Thread (through operator_cylinder.h), whose
// headers the device compiler cannot compile with MSVC and Boost 1.91, so Engine_Ext_Cylinder is not a type
// here and the extension is detected through its operator extension instead (see
// HIP_CreateExt_Cylinder())

// Closed alpha direction and r=0 axis of cylindrical meshes, see Engine_Ext_Cylinder.
// The mesh directions are (r, alpha, z) = (x, y, z) of the grid.

// Ez on the axis: one thread per z line, summing the alpha currents in the host order
__global__ void cyl_r0_ez(float* volt, const float* curr, const float* vv_R0, const float* vi_R0, HIP_GridDim N)
{
	const unsigned int z = blockIdx.x*blockDim.x + threadIdx.x;
	if (z>=N.nz)
		return;
	float v = volt[nijk(N, 2, 0, 0, z)]*vv_R0[z];
	for (unsigned int a=0; a<N.ny-1; ++a)
		v = v + vi_R0[z]*curr[nijk(N, 1, 0, a, z)];
	volt[nijk(N, 2, 0, 0, z)] = v;
}

// on the axis: no alpha voltage, the same Ez for all alpha
__global__ void cyl_r0_spread(float* volt, HIP_GridDim N)
{
	const unsigned int z = blockIdx.x*blockDim.x + threadIdx.x;
	const unsigned int a = blockIdx.y*blockDim.y + threadIdx.y;
	if (z>=N.nz || a>=N.ny)
		return;
	volt[nijk(N, 1, 0, a, z)] = 0;
	if (a>0)
		volt[nijk(N, 2, 0, a, z)] = volt[nijk(N, 2, 0, 0, z)];
}

// copy the tangential voltages from the last alpha plane to the first
__global__ void cyl_close_volt(float* volt, HIP_GridDim N)
{
	const unsigned int z = blockIdx.x*blockDim.x + threadIdx.x;
	const unsigned int r = blockIdx.y*blockDim.y + threadIdx.y;
	if (z>=N.nz || r>=N.nx)
		return;
	const unsigned int last = N.ny-2;
	volt[nijk(N, 0, r, 0, z)] = volt[nijk(N, 0, r, last, z)];
	volt[nijk(N, 2, r, 0, z)] = volt[nijk(N, 2, r, last, z)];
}

// copy the tangential currents from the first alpha plane to the last
__global__ void cyl_close_curr(float* curr, HIP_GridDim N)
{
	const unsigned int z = blockIdx.x*blockDim.x + threadIdx.x;
	const unsigned int r = blockIdx.y*blockDim.y + threadIdx.y;
	if (z+1>=N.nz || r+1>=N.nx)
		return;
	const unsigned int last = N.ny-2;
	curr[nijk(N, 0, r, last, z)] = curr[nijk(N, 0, r, 0, z)];
	curr[nijk(N, 2, r, last, z)] = curr[nijk(N, 2, r, 0, z)];
}

class HIP_Ext_Cylinder : public GPU_Extension
{
public:
	HIP_Ext_Cylinder(GPU_Backend_HIP::Impl* impl, Operator_Ext_Cylinder* op_ext);

	virtual void DoPostVoltageUpdates();
	virtual void DoPostCurrentUpdates();

protected:
	GPU_Backend_HIP::Impl* d;
	bool m_ClosedAlpha;
	bool m_R0Included;
	float *m_vv_R0, *m_vi_R0;
};

HIP_Ext_Cylinder::HIP_Ext_Cylinder(GPU_Backend_HIP::Impl* impl, Operator_Ext_Cylinder* op_ext)
{
	d = impl;
	m_ClosedAlpha = op_ext->CC_closedAlpha;
	m_R0Included = op_ext->CC_R0_included;
	m_vv_R0 = m_vi_R0 = NULL;
	if (m_ClosedAlpha && m_R0Included)
	{
		m_vv_R0 = d->Alloc<float>(d->dim.nz, op_ext->vv_R0);
		m_vi_R0 = d->Alloc<float>(d->dim.nz, op_ext->vi_R0);
	}
}

void HIP_Ext_Cylinder::DoPostVoltageUpdates()
{
	if (!m_ClosedAlpha)
		return;
	if (m_R0Included)
	{
		HIP_Launch(d, "cyl_r0_ez", cyl_r0_ez, d->dim.nz, 1, 1, d->volt, (const float*)d->curr, (const float*)m_vv_R0, (const float*)m_vi_R0, d->dim);
		HIP_Launch(d, "cyl_r0_spread", cyl_r0_spread, d->dim.nz, d->dim.ny, 1, d->volt, d->dim);
	}
	HIP_Launch(d, "cyl_close_volt", cyl_close_volt, d->dim.nz, d->dim.nx, 1, d->volt, d->dim);
}

void HIP_Ext_Cylinder::DoPostCurrentUpdates()
{
	if (!m_ClosedAlpha)
		return;
	HIP_Launch(d, "cyl_close_curr", cyl_close_curr, d->dim.nz, d->dim.nx, 1, d->curr, d->dim);
}

GPU_Extension* HIP_CreateExt_Cylinder(GPU_Backend_HIP::Impl* d, Engine_Extension* eng_ext, Engine* eng)
{
	UNUSED(eng);
	// no dynamic_cast to Engine_Ext_Cylinder here (see the includes); the engine extension of an
	// Operator_Ext_Cylinder is always one
	Operator_Ext_Cylinder* op_ext = dynamic_cast<Operator_Ext_Cylinder*>(eng_ext->GetOperatorExtension());
	if (!op_ext)
		return NULL;
	return new HIP_Ext_Cylinder(d, op_ext);
}
