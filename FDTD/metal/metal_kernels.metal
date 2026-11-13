/*
*	Copyright (C) 2026 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

// The kernels of the Metal backend: the entry points of the shared code (kernels/cell_update.h,
// kernels/ext_ops.h), which the build puts in front of this file (metal_cm_source.h, see
// FDTD/CMakeLists.txt), and the copy of the snapshots and the field energy. The structs of the
// shared code hold device pointers (MTLBuffer.gpuAddress); a kernel copies them from the constant
// buffer into its own memory before calling their functions. The node indices are 64 bit.

//! The grid of a kernel: the lines, and the first x line of the box (temporal blocking; the others 0)
struct CM_Grid
{
	uint nx, ny, nz, x0;
};

#define CM_UPDATE_KERNEL(name, update, Index) \
kernel void name(device float* X [[buffer(0)]], device const float* Y [[buffer(1)]], device const Index* index [[buffer(2)]], \
                 constant CM_Field& field [[buffer(3)]], constant CM_Grid& G [[buffer(4)]], uint3 gid [[thread_position_in_grid]]) \
{ \
	if ((gid.x>=G.nz) || (gid.y>=G.ny)) \
		return; \
	const CM_Field F = field; \
	update<Index>(X, Y, index, F, G.nx, G.ny, G.nz, G.x0 + gid.z, gid.y, gid.x); \
}

CM_UPDATE_KERNEL(cm_voltages16, cm_update_voltages, ushort)
CM_UPDATE_KERNEL(cm_voltages32, cm_update_voltages, uint)
CM_UPDATE_KERNEL(cm_currents16, cm_update_currents, ushort)
CM_UPDATE_KERNEL(cm_currents32, cm_update_currents, uint)

//! The points [r.x, r.y) of a point operation, in parallel; and in order (points sharing a field value)
#define CM_POINT_KERNELS(name, Op) \
kernel void name(constant Op& op [[buffer(0)]], constant ulong2& r [[buffer(1)]], uint gid [[thread_position_in_grid]]) \
{ \
	const ulong n = r.x + gid; \
	if (n>=r.y) \
		return; \
	const Op o = op; \
	o(n); \
} \
kernel void name##_serial(constant Op& op [[buffer(0)]], constant ulong2& r [[buffer(1)]], uint gid [[thread_position_in_grid]]) \
{ \
	if (gid!=0) \
		return; \
	const Op o = op; \
	for (ulong n=r.x; n<r.y; ++n) \
		o(n); \
}

#define CM_OP_KERNELS(kind, Op, name) CM_POINT_KERNELS(name, Op)
CM_POINT_OPS(CM_OP_KERNELS)

//! A copy of the fields for a snapshot (GPU_Backend_Metal::SnapshotFields()): a thread per value, the grid (z, y, n*nx + x)
kernel void copy_fields(const device float* volt [[buffer(0)]], const device float* curr [[buffer(1)]],
                        device float* snap_volt [[buffer(2)]], device float* snap_curr [[buffer(3)]],
                        constant CM_Grid& G [[buffer(4)]], uint3 gid [[thread_position_in_grid]])
{
	if ((gid.x>=G.nz) || (gid.y>=G.ny))
		return;
	const ulong i = ((ulong)gid.z*G.ny + gid.y)*G.nz + gid.x;
	snap_volt[i] = volt[i];
	snap_curr[i] = curr[i];
}

//! The field energy: the squared voltages and currents of the nodes below L (the lines of the fields N), summed along x in
//! float, a thread per (z, y) line; the host sums the lines in double (GPU_Backend_Metal::CalcFastEnergy())
kernel void energy_lines(const device float* volt [[buffer(0)]], const device float* curr [[buffer(1)]], device float2* partial [[buffer(2)]],
                         constant CM_Grid& N [[buffer(3)]], constant CM_Grid& L [[buffer(4)]], uint2 gid [[thread_position_in_grid]])
{
	const uint z = gid.x, y = gid.y;
	if ((z>=L.nz) || (y>=L.ny))
		return;
	const ulong sn = (ulong)N.nx*N.ny*N.nz;
	float E = 0, H = 0;
	for (uint x=0; x<L.nx; ++x)
	{
		const ulong i = ((ulong)x*N.ny + y)*N.nz + z;
		E += volt[i]*volt[i];
		E += volt[sn+i]*volt[sn+i];
		E += volt[2*sn+i]*volt[2*sn+i];
		H += curr[i]*curr[i];
		H += curr[sn+i]*curr[sn+i];
		H += curr[2*sn+i]*curr[2*sn+i];
	}
	partial[(ulong)y*L.nz + z] = float2(E, H);
}
