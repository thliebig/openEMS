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

// The memory and the launches of the device extensions of every backend on the cell model
// (FDTD/gpu_ext.h): one kernel template runs any point operation of FDTD/kernels/ext_ops.h.

#include "hip/hip_runtime.h"

#include "hip_internal.h"
#include "FDTD/kernels/ext_ops.h"

template <typename Op>
__global__ void run_points(Op op, uint64_t begin, uint64_t end)
{
	const uint64_t n = begin + (uint64_t)blockIdx.x*blockDim.x + threadIdx.x;
	if (n<end)
		op(n);
}

// points that share an edge, in their order
template <typename Op>
__global__ void run_points_serial(Op op, uint64_t begin, uint64_t end)
{
	if (blockIdx.x*blockDim.x + threadIdx.x != 0)
		return;
	for (uint64_t n=begin; n<end; ++n)
		op(n);
}

template <typename Op>
static void Run(GPU_Backend_HIP::Impl* d, const char* name, const void* op, uint64_t begin, uint64_t end, bool serial)
{
	const Op& o = *(const Op*)op;
	if (serial)
		HIP_Launch(d, name, run_points_serial<Op>, 1, 1, 1, o, begin, end);
	else
		HIP_Launch(d, name, run_points<Op>, end-begin, 1, 1, o, begin, end);
}

void* GPU_Backend_HIP::DeviceAlloc(size_t bytes, const void* host)
{
	return d->Alloc<char>(bytes, (const char*)host);
}

void GPU_Backend_HIP::DeviceDownload(void* host, const void* dev, size_t bytes)
{
	HIP_Check(hipMemcpyAsync(host, dev, bytes, hipMemcpyDeviceToHost, d->Stream()), "download");
	d->Flush();
}

float* GPU_Backend_HIP::DeviceField(bool currents)
{
	return currents ? d->curr : d->volt;
}

void GPU_Backend_HIP::RunPoints(GPU_PointOp kind, const void* op, uint64_t begin, uint64_t end, bool serial)
{
	if (end<=begin)
		return;
	switch (kind)
	{
#define RUN(K, T, name) case K: Run<T>(d, #name, op, begin, end, serial); break;
	CM_POINT_OPS(RUN)
#undef RUN
	}
}
