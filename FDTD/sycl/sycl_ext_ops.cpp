/*
*	Copyright (C) 2026 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

// The memory and the launches of the device extensions of every backend on the cell model
// (FDTD/gpu_ext.h): one kernel template runs any point operation of FDTD/kernels/ext_ops.h.

#include "sycl_internal.h"
#include "FDTD/kernels/ext_ops.h"

template <typename Op>
static void Run(GPU_Backend_SYCL::Impl* d, const char* name, const void* op, uint64_t begin, uint64_t end, bool serial, const CM_PointFilter& filter)
{
	const Op o = *(const Op*)op;
	const CM_PointFilter f = filter;
	if (serial)
	{
		// points that share an edge, in their order
		SYCL_Check(name, [&] {
			d->Queue().single_task([=]() {
				for (uint64_t n=begin; n<end; ++n)
					if (f(n))
						o(n);
			});
		});
	}
	else
		SYCL_Launch(d, name, end-begin, 1, 1, [=](size_t i, size_t, size_t) {
			const uint64_t n = begin + i;
			if (f(n))
				o(n);
		});
}

void* GPU_Backend_SYCL::DeviceAlloc(size_t bytes, const void* host)
{
	return d->Alloc<char>(bytes, (const char*)host);
}

void GPU_Backend_SYCL::DeviceDownload(void* host, const void* dev, size_t bytes)
{
	SYCL_Check("download", [&] {d->Queue().memcpy(host, dev, bytes);});
	d->Flush();
}

float* GPU_Backend_SYCL::DeviceField(bool currents)
{
	return currents ? d->curr : d->volt;
}

void GPU_Backend_SYCL::RunPoints(GPU_PointOp kind, const void* op, uint64_t begin, uint64_t end, bool serial, const CM_PointFilter& filter)
{
	if (end<=begin)
		return;
	switch (kind)
	{
#define RUN(K, T, name) case K: Run<T>(d, #name, op, begin, end, serial, filter); break;
	CM_POINT_OPS(RUN)
#undef RUN
	}
}
