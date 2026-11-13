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

#include <algorithm>
#include <cstring>
#include <stdexcept>

#include "gpu_backend_reference.h"
#include "operator.h"
#include "gpu_cell_model.h"
#include "kernels/ext_ops.h"
#include "engine_plan.h"
#include "tools/useful.h"

GPU_Backend_Reference::GPU_Backend_Reference()
{
	for (int n=0; n<3; ++n)
		numLines[n] = 0;
}

bool GPU_Backend_Reference::Init(const Operator* op, const GPU_CellModel& model)
{
	for (int n=0; n<3; ++n)
		numLines[n] = op->GetNumberOfLines(n, true);
	if ((model.numLines[0]!=numLines[0]) || (model.numLines[1]!=numLines[1]) || (model.numLines[2]!=numLines[2]))
		return false;
	m_volt.Init("gpu_volt", numLines);
	m_curr.Init("gpu_curr", numLines);

	// the "upload": a copy of the model, with its state
	m_Model.reset(new GPU_CellModel(model));
	for (int f=0; f<2; ++f)
	{
		const GPU_CellModel::Field& src = m_Model->field[f];
		CM_Field& F = m_Field[f];
		m_State[f].assign(src.numSlots*src.slotFloats, 0);
		F.combos = src.combos.data();
		F.slot = src.slot.data();
		F.state = m_State[f].data();
		F.numSlots = src.numSlots;
		F.poles = src.poles.data();
		for (int q=0; q<5; ++q)
			F.profile[q] = src.profile[q].data();
		for (int n=0; n<3; ++n)
		{
			F.gxy[n] = src.gxy[n].data();
			F.gz[n] = src.gz[n].data();
		}
	}
	return true;
}

// the kernels of the device backends (kernels/cell_update.h) on the x planes [x0, x1), one plane per item
void GPU_Backend_Reference::Update(bool voltages, unsigned int x0, unsigned int x1)
{
	const unsigned int nx = numLines[0], ny = numLines[1], nz = numLines[2];
	x1 = std::min(x1, nx);
	if (x0>=x1)
		return;
	float* V = m_volt.data();
	float* I = m_curr.data();
	const uint32_t* index = m_Model->field[voltages ? 0 : 1].index.data();
	const CM_Field& F = m_Field[voltages ? 0 : 1];
	ParallelFor(x1-x0, AvailableThreads(), [&](unsigned int, size_t begin, size_t end) {
		for (unsigned int x=x0+(unsigned int)begin; x<x0+end; ++x)
			for (unsigned int y=0; y<ny; ++y)
				for (unsigned int z=0; z<nz; ++z)
					if (voltages)
						cm_update_voltages<uint32_t>(V, I, index, F, nx, ny, nz, x, y, z);
					else
						cm_update_currents<uint32_t>(I, V, index, F, nx, ny, nz, x, y, z);
	});
}

void GPU_Backend_Reference::UpdateVoltages(const Engine_Box& box) {Update(true, box.start[0], box.stop[0]);}
void GPU_Backend_Reference::UpdateCurrents(const Engine_Box& box) {Update(false, box.start[0], box.stop[0]);}

static void CopyField(const ArrayLib::ArrayNIJK<FDTD_FLOAT>& src, ArrayLib::ArrayNIJK<FDTD_FLOAT>& dst)
{
	if (src.size()!=dst.size())
		throw std::runtime_error("GPU_Backend_Reference: host and device field size mismatch");
	std::copy_n(src.data(), src.size(), dst.data());
}

void GPU_Backend_Reference::Download(bool currents, ArrayLib::ArrayNIJK<FDTD_FLOAT>& field)
{
	CopyField(currents ? m_curr : m_volt, field);
}

void* GPU_Backend_Reference::DeviceAlloc(size_t bytes, const void* host)
{
	m_Memory.emplace_back(new char[std::max<size_t>(bytes, 1)]());
	if (host)
		std::memcpy(m_Memory.back().get(), host, bytes);
	return m_Memory.back().get();
}

void GPU_Backend_Reference::DeviceDownload(void* host, const void* dev, size_t bytes)
{
	std::memcpy(host, dev, bytes);
}

// as Engine_Interface_FDTD::CalcFastEnergy() on the host mirror, in the same order: the host mirror is stale during a batch
bool GPU_Backend_Reference::CalcFastEnergy(const unsigned int numNodes[3], double& E_energy, double& H_energy)
{
	E_energy = H_energy = 0;
	unsigned int pos[3];
	for (pos[0]=0; pos[0]<numNodes[0]; ++pos[0])
		for (pos[1]=0; pos[1]<numNodes[1]; ++pos[1])
			for (pos[2]=0; pos[2]<numNodes[2]; ++pos[2])
			{
				for (int n=0; n<3; ++n)
					E_energy += m_volt(n, pos[0], pos[1], pos[2]) * m_volt(n, pos[0], pos[1], pos[2]);
				for (int n=0; n<3; ++n)
					H_energy += m_curr(n, pos[0], pos[1], pos[2]) * m_curr(n, pos[0], pos[1], pos[2]);
			}
	return true;
}

namespace
{
//! the points [begin, end) of \a op on the host threads, or in order
template <typename Op>
void Run(const void* op, uint64_t begin, uint64_t end, bool serial)
{
	const Op& o = *(const Op*)op;
	if (serial || (end-begin<1024))
	{
		for (uint64_t n=begin; n<end; ++n)
			o(n);
		return;
	}
	ParallelFor(end-begin, AvailableThreads(), [&](unsigned int, size_t b, size_t e) {
		for (size_t n=b; n<e; ++n)
			o(begin+n);
	});
}
}

void GPU_Backend_Reference::RunPoints(GPU_PointOp kind, const void* op, uint64_t begin, uint64_t end, bool serial)
{
	switch (kind)
	{
#define RUN(K, T, name) case K: Run<T>(op, begin, end, serial); break;
	CM_POINT_OPS(RUN)
#undef RUN
	}
}
