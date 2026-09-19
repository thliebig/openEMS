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

#include <stdexcept>

#include "engine_gpu.h"
#include "gpu_backend.h"
#include "extensions/engine_extension.h"
#include "extensions/operator_extension.h"

using std::cout;
using std::cerr;
using std::endl;

//! \brief construct an Engine_GPU instance
//! it's the responsibility of the caller to free the returned pointer
Engine_GPU* Engine_GPU::New(const Operator* op, const std::string& backend, GPU_Backend* parent)
{
	Engine_GPU* e = new Engine_GPU(op, backend, parent);
	e->Init();
	return e;
}

Engine_GPU::Engine_GPU(const Operator* op, const std::string& backend, GPU_Backend* parent) : Engine(op)
{
	m_BackendName = backend;
	m_ParentBackend = parent;
	m_type = GPU;
	m_Backend = NULL;
	m_FieldsOnHost = true;
	m_SharedMemory = false;
	m_StaleVolt.stale = m_StaleCurr.stale = false;
	m_StaleVolt.full_this_batch = m_StaleCurr.full_this_batch = false;
	m_StaleVolt.full_last_batch = m_StaleCurr.full_last_batch = false;
}

Engine_GPU::~Engine_GPU()
{
	Reset();
}

void Engine_GPU::Init()
{
	// allocates the host mirror and creates the extensions
	Engine::Init();

	// a sub-grid shares the device work stream of its parent grid
	m_Backend = m_ParentBackend ? m_ParentBackend->NewSubGridBackend() : GPU_Backend::New(m_BackendName);
	cout << "Create FDTD engine (GPU, backend: " << m_Backend->GetName() << ")" << endl;
	if (!m_Backend->Init(Op))
		throw std::runtime_error("Engine_GPU::Init: GPU backend initialization failed");

	// with unified memory the host mirror is a view of the device fields
	m_SharedMemory = m_Backend->GetSharedVoltages() && m_Backend->GetSharedCurrents();
	if (m_SharedMemory)
	{
		volt_ptr->InitView("volt", {numLines[0], numLines[1], numLines[2]}, m_Backend->GetSharedVoltages());
		curr_ptr->InitView("curr", {numLines[0], numLines[1], numLines[2]}, m_Backend->GetSharedCurrents());
	}

	// the fields stay on the device only if every extension has a device implementation
	std::vector<std::string> host_exts;
	for (size_t n=0; n<m_Eng_exts.size(); ++n)
	{
		Engine_Extension* eng_ext = m_Eng_exts.at(n);
		GPU_Extension* gpu_ext = NULL;
		if (eng_ext->GetOperatorExtension()->IsGPUSave())
			gpu_ext = m_Backend->CreateExtension(eng_ext, this);
		if (gpu_ext)
			m_GPU_exts.push_back(gpu_ext);
		else
			host_exts.push_back(eng_ext->GetExtensionName());
	}

	m_FieldsOnHost = !host_exts.empty();
	if (m_FieldsOnHost)
	{
		ClearGPUExtensions();
		cout << "Engine_GPU: extensions run on the host copy of the fields (host fallback)" << endl;
		cerr << "Engine_GPU::Init: Warning: extensions without a GPU implementation run on the host, the fields are copied every timestep:" << endl;
		for (size_t n=0; n<host_exts.size(); ++n)
			cerr << "\t" << host_exts.at(n) << endl;
	}
	else
		cout << "Engine_GPU: all extensions run on the device" << endl;
}

void Engine_GPU::ForceHostFallback()
{
	if (m_FieldsOnHost)
		return;
	ClearGPUExtensions();
	m_FieldsOnHost = true;
}

void Engine_GPU::ClearGPUExtensions()
{
	for (size_t n=0; n<m_GPU_exts.size(); ++n)
		delete m_GPU_exts.at(n);
	m_GPU_exts.clear();
}

void Engine_GPU::Reset()
{
	// the device extensions refer to the backend and to the engine extensions
	ClearGPUExtensions();
	delete m_Backend;
	m_Backend = NULL;
	Engine::Reset();
}

void Engine_GPU::VoltagesToHost()
{
	if (m_SharedMemory)
		m_Backend->Synchronize();
	else
		m_Backend->DownloadVoltages(*volt_ptr);
	m_StaleVolt.stale = false;
	m_StaleVolt.lines.clear();
}

void Engine_GPU::CurrentsToHost()
{
	if (m_SharedMemory)
		m_Backend->Synchronize();
	else
		m_Backend->DownloadCurrents(*curr_ptr);
	m_StaleCurr.stale = false;
	m_StaleCurr.lines.clear();
}

FDTD_FLOAT Engine_GPU::GetVolt(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const
{
	if (m_StaleVolt.stale)
		return ReadStale(false, n, x, y, z);
	return Engine::GetVolt(n, x, y, z);
}

FDTD_FLOAT Engine_GPU::GetCurr(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const
{
	if (m_StaleCurr.stale)
		return ReadStale(true, n, x, y, z);
	return Engine::GetCurr(n, x, y, z);
}

// Reading a z-line takes about the latency of a transfer, the whole field (24 bytes
// per node) a lot longer. Probes read a few lines per batch, dumps all of them.
static const size_t MAX_STALE_LINES = 256;

FDTD_FLOAT Engine_GPU::ReadStale(bool currents, unsigned int n, unsigned int x, unsigned int y, unsigned int z) const
{
	Engine_GPU* self = const_cast<Engine_GPU*>(this);   // the host mirror is a cache of the device fields
	StaleField& field = currents ? m_StaleCurr : m_StaleVolt;
	const unsigned int key = (n*numLines[0] + x)*numLines[1] + y;
	std::unordered_map<unsigned int, std::vector<FDTD_FLOAT>>::const_iterator it = field.lines.find(key);
	if (it!=field.lines.end())
		return it->second.at(z);

	if (!field.full_last_batch && (field.lines.size()<MAX_STALE_LINES))
	{
		std::vector<FDTD_FLOAT> line(numLines[2]);
		if (m_Backend->DownloadRange(currents, (size_t)key*numLines[2], numLines[2], line.data()))
		{
			const FDTD_FLOAT value = line.at(z);
			field.lines[key].swap(line);
			return value;
		}
	}

	// many lines (or no range download): the whole field
	if (currents)
		self->CurrentsToHost();
	else
		self->VoltagesToHost();
	field.full_this_batch = true;
	return currents ? Engine::GetCurr(n, x, y, z) : Engine::GetVolt(n, x, y, z);
}

void Engine_GPU::MarkStale(StaleField& field)
{
	// a batch that read the whole field is likely followed by another one (e.g. dumps)
	field.full_last_batch = field.full_this_batch;
	field.full_this_batch = false;
	field.stale = true;
	field.lines.clear();
}

bool Engine_GPU::SnapshotFields(unsigned int slot, const FDTD_FLOAT* &volt, const FDTD_FLOAT* &curr)
{
	// only where a copy on the device is cheaper than reading the fields (shared memory)
	if (m_FieldsOnHost || !m_SharedMemory)
		return false;
	return m_Backend->SnapshotFields(slot, volt, curr);
}

void Engine_GPU::UpdateHostMirror()
{
	if (m_StaleVolt.stale)
		VoltagesToHost();
	if (m_StaleCurr.stale)
		CurrentsToHost();
}

// with shared memory the host only writes while the device is idle (after a
// Synchronize()), and the writes are visible to all work committed later
void Engine_GPU::VoltagesToDevice()
{
	if (!m_SharedMemory)
		m_Backend->UploadVoltages(*volt_ptr);
}

void Engine_GPU::CurrentsToDevice()
{
	if (!m_SharedMemory)
		m_Backend->UploadCurrents(*curr_ptr);
}

bool Engine_GPU::IterateTS(unsigned int iterTS)
{
	for (unsigned int iter=0; iter<iterTS; ++iter)
	{
		VoltageHalfStep();
		CurrentHalfStep();
		NextTimestep();
	}
	FinishBatch();
	return true;
}

// The fast path runs all extensions on the device, in the same order as the host
// extensions (see Engine).
// Host fallback: any hook may write either field, including the pre-update hooks
// (the UPML swaps its flux into the fields before the main update), so both
// fields are uploaded before and the updated one downloaded after each main update.

void Engine_GPU::VoltageHalfStep()
{
	if (!m_FieldsOnHost)
	{
		for (int n=m_GPU_exts.size()-1; n>=0; --n)
			m_GPU_exts.at(n)->DoPreVoltageUpdates();
		m_Backend->UpdateVoltages();
		for (size_t n=0; n<m_GPU_exts.size(); ++n)
			m_GPU_exts.at(n)->DoPostVoltageUpdates();
		for (size_t n=0; n<m_GPU_exts.size(); ++n)
			m_GPU_exts.at(n)->Apply2Voltages();
		return;
	}

	DoPreVoltageUpdates();
	VoltagesToDevice();
	CurrentsToDevice();
	m_Backend->UpdateVoltages();
	VoltagesToHost();
	DoPostVoltageUpdates();
	Apply2Voltages();
}

void Engine_GPU::CurrentHalfStep()
{
	if (!m_FieldsOnHost)
	{
		for (int n=m_GPU_exts.size()-1; n>=0; --n)
			m_GPU_exts.at(n)->DoPreCurrentUpdates();
		m_Backend->UpdateCurrents();
		for (size_t n=0; n<m_GPU_exts.size(); ++n)
			m_GPU_exts.at(n)->DoPostCurrentUpdates();
		for (size_t n=0; n<m_GPU_exts.size(); ++n)
			m_GPU_exts.at(n)->Apply2Current();
		return;
	}

	DoPreCurrentUpdates();
	VoltagesToDevice();
	CurrentsToDevice();
	m_Backend->UpdateCurrents();
	CurrentsToHost();
	DoPostCurrentUpdates();
	Apply2Current();
}

bool Engine_GPU::CalcFastEnergy(const unsigned int numNodes[3], double& E_energy, double& H_energy) const
{
	// in the host fallback the host mirror holds the latest fields
	if (m_FieldsOnHost)
		return false;
	return m_Backend->CalcFastEnergy(numNodes, E_energy, H_energy);
}

void Engine_GPU::NextTimestep()
{
	++numTS;
}

void Engine_GPU::FinishBatch()
{
	if (m_FieldsOnHost)
		return;   // the host mirror is up to date after every half-step
	if (m_SharedMemory)
	{
		// the host mirror is the device memory
		VoltagesToHost();
		CurrentsToHost();
	}
	else
	{
		// read on demand, see GetVolt()
		MarkStale(m_StaleVolt);
		MarkStale(m_StaleCurr);
	}
	for (size_t n=0; n<m_GPU_exts.size(); ++n)
		m_GPU_exts.at(n)->Synchronize();
}
