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
}

void Engine_GPU::CurrentsToHost()
{
	if (m_SharedMemory)
		m_Backend->Synchronize();
	else
		m_Backend->DownloadCurrents(*curr_ptr);
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

void Engine_GPU::NextTimestep()
{
	++numTS;
}

void Engine_GPU::FinishBatch()
{
	if (m_FieldsOnHost)
		return;   // the host mirror is up to date after every half-step
	// update the host mirror for the field processing
	VoltagesToHost();
	CurrentsToHost();
	for (size_t n=0; n<m_GPU_exts.size(); ++n)
		m_GPU_exts.at(n)->Synchronize();
}
