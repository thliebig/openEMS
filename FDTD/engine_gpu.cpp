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
Engine_GPU* Engine_GPU::New(const Operator_GPU* op)
{
	Engine_GPU* e = new Engine_GPU(op);
	e->Init();
	return e;
}

Engine_GPU::Engine_GPU(const Operator_GPU* op) : Engine(op)
{
	m_type = GPU;
	m_Backend = NULL;
	m_FieldsOnHost = true;
}

Engine_GPU::~Engine_GPU()
{
	Reset();
}

void Engine_GPU::Init()
{
	// allocates the host mirror and creates the extensions
	Engine::Init();

	m_Backend = GPU_Backend::New();
	cout << "Create FDTD engine (GPU, backend: " << m_Backend->GetName() << ")" << endl;
	if (!m_Backend->Init(Op))
		throw std::runtime_error("Engine_GPU::Init: GPU backend initialization failed");

	m_FieldsOnHost = false;
	for (size_t n=0; n<m_Eng_exts.size(); ++n)
	{
		if (m_Eng_exts.at(n)->GetOperatorExtension()->IsGPUSave())
			continue;
		if (!m_FieldsOnHost)
			cerr << "Engine_GPU::Init: Warning: extensions without a GPU implementation run on the host, the fields are copied every timestep:" << endl;
		cerr << "\t" << m_Eng_exts.at(n)->GetExtensionName() << endl;
		m_FieldsOnHost = true;
	}
}

void Engine_GPU::Reset()
{
	delete m_Backend;
	m_Backend = NULL;
	Engine::Reset();
}

bool Engine_GPU::IterateTS(unsigned int iterTS)
{
	ArrayLib::ArrayNIJK<FDTD_FLOAT>& volt = *volt_ptr;
	ArrayLib::ArrayNIJK<FDTD_FLOAT>& curr = *curr_ptr;

	// Host fallback: any hook may write either field, including the pre-update hooks
	// (the UPML swaps its flux into the fields before the main update), so both
	// fields are uploaded before and the updated one downloaded after each main update.
	for (unsigned int iter=0; iter<iterTS; ++iter)
	{
		//voltage updates with extensions
		DoPreVoltageUpdates();
		if (m_FieldsOnHost)
		{
			m_Backend->UploadVoltages(volt);
			m_Backend->UploadCurrents(curr);
		}
		m_Backend->UpdateVoltages();
		if (m_FieldsOnHost)
			m_Backend->DownloadVoltages(volt);
		DoPostVoltageUpdates();
		Apply2Voltages();

		//current updates with extensions
		DoPreCurrentUpdates();
		if (m_FieldsOnHost)
		{
			m_Backend->UploadVoltages(volt);
			m_Backend->UploadCurrents(curr);
		}
		m_Backend->UpdateCurrents();
		if (m_FieldsOnHost)
			m_Backend->DownloadCurrents(curr);
		DoPostCurrentUpdates();
		Apply2Current();

		++numTS;
	}

	// update the host mirror for the field processing
	if (!m_FieldsOnHost)
	{
		m_Backend->DownloadVoltages(volt);
		m_Backend->DownloadCurrents(curr);
	}
	return true;
}
