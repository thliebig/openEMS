/*
*	Copyright (C) 2010 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#include "engine.h"
#include "extensions/engine_extension.h"
#include "extensions/operator_extension.h"
#include "tools/array_ops.h"

//! \brief construct an Engine instance
//! it's the responsibility of the caller to free the returned pointer
Engine* Engine::New(const Operator* op)
{
	cout << "Create FDTD engine" << endl;
	Engine* e = new Engine(op);
	e->Init();
	return e;
}

Engine::Engine(const Operator* op)
{
	m_type = BASIC;
	numTS = 0;
	Op = op;
	for (int n=0; n<3; ++n)
		numLines[n] = Op->GetNumberOfLines(n, true);
	volt=NULL;
	curr=NULL;
}

Engine::~Engine()
{
	this->Reset();
}

void Engine::Init()
{
	numTS = 0;
	volt = Create_N_3DArray<FDTD_FLOAT>(numLines);
	curr = Create_N_3DArray<FDTD_FLOAT>(numLines);

	InitExtensions();
	SortExtensionByPriority();
}

void Engine::InitExtensions()
{
	for (size_t n=0; n<Op->GetNumberOfExtentions(); ++n)
	{
		Operator_Extension* op_ext = Op->GetExtension(n);
		Engine_Extension* eng_ext = op_ext->CreateEngineExtention();
		if (eng_ext)
		{
			eng_ext->SetEngine(this);
			m_Eng_exts.push_back(eng_ext);
		}
	}
}

void Engine::ClearExtensions()
{
	for (size_t n=0; n<m_Eng_exts.size(); ++n)
		delete m_Eng_exts.at(n);
	m_Eng_exts.clear();
}

bool CompareExtensions(Engine_Extension* i, Engine_Extension* j)
{
	return (*i<*j);
}

void Engine::SortExtensionByPriority()
{
	stable_sort(m_Eng_exts.begin(),m_Eng_exts.end(), CompareExtensions);
	reverse(m_Eng_exts.begin(),m_Eng_exts.end());

	if (g_settings.GetVerboseLevel()>1)
	{
		cout << "---  Engine::SortExtensionByPriority() ---" << endl;
		for (size_t n=0; n<m_Eng_exts.size(); ++n)
			cout << " #" << n << ": " << m_Eng_exts.at(n)->GetExtensionName() << " (" << m_Eng_exts.at(n)->GetPriority() << ")" << endl;
	}
}

void Engine::Reset()
{
	Delete_N_3DArray(volt,numLines);
	volt=NULL;
	Delete_N_3DArray(curr,numLines);
	curr=NULL;

	ClearExtensions();
}

[[gnu::target("default")]]
#include "update_voltages.inc"

[[gnu::target("avx2,fma")]]
#include "update_voltages.inc"

[[gnu::target("avx512f")]]
#include "update_voltages.inc"


[[gnu::target("default")]]
#include "update_currents.inc"

[[gnu::target("avx2,fma")]]
#include "update_currents.inc"

[[gnu::target("avx512f")]]
#include "update_currents.inc"

void Engine::DoPreVoltageUpdates()
{
	//execute extensions in reverse order -> highest priority gets access to the voltages last
	for (int n=m_Eng_exts.size()-1; n>=0; --n)
		m_Eng_exts.at(n)->DoPreVoltageUpdates();

}

void Engine::DoPostVoltageUpdates()
{
	//execute extensions in normal order -> highest priority gets access to the voltages first
	for (size_t n=0; n<m_Eng_exts.size(); ++n)
		m_Eng_exts.at(n)->DoPostVoltageUpdates();
}

void Engine::Apply2Voltages()
{
	//execute extensions in normal order -> highest priority gets access to the voltages first
	for (size_t n=0; n<m_Eng_exts.size(); ++n)
		m_Eng_exts.at(n)->Apply2Voltages();
}

void Engine::DoPreCurrentUpdates()
{
	//execute extensions in reverse order -> highest priority gets access to the currents last
	for (int n=m_Eng_exts.size()-1; n>=0; --n)
		m_Eng_exts.at(n)->DoPreCurrentUpdates();
}

void Engine::DoPostCurrentUpdates()
{
	//execute extensions in normal order -> highest priority gets access to the currents first
	for (size_t n=0; n<m_Eng_exts.size(); ++n)
		m_Eng_exts.at(n)->DoPostCurrentUpdates();
}

void Engine::Apply2Current()
{
	//execute extensions in normal order -> highest priority gets access to the currents first
	for (size_t n=0; n<m_Eng_exts.size(); ++n)
		m_Eng_exts.at(n)->Apply2Current();
}

bool Engine::IterateTS(unsigned int iterTS)
{
	for (unsigned int iter=0; iter<iterTS; ++iter)
	{
		//voltage updates with extensions
		DoPreVoltageUpdates();
		UpdateVoltages(0,numLines[0]);
		DoPostVoltageUpdates();
		Apply2Voltages();

		//current updates with extensions
		DoPreCurrentUpdates();
		UpdateCurrents(0,numLines[0]-1);
		DoPostCurrentUpdates();
		Apply2Current();

		++numTS;
	}
	return true;
}
