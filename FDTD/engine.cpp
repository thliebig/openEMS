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
	volt_ptr = NULL;
	curr_ptr = NULL;
}

Engine::~Engine()
{
	this->Reset();
}

void Engine::Init()
{
	numTS = 0;
	volt_ptr = new ArrayLib::ArrayNIJK<FDTD_FLOAT>("volt", numLines);
	curr_ptr = new ArrayLib::ArrayNIJK<FDTD_FLOAT>("curr", numLines);

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
	delete volt_ptr;
	volt_ptr = NULL;
	delete curr_ptr;
	curr_ptr = NULL;

	ClearExtensions();
}

void Engine::UpdateVoltages(unsigned int startX, unsigned int numX)
{
	ArrayLib::ArrayNIJK<FDTD_FLOAT>& volt = *volt_ptr;
	ArrayLib::ArrayNIJK<FDTD_FLOAT>& curr = *curr_ptr;
	ArrayLib::ArrayNIJK<FDTD_FLOAT>& vv = *Op->vv_ptr;
	ArrayLib::ArrayNIJK<FDTD_FLOAT>& vi = *Op->vi_ptr;

	unsigned int pos[3];
	bool shift[3];

	pos[0] = startX;
	//voltage updates
	for (unsigned int posX=0; posX<numX; ++posX)
	{
		shift[0]=pos[0];
		for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
		{
			shift[1]=pos[1];
			for (pos[2]=0; pos[2]<numLines[2]; ++pos[2])
			{
				shift[2]=pos[2];
				//do the updates here
				//for x
				volt[0][pos[0]][pos[1]][pos[2]] *=
				    vv[0][pos[0]][pos[1]][pos[2]];
				volt[0][pos[0]][pos[1]][pos[2]] +=
				    vi[0][pos[0]][pos[1]][pos[2]] * (
				        curr[2][pos[0]][pos[1]         ][pos[2]         ] -
				        curr[2][pos[0]][pos[1]-shift[1]][pos[2]         ] -
				        curr[1][pos[0]][pos[1]         ][pos[2]         ] +
				        curr[1][pos[0]][pos[1]         ][pos[2]-shift[2]]
				    );

				//for y
				volt[1][pos[0]][pos[1]][pos[2]] *=
				    vv[1][pos[0]][pos[1]][pos[2]];
				volt[1][pos[0]][pos[1]][pos[2]] +=
				    vi[1][pos[0]][pos[1]][pos[2]] * (
				        curr[0][pos[0]         ][pos[1]][pos[2]         ] -
				        curr[0][pos[0]         ][pos[1]][pos[2]-shift[2]] -
				        curr[2][pos[0]         ][pos[1]][pos[2]         ] +
				        curr[2][pos[0]-shift[0]][pos[1]][pos[2]         ]
				    );

				//for z
				volt[2][pos[0]][pos[1]][pos[2]] *=
				    vv[2][pos[0]][pos[1]][pos[2]];
				volt[2][pos[0]][pos[1]][pos[2]] +=
				    vi[2][pos[0]][pos[1]][pos[2]] * (
				        curr[1][pos[0]         ][pos[1]         ][pos[2]] -
				        curr[1][pos[0]-shift[0]][pos[1]         ][pos[2]] -
				        curr[0][pos[0]         ][pos[1]         ][pos[2]] +
				        curr[0][pos[0]         ][pos[1]-shift[1]][pos[2]]
				    );
			}
		}
		++pos[0];
	}
}

void Engine::UpdateCurrents(unsigned int startX, unsigned int numX)
{
	ArrayLib::ArrayNIJK<FDTD_FLOAT>& curr = *curr_ptr;
	ArrayLib::ArrayNIJK<FDTD_FLOAT>& volt = *volt_ptr;
	ArrayLib::ArrayNIJK<FDTD_FLOAT>& ii = *Op->ii_ptr;
	ArrayLib::ArrayNIJK<FDTD_FLOAT>& iv = *Op->iv_ptr;

	unsigned int pos[3];
	pos[0] = startX;
	for (unsigned int posX=0; posX<numX; ++posX)
	{
		for (pos[1]=0; pos[1]<numLines[1]-1; ++pos[1])
		{
			for (pos[2]=0; pos[2]<numLines[2]-1; ++pos[2])
			{
				//do the updates here
				//for x
				curr[0][pos[0]][pos[1]][pos[2]] *=
				    ii[0][pos[0]][pos[1]][pos[2]];
				curr[0][pos[0]][pos[1]][pos[2]] +=
				    iv[0][pos[0]][pos[1]][pos[2]] * (
				        volt[2][pos[0]][pos[1]  ][pos[2]  ] -
				        volt[2][pos[0]][pos[1]+1][pos[2]  ] -
				        volt[1][pos[0]][pos[1]  ][pos[2]  ] +
				        volt[1][pos[0]][pos[1]  ][pos[2]+1]
				    );

				//for y
				curr[1][pos[0]][pos[1]][pos[2]] *=
				    ii[1][pos[0]][pos[1]][pos[2]];
				curr[1][pos[0]][pos[1]][pos[2]] +=
				    iv[1][pos[0]][pos[1]][pos[2]] * (
				        volt[0][pos[0]  ][pos[1]][pos[2]  ] -
				        volt[0][pos[0]  ][pos[1]][pos[2]+1] -
				        volt[2][pos[0]  ][pos[1]][pos[2]  ] +
				        volt[2][pos[0]+1][pos[1]][pos[2]  ]
				    );

				//for z
				curr[2][pos[0]][pos[1]][pos[2]] *=
				    ii[2][pos[0]][pos[1]][pos[2]];
				curr[2][pos[0]][pos[1]][pos[2]] +=
				    iv[2][pos[0]][pos[1]][pos[2]] * (
				        volt[1][pos[0]  ][pos[1]  ][pos[2]] -
				        volt[1][pos[0]+1][pos[1]  ][pos[2]] -
				        volt[0][pos[0]  ][pos[1]  ][pos[2]] +
				        volt[0][pos[0]  ][pos[1]+1][pos[2]]
				    );
			}
		}
		++pos[0];
	}
}

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
