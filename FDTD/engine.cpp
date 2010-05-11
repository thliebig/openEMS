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
#include "engine_extension.h"
#include "operator_extension.h"
#include "tools/array_ops.h"

//! \brief construct an Engine instance
//! it's the responsibility of the caller to free the returned pointer
Engine* Engine::New(const Operator* op)
{
	Engine* e = new Engine(op);
	e->Init();
	return e;
}

Engine::Engine(const Operator* op)
{
	Op = op;
	for (int n=0;n<3;++n)
		numLines[n] = Op->GetOriginalNumLines(n);
	volt=NULL;
	curr=NULL;
}

Engine::~Engine()
{
	this->Reset();
}

void Engine::Init()
{
	Reset();
	numTS = 0;
	volt = Create_N_3DArray(numLines);
	curr = Create_N_3DArray(numLines);

	file_et.open( "et" );

	InitExtensions();
}

void Engine::InitExtensions()
{
	for (size_t n=0;n<Op->GetNumberOfExtentions();++n)
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

void Engine::Reset()
{
	Delete_N_3DArray(volt,numLines);
	volt=NULL;
	Delete_N_3DArray(curr,numLines);
	curr=NULL;

	file_et.close();

	for (size_t n=0;n<m_Eng_exts.size();++n)
		delete m_Eng_exts.at(n);
	m_Eng_exts.clear();
}

void Engine::UpdateVoltages(unsigned int startX, unsigned int numX)
{
	unsigned int pos[3];
	bool shift[3];

	pos[0] = startX;
	//voltage updates
	for (unsigned int posX=0;posX<numX;++posX)
	{
		shift[0]=pos[0];
		for (pos[1]=0;pos[1]<numLines[1];++pos[1])
		{
			shift[1]=pos[1];
			for (pos[2]=0;pos[2]<numLines[2];++pos[2])
			{
				shift[2]=pos[2];
				//do the updates here
				//for x
				volt[0][pos[0]][pos[1]][pos[2]] *= Op->vv[0][pos[0]][pos[1]][pos[2]];
				volt[0][pos[0]][pos[1]][pos[2]] += Op->vi[0][pos[0]][pos[1]][pos[2]] * ( curr[2][pos[0]][pos[1]][pos[2]] - curr[2][pos[0]][pos[1]-shift[1]][pos[2]] - curr[1][pos[0]][pos[1]][pos[2]] + curr[1][pos[0]][pos[1]][pos[2]-shift[2]]);

				//for y
				volt[1][pos[0]][pos[1]][pos[2]] *= Op->vv[1][pos[0]][pos[1]][pos[2]];
				volt[1][pos[0]][pos[1]][pos[2]] += Op->vi[1][pos[0]][pos[1]][pos[2]] * ( curr[0][pos[0]][pos[1]][pos[2]] - curr[0][pos[0]][pos[1]][pos[2]-shift[2]] - curr[2][pos[0]][pos[1]][pos[2]] + curr[2][pos[0]-shift[0]][pos[1]][pos[2]]);

				//for z
				volt[2][pos[0]][pos[1]][pos[2]] *= Op->vv[2][pos[0]][pos[1]][pos[2]];
				volt[2][pos[0]][pos[1]][pos[2]] += Op->vi[2][pos[0]][pos[1]][pos[2]] * ( curr[1][pos[0]][pos[1]][pos[2]] - curr[1][pos[0]-shift[0]][pos[1]][pos[2]] - curr[0][pos[0]][pos[1]][pos[2]] + curr[0][pos[0]][pos[1]-shift[1]][pos[2]]);
			}
		}
		++pos[0];
	}
}

void Engine::ApplyVoltageExcite()
{
	int exc_pos;
	//soft voltage excitation here (E-field excite)
	for (unsigned int n=0;n<Op->Exc->E_Count;++n)
	{
		exc_pos = (int)numTS - (int)Op->Exc->E_delay[n];
		exc_pos *= (exc_pos>0 && exc_pos<=(int)Op->Exc->Length);
//			if (n==0) cerr << numTS << " => " << Op->ExciteSignal[exc_pos] << endl;
		GetVolt(Op->Exc->E_dir[n],Op->Exc->E_index[0][n],Op->Exc->E_index[1][n],Op->Exc->E_index[2][n]) += Op->Exc->E_amp[n]*Op->Exc->Signal_volt[exc_pos];
	}

	// write the first excitation into the file "et1"
	if (numTS < Op->Exc->Length)
		file_et << numTS * Op->GetTimestep() << "\t" << Op->Exc->Signal_volt[numTS] << "\n"; // do not use std::endl here, because it will do an implicit flush
	else
		file_et << numTS * Op->GetTimestep() << "\t0" << "\n"; // do not use std::endl here, because it will do an implicit flush
}

void Engine::UpdateCurrents(unsigned int startX, unsigned int numX)
{
	unsigned int pos[3];
	pos[0] = startX;
	for (unsigned int posX=0;posX<numX;++posX)
	{
		for (pos[1]=0;pos[1]<numLines[1]-1;++pos[1])
		{
			for (pos[2]=0;pos[2]<numLines[2]-1;++pos[2])
			{
				//do the updates here
				//for x
				curr[0][pos[0]][pos[1]][pos[2]] *= Op->ii[0][pos[0]][pos[1]][pos[2]];
				curr[0][pos[0]][pos[1]][pos[2]] += Op->iv[0][pos[0]][pos[1]][pos[2]] * ( volt[2][pos[0]][pos[1]][pos[2]] - volt[2][pos[0]][pos[1]+1][pos[2]] - volt[1][pos[0]][pos[1]][pos[2]] + volt[1][pos[0]][pos[1]][pos[2]+1]);

				//for y
				curr[1][pos[0]][pos[1]][pos[2]] *= Op->ii[1][pos[0]][pos[1]][pos[2]];
				curr[1][pos[0]][pos[1]][pos[2]] += Op->iv[1][pos[0]][pos[1]][pos[2]] * ( volt[0][pos[0]][pos[1]][pos[2]] - volt[0][pos[0]][pos[1]][pos[2]+1] - volt[2][pos[0]][pos[1]][pos[2]] + volt[2][pos[0]+1][pos[1]][pos[2]]);

				//for z
				curr[2][pos[0]][pos[1]][pos[2]] *= Op->ii[2][pos[0]][pos[1]][pos[2]];
				curr[2][pos[0]][pos[1]][pos[2]] += Op->iv[2][pos[0]][pos[1]][pos[2]] * ( volt[1][pos[0]][pos[1]][pos[2]] - volt[1][pos[0]+1][pos[1]][pos[2]] - volt[0][pos[0]][pos[1]][pos[2]] + volt[0][pos[0]][pos[1]+1][pos[2]]);
			}
		}
		++pos[0];
	}
}

void Engine::ApplyCurrentExcite()
{
	int exc_pos;
	//soft current excitation here (H-field excite)
	for (unsigned int n=0;n<Op->Exc->Curr_Count;++n)
	{
		exc_pos = (int)numTS - (int)Op->Exc->Curr_delay[n];
		exc_pos *= (exc_pos>0 && exc_pos<=(int)Op->Exc->Length);
//			if (n==0) cerr << numTS << " => " << Op->ExciteSignal[exc_pos] << endl;
		GetCurr(Op->Exc->Curr_dir[n],Op->Exc->Curr_index[0][n],Op->Exc->Curr_index[1][n],Op->Exc->Curr_index[2][n]) += Op->Exc->Curr_amp[n]*Op->Exc->Signal_curr[exc_pos];
	}
}

bool Engine::IterateTS(unsigned int iterTS)
{
	for (unsigned int iter=0;iter<iterTS;++iter)
	{
		//voltage updates with extensions
		for (size_t n=0;n<m_Eng_exts.size();++n)
			m_Eng_exts.at(n)->DoPreVoltageUpdates();

		UpdateVoltages(0,numLines[0]);

		for (size_t n=0;n<m_Eng_exts.size();++n)
			m_Eng_exts.at(n)->DoPostVoltageUpdates();
		for (size_t n=0;n<m_Eng_exts.size();++n)
			m_Eng_exts.at(n)->Apply2Voltages();

		ApplyVoltageExcite();

		//current updates with extensions
		for (size_t n=0;n<m_Eng_exts.size();++n)
			m_Eng_exts.at(n)->DoPreCurrentUpdates();

		UpdateCurrents(0,numLines[0]-1);

		for (size_t n=0;n<m_Eng_exts.size();++n)
			m_Eng_exts.at(n)->DoPostCurrentUpdates();
		for (size_t n=0;n<m_Eng_exts.size();++n)
			m_Eng_exts.at(n)->Apply2Current();

		ApplyCurrentExcite();

		++numTS;
	}
	return true;
}
