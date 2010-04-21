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
	{
		numLines[n] = Op->GetNumberOfLines(n);
	}
}

Engine::~Engine()
{
	this->Reset();
}

void Engine::Init()
{
	numTS = 0;
	volt = Create_N_3DArray(numLines);
	curr = Create_N_3DArray(numLines);

	file_et1.open( "et1" );
}

void Engine::Reset()
{
	Delete_N_3DArray(volt,numLines);
	volt=NULL;
	Delete_N_3DArray(curr,numLines);
	curr=NULL;

	file_et1.close();
}

void Engine::UpdateVoltages()
{
	unsigned int pos[3];
	bool shift[3];

	//voltage updates
	for (pos[0]=0;pos[0]<numLines[0];++pos[0])
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
	}
}

void Engine::ApplyVoltageExcite()
{
	int exc_pos;
	//soft voltage excitation here (E-field excite)
	for (unsigned int n=0;n<Op->E_Exc_Count;++n)
	{
		exc_pos = (int)numTS - (int)Op->E_Exc_delay[n];
		exc_pos *= (exc_pos>0 && exc_pos<=(int)Op->ExciteLength);
//			if (n==0) cerr << numTS << " => " << Op->ExciteSignal[exc_pos] << endl;
		volt[Op->E_Exc_dir[n]][Op->E_Exc_index[0][n]][Op->E_Exc_index[1][n]][Op->E_Exc_index[2][n]] += Op->E_Exc_amp[n]*Op->ExciteSignal[exc_pos];
	}

	// write the first excitation into the file "et1"
	if (Op->E_Exc_Count >= 1) {
		exc_pos = (int)numTS - (int)Op->E_Exc_delay[0];
		exc_pos *= (exc_pos>0 && exc_pos<=(int)Op->ExciteLength);
		file_et1 << numTS * Op->GetTimestep() << "\t" << Op->E_Exc_amp[0]*Op->ExciteSignal[exc_pos] << "\n"; // do not use std::endl here, because it will do an implicit flush
	}
}

void Engine::UpdateCurrents()
{
	unsigned int pos[3];
	for (pos[0]=0;pos[0]<numLines[0]-1;++pos[0])
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
	}
}

void Engine::ApplyCurrentExcite()
{
}

bool Engine::IterateTS(unsigned int iterTS)
{
	for (unsigned int iter=0;iter<iterTS;++iter)
	{
		UpdateVoltages();
		ApplyVoltageExcite();
		UpdateCurrents();
		ApplyCurrentExcite();
		++numTS;
	}
	return true;
}
