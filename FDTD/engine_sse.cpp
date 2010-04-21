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

#include "engine_sse.h"
#include "tools/array_ops.h"

//! \brief construct an Engine_sse instance
//! it's the responsibility of the caller to free the returned pointer
Engine_sse* Engine_sse::New(const Operator* op)
{
	Engine_sse* e = new Engine_sse(op);
	e->Init();
	return e;
}

Engine_sse::Engine_sse(const Operator* op) : Engine(op)
{
	Op = op;
	for (int n=0;n<3;++n)
	{
		numLines[n] = Op->GetNumberOfLines(n);
	}
}

Engine_sse::~Engine_sse()
{
	this->Reset();
}

void Engine_sse::Init()
{
	numTS = 0;
	volt_ = Create_N_3DArray_v4sf(numLines);
	curr = Create_N_3DArray(numLines);
}

void Engine_sse::Reset()
{
	Delete_N_3DArray_v4sf(volt_,numLines);
	volt=NULL;
	Delete_N_3DArray(curr,numLines);
	curr=NULL;
}

void Engine_sse::UpdateVoltages()
{
	unsigned int pos[4];
	bool shift[3];

	//voltage updates
	for (pos[0]=0;pos[0]<numLines[0];++pos[0])
	{
		shift[0]=pos[0];
		for (pos[1]=0;pos[1]<numLines[1];++pos[1])
		{
			shift[1]=pos[1];
			for (pos[2]=0;pos[2]<numLines[2]/4;++pos[2])
			{
				//do the updates here
				//for x
				volt_[0][pos[0]][pos[1]][pos[2]].v *= Op->vv_[0][pos[0]][pos[1]][pos[2]].v;
				volt_[1][pos[0]][pos[1]][pos[2]].v *= Op->vv_[1][pos[0]][pos[1]][pos[2]].v;
				volt_[2][pos[0]][pos[1]][pos[2]].v *= Op->vv_[2][pos[0]][pos[1]][pos[2]].v;


				for (pos[3]=0;pos[3]<4;++pos[3]) {
					shift[2]=pos[2]+pos[3];
					volt_[0][pos[0]][pos[1]][pos[2]].f[pos[3]] += Op->vi_[0][pos[0]][pos[1]][pos[2]].f[pos[3]] * ( curr[2][pos[0]][pos[1]][pos[2]] - curr[2][pos[0]][pos[1]-shift[1]][pos[2]] - curr[1][pos[0]][pos[1]][pos[2]] + curr[1][pos[0]][pos[1]][pos[2]-shift[2]]);
					volt_[1][pos[0]][pos[1]][pos[2]].f[pos[3]] += Op->vi_[1][pos[0]][pos[1]][pos[2]].f[pos[3]] * ( curr[0][pos[0]][pos[1]][pos[2]] - curr[0][pos[0]][pos[1]][pos[2]-shift[2]] - curr[2][pos[0]][pos[1]][pos[2]] + curr[2][pos[0]-shift[0]][pos[1]][pos[2]]);
					volt_[2][pos[0]][pos[1]][pos[2]].f[pos[3]] += Op->vi_[2][pos[0]][pos[1]][pos[2]].f[pos[3]] * ( curr[1][pos[0]][pos[1]][pos[2]] - curr[1][pos[0]-shift[0]][pos[1]][pos[2]] - curr[0][pos[0]][pos[1]][pos[2]] + curr[0][pos[0]][pos[1]-shift[1]][pos[2]]);
				}
			}
		}
	}
}

void Engine_sse::ApplyVoltageExcite()
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
}

void Engine_sse::UpdateCurrents()
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

void Engine_sse::ApplyCurrentExcite()
{
}

bool Engine_sse::IterateTS(unsigned int iterTS)
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
