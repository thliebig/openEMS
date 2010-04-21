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

//! \brief construct an Engine_sse instance
//! it's the responsibility of the caller to free the returned pointer
Engine_sse* Engine_sse::New(const Operator_sse* op)
{
	Engine_sse* e = new Engine_sse(op);
	e->Init();
	return e;
}

Engine_sse::Engine_sse(const Operator_sse* op) : Engine(op)
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
	curr_ = Create_N_3DArray_v4sf(numLines);
	volt = 0; // not used
	curr = 0; // not used
}

void Engine_sse::Reset()
{
	Delete_N_3DArray_v4sf(volt_,numLines);
	volt_ = 0;
	Delete_N_3DArray_v4sf(curr_,numLines);
	curr_ = 0;
}

void Engine_sse::UpdateVoltages()
{
	unsigned int pos[3];
	bool shift[2];
	f4vector temp;

	for (pos[0]=0;pos[0]<numLines[0];++pos[0])
	{
		shift[0]=pos[0];
		for (pos[1]=0;pos[1]<numLines[1];++pos[1])
		{
			shift[1]=pos[1];
			for (pos[2]=0;pos[2]<numLines[2]/4;++pos[2])
			{
				// x-polarization
				temp.f[0] = curr_[1][pos[0]][pos[1]][pos[2]-(bool)pos[2]].f[3];
				temp.f[1] = curr_[1][pos[0]][pos[1]][pos[2]].f[0];
				temp.f[2] = curr_[1][pos[0]][pos[1]][pos[2]].f[1];
				temp.f[3] = curr_[1][pos[0]][pos[1]][pos[2]].f[2];
				volt_[0][pos[0]][pos[1]][pos[2]].v *= Op->vv_[0][pos[0]][pos[1]][pos[2]].v;
				volt_[0][pos[0]][pos[1]][pos[2]].v += Op->vi_[0][pos[0]][pos[1]][pos[2]].v * ( curr_[2][pos[0]][pos[1]][pos[2]].v - curr_[2][pos[0]][pos[1]-shift[1]][pos[2]].v - curr_[1][pos[0]][pos[1]][pos[2]].v + temp.v );

				// y-polarization
				temp.f[0] = curr_[0][pos[0]][pos[1]][pos[2]-(bool)pos[2]].f[3];
				temp.f[1] = curr_[0][pos[0]][pos[1]][pos[2]].f[0];
				temp.f[2] = curr_[0][pos[0]][pos[1]][pos[2]].f[1];
				temp.f[3] = curr_[0][pos[0]][pos[1]][pos[2]].f[2];
				volt_[1][pos[0]][pos[1]][pos[2]].v *= Op->vv_[1][pos[0]][pos[1]][pos[2]].v;
				volt_[1][pos[0]][pos[1]][pos[2]].v += Op->vi_[1][pos[0]][pos[1]][pos[2]].v * ( curr_[0][pos[0]][pos[1]][pos[2]].v - temp.v - curr_[2][pos[0]][pos[1]][pos[2]].v + curr_[2][pos[0]-shift[0]][pos[1]][pos[2]].v);

				// z-polarization
				volt_[2][pos[0]][pos[1]][pos[2]].v *= Op->vv_[2][pos[0]][pos[1]][pos[2]].v;
				volt_[2][pos[0]][pos[1]][pos[2]].v += Op->vi_[2][pos[0]][pos[1]][pos[2]].v * ( curr_[1][pos[0]][pos[1]][pos[2]].v - curr_[1][pos[0]-shift[0]][pos[1]][pos[2]].v - curr_[0][pos[0]][pos[1]][pos[2]].v + curr_[0][pos[0]][pos[1]-shift[1]][pos[2]].v);
			}
		}
	}
}

void Engine_sse::ApplyVoltageExcite()
{
	int exc_pos;
	unsigned int pos;
	//soft voltage excitation here (E-field excite)
	for (unsigned int n=0;n<Op->E_Exc_Count;++n)
	{
		exc_pos = (int)numTS - (int)Op->E_Exc_delay[n];
		exc_pos *= (exc_pos>0 && exc_pos<=(int)Op->ExciteLength);
		pos = Op->E_Exc_index[2][n];
		volt_[Op->E_Exc_dir[n]][Op->E_Exc_index[0][n]][Op->E_Exc_index[1][n]][pos/4].f[pos%4] += Op->E_Exc_amp[n]*Op->ExciteSignal[exc_pos];
	}
}

void Engine_sse::UpdateCurrents()
{
	unsigned int pos[5];
	f4vector temp;

	for (pos[0]=0;pos[0]<numLines[0]-1;++pos[0])
	{
		for (pos[1]=0;pos[1]<numLines[1]-1;++pos[1])
		{
			for (pos[2]=0;pos[2]<numLines[2]/4;++pos[2]) // FIXME is this correct?
			{
				// x-pol
				temp.f[0] = volt_[1][pos[0]][pos[1]][pos[2]].f[1];
				temp.f[1] = volt_[1][pos[0]][pos[1]][pos[2]].f[2];
				temp.f[2] = volt_[1][pos[0]][pos[1]][pos[2]].f[3];
				temp.f[3] = volt_[1][pos[0]][pos[1]][pos[2]+1].f[0]; // FIXME outside sim area
				curr_[0][pos[0]][pos[1]][pos[2]].v *= Op->ii_[0][pos[0]][pos[1]][pos[2]].v;
				curr_[0][pos[0]][pos[1]][pos[2]].v += Op->iv_[0][pos[0]][pos[1]][pos[2]].v * ( volt_[2][pos[0]][pos[1]][pos[2]].v - volt_[2][pos[0]][pos[1]+1][pos[2]].v - volt_[1][pos[0]][pos[1]][pos[2]].v + temp.v);

				// y-pol
				temp.f[0] = volt_[0][pos[0]][pos[1]][pos[2]].f[1];
				temp.f[1] = volt_[0][pos[0]][pos[1]][pos[2]].f[2];
				temp.f[2] = volt_[0][pos[0]][pos[1]][pos[2]].f[3];
				temp.f[3] = volt_[0][pos[0]][pos[1]][pos[2]+1].f[0]; // FIXME outside sim area
				curr_[1][pos[0]][pos[1]][pos[2]].v *= Op->ii_[1][pos[0]][pos[1]][pos[2]].v;
				curr_[1][pos[0]][pos[1]][pos[2]].v += Op->iv_[1][pos[0]][pos[1]][pos[2]].v * ( volt_[0][pos[0]][pos[1]][pos[2]].v - temp.v - volt_[2][pos[0]][pos[1]][pos[2]].v + volt_[2][pos[0]+1][pos[1]][pos[2]].v);

				// z-pol
				curr_[2][pos[0]][pos[1]][pos[2]].v *= Op->ii_[2][pos[0]][pos[1]][pos[2]].v;
				curr_[2][pos[0]][pos[1]][pos[2]].v += Op->iv_[2][pos[0]][pos[1]][pos[2]].v * ( volt_[1][pos[0]][pos[1]][pos[2]].v - volt_[1][pos[0]+1][pos[1]][pos[2]].v - volt_[0][pos[0]][pos[1]][pos[2]].v + volt_[0][pos[0]][pos[1]+1][pos[2]].v);
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
