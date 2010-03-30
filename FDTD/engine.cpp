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
Engine* Engine::createEngine(Operator* op)
{
	Engine* e = new Engine(op);
	e->Init();
	return e;
}

Engine::Engine(Operator* op)
{
	Op = op;
}

Engine::~Engine()
{
	Reset();
}

void Engine::Init()
{
	numTS = 0;
	volt = Create_N_3DArray(Op->numLines);
	curr = Create_N_3DArray(Op->numLines);
}

void Engine::Reset()
{
	Delete_N_3DArray(volt,Op->numLines);
	volt=NULL;
	Delete_N_3DArray(curr,Op->numLines);
	curr=NULL;
}

bool Engine::IterateTS(unsigned int iterTS)
{
	unsigned int pos[3];
	int exc_pos;
	bool shift[3];

	for (unsigned int iter=0;iter<iterTS;++iter)
	{
		//voltage updates
		for (pos[0]=0;pos[0]<Op->numLines[0];++pos[0])
		{
			shift[0]=pos[0];
			for (pos[1]=0;pos[1]<Op->numLines[1];++pos[1])
			{
				shift[1]=pos[1];
				for (pos[2]=0;pos[2]<Op->numLines[2];++pos[2])
				{
					shift[2]=pos[2];
					//do the updates here
					//for x
					volt[0][pos[0]][pos[1]][pos[2]] *= Op->vv[0][pos[0]][pos[1]][pos[2]];
					volt[0][pos[0]][pos[1]][pos[2]] += Op->vi[0][pos[0]][pos[1]][pos[2]] * ( curr[2][pos[0]][pos[1]][pos[2]] - curr[2][pos[0]][pos[1]-shift[1]][pos[2]] - curr[1][pos[0]][pos[1]][pos[2]] + curr[1][pos[0]][pos[1]][pos[2]-shift[2]]);

					//for y
					volt[1][pos[0]][pos[1]][pos[2]] *= Op->vv[1][pos[0]][pos[1]][pos[2]];
					volt[1][pos[0]][pos[1]][pos[2]] += Op->vi[1][pos[0]][pos[1]][pos[2]] * ( curr[0][pos[0]][pos[1]][pos[2]] - curr[0][pos[0]][pos[1]][pos[2]-shift[2]] - curr[2][pos[0]][pos[1]][pos[2]] + curr[2][pos[0]-shift[0]][pos[1]][pos[2]]);

					//for x
					volt[2][pos[0]][pos[1]][pos[2]] *= Op->vv[2][pos[0]][pos[1]][pos[2]];
					volt[2][pos[0]][pos[1]][pos[2]] += Op->vi[2][pos[0]][pos[1]][pos[2]] * ( curr[1][pos[0]][pos[1]][pos[2]] - curr[1][pos[0]-shift[0]][pos[1]][pos[2]] - curr[0][pos[0]][pos[1]][pos[2]] + curr[0][pos[0]][pos[1]-shift[1]][pos[2]]);
				}
			}
		}

		//soft voltage excitation here (E-field excite)
		for (unsigned int n=0;n<Op->E_Exc_Count;++n)
		{
			exc_pos = (int)numTS - (int)Op->E_Exc_delay[n];
			exc_pos *= (exc_pos>0 && exc_pos<=(int)Op->ExciteLength);
//			if (n==0) cerr << numTS << " => " << Op->ExciteSignal[exc_pos] << endl;
			volt[Op->E_Exc_dir[n]][Op->E_Exc_index[0][n]][Op->E_Exc_index[1][n]][Op->E_Exc_index[2][n]] += Op->E_Exc_amp[n]*Op->ExciteSignal[exc_pos];
		}

		//current updates
		for (pos[0]=0;pos[0]<Op->numLines[0]-1;++pos[0])
		{
			for (pos[1]=0;pos[1]<Op->numLines[1]-1;++pos[1])
			{
				for (pos[2]=0;pos[2]<Op->numLines[2]-1;++pos[2])
				{
					//do the updates here
					//for x
					curr[0][pos[0]][pos[1]][pos[2]] *= Op->ii[0][pos[0]][pos[1]][pos[2]];
					curr[0][pos[0]][pos[1]][pos[2]] += Op->iv[0][pos[0]][pos[1]][pos[2]] * ( volt[2][pos[0]][pos[1]][pos[2]] - volt[2][pos[0]][pos[1]+1][pos[2]] - volt[1][pos[0]][pos[1]][pos[2]] + volt[1][pos[0]][pos[1]][pos[2]+1]);

					//for y
					curr[1][pos[0]][pos[1]][pos[2]] *= Op->ii[1][pos[0]][pos[1]][pos[2]];
					curr[1][pos[0]][pos[1]][pos[2]] += Op->iv[1][pos[0]][pos[1]][pos[2]] * ( volt[0][pos[0]][pos[1]][pos[2]] - volt[0][pos[0]][pos[1]][pos[2]+1] - volt[2][pos[0]][pos[1]][pos[2]] + volt[2][pos[0]+1][pos[1]][pos[2]]);

					//for x
					curr[2][pos[0]][pos[1]][pos[2]] *= Op->ii[2][pos[0]][pos[1]][pos[2]];
					curr[2][pos[0]][pos[1]][pos[2]] += Op->iv[2][pos[0]][pos[1]][pos[2]] * ( volt[1][pos[0]][pos[1]][pos[2]] - volt[1][pos[0]+1][pos[1]][pos[2]] - volt[0][pos[0]][pos[1]][pos[2]] + volt[0][pos[0]][pos[1]+1][pos[2]]);
				}
			}
		}

		//soft current excitation here (H-field excite)
		++numTS;
	}
	return true;
}
