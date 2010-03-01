#include "engine.h"
#include "tools/AdrOp.h"
#include "tools/array_ops.h"

Engine::Engine(Operator* op)
{
	Op = op;
	Init();
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

	for (unsigned int iter=0;iter<iterTS;++iter)
	{
		//voltage updates
		for (pos[2]=1;pos[2]<Op->numLines[2]-1;++pos[2])
		{
			for (pos[1]=1;pos[1]<Op->numLines[1]-1;++pos[1])
			{
				for (pos[0]=1;pos[0]<Op->numLines[0]-1;++pos[0])
				{
					//do the updates here
					//for x
					volt[0][pos[0]][pos[1]][pos[2]] *= Op->vv[0][pos[0]][pos[1]][pos[2]];
					volt[0][pos[0]][pos[1]][pos[2]] += Op->vi[0][pos[0]][pos[1]][pos[2]] * ( curr[2][pos[0]][pos[1]][pos[2]] - curr[2][pos[0]][pos[1]-1][pos[2]] - curr[1][pos[0]][pos[1]][pos[2]] + curr[1][pos[0]][pos[1]][pos[2]-1]);

					//for y
					volt[1][pos[0]][pos[1]][pos[2]] *= Op->vv[1][pos[0]][pos[1]][pos[2]];
					volt[1][pos[0]][pos[1]][pos[2]] += Op->vi[1][pos[0]][pos[1]][pos[2]] * ( curr[0][pos[0]][pos[1]][pos[2]] - curr[0][pos[0]][pos[1]][pos[2]-1] - curr[2][pos[0]][pos[1]][pos[2]] + curr[2][pos[0]-1][pos[1]][pos[2]]);

					//for x
					volt[2][pos[0]][pos[1]][pos[2]] *= Op->vv[2][pos[0]][pos[1]][pos[2]];
					volt[2][pos[0]][pos[1]][pos[2]] += Op->vi[2][pos[0]][pos[1]][pos[2]] * ( curr[1][pos[0]][pos[1]][pos[2]] - curr[1][pos[0]-1][pos[1]][pos[2]] - curr[0][pos[0]][pos[1]][pos[2]] + curr[0][pos[0]][pos[1]-1][pos[2]]);
				}
			}
		}

		//soft voltage excitation here (E-field excite)
		for (unsigned int n=0;n<Op->E_Ex_Count;++n)
		{
			exc_pos = (int)numTS - (int)Op->E_Ex_delay[n];
			exc_pos*= (exc_pos>0 && exc_pos<Op->ExciteLength);
			volt[0][Op->E_Ex_index[0][n]][Op->E_Ex_index[1][n]][Op->E_Ex_index[2][n]] += Op->E_Ex_amp[0][n]*Op->ExciteSignal[exc_pos];
			volt[1][Op->E_Ex_index[0][n]][Op->E_Ex_index[1][n]][Op->E_Ex_index[2][n]] += Op->E_Ex_amp[1][n]*Op->ExciteSignal[exc_pos];
			volt[2][Op->E_Ex_index[0][n]][Op->E_Ex_index[1][n]][Op->E_Ex_index[2][n]] += Op->E_Ex_amp[2][n]*Op->ExciteSignal[exc_pos];
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
					curr[0][pos[0]][pos[1]][pos[2]] -= Op->iv[0][pos[0]][pos[1]][pos[2]] * ( volt[2][pos[0]][pos[1]+1][pos[2]] - volt[2][pos[0]][pos[1]][pos[2]] - volt[1][pos[0]][pos[1]][pos[2]+1] + volt[1][pos[0]][pos[1]][pos[2]]);

					//for y
					curr[1][pos[0]][pos[1]][pos[2]] *= Op->ii[1][pos[0]][pos[1]][pos[2]];
					curr[1][pos[0]][pos[1]][pos[2]] -= Op->iv[1][pos[0]][pos[1]][pos[2]] * ( volt[0][pos[0]][pos[1]][pos[2]+1] - volt[0][pos[0]][pos[1]][pos[2]] - volt[2][pos[0]+1][pos[1]][pos[2]] + volt[2][pos[0]][pos[1]][pos[2]]);

					//for x
					curr[2][pos[0]][pos[1]][pos[2]] *= Op->ii[2][pos[0]][pos[1]][pos[2]];
					curr[2][pos[0]][pos[1]][pos[2]] -= Op->iv[2][pos[0]][pos[1]][pos[2]] * ( volt[1][pos[0]+1][pos[1]][pos[2]] - volt[1][pos[0]][pos[1]][pos[2]] - volt[0][pos[0]][pos[1]+1][pos[2]] + volt[0][pos[0]][pos[1]][pos[2]]);
				}
			}
		}

		//soft current excitation here (H-field excite)
		++numTS;
	}
	return true;
}
