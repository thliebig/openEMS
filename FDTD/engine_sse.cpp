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

#ifndef SSE_CORRECT_DENORMALS
#include <xmmintrin.h>
#endif

#include "engine_sse.h"

//! \brief construct an Engine_sse instance
//! it's the responsibility of the caller to free the returned pointer
Engine_sse* Engine_sse::New(const Operator_sse* op)
{
	cout << "Create FDTD engine (SSE)" << endl;
	Engine_sse* e = new Engine_sse(op);
	e->Init();
	return e;
}

Engine_sse::Engine_sse(const Operator_sse* op) : Engine(op)
{
	m_type = SSE;
	Op = op;
	f4_volt = 0;
	f4_curr = 0;
	numVectors =  ceil((double)numLines[2]/4.0);

	// speed up the calculation of denormal floating point values (flush-to-zero)
#ifndef SSE_CORRECT_DENORMALS
	unsigned int oldMXCSR = _mm_getcsr(); //read the old MXCSR setting
	unsigned int newMXCSR = oldMXCSR | 0x8040; // set DAZ and FZ bits
	_mm_setcsr( newMXCSR ); //write the new MXCSR setting to the MXCSR
#endif
}

Engine_sse::~Engine_sse()
{
	//_mm_setcsr( oldMXCSR ); // restore old setting
	Reset();
}

void Engine_sse::Init()
{
	Engine::Init();

	Delete_N_3DArray(volt,numLines);
	volt=NULL; // not used
	Delete_N_3DArray(curr,numLines);
	curr=NULL; // not used

	f4_volt = Create_N_3DArray_v4sf(numLines);
	f4_curr = Create_N_3DArray_v4sf(numLines);
}

void Engine_sse::Reset()
{
	Engine::Reset();
	Delete_N_3DArray_v4sf(f4_volt,numLines);
	f4_volt = 0;
	Delete_N_3DArray_v4sf(f4_curr,numLines);
	f4_curr = 0;
}

void Engine_sse::UpdateVoltages(unsigned int startX, unsigned int numX)
{
	unsigned int pos[3];
	bool shift[2];
	f4vector temp;

	pos[0] = startX;
	for (unsigned int posX=0; posX<numX; ++posX)
	{
		shift[0]=pos[0];
		for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
		{
			shift[1]=pos[1];
			for (pos[2]=1; pos[2]<numVectors; ++pos[2])
			{
				// x-polarization
				f4_volt[0][pos[0]][pos[1]][pos[2]].v *= Op->f4_vv[0][pos[0]][pos[1]][pos[2]].v;
				f4_volt[0][pos[0]][pos[1]][pos[2]].v += Op->f4_vi[0][pos[0]][pos[1]][pos[2]].v * ( f4_curr[2][pos[0]][pos[1]][pos[2]].v - f4_curr[2][pos[0]][pos[1]-shift[1]][pos[2]].v - f4_curr[1][pos[0]][pos[1]][pos[2]].v + f4_curr[1][pos[0]][pos[1]][pos[2]-1].v );

				// y-polarization
				f4_volt[1][pos[0]][pos[1]][pos[2]].v *= Op->f4_vv[1][pos[0]][pos[1]][pos[2]].v;
				f4_volt[1][pos[0]][pos[1]][pos[2]].v += Op->f4_vi[1][pos[0]][pos[1]][pos[2]].v * ( f4_curr[0][pos[0]][pos[1]][pos[2]].v - f4_curr[0][pos[0]][pos[1]][pos[2]-1].v - f4_curr[2][pos[0]][pos[1]][pos[2]].v + f4_curr[2][pos[0]-shift[0]][pos[1]][pos[2]].v);

				// z-polarization
				f4_volt[2][pos[0]][pos[1]][pos[2]].v *= Op->f4_vv[2][pos[0]][pos[1]][pos[2]].v;
				f4_volt[2][pos[0]][pos[1]][pos[2]].v += Op->f4_vi[2][pos[0]][pos[1]][pos[2]].v * ( f4_curr[1][pos[0]][pos[1]][pos[2]].v - f4_curr[1][pos[0]-shift[0]][pos[1]][pos[2]].v - f4_curr[0][pos[0]][pos[1]][pos[2]].v + f4_curr[0][pos[0]][pos[1]-shift[1]][pos[2]].v);
			}

			// for pos[2] = 0
			// x-polarization
			temp.f[0] = 0;
			temp.f[1] = f4_curr[1][pos[0]][pos[1]][numVectors-1].f[0];
			temp.f[2] = f4_curr[1][pos[0]][pos[1]][numVectors-1].f[1];
			temp.f[3] = f4_curr[1][pos[0]][pos[1]][numVectors-1].f[2];
			f4_volt[0][pos[0]][pos[1]][0].v *= Op->f4_vv[0][pos[0]][pos[1]][0].v;
			f4_volt[0][pos[0]][pos[1]][0].v += Op->f4_vi[0][pos[0]][pos[1]][0].v * ( f4_curr[2][pos[0]][pos[1]][0].v - f4_curr[2][pos[0]][pos[1]-shift[1]][0].v - f4_curr[1][pos[0]][pos[1]][0].v + temp.v );

			// y-polarization
			temp.f[0] = 0;
			temp.f[1] = f4_curr[0][pos[0]][pos[1]][numVectors-1].f[0];
			temp.f[2] = f4_curr[0][pos[0]][pos[1]][numVectors-1].f[1];
			temp.f[3] = f4_curr[0][pos[0]][pos[1]][numVectors-1].f[2];
			f4_volt[1][pos[0]][pos[1]][0].v *= Op->f4_vv[1][pos[0]][pos[1]][0].v;
			f4_volt[1][pos[0]][pos[1]][0].v += Op->f4_vi[1][pos[0]][pos[1]][0].v * ( f4_curr[0][pos[0]][pos[1]][0].v - temp.v - f4_curr[2][pos[0]][pos[1]][0].v + f4_curr[2][pos[0]-shift[0]][pos[1]][0].v);

			// z-polarization
			f4_volt[2][pos[0]][pos[1]][0].v *= Op->f4_vv[2][pos[0]][pos[1]][0].v;
			f4_volt[2][pos[0]][pos[1]][0].v += Op->f4_vi[2][pos[0]][pos[1]][0].v * ( f4_curr[1][pos[0]][pos[1]][0].v - f4_curr[1][pos[0]-shift[0]][pos[1]][0].v - f4_curr[0][pos[0]][pos[1]][0].v + f4_curr[0][pos[0]][pos[1]-shift[1]][0].v);
		}
		++pos[0];
	}
}

void Engine_sse::UpdateCurrents(unsigned int startX, unsigned int numX)
{
	unsigned int pos[5];
	f4vector temp;

	pos[0] = startX;
	for (unsigned int posX=0; posX<numX; ++posX)
	{
		for (pos[1]=0; pos[1]<numLines[1]-1; ++pos[1])
		{
			for (pos[2]=0; pos[2]<numVectors-1; ++pos[2])
			{
				// x-pol
				f4_curr[0][pos[0]][pos[1]][pos[2]].v *= Op->f4_ii[0][pos[0]][pos[1]][pos[2]].v;
				f4_curr[0][pos[0]][pos[1]][pos[2]].v += Op->f4_iv[0][pos[0]][pos[1]][pos[2]].v * ( f4_volt[2][pos[0]][pos[1]][pos[2]].v - f4_volt[2][pos[0]][pos[1]+1][pos[2]].v - f4_volt[1][pos[0]][pos[1]][pos[2]].v + f4_volt[1][pos[0]][pos[1]][pos[2]+1].v);

				// y-pol
				f4_curr[1][pos[0]][pos[1]][pos[2]].v *= Op->f4_ii[1][pos[0]][pos[1]][pos[2]].v;
				f4_curr[1][pos[0]][pos[1]][pos[2]].v += Op->f4_iv[1][pos[0]][pos[1]][pos[2]].v * ( f4_volt[0][pos[0]][pos[1]][pos[2]].v - f4_volt[0][pos[0]][pos[1]][pos[2]+1].v - f4_volt[2][pos[0]][pos[1]][pos[2]].v + f4_volt[2][pos[0]+1][pos[1]][pos[2]].v);

				// z-pol
				f4_curr[2][pos[0]][pos[1]][pos[2]].v *= Op->f4_ii[2][pos[0]][pos[1]][pos[2]].v;
				f4_curr[2][pos[0]][pos[1]][pos[2]].v += Op->f4_iv[2][pos[0]][pos[1]][pos[2]].v * ( f4_volt[1][pos[0]][pos[1]][pos[2]].v - f4_volt[1][pos[0]+1][pos[1]][pos[2]].v - f4_volt[0][pos[0]][pos[1]][pos[2]].v + f4_volt[0][pos[0]][pos[1]+1][pos[2]].v);
			}

			// for pos[2] = numVectors-1
			// x-pol
			temp.f[0] = f4_volt[1][pos[0]][pos[1]][0].f[1];
			temp.f[1] = f4_volt[1][pos[0]][pos[1]][0].f[2];
			temp.f[2] = f4_volt[1][pos[0]][pos[1]][0].f[3];
			temp.f[3] = 0;
			f4_curr[0][pos[0]][pos[1]][numVectors-1].v *= Op->f4_ii[0][pos[0]][pos[1]][numVectors-1].v;
			f4_curr[0][pos[0]][pos[1]][numVectors-1].v += Op->f4_iv[0][pos[0]][pos[1]][numVectors-1].v * ( f4_volt[2][pos[0]][pos[1]][numVectors-1].v - f4_volt[2][pos[0]][pos[1]+1][numVectors-1].v - f4_volt[1][pos[0]][pos[1]][numVectors-1].v + temp.v);

			// y-pol
			temp.f[0] = f4_volt[0][pos[0]][pos[1]][0].f[1];
			temp.f[1] = f4_volt[0][pos[0]][pos[1]][0].f[2];
			temp.f[2] = f4_volt[0][pos[0]][pos[1]][0].f[3];
			temp.f[3] = 0;
			f4_curr[1][pos[0]][pos[1]][numVectors-1].v *= Op->f4_ii[1][pos[0]][pos[1]][numVectors-1].v;
			f4_curr[1][pos[0]][pos[1]][numVectors-1].v += Op->f4_iv[1][pos[0]][pos[1]][numVectors-1].v * ( f4_volt[0][pos[0]][pos[1]][numVectors-1].v - temp.v - f4_volt[2][pos[0]][pos[1]][numVectors-1].v + f4_volt[2][pos[0]+1][pos[1]][numVectors-1].v);

			// z-pol
			f4_curr[2][pos[0]][pos[1]][numVectors-1].v *= Op->f4_ii[2][pos[0]][pos[1]][numVectors-1].v;
			f4_curr[2][pos[0]][pos[1]][numVectors-1].v += Op->f4_iv[2][pos[0]][pos[1]][numVectors-1].v * ( f4_volt[1][pos[0]][pos[1]][numVectors-1].v - f4_volt[1][pos[0]+1][pos[1]][numVectors-1].v - f4_volt[0][pos[0]][pos[1]][numVectors-1].v + f4_volt[0][pos[0]][pos[1]+1][numVectors-1].v);
		}
		++pos[0];
	}
}
