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

#include "engine_sse_compressed.h"
#ifdef __SSE2__
#include <emmintrin.h>
#endif

Engine_SSE_Compressed* Engine_SSE_Compressed::New(const Operator_SSE_Compressed* op)
{
	cout << "Create FDTD engine (compressed SSE)" << endl;
	Engine_SSE_Compressed* e = new Engine_SSE_Compressed(op);
	e->Init();
	return e;
}

Engine_SSE_Compressed::Engine_SSE_Compressed(const Operator_SSE_Compressed* op) : Engine_sse(op)
{
	Op = op;
}

Engine_SSE_Compressed::~Engine_SSE_Compressed()
{
}

void Engine_SSE_Compressed::UpdateVoltages(unsigned int startX, unsigned int numX)
{
	unsigned int pos[3];
	bool shift[2];
	f4vector temp;

	pos[0] = startX;
	unsigned int index=0;
	for (unsigned int posX=0; posX<numX; ++posX)
	{
		shift[0]=pos[0];
		for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
		{
			shift[1]=pos[1];
			for (pos[2]=1; pos[2]<numVectors; ++pos[2])
			{
				index = Op->m_Op_index[pos[0]][pos[1]][pos[2]];
				// x-polarization
				f4_volt[0][pos[0]][pos[1]][pos[2]].v *= Op->f4_vv_Compressed[0][index].v;
				f4_volt[0][pos[0]][pos[1]][pos[2]].v += Op->f4_vi_Compressed[0][index].v * ( f4_curr[2][pos[0]][pos[1]][pos[2]].v - f4_curr[2][pos[0]][pos[1]-shift[1]][pos[2]].v - f4_curr[1][pos[0]][pos[1]][pos[2]].v + f4_curr[1][pos[0]][pos[1]][pos[2]-1].v );

				// y-polarization
				f4_volt[1][pos[0]][pos[1]][pos[2]].v *= Op->f4_vv_Compressed[1][index].v;
				f4_volt[1][pos[0]][pos[1]][pos[2]].v += Op->f4_vi_Compressed[1][index].v * ( f4_curr[0][pos[0]][pos[1]][pos[2]].v - f4_curr[0][pos[0]][pos[1]][pos[2]-1].v - f4_curr[2][pos[0]][pos[1]][pos[2]].v + f4_curr[2][pos[0]-shift[0]][pos[1]][pos[2]].v);

				// z-polarization
				f4_volt[2][pos[0]][pos[1]][pos[2]].v *= Op->f4_vv_Compressed[2][index].v;
				f4_volt[2][pos[0]][pos[1]][pos[2]].v += Op->f4_vi_Compressed[2][index].v * ( f4_curr[1][pos[0]][pos[1]][pos[2]].v - f4_curr[1][pos[0]-shift[0]][pos[1]][pos[2]].v - f4_curr[0][pos[0]][pos[1]][pos[2]].v + f4_curr[0][pos[0]][pos[1]-shift[1]][pos[2]].v);
			}

			// for pos[2] = 0
			// x-polarization
			index = Op->m_Op_index[pos[0]][pos[1]][0];
#ifdef __SSE2__
			temp.v = (__m128)_mm_slli_si128( (__m128i)f4_curr[1][pos[0]][pos[1]][numVectors-1].v, 4 );
#else
			temp.f[0] = 0;
			temp.f[1] = f4_curr[1][pos[0]][pos[1]][numVectors-1].f[0];
			temp.f[2] = f4_curr[1][pos[0]][pos[1]][numVectors-1].f[1];
			temp.f[3] = f4_curr[1][pos[0]][pos[1]][numVectors-1].f[2];
#endif
			f4_volt[0][pos[0]][pos[1]][0].v *= Op->f4_vv_Compressed[0][index].v;
			f4_volt[0][pos[0]][pos[1]][0].v += Op->f4_vi_Compressed[0][index].v * ( f4_curr[2][pos[0]][pos[1]][0].v - f4_curr[2][pos[0]][pos[1]-shift[1]][0].v - f4_curr[1][pos[0]][pos[1]][0].v + temp.v );

			// y-polarization
#ifdef __SSE2__
			temp.v = (__m128)_mm_slli_si128( (__m128i)f4_curr[0][pos[0]][pos[1]][numVectors-1].v, 4 );
#else
			temp.f[0] = 0;
			temp.f[1] = f4_curr[0][pos[0]][pos[1]][numVectors-1].f[0];
			temp.f[2] = f4_curr[0][pos[0]][pos[1]][numVectors-1].f[1];
			temp.f[3] = f4_curr[0][pos[0]][pos[1]][numVectors-1].f[2];
#endif
			f4_volt[1][pos[0]][pos[1]][0].v *= Op->f4_vv_Compressed[1][index].v;
			f4_volt[1][pos[0]][pos[1]][0].v += Op->f4_vi_Compressed[1][index].v * ( f4_curr[0][pos[0]][pos[1]][0].v - temp.v - f4_curr[2][pos[0]][pos[1]][0].v + f4_curr[2][pos[0]-shift[0]][pos[1]][0].v);

			// z-polarization
			f4_volt[2][pos[0]][pos[1]][0].v *= Op->f4_vv_Compressed[2][index].v;
			f4_volt[2][pos[0]][pos[1]][0].v += Op->f4_vi_Compressed[2][index].v * ( f4_curr[1][pos[0]][pos[1]][0].v - f4_curr[1][pos[0]-shift[0]][pos[1]][0].v - f4_curr[0][pos[0]][pos[1]][0].v + f4_curr[0][pos[0]][pos[1]-shift[1]][0].v);
		}
		++pos[0];
	}
}

void Engine_SSE_Compressed::UpdateCurrents(unsigned int startX, unsigned int numX)
{
	unsigned int pos[3];
	f4vector temp;

	pos[0] = startX;
	unsigned int index;
	for (unsigned int posX=0; posX<numX; ++posX)
	{
		for (pos[1]=0; pos[1]<numLines[1]-1; ++pos[1])
		{
			for (pos[2]=0; pos[2]<numVectors-1; ++pos[2])
			{
				index = Op->m_Op_index[pos[0]][pos[1]][pos[2]];
				// x-pol
				f4_curr[0][pos[0]][pos[1]][pos[2]].v *= Op->f4_ii_Compressed[0][index].v;
				f4_curr[0][pos[0]][pos[1]][pos[2]].v += Op->f4_iv_Compressed[0][index].v * ( f4_volt[2][pos[0]][pos[1]][pos[2]].v - f4_volt[2][pos[0]][pos[1]+1][pos[2]].v - f4_volt[1][pos[0]][pos[1]][pos[2]].v + f4_volt[1][pos[0]][pos[1]][pos[2]+1].v);

				// y-pol
				f4_curr[1][pos[0]][pos[1]][pos[2]].v *= Op->f4_ii_Compressed[1][index].v;
				f4_curr[1][pos[0]][pos[1]][pos[2]].v += Op->f4_iv_Compressed[1][index].v * ( f4_volt[0][pos[0]][pos[1]][pos[2]].v - f4_volt[0][pos[0]][pos[1]][pos[2]+1].v - f4_volt[2][pos[0]][pos[1]][pos[2]].v + f4_volt[2][pos[0]+1][pos[1]][pos[2]].v);

				// z-pol
				f4_curr[2][pos[0]][pos[1]][pos[2]].v *= Op->f4_ii_Compressed[2][index].v;
				f4_curr[2][pos[0]][pos[1]][pos[2]].v += Op->f4_iv_Compressed[2][index].v * ( f4_volt[1][pos[0]][pos[1]][pos[2]].v - f4_volt[1][pos[0]+1][pos[1]][pos[2]].v - f4_volt[0][pos[0]][pos[1]][pos[2]].v + f4_volt[0][pos[0]][pos[1]+1][pos[2]].v);
			}

			index = Op->m_Op_index[pos[0]][pos[1]][numVectors-1];
			// for pos[2] = numVectors-1
			// x-pol
#ifdef __SSE2__
			temp.v = (__m128)_mm_srli_si128( (__m128i)f4_volt[1][pos[0]][pos[1]][0].v, 4 );
#else
			temp.f[0] = f4_volt[1][pos[0]][pos[1]][0].f[1];
			temp.f[1] = f4_volt[1][pos[0]][pos[1]][0].f[2];
			temp.f[2] = f4_volt[1][pos[0]][pos[1]][0].f[3];
			temp.f[3] = 0;
#endif
			f4_curr[0][pos[0]][pos[1]][numVectors-1].v *= Op->f4_ii_Compressed[0][index].v;
			f4_curr[0][pos[0]][pos[1]][numVectors-1].v += Op->f4_iv_Compressed[0][index].v * ( f4_volt[2][pos[0]][pos[1]][numVectors-1].v - f4_volt[2][pos[0]][pos[1]+1][numVectors-1].v - f4_volt[1][pos[0]][pos[1]][numVectors-1].v + temp.v);

			// y-pol
#ifdef __SSE2__
			temp.v = (__m128)_mm_srli_si128( (__m128i)f4_volt[0][pos[0]][pos[1]][0].v, 4 );
#else
			temp.f[0] = f4_volt[0][pos[0]][pos[1]][0].f[1];
			temp.f[1] = f4_volt[0][pos[0]][pos[1]][0].f[2];
			temp.f[2] = f4_volt[0][pos[0]][pos[1]][0].f[3];
			temp.f[3] = 0;
#endif
			f4_curr[1][pos[0]][pos[1]][numVectors-1].v *= Op->f4_ii_Compressed[1][index].v;
			f4_curr[1][pos[0]][pos[1]][numVectors-1].v += Op->f4_iv_Compressed[1][index].v * ( f4_volt[0][pos[0]][pos[1]][numVectors-1].v - temp.v - f4_volt[2][pos[0]][pos[1]][numVectors-1].v + f4_volt[2][pos[0]+1][pos[1]][numVectors-1].v);

			// z-pol
			f4_curr[2][pos[0]][pos[1]][numVectors-1].v *= Op->f4_ii_Compressed[2][index].v;
			f4_curr[2][pos[0]][pos[1]][numVectors-1].v += Op->f4_iv_Compressed[2][index].v * ( f4_volt[1][pos[0]][pos[1]][numVectors-1].v - f4_volt[1][pos[0]+1][pos[1]][numVectors-1].v - f4_volt[0][pos[0]][pos[1]][numVectors-1].v + f4_volt[0][pos[0]][pos[1]+1][numVectors-1].v);
		}
		++pos[0];
	}
}
