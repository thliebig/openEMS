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

// enable AltiVec mode on ppc64le in GCC's emmintrin.h
#include <boost/predef.h>
#if BOOST_ARCH_PPC_64 && BOOST_ENDIAN_LITTLE_BYTE
#define NO_WARN_X86_INTRINSICS
#include <emmintrin.h>
#endif

#ifdef __SSE2__
#include <emmintrin.h>
#endif

using std::cout;
using std::endl;

//! \brief construct an Engine_sse instance
//! it's the responsibility of the caller to free the returned pointer
Engine_sse* Engine_sse::New(const Operator* op, unsigned int numThreads)
{
	cout << "Create FDTD engine (" << GetTypeName(SSE) << ")" << endl;
	Engine_sse* e = new Engine_sse(op);
	e->SetNumThreads(numThreads);
	e->Init();
	return e;
}

Engine_sse::Engine_sse(const Operator* op) : Engine(op)
{
	m_type = SSE;
	f4_volt_ptr = NULL;
	f4_curr_ptr = NULL;
	numVectors =  ceil((double)numLines[2]/4.0);
}

Engine_sse::~Engine_sse()
{
	Reset();
}

void Engine_sse::Init()
{
	Engine::Init();

	// This engine uses its own SIMD arrays to represent E&M fields, so
	// free the arrays from the base class.
	delete volt_ptr;
	volt_ptr = NULL;
	delete curr_ptr;
	curr_ptr = NULL;

	f4_volt_ptr = new ArrayLib::ArrayENG<f4vector>(
		"f4_volt", {numLines[0], numLines[1], numVectors}
	);
	f4_curr_ptr = new ArrayLib::ArrayENG<f4vector>(
		"f4_curr", {numLines[0], numLines[1], numVectors}
	);
}

void Engine_sse::Reset()
{
	Engine::Reset();
	delete f4_volt_ptr;
	f4_volt_ptr = NULL;
	delete f4_curr_ptr;
	f4_curr_ptr = NULL;
}

void Engine_sse::UpdateVoltages(unsigned int startX, unsigned int numX)
{
	f4vector* f4_volt = f4_volt_ptr->data();
	f4vector* f4_curr = f4_curr_ptr->data();
	const Coeffs_SSE_Compressed& coeffs = Op->GetCoeffsAs<Coeffs_SSE_Compressed>();

	unsigned int pos[3];
	size_t v_pos;
	size_t i_pos;
	size_t i_pos_z_start, i_pos_z_end;
	size_t i_shift[3]; // shift of current index by x/y/z
	f4vector temp;

	// array index shift for x/y/z field component
	size_t v_N_shift  = f4_volt_ptr->stride(0);
	size_t i_N_shift  = f4_curr_ptr->stride(0);

	// index shift for a z-position
	i_shift[2] = f4_curr_ptr->stride(3);

	pos[0] = startX;
	unsigned int index=0;
	for (unsigned int posX=0; posX<numX; ++posX)
	{
		// index shift for a x-position
		i_shift[0] = (pos[0]>0)*f4_curr_ptr->stride(1);
		for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
		{
			// index shift for a y-position
			i_shift[1] = (pos[1]>0)*f4_curr_ptr->stride(2);
			for (pos[2]=1; pos[2]<numVectors; ++pos[2])
			{
				index = coeffs.index(pos[0], pos[1], pos[2]);
				i_pos = f4_curr_ptr->linearIndex({0, pos[0], pos[1], pos[2]});
				v_pos = f4_volt_ptr->linearIndex({0, pos[0], pos[1], pos[2]});

				// x-polarization
				f4_volt[v_pos].v *= coeffs.f4_vv_Compressed[0][index].v;
				f4_volt[v_pos].v +=
				    coeffs.f4_vi_Compressed[0][index].v * (
				        f4_curr[i_pos + 2*i_N_shift             ].v -
				        f4_curr[i_pos + 2*i_N_shift - i_shift[1]].v -
				        f4_curr[i_pos +   i_N_shift             ].v +
				        f4_curr[i_pos +   i_N_shift - i_shift[2]].v
				    );

				// y-polarization
				v_pos += v_N_shift; // go to next component
				f4_volt[v_pos].v *= coeffs.f4_vv_Compressed[1][index].v;
				f4_volt[v_pos].v +=
				    coeffs.f4_vi_Compressed[1][index].v * (
				        f4_curr[i_pos                           ].v -
				        f4_curr[i_pos               - i_shift[2]].v -
				        f4_curr[i_pos + 2*i_N_shift             ].v +
				        f4_curr[i_pos + 2*i_N_shift - i_shift[0]].v
				    );

				// z-polarization
				v_pos += v_N_shift; // go to next component
				f4_volt[v_pos].v *= coeffs.f4_vv_Compressed[2][index].v;
				f4_volt[v_pos].v +=
				    coeffs.f4_vi_Compressed[2][index].v * (
				        f4_curr[i_pos + i_N_shift             ].v -
				        f4_curr[i_pos + i_N_shift - i_shift[0]].v -
				        f4_curr[i_pos                         ].v +
				        f4_curr[i_pos             - i_shift[1]].v
				    );
			}

			i_pos_z_start = f4_curr_ptr->linearIndex({0, pos[0], pos[1], 0});
			i_pos_z_end   = f4_curr_ptr->linearIndex({0, pos[0], pos[1], numVectors-1});

			// for pos[2] = 0
			// x-polarization
			index = coeffs.index(pos[0], pos[1], 0);
#ifdef __SSE2__
			temp.v = (__m128)_mm_slli_si128(
			             (__m128i)f4_curr[i_pos_z_end + i_N_shift].v, 4
			         );
#else
			temp.f[0] = 0;
			temp.f[1] = f4_curr[i_pos_z_end + i_N_shift].f[0];
			temp.f[2] = f4_curr[i_pos_z_end + i_N_shift].f[1];
			temp.f[3] = f4_curr[i_pos_z_end + i_N_shift].f[2];
#endif
			v_pos = f4_volt_ptr->linearIndex({0, pos[0], pos[1], 0});
			f4_volt[v_pos].v *= coeffs.f4_vv_Compressed[0][index].v;
			f4_volt[v_pos].v +=
			    coeffs.f4_vi_Compressed[0][index].v * (
			        f4_curr[i_pos_z_start + i_N_shift*2             ].v -
			        f4_curr[i_pos_z_start + i_N_shift*2 - i_shift[1]].v -
			        f4_curr[i_pos_z_start + i_N_shift  ].v +
			        temp.v
			    );

			// y-polarization
#ifdef __SSE2__
			temp.v = (__m128)_mm_slli_si128(
			             (__m128i)f4_curr[i_pos_z_end].v, 4
			         );
#else
			temp.f[0] = 0;
			temp.f[1] = f4_curr[i_pos_z_end].f[0];
			temp.f[2] = f4_curr[i_pos_z_end].f[1];
			temp.f[3] = f4_curr[i_pos_z_end].f[2];
#endif
			v_pos += v_N_shift; // go to next component
			f4_volt[v_pos].v *= coeffs.f4_vv_Compressed[1][index].v;
			f4_volt[v_pos].v +=
			    coeffs.f4_vi_Compressed[1][index].v * (
			        f4_curr[i_pos_z_start].v -
			        temp.v -
			        f4_curr[i_pos_z_start + i_N_shift*2             ].v +
			        f4_curr[i_pos_z_start + i_N_shift*2 - i_shift[0]].v
			    );

			// z-polarization
			v_pos += v_N_shift; // go to next component
			f4_volt[v_pos].v *= coeffs.f4_vv_Compressed[2][index].v;
			f4_volt[v_pos].v +=
			    coeffs.f4_vi_Compressed[2][index].v * (
			        f4_curr[i_pos_z_start + i_N_shift             ].v -
			        f4_curr[i_pos_z_start + i_N_shift - i_shift[0]].v -
			        f4_curr[i_pos_z_start                         ].v +
			        f4_curr[i_pos_z_start             - i_shift[1]].v
			    );

		}
		++pos[0];
	}
}

void Engine_sse::UpdateCurrents(unsigned int startX, unsigned int numX)
{
	f4vector* f4_volt = f4_volt_ptr->data();
	f4vector* f4_curr = f4_curr_ptr->data();
	const Coeffs_SSE_Compressed& coeffs = Op->GetCoeffsAs<Coeffs_SSE_Compressed>();

	size_t v_pos;
	size_t i_pos;
	size_t v_shift[3]; // shift of voltages index by x/y/z
	size_t v_pos_z_start, v_pos_z_end;

	unsigned int pos[3];
	f4vector temp;

	// array index shift for x/y/z field component
	size_t v_N_shift  = f4_volt_ptr->stride(0);
	size_t i_N_shift  = f4_curr_ptr->stride(0);

	// index shift for a x-, y- and z-position
	for (unsigned int n=0;n<3;++n)
		v_shift[n]= f4_volt_ptr->stride(n+1);

	pos[0] = startX;
	unsigned int index;
	for (unsigned int posX=0; posX<numX; ++posX)
	{
		for (pos[1]=0; pos[1]<numLines[1]-1; ++pos[1])
		{
			for (pos[2]=0; pos[2]<numVectors-1; ++pos[2])
			{
				index = coeffs.index(pos[0], pos[1], pos[2]);
				i_pos = f4_curr_ptr->linearIndex({0, pos[0], pos[1], pos[2]});
				v_pos = f4_volt_ptr->linearIndex({0, pos[0], pos[1], pos[2]});

				// x-pol
				f4_curr[i_pos].v *= coeffs.f4_ii_Compressed[0][index].v;
				f4_curr[i_pos].v +=
				    coeffs.f4_iv_Compressed[0][index].v * (
				        f4_volt[v_pos + 2*v_N_shift             ].v -
				        f4_volt[v_pos + 2*v_N_shift + v_shift[1]].v -
				        f4_volt[v_pos +   v_N_shift             ].v +
				        f4_volt[v_pos +   v_N_shift + v_shift[2]].v
				    );

				// y-pol
				i_pos += i_N_shift;
				f4_curr[i_pos].v *= coeffs.f4_ii_Compressed[1][index].v;
				f4_curr[i_pos].v +=
				    coeffs.f4_iv_Compressed[1][index].v * (
				        f4_volt[v_pos                           ].v -
				        f4_volt[v_pos +               v_shift[2]].v -
				        f4_volt[v_pos + 2*v_N_shift             ].v +
				        f4_volt[v_pos + 2*v_N_shift + v_shift[0]].v
				    );

				// z-pol
				i_pos += i_N_shift;
				f4_curr[i_pos].v *= coeffs.f4_ii_Compressed[2][index].v;
				f4_curr[i_pos].v +=
				    coeffs.f4_iv_Compressed[2][index].v * (
				        f4_volt[v_pos +   v_N_shift             ].v -
				        f4_volt[v_pos +   v_N_shift + v_shift[0]].v -
				        f4_volt[v_pos                           ].v +
				        f4_volt[v_pos +               v_shift[1]].v
				    );
			}

			v_pos_z_start = f4_volt_ptr->linearIndex({0, pos[0], pos[1], 0});
			v_pos_z_end   = f4_volt_ptr->linearIndex({0, pos[0], pos[1], numVectors-1});

			index = coeffs.index(pos[0], pos[1], numVectors-1);
			// for pos[2] = numVectors-1
			// x-pol
#ifdef __SSE2__
			temp.v = (__m128)_mm_srli_si128(
			             (__m128i)f4_volt[v_pos_z_start + v_N_shift].v, 4
			         );
#else
			temp.f[0] = f4_volt[v_pos_z_start + v_N_shift].f[1];
			temp.f[1] = f4_volt[v_pos_z_start + v_N_shift].f[2];
			temp.f[2] = f4_volt[v_pos_z_start + v_N_shift].f[3];
			temp.f[3] = 0;
#endif
			i_pos = f4_curr_ptr->linearIndex({0, pos[0], pos[1], numVectors-1});
			f4_curr[i_pos].v *=
			    coeffs.f4_ii_Compressed[0][index].v;
			f4_curr[i_pos].v +=
			    coeffs.f4_iv_Compressed[0][index].v * (
			        f4_volt[v_pos_z_end + 2*v_N_shift             ].v -
			        f4_volt[v_pos_z_end + 2*v_N_shift + v_shift[1]].v -
			        f4_volt[v_pos_z_end + v_N_shift               ].v +
			        temp.v
			    );

			// y-pol
#ifdef __SSE2__
			temp.v = (__m128)_mm_srli_si128(
			             (__m128i)f4_volt[v_pos_z_start].v, 4
			         );
#else
			temp.f[0] = f4_volt[v_pos_z_start].f[1];
			temp.f[1] = f4_volt[v_pos_z_start].f[2];
			temp.f[2] = f4_volt[v_pos_z_start].f[3];
			temp.f[3] = 0;
#endif
			i_pos += i_N_shift;
			f4_curr[i_pos].v *=
			    coeffs.f4_ii_Compressed[1][index].v;
			f4_curr[i_pos].v +=
			    coeffs.f4_iv_Compressed[1][index].v * (
			        f4_volt[v_pos_z_end].v -
			        temp.v -
			        f4_volt[v_pos_z_end + 2*v_N_shift             ].v +
			        f4_volt[v_pos_z_end + 2*v_N_shift + v_shift[0]].v
			    );

			// z-pol
			i_pos += i_N_shift;
			f4_curr[i_pos].v *=
			    coeffs.f4_ii_Compressed[2][index].v;
			f4_curr[i_pos].v +=
			    coeffs.f4_iv_Compressed[2][index].v * (
			        f4_volt[v_pos_z_end + v_N_shift             ].v -
			        f4_volt[v_pos_z_end + v_N_shift + v_shift[0]].v -
			        f4_volt[v_pos_z_end                         ].v +
			        f4_volt[v_pos_z_end             + v_shift[1]].v
			    );
		}
		++pos[0];
	}
}
