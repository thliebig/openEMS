/*
*	Copyright (C) 2010 Sebastian Held (Sebastian.Held@gmx.de)
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

#include "operator_sse.h"
#include "tools/array_ops.h"

Operator_sse* Operator_sse::New()
{
	Operator_sse* op = new Operator_sse();
	op->Init();
	return op;
}

Operator_sse::Operator_sse() : Operator()
{
}

Operator_sse::~Operator_sse()
{
	Reset();
}

void Operator_sse::Init()
{
	Operator::Init();
	vv_ = 0;
	vi_ = 0;
	iv_ = 0;
	ii_ = 0;
}

void Operator_sse::Reset()
{
	Delete_N_3DArray_v4sf(vv_,numLines);
	Delete_N_3DArray_v4sf(vi_,numLines);
	Delete_N_3DArray_v4sf(iv_,numLines);
	Delete_N_3DArray_v4sf(ii_,numLines);
	Operator::Reset();
	Init(); // FIXME this calls Operator::Init() twice...
}

void Operator_sse::InitOperator()
{
	Operator::InitOperator();
	Delete_N_3DArray_v4sf(vv_,numLines);
	Delete_N_3DArray_v4sf(vi_,numLines);
	Delete_N_3DArray_v4sf(iv_,numLines);
	Delete_N_3DArray_v4sf(ii_,numLines);
	vv_ = Create_N_3DArray_v4sf(numLines);
	vi_ = Create_N_3DArray_v4sf(numLines);
	iv_ = Create_N_3DArray_v4sf(numLines);
	ii_ = Create_N_3DArray_v4sf(numLines);
}

int Operator_sse::CalcECOperator()
{
	Operator::CalcECOperator();

	// copy operator to aligned memory
	// FIXME this is really inefficient!
	unsigned int pos[3];
	for (int n=0;n<3;++n)
	{
		for (pos[0]=0;pos[0]<numLines[0];++pos[0])
		{
			for (pos[1]=0;pos[1]<numLines[1];++pos[1])
			{
				for (pos[2]=0;pos[2]<numLines[2];++pos[2])
				{
					vv_[n][pos[0]][pos[1]][pos[2]/4].f[pos[2]%4] = vv[n][pos[0]][pos[1]][pos[2]];
					vi_[n][pos[0]][pos[1]][pos[2]/4].f[pos[2]%4] = vi[n][pos[0]][pos[1]][pos[2]];
					iv_[n][pos[0]][pos[1]][pos[2]/4].f[pos[2]%4] = iv[n][pos[0]][pos[1]][pos[2]];
					ii_[n][pos[0]][pos[1]][pos[2]/4].f[pos[2]%4] = ii[n][pos[0]][pos[1]][pos[2]];
				}
			}
		}
	}

	return 0;
}

void Operator_sse::ApplyMagneticBC(bool* dirs)
{
	Operator::ApplyMagneticBC(dirs);

	// copy operator to aligned memory
	// FIXME this is really inefficient!
	unsigned int pos[3];
	for (int n=0;n<3;++n)
	{
		for (pos[0]=0;pos[0]<numLines[0];++pos[0])
		{
			for (pos[1]=0;pos[1]<numLines[1];++pos[1])
			{
				for (pos[2]=0;pos[2]<numLines[2];++pos[2])
				{
					vv_[n][pos[0]][pos[1]][pos[2]/4].f[pos[2]%4] = vv[n][pos[0]][pos[1]][pos[2]];
					vi_[n][pos[0]][pos[1]][pos[2]/4].f[pos[2]%4] = vi[n][pos[0]][pos[1]][pos[2]];
					iv_[n][pos[0]][pos[1]][pos[2]/4].f[pos[2]%4] = iv[n][pos[0]][pos[1]][pos[2]];
					ii_[n][pos[0]][pos[1]][pos[2]/4].f[pos[2]%4] = ii[n][pos[0]][pos[1]][pos[2]];
				}
			}
		}
	}
}
