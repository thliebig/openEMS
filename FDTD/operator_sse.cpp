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

#include "engine_sse.h"
#include "operator_sse.h"
#include "tools/array_ops.h"
#include "processfields.h"

Operator_sse* Operator_sse::New()
{
	cout << "Create FDTD operator (SSE)" << endl;
	Operator_sse* op = new Operator_sse();
	op->Init();
	return op;
}

Operator_sse::Operator_sse() : Operator()
{
	f4_vv = 0;
	f4_vi = 0;
	f4_iv = 0;
	f4_ii = 0;
}

Operator_sse::~Operator_sse()
{
	Reset();
}

Engine* Operator_sse::CreateEngine() const
{
	//!create a special sse-engine
	Engine_sse* eng = Engine_sse::New(this);
	return eng;
}

void Operator_sse::Init()
{
	Operator::Init();
	f4_vv = 0;
	f4_vi = 0;
	f4_iv = 0;
	f4_ii = 0;
}

void Operator_sse::Reset()
{
	Delete_N_3DArray_v4sf(f4_vv,numLines);
	Delete_N_3DArray_v4sf(f4_vi,numLines);
	Delete_N_3DArray_v4sf(f4_iv,numLines);
	Delete_N_3DArray_v4sf(f4_ii,numLines);
	f4_vv = 0;
	f4_vi = 0;
	f4_iv = 0;
	f4_ii = 0;
	Operator::Reset();
}

void Operator_sse::InitOperator()
{
	Delete_N_3DArray_v4sf(f4_vv,numLines);
	Delete_N_3DArray_v4sf(f4_vi,numLines);
	Delete_N_3DArray_v4sf(f4_iv,numLines);
	Delete_N_3DArray_v4sf(f4_ii,numLines);
	f4_vv = Create_N_3DArray_v4sf(numLines);
	f4_vi = Create_N_3DArray_v4sf(numLines);
	f4_iv = Create_N_3DArray_v4sf(numLines);
	f4_ii = Create_N_3DArray_v4sf(numLines);

	numVectors =  ceil((double)numLines[2]/4.0);
}

void Operator_sse::DumpOperator2File(string filename)
{
#ifdef OUTPUT_IN_DRAWINGUNITS
	double discLines_scaling = 1;
#else
	double discLines_scaling = GetGridDelta();
#endif

	ofstream file(filename.c_str(),ios_base::out);
	if (file.is_open()==false)
	{
		cerr << "Operator_sse::DumpOperator2File: Can't open file: " << filename << endl;
		return;
	}

	cout << "Dumping FDTD operator information to vtk file: " << filename << " ..." << flush ;


	FDTD_FLOAT**** exc = Create_N_3DArray<FDTD_FLOAT>(numLines);
	if (Exc) {
		for (unsigned int n=0;n<Exc->E_Count;++n)
			exc[Exc->E_dir[n]][Exc->E_index[0][n]][Exc->E_index[1][n]][Exc->E_index[2][n]] = Exc->E_amp[n];
	}

	vv = Create_N_3DArray<FDTD_FLOAT>(numLines);
	vi = Create_N_3DArray<FDTD_FLOAT>(numLines);
	iv = Create_N_3DArray<FDTD_FLOAT>(numLines);
	ii = Create_N_3DArray<FDTD_FLOAT>(numLines);

	unsigned int pos[3];
	for (pos[0]=0;pos[0]<numLines[0];++pos[0])
	{
		for (pos[1]=0;pos[1]<numLines[1];++pos[1])
		{
			for (pos[2]=0;pos[2]<numLines[2];++pos[2])
			{
				for (int n=0;n<3;++n)
				{
					vv[n][pos[0]][pos[1]][pos[2]] = GetVV(n,pos[0],pos[1],pos[2]);
					vi[n][pos[0]][pos[1]][pos[2]] = GetVI(n,pos[0],pos[1],pos[2]);
					ii[n][pos[0]][pos[1]][pos[2]] = GetII(n,pos[0],pos[1],pos[2]);
					iv[n][pos[0]][pos[1]][pos[2]] = GetIV(n,pos[0],pos[1],pos[2]);
				}
			}
		}
	}

	string names[] = {"vv", "vi", "iv" , "ii", "exc"};
	FDTD_FLOAT**** array[] = {vv,vi,iv,ii,exc};

	ProcessFields::DumpMultiVectorArray2VTK(file, names , array , 5, discLines, numLines,  6, "Operator dump" , (ProcessFields::MeshType)m_MeshType, discLines_scaling);

	Delete_N_3DArray(exc,numLines);
	Delete_N_3DArray(vv,numLines);vv=NULL;
	Delete_N_3DArray(vi,numLines);vi=NULL;
	Delete_N_3DArray(iv,numLines);iv=NULL;
	Delete_N_3DArray(ii,numLines);ii=NULL;

	file.close();

	cout << " done!" << endl;
}

