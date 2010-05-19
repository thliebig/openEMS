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

#include "operator_sse_compressed.h"
#include "engine_sse_compressed.h"
#include "engine_sse.h"
#include "tools/array_ops.h"

#include <list>

Operator_SSE_Compressed* Operator_SSE_Compressed::New()
{
	cout << "Create FDTD operator (compressed SSE)" << endl;
	Operator_SSE_Compressed* op = new Operator_SSE_Compressed();
	op->Init();
	return op;
}

Operator_SSE_Compressed::Operator_SSE_Compressed() : Operator_sse()
{
	m_Op_index = NULL;
	m_Use_Compression = false;
	m_Max_Compression = false;

	m_Compression_Threshold = 0.8;

	for (int n=0;n<3;++n)
	{
		f4_vv_Compressed[n]=NULL;
		f4_vi_Compressed[n]=NULL;
		f4_ii_Compressed[n]=NULL;
		f4_iv_Compressed[n]=NULL;
	}

	m_max_fifo = 1000;
}

Operator_SSE_Compressed::~Operator_SSE_Compressed()
{
	Reset();
}

Engine* Operator_SSE_Compressed::CreateEngine() const
{
	if (m_Use_Compression==false)
	{
		//!create a default sse-engine
		Engine_sse* eng = Engine_sse::New(this);
		return eng;
	}
	Engine_SSE_Compressed* eng = Engine_SSE_Compressed::New(this);
	return eng;
}

int Operator_SSE_Compressed::CalcECOperator()
{
	Operator_sse::CalcECOperator();
	m_Use_Compression = CompressOperator();

	return 0;
}

void Operator_SSE_Compressed::Init()
{
	Operator_sse::Init();
	m_Op_index = NULL;
	for (int n=0;n<3;++n)
	{
		f4_vv_Compressed[n]=NULL;
		f4_vi_Compressed[n]=NULL;
		f4_ii_Compressed[n]=NULL;
		f4_iv_Compressed[n]=NULL;
	}
}

void Operator_SSE_Compressed::Reset()
{
	Operator_sse::Reset();

	if (m_Op_index)
	{
		unsigned int pos[3];
		for (pos[0]=0;pos[0]<numLines[0];++pos[0])
		{
			for (pos[1]=0;pos[1]<numLines[1];++pos[1])
			{
				delete[] m_Op_index[pos[0]][pos[1]];
			}
			delete[] m_Op_index[pos[0]];
		}
		delete[] m_Op_index;
		m_Op_index = NULL;
	}

	for (int n=0;n<3;++n)
	{
		Delete1DArray_v4sf(f4_vv_Compressed[n]);
		f4_vv_Compressed[n]=NULL;
		Delete1DArray_v4sf(f4_vi_Compressed[n]);
		f4_vi_Compressed[n]=NULL;
		Delete1DArray_v4sf(f4_ii_Compressed[n]);
		f4_ii_Compressed[n]=NULL;
		Delete1DArray_v4sf(f4_iv_Compressed[n]);
		f4_iv_Compressed[n]=NULL;
	}
}

void Operator_SSE_Compressed::InitOperator()
{
	Operator_sse::InitOperator();

	unsigned int pos[3];
	m_Op_index = new unsigned int**[numLines[0]];
	for (pos[0]=0;pos[0]<numLines[0];++pos[0])
	{
		m_Op_index[pos[0]] = new unsigned int*[numLines[1]];
		for (pos[1]=0;pos[1]<numLines[1];++pos[1])
		{
			m_Op_index[pos[0]][pos[1]] = new unsigned int[numVectors];
			for (pos[2]=0;pos[2]<numVectors;++pos[2])
			{
				m_Op_index[pos[0]][pos[1]][pos[2]] = 0;
			}
		}
	}
}

void Operator_SSE_Compressed::ShowStat() const
{
	Operator_sse::ShowStat();

	cout << "SSE compression enabled\t: " << m_Use_Compression << endl;
	cout << "Unique SSE operator\t: " << m_Op_Count << endl;
	cout << "Optimal compression\t: " << m_Max_Compression << endl;
	cout << "-----------------------------------" << endl;
}


bool Operator_SSE_Compressed::CompareOperators(unsigned int pos1[3], unsigned int pos2[3])
{
//	cerr << pos1[0] << " " << pos1[1] << " " << pos1[2] << endl;
	for (int n=0;n<3;++n)
	{
		for (int m=0;m<4;++m)
		{
			if (f4_vv[n][pos1[0]][pos1[1]][pos1[2]].f[m] != f4_vv[n][pos2[0]][pos2[1]][pos2[2]].f[m])	return false;
			if (f4_vi[n][pos1[0]][pos1[1]][pos1[2]].f[m] != f4_vi[n][pos2[0]][pos2[1]][pos2[2]].f[m])	return false;
			if (f4_iv[n][pos1[0]][pos1[1]][pos1[2]].f[m] != f4_iv[n][pos2[0]][pos2[1]][pos2[2]].f[m])	return false;
			if (f4_ii[n][pos1[0]][pos1[1]][pos1[2]].f[m] != f4_ii[n][pos2[0]][pos2[1]][pos2[2]].f[m])	return false;
		}
	}
	return true;
}

bool Operator_SSE_Compressed::CompressOperator()
{
	cout << "Compressing the FDTD operator... this my take a while..." << endl;

	m_Max_Compression = true;

	list<unsigned int> fifo;
	vector<unsigned int> index_list[3];

	bool found;
	unsigned int pos[3];
	unsigned int index_pos[3];
	list<unsigned int>::iterator it;
	for (pos[0]=0;pos[0]<numLines[0];++pos[0])
	{
		for (pos[1]=0;pos[1]<numLines[1];++pos[1])
		{
			for (pos[2]=0;pos[2]<numVectors;++pos[2])
			{
				found = false;
				for ( it=fifo.begin() ; it != fifo.end(); it++ )
				{
					index_pos[0] = index_list[0].at(*it);
					index_pos[1] = index_list[1].at(*it);
					index_pos[2] = index_list[2].at(*it);
					found = CompareOperators(index_pos, pos);
					if (found)
					{
						m_Op_index[pos[0]][pos[1]][pos[2]] = *it;
						fifo.erase(it);
						fifo.push_front(*it);	//push already existing value to the front
						it = fifo.end();
						++it;
						break;
					}
				}
				if (found==false)
				{
					fifo.push_front(index_list[0].size());
					m_Op_index[pos[0]][pos[1]][pos[2]] = index_list[0].size();
					index_list[0].push_back(pos[0]);
					index_list[1].push_back(pos[1]);
					index_list[2].push_back(pos[2]);
				}
				if (fifo.size()>m_max_fifo)
				{
					fifo.pop_back();
					m_Max_Compression = false;
				}
			}
		}
	}

	m_Op_Count = index_list[0].size();
	if ( ((double)m_Op_Count/(double)GetNumberCells()) > m_Compression_Threshold)
	{
		cerr << "Operator_SSE_Compressed::CompressOperator: Warning: Compression unsuccessful, ratio is " << (double)m_Op_Count/(double)GetNumberCells() << "! Fallback to conventional sse-engine..." << endl;
		return false;
	}
//	cerr << "total found: " << index_list[0].size() << endl;

	for (int n=0;n<3;++n)
	{
		f4_vv_Compressed[n] = Create1DArray_v4sf(m_Op_Count);
		f4_vi_Compressed[n] = Create1DArray_v4sf(m_Op_Count);
		f4_ii_Compressed[n] = Create1DArray_v4sf(m_Op_Count);
		f4_iv_Compressed[n] = Create1DArray_v4sf(m_Op_Count);

		for (unsigned int m=0;m<m_Op_Count;++m)
		{
			for (unsigned int v=0;v<4;++v)
			{
				f4_vv_Compressed[n][m].f[v] = f4_vv[n][index_list[0].at(m)][index_list[1].at(m)][index_list[2].at(m)].f[v];
				f4_vi_Compressed[n][m].f[v] = f4_vi[n][index_list[0].at(m)][index_list[1].at(m)][index_list[2].at(m)].f[v];
				f4_ii_Compressed[n][m].f[v] = f4_ii[n][index_list[0].at(m)][index_list[1].at(m)][index_list[2].at(m)].f[v];
				f4_iv_Compressed[n][m].f[v] = f4_iv[n][index_list[0].at(m)][index_list[1].at(m)][index_list[2].at(m)].f[v];
			}
		}
	}

	Delete_N_3DArray_v4sf(f4_vv,numLines);
	Delete_N_3DArray_v4sf(f4_vi,numLines);
	Delete_N_3DArray_v4sf(f4_iv,numLines);
	Delete_N_3DArray_v4sf(f4_ii,numLines);
	f4_vv = 0;
	f4_vi = 0;
	f4_iv = 0;
	f4_ii = 0;

	return true;
}
