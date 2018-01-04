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

#include <map>
#include <cstring>

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
}

Operator_SSE_Compressed::~Operator_SSE_Compressed()
{
	Delete();
}

Engine* Operator_SSE_Compressed::CreateEngine()
{
	if (!m_Use_Compression)
		//! create a default sse-engine
		m_Engine = Engine_sse::New(this);
  else
	  m_Engine = Engine_SSE_Compressed::New(this);
	return m_Engine;
}

int Operator_SSE_Compressed::CalcECOperator( DebugFlags debugFlags )
{
	int ErrCode = Operator_sse::CalcECOperator( debugFlags );
	m_Use_Compression = false;
	m_Use_Compression = CompressOperator();

	return ErrCode;
}

void Operator_SSE_Compressed::Init()
{
	Operator_sse::Init();
	m_Use_Compression = false;
	m_Op_index = NULL;
}

void Operator_SSE_Compressed::Delete()
{
	if (m_Op_index)
	{
		Delete3DArray<unsigned int>( m_Op_index, numLines );
		m_Op_index = 0;
	}

	m_Use_Compression = false;
	for (int n=0; n<3; n++)
	{
		f4_vv_Compressed[n].clear();
		f4_vi_Compressed[n].clear();
		f4_iv_Compressed[n].clear();
		f4_ii_Compressed[n].clear();
	}
}

void Operator_SSE_Compressed::Reset()
{
	Delete();
	Operator_sse::Reset();
}

void Operator_SSE_Compressed::InitOperator()
{
	//cleanup compression
	m_Use_Compression = false;
	for (int n=0; n<3; n++)
	{
		f4_vv_Compressed[n].clear();
		f4_vi_Compressed[n].clear();
		f4_iv_Compressed[n].clear();
		f4_ii_Compressed[n].clear();
	}

	Operator_sse::InitOperator();
	m_Op_index = Create3DArray<unsigned int>( numLines );
}

void Operator_SSE_Compressed::ShowStat() const
{
	Operator_sse::ShowStat();

	cout << "SSE compression enabled\t: " << (m_Use_Compression?"yes":"no") << endl;
	cout << "Unique SSE operators\t: " << f4_vv_Compressed->size() << endl;
	cout << "-----------------------------------" << endl;
}

bool Operator_SSE_Compressed::CompressOperator()
{
	if (g_settings.GetVerboseLevel()>0)
		cout << "Compressing the FDTD operator... this may take a while..." << endl;

	map<SSE_coeff,unsigned int> lookUpMap;

	unsigned int pos[3];
	for (pos[0]=0; pos[0]<numLines[0]; ++pos[0])
	{
		for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
		{
			for (pos[2]=0; pos[2]<numVectors; ++pos[2])
			{
				f4vector vv[3] = { f4_vv[0][pos[0]][pos[1]][pos[2]], f4_vv[1][pos[0]][pos[1]][pos[2]], f4_vv[2][pos[0]][pos[1]][pos[2]] };
				f4vector vi[3] = { f4_vi[0][pos[0]][pos[1]][pos[2]], f4_vi[1][pos[0]][pos[1]][pos[2]], f4_vi[2][pos[0]][pos[1]][pos[2]] };
				f4vector iv[3] = { f4_iv[0][pos[0]][pos[1]][pos[2]], f4_iv[1][pos[0]][pos[1]][pos[2]], f4_iv[2][pos[0]][pos[1]][pos[2]] };
				f4vector ii[3] = { f4_ii[0][pos[0]][pos[1]][pos[2]], f4_ii[1][pos[0]][pos[1]][pos[2]], f4_ii[2][pos[0]][pos[1]][pos[2]] };
				SSE_coeff c( vv, vi, iv, ii );

				map<SSE_coeff,unsigned int>::iterator it;
				it = lookUpMap.find(c);
				if (it == lookUpMap.end())
				{
					// not found -> insert
					unsigned int index = f4_vv_Compressed[0].size();
					for (int n=0; n<3; n++)
					{
						f4_vv_Compressed[n].push_back( vv[n] );
						f4_vi_Compressed[n].push_back( vi[n] );
						f4_iv_Compressed[n].push_back( iv[n] );
						f4_ii_Compressed[n].push_back( ii[n] );
					}
					lookUpMap[c] = index;
					m_Op_index[pos[0]][pos[1]][pos[2]] = index;
				}
				else
				{
					// this operator is already in the list
					unsigned int index = (*it).second;
					m_Op_index[pos[0]][pos[1]][pos[2]] = index;
				}
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




// ----------------------------------------------------------------------------

SSE_coeff::SSE_coeff( f4vector vv[3], f4vector vi[3], f4vector iv[3], f4vector ii[3] )
{
	for (int n=0; n<3; n++)
	{
		m_vv[n] = vv[n];
		m_vi[n] = vi[n];
		m_iv[n] = iv[n];
		m_ii[n] = ii[n];
	}
}

bool SSE_coeff::operator==( const SSE_coeff& other ) const
{
	for (int n=0; n<3; n++)
	{
		if (memcmp( &(m_vv[n]), &(other.m_vv[n]), sizeof(f4vector) ) != 0) return false;
		if (memcmp( &(m_vi[n]), &(other.m_vi[n]), sizeof(f4vector) ) != 0) return false;
		if (memcmp( &(m_iv[n]), &(other.m_iv[n]), sizeof(f4vector) ) != 0) return false;
		if (memcmp( &(m_ii[n]), &(other.m_ii[n]), sizeof(f4vector) ) != 0) return false;
	}
	return true;
}
bool SSE_coeff::operator!=( const SSE_coeff& other ) const
{
	return !(*this == other);
}
bool SSE_coeff::operator<( const SSE_coeff& other ) const
{
	for (int n=0; n<3; n++)
	{
		for (int c=0; c<4; c++)
		{
			if (m_vv[n].f[c] > other.m_vv[n].f[c]) return false;
			if (m_vv[n].f[c] < other.m_vv[n].f[c]) return true;
			if (m_vi[n].f[c] > other.m_vi[n].f[c]) return false;
			if (m_vi[n].f[c] < other.m_vi[n].f[c]) return true;
			if (m_iv[n].f[c] > other.m_iv[n].f[c]) return false;
			if (m_iv[n].f[c] < other.m_iv[n].f[c]) return true;
			if (m_ii[n].f[c] > other.m_ii[n].f[c]) return false;
			if (m_ii[n].f[c] < other.m_ii[n].f[c]) return true;
		}
	}
	return false;
}

void SSE_coeff::print( ostream& stream ) const
{
	stream << "SSE_coeff: (" << endl;
	for (int n=0; n<3; n++)
	{
		stream << "n=" << n << ":" << endl;
		stream << "vv=";
		for (int c=0; c<4; c++)
			stream << m_vv[n].f[c] << " ";
		stream << endl << "vi=";
		for (int c=0; c<4; c++)
			stream << m_vi[n].f[c] << " ";
		stream << endl << "iv=";
		for (int c=0; c<4; c++)
			stream << m_iv[n].f[c] << " ";
		stream << endl << "ii=";
		for (int c=0; c<4; c++)
			stream << m_ii[n].f[c] << " ";
		stream << endl;
	}
	stream << ")" << endl;
}
