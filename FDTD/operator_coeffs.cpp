/*
*	Copyright (C) 2010-2026 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#include <cmath>
#include <cstring>
#include <iostream>
#include <map>
#include <stdexcept>

#include "operator_coeffs.h"
#include "tools/global.h"

using std::cout;
using std::endl;

// ----------------------------------------------------------------------------

Coeffs_Basic::Coeffs_Basic()
{
	vv = NULL;
	vi = NULL;
	iv = NULL;
	ii = NULL;
}

Coeffs_Basic::~Coeffs_Basic()
{
	Free();
}

void Coeffs_Basic::Init(const unsigned int numLines[3])
{
	Coeffs_Basic::Free();
	unsigned int lines[3] = {numLines[0], numLines[1], numLines[2]};
	vv = new ArrayLib::ArrayNIJK<FDTD_FLOAT>("vv", lines);
	vi = new ArrayLib::ArrayNIJK<FDTD_FLOAT>("vi", lines);
	iv = new ArrayLib::ArrayNIJK<FDTD_FLOAT>("iv", lines);
	ii = new ArrayLib::ArrayNIJK<FDTD_FLOAT>("ii", lines);
}

void Coeffs_Basic::Free()
{
	delete vv;
	delete vi;
	delete iv;
	delete ii;
	vv = NULL;
	vi = NULL;
	iv = NULL;
	ii = NULL;
}

// ----------------------------------------------------------------------------

//! the 24 f4vector coefficients of one node, the key for finding identical ones (Coeffs_SSE_Compressed::Finalize)
class SSE_coeff
{
public:
	SSE_coeff( f4vector vv[3], f4vector vi[3], f4vector iv[3], f4vector ii[3] );
	bool operator==( const SSE_coeff& ) const;
	bool operator!=( const SSE_coeff& ) const;
	bool operator<( const SSE_coeff& ) const;
	void print( std::ostream& stream ) const;
protected:
	f4vector m_vv[3];
	f4vector m_vi[3];
	f4vector m_iv[3];
	f4vector m_ii[3];
};

Coeffs_SSE_Compressed::Coeffs_SSE_Compressed() : Coeffs_Basic()
{
	numVectors = 0;
	use_compression = false;
	m_numLines[0] = m_numLines[1] = m_numLines[2] = 0;
}

Coeffs_SSE_Compressed::~Coeffs_SSE_Compressed()
{
	ClearCompressed();
}

void Coeffs_SSE_Compressed::ClearCompressed()
{
	use_compression = false;
	for (int n=0; n<3; n++)
	{
		f4_vv_Compressed[n].clear();
		f4_vi_Compressed[n].clear();
		f4_iv_Compressed[n].clear();
		f4_ii_Compressed[n].clear();
	}
}

void Coeffs_SSE_Compressed::Init(const unsigned int numLines[3])
{
	ClearCompressed();
	index.Reset();
	Coeffs_Basic::Init(numLines);
	for (int n=0; n<3; ++n)
		m_numLines[n] = numLines[n];
	numVectors = ceil((double)numLines[2]/4.0);
}

void Coeffs_SSE_Compressed::Free()
{
	ClearCompressed();
	index.Reset();
	Coeffs_Basic::Free();
}

void Coeffs_SSE_Compressed::ThrowCompressed()
{
	throw std::logic_error("Coeffs_SSE_Compressed: the coefficients are compressed, a Set function would change all nodes sharing the set");
}

void Coeffs_SSE_Compressed::ShowStat(std::ostream& ostr) const
{
	ostr << "SSE compression enabled\t: " << (use_compression?"yes":"no") << endl;
	ostr << "Unique SSE operators\t: " << f4_vv_Compressed->size() << endl;
	ostr << "-----------------------------------" << endl;
}

void Coeffs_SSE_Compressed::Finalize()
{
	use_compression = false;

	if (g_settings.GetVerboseLevel()>0)
		cout << "Compressing the FDTD operator... this may take a while..." << endl;

	unsigned int idx_lines[3] = {m_numLines[0], m_numLines[1], numVectors};
	index.Init("Op_index", idx_lines);

	// the coefficients of the four z-lines z+k*numVectors in lane k, lanes beyond the last line 0
	ArrayLib::ArrayNIJK<FDTD_FLOAT>* full[4] = {vv, vi, iv, ii};
	auto gather = [&](int c, int n, const unsigned int pos[3])
	{
		f4vector f;
		for (unsigned int k=0; k<4; ++k)
		{
			unsigned int z = pos[2] + k*numVectors;
			f.f[k] = (z<m_numLines[2]) ? (*full[c])(n, pos[0], pos[1], z) : 0;
		}
		return f;
	};

	std::map<SSE_coeff,unsigned int> lookUpMap;

	unsigned int pos[3];
	for (pos[0]=0; pos[0]<m_numLines[0]; ++pos[0])
	{
		for (pos[1]=0; pos[1]<m_numLines[1]; ++pos[1])
		{
			for (pos[2]=0; pos[2]<numVectors; ++pos[2])
			{
				f4vector vv4[3] = { gather(0, 0, pos), gather(0, 1, pos), gather(0, 2, pos) };
				f4vector vi4[3] = { gather(1, 0, pos), gather(1, 1, pos), gather(1, 2, pos) };
				f4vector iv4[3] = { gather(2, 0, pos), gather(2, 1, pos), gather(2, 2, pos) };
				f4vector ii4[3] = { gather(3, 0, pos), gather(3, 1, pos), gather(3, 2, pos) };
				SSE_coeff c( vv4, vi4, iv4, ii4 );

				std::map<SSE_coeff,unsigned int>::iterator it;
				it = lookUpMap.find(c);
				if (it == lookUpMap.end())
				{
					// not found -> insert
					unsigned int idx = f4_vv_Compressed[0].size();
					for (int n=0; n<3; n++)
					{
						f4_vv_Compressed[n].push_back( vv4[n] );
						f4_vi_Compressed[n].push_back( vi4[n] );
						f4_iv_Compressed[n].push_back( iv4[n] );
						f4_ii_Compressed[n].push_back( ii4[n] );
					}
					lookUpMap[c] = idx;
					index(pos[0], pos[1], pos[2]) = idx;
				}
				else
				{
					// this operator is already in the list
					index(pos[0], pos[1], pos[2]) = (*it).second;
				}
			}
		}
	}

	Coeffs_Basic::Free();
	use_compression = true;
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

void SSE_coeff::print( std::ostream& stream ) const
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
