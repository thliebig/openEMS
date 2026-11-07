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

#include <algorithm>
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

Coeffs_EC::Coeffs_EC()
{
	m_Lines[0] = m_Lines[1] = m_Lines[2] = 0;
	m_dT = 0;
	for (int e=0; e<4; ++e)
		for (int n=0; n<3; ++n)
			m_EC[e][n] = NULL;
}

Coeffs_EC::~Coeffs_EC()
{
	Free();
}

void Coeffs_EC::Init(const unsigned int numLines[3])
{
	Free();
	for (int n=0; n<3; ++n)
		m_Lines[n] = numLines[n];
	m_Flags.assign((size_t)numLines[0]*numLines[1]*numLines[2], 0);
}

void Coeffs_EC::Free()
{
	std::vector<uint16_t>().swap(m_Flags);
	std::vector<Box>().swap(m_Boxes);
	std::unordered_map<size_t, size_t>().swap(m_MapIndex);
	std::vector<std::array<float,12> >().swap(m_MapValues);
	for (int e=0; e<4; ++e)
		for (int n=0; n<3; ++n)
			m_EC[e][n] = NULL;
}

void Coeffs_EC::ShowStat(std::ostream& ostr) const
{
	size_t box = 0;
	for (const Box& b : m_Boxes)
		box += b.values.size()/12;
	ostr << "Coefficients\t\t: computed from the circuit values; stored for " << box << " nodes in boxes (UPML), " << m_MapValues.size() << " in a map" << endl;
}

void Coeffs_EC::SetEC(FDTD_FLOAT* const C[3], FDTD_FLOAT* const G[3], FDTD_FLOAT* const L[3], FDTD_FLOAT* const R[3], double dT)
{
	for (int n=0; n<3; ++n)
	{
		m_EC[0][n] = C[n];
		m_EC[1][n] = G[n];
		m_EC[2][n] = L[n];
		m_EC[3][n] = R[n];
	}
	m_dT = dT;
}

void Coeffs_EC::ReserveBox(const unsigned int start[3], const unsigned int lines[3])
{
	// the box index has 3 bits in the flags: beyond MAX_BOXES the nodes go to the map (a UPML has at most 6 regions)
	if (m_Boxes.size()>=MAX_BOXES)
		return;
	Box b;
	size_t nodes = 1;
	for (int n=0; n<3; ++n)
	{
		if (start[n]+lines[n]>m_Lines[n])
			throw std::out_of_range("Coeffs_EC::ReserveBox: the box exceeds the mesh");
		b.start[n] = start[n];
		b.lines[n] = lines[n];
		nodes *= lines[n];
	}
	if (nodes==0)
		return;
	b.values.resize(12*nodes);
	m_Boxes.push_back(std::move(b));
}

const float* Coeffs_EC::Own(size_t i, unsigned int x, unsigned int y, unsigned int z) const
{
	return const_cast<Coeffs_EC*>(this)->Own(i, x, y, z);
}

float* Coeffs_EC::Own(size_t i, unsigned int x, unsigned int y, unsigned int z)
{
	const unsigned int box = m_Flags[i]>>BOX_SHIFT;
	if (box==0)
		return m_MapValues[m_MapIndex.at(i)].data();
	Box& b = m_Boxes[box-1];
	const size_t l = (x-b.start[0]) + (size_t)b.lines[0]*((y-b.start[1]) + (size_t)b.lines[1]*(z-b.start[2]));
	return &b.values[12*l];
}

float* Coeffs_EC::MakeOwn(size_t i, unsigned int x, unsigned int y, unsigned int z)
{
	float values[12];
	for (int k=0; k<12; ++k)
		values[k] = (m_Flags[i] & (1u<<k)) ? 0 : Compute(k, i);
	const unsigned int p[3] = {x, y, z};
	for (size_t n=0; n<m_Boxes.size(); ++n)
	{
		const Box& b = m_Boxes[n];
		bool inside = true;
		for (int a=0; a<3 && inside; ++a)
			inside = (p[a]-b.start[a])<b.lines[a];   // wraps below the start
		if (inside)
		{
			m_Flags[i] = OWN | (uint16_t)((n+1)<<BOX_SHIFT);
			float* own = Own(i, x, y, z);
			std::copy(values, values+12, own);
			return own;
		}
	}
	// the map is not locked: a single thread only (the class doc)
	m_MapIndex[i] = m_MapValues.size();
	m_MapValues.emplace_back();
	std::copy(values, values+12, m_MapValues.back().begin());
	m_Flags[i] = OWN;
	return m_MapValues.back().data();
}

void Coeffs_EC::GetPair(int f, unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT& self, FDTD_FLOAT& curl) const
{
	const size_t i = Index(x, y, z);
	const uint16_t flags = m_Flags[i];
	const int k = 4*n + 2*f;
	if (flags & OWN)
	{
		const float* own = Own(i, x, y, z);
		self = own[k];
		curl = own[k+1];
		return;
	}
	Coeffs_FromEC(m_EC[2*f][n][i], m_EC[2*f+1][n][i], m_dT, self, curl);
	if (flags & (1u<<k))
		self = 0;
	if (flags & (1u<<(k+1)))
		curl = 0;
}

void Coeffs_EC::Set(int k, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value)
{
	const size_t i = Index(x, y, z);
	const uint16_t f = m_Flags[i];
	if (f & OWN)
		Own(i, x, y, z)[k] = value;
	else if (value==0)
		m_Flags[i] = f | (1u<<k);
	else
		MakeOwn(i, x, y, z)[k] = value;
}

void Coeffs_EC::Reset(int n, unsigned int x, unsigned int y, unsigned int z)
{
	const size_t i = Index(x, y, z);
	if (m_Flags[i] & OWN)
	{
		float* own = Own(i, x, y, z);
		for (int c=0; c<4; ++c)
			own[4*n+c] = Compute(4*n+c, i);
	}
	else
		m_Flags[i] &= ~(uint16_t)(0xF<<(4*n));
}

void Coeffs_EC::Freeze(unsigned int x, unsigned int y, unsigned int z)
{
	const size_t i = Index(x, y, z);
	if (!(m_Flags[i] & OWN))
		MakeOwn(i, x, y, z);
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
