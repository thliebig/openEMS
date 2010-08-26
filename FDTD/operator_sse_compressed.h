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

#ifndef OPERATOR_SSE_COMPRESSED_H
#define OPERATOR_SSE_COMPRESSED_H

#include "operator_sse.h"

class SSE_coeff
{
public:
        SSE_coeff( f4vector vv[3], f4vector vi[3], f4vector iv[3], f4vector ii[3] );
        bool operator==( const SSE_coeff& ) const;
        bool operator!=( const SSE_coeff& ) const;
        bool operator<( const SSE_coeff& ) const;
        void print( ostream& stream ) const;
protected:
        f4vector m_vv[3];
        f4vector m_vi[3];
        f4vector m_iv[3];
        f4vector m_ii[3];
};

class Operator_SSE_Compressed : public Operator_sse
{
public:
	//! Create a new operator
	static Operator_SSE_Compressed* New();
	virtual ~Operator_SSE_Compressed();

	virtual Engine* CreateEngine() const;

	virtual int CalcECOperator();

        inline virtual FDTD_FLOAT& GetVV( unsigned int n, unsigned int x, unsigned int y, unsigned int z ) { if (m_Use_Compression) return f4_vv_Compressed[n][m_Op_index[x][y][z%numVectors]].f[z/numVectors]; else return Operator_sse::GetVV(n,x,y,z);}
        inline virtual FDTD_FLOAT& GetVI( unsigned int n, unsigned int x, unsigned int y, unsigned int z ) { if (m_Use_Compression) return f4_vi_Compressed[n][m_Op_index[x][y][z%numVectors]].f[z/numVectors]; else return Operator_sse::GetVI(n,x,y,z);}

        inline virtual FDTD_FLOAT& GetII( unsigned int n, unsigned int x, unsigned int y, unsigned int z ) { if (m_Use_Compression) return f4_ii_Compressed[n][m_Op_index[x][y][z%numVectors]].f[z/numVectors]; else return Operator_sse::GetII(n,x,y,z);}
        inline virtual FDTD_FLOAT& GetIV( unsigned int n, unsigned int x, unsigned int y, unsigned int z ) { if (m_Use_Compression) return f4_iv_Compressed[n][m_Op_index[x][y][z%numVectors]].f[z/numVectors]; else return Operator_sse::GetIV(n,x,y,z);}

	virtual void ShowStat() const;

	bool CompressOperator();

protected:
	Operator_SSE_Compressed();

	bool m_Use_Compression;

	virtual void Init();
	virtual void Reset();
	virtual void InitOperator();

	bool CompareOperators(unsigned int pos1[3], unsigned int pos2[3]);

	// engine needs access
public:
	unsigned int*** m_Op_index;
        vector<f4vector> f4_vv_Compressed[3]; //!< coefficient: calc new voltage from old voltage
        vector<f4vector> f4_vi_Compressed[3]; //!< coefficient: calc new voltage from old current
        vector<f4vector> f4_iv_Compressed[3]; //!< coefficient: calc new current from old current
        vector<f4vector> f4_ii_Compressed[3]; //!< coefficient: calc new current from old voltage

};

#endif // OPERATOR_SSE_Compressed_H
