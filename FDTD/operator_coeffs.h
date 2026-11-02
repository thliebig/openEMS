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

#ifndef OPERATOR_COEFFS_H
#define OPERATOR_COEFFS_H

#include <ostream>
#include <vector>

#include "tools/global.h"
#include "tools/constants.h"
#include "tools/array_ops.h"
#include "tools/aligned_allocator.h"
#include "tools/arraylib/array_nijk.h"
#include "tools/arraylib/array_ijk.h"

//! The update coefficients of the operator: vv, vi (voltage update) and ii, iv (current update) per node and direction
/*!
  The operator builds them through the Set functions, whatever the layout; an
  engine reads the layout it was written for directly, see Operator::GetCoeffs().
  */
class Operator_Coeffs
{
public:
	virtual ~Operator_Coeffs() {}

	//! Allocate the coefficients for \a numLines mesh lines, dropping any old ones
	virtual void Init(const unsigned int numLines[3]) = 0;
	//! Free the coefficients; the Get and Set functions are invalid after this
	virtual void Free() = 0;
	//! Called once the operator and its extensions are built
	virtual void Finalize() {}
	virtual void ShowStat(std::ostream& ostr) const {(void)ostr;}

	virtual FDTD_FLOAT GetVV(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const = 0;
	virtual FDTD_FLOAT GetVI(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const = 0;
	virtual FDTD_FLOAT GetII(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const = 0;
	virtual FDTD_FLOAT GetIV(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const = 0;

	virtual void SetVV(unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value) = 0;
	virtual void SetVI(unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value) = 0;
	virtual void SetII(unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value) = 0;
	virtual void SetIV(unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value) = 0;
};

//! One float per coefficient, node and direction (Engine, and the host side of Engine_GPU)
class Coeffs_Basic : public Operator_Coeffs
{
public:
	Coeffs_Basic();
	virtual ~Coeffs_Basic();

	virtual void Init(const unsigned int numLines[3]);
	virtual void Free();

	inline virtual FDTD_FLOAT GetVV(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const {return (*vv)(n, x, y, z);}
	inline virtual FDTD_FLOAT GetVI(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const {return (*vi)(n, x, y, z);}
	inline virtual FDTD_FLOAT GetII(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const {return (*ii)(n, x, y, z);}
	inline virtual FDTD_FLOAT GetIV(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const {return (*iv)(n, x, y, z);}

	inline virtual void SetVV(unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value) {(*vv)(n, x, y, z) = value;}
	inline virtual void SetVI(unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value) {(*vi)(n, x, y, z) = value;}
	inline virtual void SetII(unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value) {(*ii)(n, x, y, z) = value;}
	inline virtual void SetIV(unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value) {(*iv)(n, x, y, z) = value;}

	ArrayLib::ArrayNIJK<FDTD_FLOAT>* vv; //!< calc new voltage from old voltage
	ArrayLib::ArrayNIJK<FDTD_FLOAT>* vi; //!< calc new voltage from old current
	ArrayLib::ArrayNIJK<FDTD_FLOAT>* iv; //!< calc new current from old voltage
	ArrayLib::ArrayNIJK<FDTD_FLOAT>* ii; //!< calc new current from old current
};

//! The coefficients of Engine_sse: the distinct sets of four z-lines (one f4vector per coefficient and direction) plus an index per node
/*!
  The operator builds into the Coeffs_Basic arrays. Finalize() gathers the
  z-lines z, z+numVectors, z+2*numVectors, z+3*numVectors into the four lanes of
  one f4vector (the field layout of Engine_sse, lanes beyond the last line 0),
  keeps each distinct set of the 24 f4vectors of a node once and frees the basic
  arrays. After that the Get functions read the compressed sets; the Set
  functions throw, as a set is shared by many nodes.
  */
class Coeffs_SSE_Compressed : public Coeffs_Basic
{
public:
	Coeffs_SSE_Compressed();
	virtual ~Coeffs_SSE_Compressed();

	virtual void Init(const unsigned int numLines[3]);
	virtual void Free();
	//! Compress the coefficients and free the basic arrays
	virtual void Finalize();
	virtual void ShowStat(std::ostream& ostr) const;

	bool IsCompressed() const {return use_compression;}

	inline virtual FDTD_FLOAT GetVV(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const {if (use_compression) return f4_vv_Compressed[n][index(x, y, z%numVectors)].f[z/numVectors]; return Coeffs_Basic::GetVV(n, x, y, z);}
	inline virtual FDTD_FLOAT GetVI(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const {if (use_compression) return f4_vi_Compressed[n][index(x, y, z%numVectors)].f[z/numVectors]; return Coeffs_Basic::GetVI(n, x, y, z);}
	inline virtual FDTD_FLOAT GetII(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const {if (use_compression) return f4_ii_Compressed[n][index(x, y, z%numVectors)].f[z/numVectors]; return Coeffs_Basic::GetII(n, x, y, z);}
	inline virtual FDTD_FLOAT GetIV(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const {if (use_compression) return f4_iv_Compressed[n][index(x, y, z%numVectors)].f[z/numVectors]; return Coeffs_Basic::GetIV(n, x, y, z);}

	inline virtual void SetVV(unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value) {if (use_compression) ThrowCompressed(); Coeffs_Basic::SetVV(n, x, y, z, value);}
	inline virtual void SetVI(unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value) {if (use_compression) ThrowCompressed(); Coeffs_Basic::SetVI(n, x, y, z, value);}
	inline virtual void SetII(unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value) {if (use_compression) ThrowCompressed(); Coeffs_Basic::SetII(n, x, y, z, value);}
	inline virtual void SetIV(unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value) {if (use_compression) ThrowCompressed(); Coeffs_Basic::SetIV(n, x, y, z, value);}

	unsigned int numVectors;   //!< f4vectors per z-line, ceil(numLines[2]/4)
	bool use_compression;
	ArrayLib::ArrayIJK<unsigned int> index;   //!< the set of each node, x, y and the f4vector along z
	std::vector<f4vector,aligned_allocator<f4vector> > f4_vv_Compressed[3];
	std::vector<f4vector,aligned_allocator<f4vector> > f4_vi_Compressed[3];
	std::vector<f4vector,aligned_allocator<f4vector> > f4_iv_Compressed[3];
	std::vector<f4vector,aligned_allocator<f4vector> > f4_ii_Compressed[3];

protected:
	unsigned int m_numLines[3];
	void ClearCompressed();
	[[noreturn]] static void ThrowCompressed();
};

#endif // OPERATOR_COEFFS_H
