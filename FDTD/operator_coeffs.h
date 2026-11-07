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

#include <array>
#include <cstdint>
#include <ostream>
#include <unordered_map>
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

	//! The voltage (\a f 0: vv, vi) or current (1: ii, iv) coefficients of direction \a n at a node, in one go
	virtual void GetPair(int f, unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT& self, FDTD_FLOAT& curl) const
	{
		self = (f==0) ? GetVV(n, x, y, z) : GetII(n, x, y, z);
		curl = (f==0) ? GetVI(n, x, y, z) : GetIV(n, x, y, z);
	}

	//! Hint: the Set functions will give most nodes of the box from \a start, \a lines lines per direction, values of their own (the UPML)
	virtual void ReserveBox(const unsigned int start[3], const unsigned int lines[3]) {(void)start; (void)lines;}
};

//! The coefficients of one edge from its circuit values: (vv, vi) from (C, G), or (ii, iv) from (L, R); both 0 without C or L
inline void Coeffs_FromEC(double C, double G, double dT, FDTD_FLOAT& self, FDTD_FLOAT& curl)
{
	if (C>0)
	{
		self = (1.0-dT*G/2.0/C)/(1.0+dT*G/2.0/C);
		curl = (dT/C)/(1.0+dT*G/2.0/C);
	}
	else
	{
		self = 0;
		curl = 0;
	}
}

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

//! The coefficients computed from the circuit values when read (Engine_CPU, Engine_GPU)
/*!
  These engines read every coefficient once to build their own form of them
  (the cell runs, the device arrays). Instead of an array per coefficient (48
  bytes per node) the coefficients are computed from the circuit values C, G,
  L, R of the operator and the timestep when read, as
  Operator::Calc_ECOperatorPos() does; the operator keeps the circuit values
  until the engine is built.

  What is set afterwards is kept per node: a bit per coefficient set to 0 (the
  PEC, PMC and boundaries), and all 12 coefficients of a node with any other
  value, in a box reserved for them (ReserveBox(), the UPML) or else in a map.

  The circuit values of an edge must not change after this unless
  Reset() follows (Operator::Calc_ECOperatorPos()), or Freeze() came first.
  The Set functions of different nodes may run on several threads only while
  all of them are in reserved boxes (the UPML build): a node that goes to the
  map changes the map, which is not locked. The Get functions may run on
  several threads while nothing is set.
  */
class Coeffs_EC : public Operator_Coeffs
{
public:
	Coeffs_EC();
	virtual ~Coeffs_EC();

	virtual void Init(const unsigned int numLines[3]);
	virtual void Free();
	virtual void ShowStat(std::ostream& ostr) const;
	virtual void ReserveBox(const unsigned int start[3], const unsigned int lines[3]);

	//! Compute the coefficients from these circuit values (per direction, indexed z + Nz*(y + Ny*x), the order of MainOp) and the timestep \a dT
	void SetEC(FDTD_FLOAT* const C[3], FDTD_FLOAT* const G[3], FDTD_FLOAT* const L[3], FDTD_FLOAT* const R[3], double dT);
	//! The coefficients of the edge \a n at the node from the circuit values again (Operator::Calc_ECOperatorPos())
	void Reset(int n, unsigned int x, unsigned int y, unsigned int z);
	//! Keep the coefficients of the node as they are, before its circuit values change
	void Freeze(unsigned int x, unsigned int y, unsigned int z);
	virtual void GetPair(int f, unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT& self, FDTD_FLOAT& curl) const;

	inline virtual FDTD_FLOAT GetVV(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const {return Get(4*n+VV, x, y, z);}
	inline virtual FDTD_FLOAT GetVI(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const {return Get(4*n+VI, x, y, z);}
	inline virtual FDTD_FLOAT GetII(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const {return Get(4*n+II, x, y, z);}
	inline virtual FDTD_FLOAT GetIV(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const {return Get(4*n+IV, x, y, z);}

	inline virtual void SetVV(unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value) {Set(4*n+VV, x, y, z, value);}
	inline virtual void SetVI(unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value) {Set(4*n+VI, x, y, z, value);}
	inline virtual void SetII(unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value) {Set(4*n+II, x, y, z, value);}
	inline virtual void SetIV(unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value) {Set(4*n+IV, x, y, z, value);}

protected:
	//! the coefficient k = 4*n + c of a node
	enum Coeff {VV=0, VI=1, II=2, IV=3};
	//! the flag of a node with coefficients of its own, the bits below it are the coefficients set to 0
	static const uint16_t OWN = 1<<12;
	//! the bits above OWN: the box of a node of its own plus one, 0 for the map
	static const int BOX_SHIFT = 13;
	static const size_t MAX_BOXES = 7;

	struct Box
	{
		unsigned int start[3], lines[3];
		std::vector<float> values;   //!< 12 per node
	};

	unsigned int m_Lines[3];
	double m_dT;
	FDTD_FLOAT* m_EC[4][3];        //!< C, G, L, R
	std::vector<uint16_t> m_Flags;  //!< per node
	std::vector<Box> m_Boxes;
	std::unordered_map<size_t, size_t> m_MapIndex;   //!< node -> m_MapValues
	std::vector<std::array<float,12> > m_MapValues;

	size_t Index(unsigned int x, unsigned int y, unsigned int z) const {return z + (size_t)m_Lines[2]*(y + (size_t)m_Lines[1]*x);}
	//! the coefficient k from the circuit values
	FDTD_FLOAT Compute(int k, size_t i) const
	{
		const int n = k>>2, c = k&3;
		FDTD_FLOAT self, curl;
		if (c<II)
			Coeffs_FromEC(m_EC[0][n][i], m_EC[1][n][i], m_dT, self, curl);
		else
			Coeffs_FromEC(m_EC[2][n][i], m_EC[3][n][i], m_dT, self, curl);
		return (c&1) ? curl : self;
	}
	//! the coefficients of a node of its own (OWN)
	const float* Own(size_t i, unsigned int x, unsigned int y, unsigned int z) const;
	float* Own(size_t i, unsigned int x, unsigned int y, unsigned int z);
	//! make the node one of its own, with its current coefficients
	float* MakeOwn(size_t i, unsigned int x, unsigned int y, unsigned int z);

	inline FDTD_FLOAT Get(int k, unsigned int x, unsigned int y, unsigned int z) const
	{
		const size_t i = Index(x, y, z);
		const uint16_t f = m_Flags[i];
		if (f & OWN)
			return Own(i, x, y, z)[k];
		if (f & (1u<<k))
			return 0;
		return Compute(k, i);
	}
	void Set(int k, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value);
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
