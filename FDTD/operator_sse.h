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

#ifndef OPERATOR_SSE_H
#define OPERATOR_SSE_H

#include "operator.h"
#include "tools/array_ops.h"

class Operator_sse : public Operator
{
	friend class Engine_Interface_SSE_FDTD;
public:
	//! Create a new operator
	static Operator_sse* New();
	virtual ~Operator_sse();

	virtual Engine* CreateEngine();

	inline virtual FDTD_FLOAT GetVV(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const
	{
		ArrayLib::ArrayNIJK<f4vector>& f4_vv = *f4_vv_ptr;
		return f4_vv[n][x][y][z%numVectors].f[z/numVectors];
	}

	inline virtual FDTD_FLOAT GetVI(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const
	{
		ArrayLib::ArrayNIJK<f4vector>& f4_vi = *f4_vi_ptr;
		return f4_vi[n][x][y][z%numVectors].f[z/numVectors];
	}

	inline virtual FDTD_FLOAT GetII(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const
	{
		ArrayLib::ArrayNIJK<f4vector>& f4_ii = *f4_ii_ptr;
		return f4_ii[n][x][y][z%numVectors].f[z/numVectors];
	}

	inline virtual FDTD_FLOAT GetIV(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const
	{
		ArrayLib::ArrayNIJK<f4vector>& f4_iv = *f4_iv_ptr;
		return f4_iv[n][x][y][z%numVectors].f[z/numVectors];
	}

	inline virtual void SetVV(unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value)
	{
		ArrayLib::ArrayNIJK<f4vector>& f4_vv = *f4_vv_ptr;
		f4_vv[n][x][y][z%numVectors].f[z/numVectors] = value;
	}

	inline virtual void SetVI(unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value)
	{
		ArrayLib::ArrayNIJK<f4vector>& f4_vi = *f4_vi_ptr;
		f4_vi[n][x][y][z%numVectors].f[z/numVectors] = value;
	}

	inline virtual void SetII(unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value)
	{
		ArrayLib::ArrayNIJK<f4vector>& f4_ii = *f4_ii_ptr;
		f4_ii[n][x][y][z%numVectors].f[z/numVectors] = value;
	}

	inline virtual void SetIV(unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value)
	{
		ArrayLib::ArrayNIJK<f4vector>& f4_iv = *f4_iv_ptr;
		f4_iv[n][x][y][z%numVectors].f[z/numVectors] = value;
	}

protected:
	//! use New() for creating a new Operator
	Operator_sse();

	virtual void Init();
	void Delete();
	virtual void Reset();
	virtual void InitOperator();

	unsigned int numVectors;

	// engine/post-proc needs access
public:
	ArrayLib::ArrayNIJK<f4vector>* f4_vv_ptr; //calc new voltage from old voltage
	ArrayLib::ArrayNIJK<f4vector>* f4_vi_ptr; //calc new voltage from old current
	ArrayLib::ArrayNIJK<f4vector>* f4_iv_ptr; //calc new current from old current
	ArrayLib::ArrayNIJK<f4vector>* f4_ii_ptr; //calc new current from old voltage
};

#endif // OPERATOR_SSE_H
