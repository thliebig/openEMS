/*
*	Copyright (C) 2010,2019 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#ifndef ARRAY_OPS_H
#define ARRAY_OPS_H

#ifdef __SIZEOF_FLOAT__
#if __SIZEOF_FLOAT__ != 4
#error wrong size of float
#endif
#endif

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <math.h>
#include "constants.h"

#define F4VECTOR_SIZE 16 // sizeof(typeid(f4vector))

#ifdef __GNUC__ // GCC
typedef float v4sf __attribute__ ((vector_size (F4VECTOR_SIZE))); // vector of four single floats
union f4vector
{
	v4sf v;
	float f[4];
};
#else // MSVC
#include <emmintrin.h>
union f4vector
{
	__m128 v;
	float f[4];
};
inline __m128 operator + (__m128 a, __m128 b) {return _mm_add_ps(a, b);}
inline __m128 operator - (__m128 a, __m128 b) {return _mm_sub_ps(a, b);}
inline __m128 operator * (__m128 a, __m128 b) {return _mm_mul_ps(a, b);}
inline __m128 operator / (__m128 a, __m128 b) {return _mm_div_ps(a, b);}

inline __m128 & operator += (__m128 & a, __m128 b){a = a + b; return a;}
inline __m128 & operator -= (__m128 & a, __m128 b){a = a - b; return a;}
inline __m128 & operator *= (__m128 & a, __m128 b){a = a * b; return a;}
inline __m128 & operator /= (__m128 & a, __m128 b){a = a / b; return a;}
#endif

void Delete1DArray_v4sf(f4vector* array);
f4vector* Create1DArray_v4sf(const unsigned int numLines);

// *************************************************************************************
// templates
// *************************************************************************************

#endif // ARRAY_OPS_H
