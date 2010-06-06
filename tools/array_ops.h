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

using namespace std;

typedef float v4sf __attribute__ ((vector_size (16))); // vector of four single floats

union f4vector
{
  v4sf v;
  float f[4];
};

void Delete1DArray_v4sf(f4vector* array);
void Delete3DArray_v4sf(f4vector*** array, const unsigned int* numLines);
void Delete_N_3DArray_v4sf(f4vector**** array, const unsigned int* numLines);
f4vector* Create1DArray_v4sf(const unsigned int numLines);
f4vector*** Create3DArray_v4sf(const unsigned int* numLines);
f4vector**** Create_N_3DArray_v4sf(const unsigned int* numLines);

// *************************************************************************************
// templates
// *************************************************************************************
template <typename T>
T** Create2DArray(const unsigned int* numLines)
{
	T** array=NULL;
	unsigned int pos[3];
	array = new T*[numLines[0]];
	for (pos[0]=0;pos[0]<numLines[0];++pos[0])
	{
		array[pos[0]] = new T[numLines[1]];
		for (pos[1]=0;pos[1]<numLines[1];++pos[1])
		{
			array[pos[0]][pos[1]] = 0;
		}
	}
	return array;
}

template <typename T>
void Delete2DArray(T** array, const unsigned int* numLines)
{
	if (array==NULL) return;
	unsigned int pos[3];
	for (pos[0]=0;pos[0]<numLines[0];++pos[0])
	{
		delete[] array[pos[0]];
	}
	delete[] array;
}

template <typename T>
inline T& Access_N_3DArray(T**** array, unsigned int n, unsigned int* pos)
{
	return array[n][pos[0]][pos[1]][pos[2]];
}

template <typename T>
inline T& Access_N_3DArray(T**** array, unsigned int n, unsigned int x, unsigned int y, unsigned int z )
{
	return array[n][x][y][z];
}

template <typename T>
T*** Create3DArray(const unsigned int* numLines)
{
	T*** array=NULL;
	unsigned int pos[3];
	array = new T**[numLines[0]];
	for (pos[0]=0;pos[0]<numLines[0];++pos[0])
	{
		array[pos[0]] = new T*[numLines[1]];
		for (pos[1]=0;pos[1]<numLines[1];++pos[1])
		{
			array[pos[0]][pos[1]] = new T[numLines[2]];
			for (pos[2]=0;pos[2]<numLines[2];++pos[2])
			{
				array[pos[0]][pos[1]][pos[2]] = 0;
			}
		}
	}
	return array;
}

template <typename T>
T**** Create_N_3DArray(const unsigned int* numLines)
{
	T**** array=NULL;
	array = new T***[3];
	for (int n=0;n<3;++n)
	{
		array[n]=Create3DArray<T>( numLines );
	}
	return array;
}

template <typename T>
void Delete3DArray(T*** array, const unsigned int* numLines)
{
	if (!array) return;
	unsigned int pos[3];
	for (pos[0]=0;pos[0]<numLines[0];++pos[0])
	{
		for (pos[1]=0;pos[1]<numLines[1];++pos[1])
		{
			delete[] array[pos[0]][pos[1]];
		}
		delete[] array[pos[0]];
	}
	delete[] array;
}

template <typename T>
void Delete_N_3DArray(T**** array, const unsigned int* numLines)
{
	if (!array) return;
	for (int n=0;n<3;++n)
	{
		Delete3DArray<T>(array[n],numLines);
	}
	delete[] array;
}

template <typename T>
void Dump_N_3DArray2File(ostream &file, T**** array, const unsigned int* numLines)
{
	unsigned int pos[3];
	for (pos[0]=0;pos[0]<numLines[0];++pos[0])
	{
		for (pos[1]=0;pos[1]<numLines[1];++pos[1])
		{
			for (pos[2]=0;pos[2]<numLines[2];++pos[2])
			{
				file << pos[0] << "\t" << pos[1] << "\t" << pos[2];
				for (int n=0;n<3;++n)
					file << "\t" << (float)array[n][pos[0]][pos[1]][pos[2]];
				file << endl;
			}
		}
	}
}

#endif // ARRAY_OPS_H
