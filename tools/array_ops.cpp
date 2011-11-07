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

#include "array_ops.h"
#include <ostream>

#ifdef WIN32
#define __MSVCRT_VERSION__ 0x0700
#include <malloc.h>
//(void**)&array, 16, sizeof(typeof(f4vector**))*numLines[0]
#define MEMALIGN( array, alignment, size ) !(*array = _aligned_malloc( size, alignment ))
#define FREE( array ) _aligned_free( array )
#else
#define MEMALIGN( array, alignment, size ) posix_memalign( array, alignment, size )
#define FREE( array ) free( array )
#endif

void Delete1DArray_v4sf(f4vector* array)
{
	if (array==NULL) return;
	FREE( array );
}


void Delete3DArray_v4sf(f4vector*** array, const unsigned int* numLines)
{
	if (array==NULL) return;
	unsigned int pos[3];
	for (pos[0]=0; pos[0]<numLines[0]; ++pos[0])
	{
		for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
		{
			FREE( array[pos[0]][pos[1]] );
			//delete[] array[pos[0]][pos[1]];
		}
		FREE( array[pos[0]] );
		//delete[] array[pos[0]];
	}
	FREE( array );
	//delete[] array;
}

void Delete_N_3DArray_v4sf(f4vector**** array, const unsigned int* numLines)
{
	if (array==NULL) return;
	for (int n=0; n<3; ++n)
	{
		Delete3DArray_v4sf(array[n],numLines);
	}
	FREE( array );
	//delete[] array;
}

f4vector* Create1DArray_v4sf(const unsigned int numLines)
{
	f4vector* array=NULL;
	if (MEMALIGN( (void**)&array, 16, sizeof(typeof(f4vector))*numLines ))
	{
		cerr << "cannot allocate aligned memory" << endl;
		exit(3);
	}
	for (unsigned int pos=0; pos<numLines; ++pos)
	{
		array[pos].f[0] = 0;
		array[pos].f[1] = 0;
		array[pos].f[2] = 0;
		array[pos].f[3] = 0;
	}
	return array;
}

//! \brief this function allocates a 3D array, which is aligned to 16 byte
f4vector*** Create3DArray_v4sf(const unsigned int* numLines)
{
	unsigned int numZ = ceil((double)numLines[2]/4.0);

	f4vector*** array=NULL;
	unsigned int pos[3];
	if (MEMALIGN( (void**)&array, 16, sizeof(typeof(f4vector**))*numLines[0] ))
	{
		cerr << "cannot allocate aligned memory" << endl;
		exit(3);
	}
	//array = new f4vector**[numLines[0]];
	for (pos[0]=0; pos[0]<numLines[0]; ++pos[0])
	{
		if (MEMALIGN( (void**)&array[pos[0]], 16, sizeof(typeof(f4vector*))*numLines[1] ))
		{
			cerr << "cannot allocate aligned memory" << endl;
			exit(3);
		}
		//array[pos[0]] = new f4vector*[numLines[1]];
		for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
		{
			if (MEMALIGN( (void**)&array[pos[0]][pos[1]], 16, sizeof(typeof(f4vector))*numZ ))
			{
				cerr << "cannot allocate aligned memory" << endl;
				exit(3);
			}
			//array[pos[0]][pos[1]] = new f4vector[numZ];
			for (pos[2]=0; pos[2]<numZ; ++pos[2])
			{
				array[pos[0]][pos[1]][pos[2]].f[0] = 0;
				array[pos[0]][pos[1]][pos[2]].f[1] = 0;
				array[pos[0]][pos[1]][pos[2]].f[2] = 0;
				array[pos[0]][pos[1]][pos[2]].f[3] = 0;
			}
		}
	}
	return array;
}

f4vector**** Create_N_3DArray_v4sf(const unsigned int* numLines)
{
	f4vector**** array=NULL;
	if (MEMALIGN( (void**)&array, 16, sizeof(typeof(f4vector***))*3 ))
	{
		cerr << "cannot allocate aligned memory" << endl;
		exit(3);
	}
	//array = new f4vector***[3];
	for (int n=0; n<3; ++n)
	{
		array[n]=Create3DArray_v4sf(numLines);
	}
	return array;
}

