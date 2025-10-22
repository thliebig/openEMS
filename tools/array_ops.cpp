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

using namespace std;

#ifdef WIN32
#include <malloc.h>
#define MEMALIGN( array, alignment, size ) !(*array = _mm_malloc( size, alignment ))
#define FREE( array ) _mm_free( array )
#else
#define MEMALIGN( array, alignment, size ) posix_memalign( array, alignment, size )
#define FREE( array ) free( array )
#endif

void Delete1DArray_v4sf(f4vector* array)
{
	if (array==NULL) return;
	FREE( array );
}

f4vector* Create1DArray_v4sf(const unsigned int numLines)
{
	f4vector* array=NULL;
	if (MEMALIGN( (void**)&array, 16, F4VECTOR_SIZE*numLines ))
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
