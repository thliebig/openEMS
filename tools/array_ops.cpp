#include "array_ops.h"
#include <ostream>

FDTD_FLOAT*** Create3DArray(unsigned int* numLines)
{
	FDTD_FLOAT*** array;
	unsigned int pos[3];
	array = new FDTD_FLOAT**[numLines[0]];
	for (pos[0]=0;pos[0]<numLines[0];++pos[0])
	{
		array[pos[0]] = new FDTD_FLOAT*[numLines[1]];
		for (pos[1]=0;pos[1]<numLines[1];++pos[1])
		{
			array[pos[0]][pos[1]] = new FDTD_FLOAT[numLines[2]];
			for (pos[2]=0;pos[2]<numLines[2];++pos[2])
			{
				array[pos[0]][pos[1]][pos[2]] = 0;
			}
		}
	}
	return array;
}

void Delete3DArray(FDTD_FLOAT*** array, unsigned int* numLines)
{
	if (array==NULL) return;
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

FDTD_FLOAT**** Create_N_3DArray(unsigned int* numLines)
{
	FDTD_FLOAT**** array;
	array = new FDTD_FLOAT***[3];
	for (int n=0;n<3;++n)
	{
		array[n]=Create3DArray(numLines);
	}
	return array;
}

void Delete_N_3DArray(FDTD_FLOAT**** array, unsigned int* numLines)
{
	if (array==NULL) return;
	unsigned int pos[3];
	for (int n=0;n<3;++n)
	{
		Delete3DArray(array[n],numLines);
	}
	delete[] array;
}

void Dump_N_3DArray2File(ostream &file, FDTD_FLOAT**** array, unsigned int* numLines)
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
					file << "\t" << array[n][pos[0]][pos[1]][pos[2]];
				file << endl;
			}
		}
	}
}
