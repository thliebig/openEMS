#include "processfields.h"

ProcessFields::ProcessFields(Operator* op, Engine* eng) : Processing(op, eng)
{
	DumpMode=0;
	DumpType = 0;

	for (int n=0;n<3;++n)
	{
		discDLines[n]=NULL;
	}
}

ProcessFields::~ProcessFields()
{
	for (int n=0;n<3;++n)
	{
		delete[] discDLines[n];
		discDLines[n]=NULL;
	}
}

void ProcessFields::DefineStartStopCoord(double* dstart, double* dstop)
{
	if (Op->SnapToMesh(dstart,start,true)==false) cerr << "ProcessFields::DefineStartStopCoord: Snapping error in mesh, check start value!!" << endl;
	if (Op->SnapToMesh(dstop,stop,true)==false) cerr << "ProcessFields::DefineStartStopCoord: Snapping error in mesh, check stop value!!" << endl;

	for (int n=0;n<3;++n)
	{
//		cerr << "start " << start[n] << "stop " << stop[n];
		if (start[n]>stop[n])
		{
			unsigned int help = start[n];
			start[n]=stop[n];
			stop[n]=help;
		}
		++stop[n];
		numDLines[n]=stop[n]-start[n];
//		cerr << " number of lines " << numDLines[n] << endl;
		delete[] discDLines[n];
		discDLines[n] = new double[numDLines[n]];
		for (unsigned int i=0;i<numDLines[n];++i)
		{
			discDLines[n][i] = 0.5*(Op->discLines[n][start[n]+i+1] +  Op->discLines[n][start[n]+i]);
//			cerr << n << " : " << discDLines[n][i] << endl;
		}
	}
}

bool ProcessFields::DumpFieldArray2VTK(ofstream &file, string name, FDTD_FLOAT**** array, double** discLines, unsigned int* numLines)
{
	file << "# vtk DataFile Version 2.0" << endl;
	file << "Rectilinear Grid openEMS_ProcessFields" << endl;
	file << "ASCII" << endl;
	file << "DATASET RECTILINEAR_GRID " << endl;
	file << "DIMENSIONS " << numLines[0] << " " << numLines[1] << " " << numLines[2] << endl;
	file << "X_COORDINATES " << numLines[0] << " float" << endl;
	for (unsigned int i=0;i<numLines[0];++i)
		file << discLines[0][i] << " ";
	file << endl;
	file << "Y_COORDINATES " << numLines[1] << " float" << endl;
	for (unsigned int i=0;i<numLines[1];++i)
		file << discLines[1][i] << " ";
	file << endl;
	file << "Z_COORDINATES " << numLines[2] << " float" << endl;
	for (unsigned int i=0;i<numLines[2];++i)
		file << discLines[2][i] << " ";
	file << endl << endl;

	file << "POINT_DATA " << numLines[0]*numLines[1]*numLines[2] << endl;
	file << "VECTORS " << name << " float " << endl;

	unsigned int pos[3];
	for (pos[2]=0;pos[2]<numLines[2];++pos[2])
	{
		for (pos[1]=0;pos[1]<numLines[1];++pos[1])
		{
			for (pos[0]=0;pos[0]<numLines[0];++pos[0])
			{
				//in x
				file << array[0][pos[0]][pos[1]][pos[2]] << " ";
				//in y
				file << array[1][pos[0]][pos[1]][pos[2]] << " ";
				//in z
				file << array[2][pos[0]][pos[1]][pos[2]] << endl;
			}
		}
	}
}
