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

#include "processfields.h"

ProcessFields::ProcessFields(Operator* op, Engine* eng) : Processing(op, eng)
{
	DumpMode=0;
	DumpType = 0;
//	SetSubSampling(1);

	for (int n=0;n<3;++n)
	{
		numLines[n]=0;
		discDLines[n]=NULL;
		discLines[n]=NULL;
	}
}

ProcessFields::~ProcessFields()
{
	for (int n=0;n<3;++n)
	{
		delete[] discDLines[n];
		discDLines[n]=NULL;
		delete[] discLines[n];
		discLines[n]=NULL;
	}
}

void ProcessFields::DefineStartStopCoord(double* dstart, double* dstop)
{
	if (DumpMode==0)
	{
		if (Op->SnapToMesh(dstart,start)==false) cerr << "ProcessFields::DefineStartStopCoord: Warning: Snapping problem, check start value!!" << endl;
		if (Op->SnapToMesh(dstop,stop)==false) cerr << "ProcessFields::DefineStartStopCoord: Warning: Snapping problem, check stop value!!" << endl;
		//create dual mesh
		for (int n=0;n<3;++n)
		{
	//		cerr << "start " << start[n] << "stop " << stop[n];
			if (start[n]>stop[n])
			{
				unsigned int help = start[n];
				start[n]=stop[n];
				stop[n]=help;
			}
			numLines[n]=stop[n]-start[n]+1;
	//		cerr << " number of lines " << numDLines[n] << endl;
			delete[] discLines[n];
			discLines[n] = new double[numLines[n]];
			for (unsigned int i=0;i<numLines[n];++i)
			{
				discLines[n][i] = Op->discLines[n][start[n]+i];
	//			cerr << n << " : " << discDLines[n][i] << endl;
			}
		}
	}
	else if (DumpMode==2)
	{
		if (Op->SnapToMesh(dstart,start,true)==false) cerr << "ProcessFields::DefineStartStopCoord: Warning: Snapping problem, check start value!!" << endl;
		if (Op->SnapToMesh(dstop,stop,true)==false) cerr << "ProcessFields::DefineStartStopCoord: Warning: Snapping problem, check stop value!!" << endl;

		//create dual mesh
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
}


double ProcessFields::CalcTotalEnergy()
{
	FDTD_FLOAT**** volt = Eng->GetVoltages();
	FDTD_FLOAT**** curr = Eng->GetCurrents();
	double energy=0;
	if (Eng==NULL) return 0.0;
	unsigned int pos[3];
	for (pos[0]=0;pos[0]<Op->numLines[0];++pos[0])
	{
		for (pos[1]=0;pos[1]<Op->numLines[1];++pos[1])
		{
			for (pos[2]=0;pos[2]<Op->numLines[2];++pos[2])
			{
				energy+=fabs(volt[0][pos[0]][pos[1]][pos[2]] * curr[1][pos[0]][pos[1]][pos[2]]);
				energy+=fabs(volt[0][pos[0]][pos[1]][pos[2]] * curr[2][pos[0]][pos[1]][pos[2]]);
				energy+=fabs(volt[1][pos[0]][pos[1]][pos[2]] * curr[0][pos[0]][pos[1]][pos[2]]);
				energy+=fabs(volt[1][pos[0]][pos[1]][pos[2]] * curr[2][pos[0]][pos[1]][pos[2]]);
				energy+=fabs(volt[2][pos[0]][pos[1]][pos[2]] * curr[0][pos[0]][pos[1]][pos[2]]);
				energy+=fabs(volt[2][pos[0]][pos[1]][pos[2]] * curr[1][pos[0]][pos[1]][pos[2]]);
			}
		}
	}
	return energy*0.5;
}

//void  ProcessFields::SetSubSampling(unsigned int subSampleRate, int dir)
//{
//	if (dir>2) return;
//	if (dir<0)
//	{
//		subSample[0]=subSampleRate;
//		subSample[1]=subSampleRate;
//		subSample[2]=subSampleRate;
//	}
//	else subSample[dir]=subSampleRate;
//}

void ProcessFields::WriteVTKHeader(ofstream &file, double** discLines, unsigned int* numLines)
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
}

void ProcessFields::WriteVTKVectorArray(ofstream &file, string name, FDTD_FLOAT**** array, unsigned int* numLines)
{
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


bool ProcessFields::DumpVectorArray2VTK(ofstream &file, string name, FDTD_FLOAT**** array, double** discLines, unsigned int* numLines)
{
	WriteVTKHeader(file, discLines, numLines);
	WriteVTKVectorArray(file, name, array, numLines);
	return true;
}

bool ProcessFields::DumpMultiVectorArray2VTK(ofstream &file, string names[], FDTD_FLOAT**** array[], unsigned int numFields, double** discLines, unsigned int* numLines)
{
	WriteVTKHeader(file, discLines, numLines);
	for (unsigned int n=0;n<numFields;++n)
	{
		WriteVTKVectorArray(file, names[n], array[n], numLines);
		file << endl;
	}
	return true;
}

void ProcessFields::WriteVTKScalarArray(ofstream &file, string name, FDTD_FLOAT*** array, unsigned int* numLines)
{
	file << "SCALARS " << name << " float " << 1 << endl;
	file << "LOOKUP_TABLE default" << endl;
	unsigned int pos[3];
	int count=0;
	for (pos[2]=0;pos[2]<numLines[2];++pos[2])
	{
		for (pos[1]=0;pos[1]<numLines[1];++pos[1])
		{
			for (pos[0]=0;pos[0]<numLines[0];++pos[0])
			{
				file << array[pos[0]][pos[1]][pos[2]] << " ";
				++count;
				if (count%10==0)
					file << endl;
			}
		}
	}
}

bool ProcessFields::DumpScalarArray2VTK(ofstream &file, string name, FDTD_FLOAT*** array, double** discLines, unsigned int* numLines)
{
	WriteVTKHeader(file, discLines, numLines);
	WriteVTKScalarArray(file, name, array, numLines);
	return true;
}

bool ProcessFields::DumpMultiScalarArray2VTK(ofstream &file, string names[], FDTD_FLOAT*** array[], unsigned int numFields, double** discLines, unsigned int* numLines)
{
	WriteVTKHeader(file, discLines, numLines);
	for (unsigned int n=0;n<numFields;++n)
	{
		WriteVTKScalarArray(file, names[n], array[n], numLines);
		file << endl;
	}
	return true;
}


