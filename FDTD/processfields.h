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

#ifndef PROCESSFIELDS_H
#define PROCESSFIELDS_H

#include "processing.h"
#include "../tools/array_ops.h"

class ProcessFields : public Processing
{
public:
	virtual ~ProcessFields();

	virtual void DefineStartStopCoord(double* dstart, double* dstop);

//	virtual void SetSubSampling(unsigned int subSampleRate, int dir=-1);

	//! Used file pattern e.g. pattern="tmp/efield_" --> "tmp/efield_000045.vtk" for timestep 45 or "tmp/efield_2.40000e9.vtk" for 2.4GHz E-field dump.
	void SetFilePattern(string fp) {filePattern=fp;}

	//! Define the Dump-Mode
	void SetDumpMode(int mode) {DumpMode=mode;}
	//! This methode will dump all fields in the center of a main cell (dual-node) using 4 E-field and 2 H-fields per direction. (default)
	void SetDumpMode2Cell() {DumpMode=2;}

	//! Set dump type: 0 for E-fields, 1 for H-fields, 2 for D-fields, 3 for B-fields, 4 for J-fields, etc...
	void SetDumpType(int type) {DumpType=type;}

	static bool DumpVectorArray2VTK(ofstream &file, string name, FDTD_FLOAT**** array, double** discLines, unsigned int* numLines);
	static bool DumpMultiVectorArray2VTK(ofstream &file, string names[], FDTD_FLOAT**** array[], unsigned int numFields, double** discLines, unsigned int* numLines);
	static bool DumpScalarArray2VTK(ofstream &file, string name, FDTD_FLOAT*** array, double** discLines, unsigned int* numLines);
	static bool DumpMultiScalarArray2VTK(ofstream &file, string names[], FDTD_FLOAT*** array[], unsigned int numFields, double** discLines, unsigned int* numLines);

//	virtual void Process();
protected:
	ProcessFields(Operator* op, Engine* eng);

	static void WriteVTKHeader(ofstream &file, double** discLines, unsigned int* numLines);
	static void WriteVTKVectorArray(ofstream &file, string name, FDTD_FLOAT**** array, unsigned int* numLines);
	static void WriteVTKScalarArray(ofstream &file, string name, FDTD_FLOAT*** array, unsigned int* numLines);

	int DumpMode;
	int DumpType;
	string filePattern;

//	unsigned int subSample[3];

	//! dump mesh
	unsigned int numLines[3];
	double* discLines[3];
	//! dual dump mesh
	unsigned int numDLines[3];
	double* discDLines[3];
};

#endif // PROCESSFIELDS_H
