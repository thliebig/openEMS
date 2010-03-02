#ifndef PROCESSFIELDS_H
#define PROCESSFIELDS_H

#include "processing.h"
#include "../tools/array_ops.h"

class ProcessFields : public Processing
{
public:
	virtual ~ProcessFields();

	virtual void DefineStartStopCoord(double* dstart, double* dstop);

	//! Used file pattern e.g. pattern="tmp/efield_" --> "tmp/efield_000045.vtk" for timestep 45 or "tmp/efield_2.40000e9.vtk" for 2.4GHz E-field dump.
	void SetFilePattern(string fp) {filePattern=fp;}

	//! This methode will dump all fields in the center of a main cell (dual-node) using 4 E-field and 2 H-fields per direction. (default)
	void SetDumpMode2Cell() {DumpMode=0;}

	//! Set dump type: 0 for E-fields, 1 for H-fields, 2 for D-fields, 3 for B-fields, 4 for J-fields, etc...
	void SetDumpType(int type) {DumpType=type;}

//	virtual void Process();
protected:
	ProcessFields(Operator* op, Engine* eng);

	bool DumpFieldArray2VTK(ofstream &file, string name, FDTD_FLOAT**** array, double** discLines, unsigned int* numLines);

	int DumpMode;
	int DumpType;
	string filePattern;

	unsigned int numDLines[3];
	double* discDLines[3];
};

#endif // PROCESSFIELDS_H
