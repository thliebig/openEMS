#include "processing.h"

Processing::Processing(Operator* op, Engine* eng)
{
	Op=op;
	Eng=eng;
}

Processing::~Processing()
{
}


void Processing::DefineStartStopCoord(double* dstart, double* dstop)
{
	if (Op->SnapToMesh(dstart,start)==false) cerr << "Processing::DefineStartStopCoord: Snapping error in mesh, check start value!!" << endl;
	if (Op->SnapToMesh(dstop,stop)==false) cerr << "Processing::DefineStartStopCoord: Snapping error in mesh, check start value!!" << endl;
}

