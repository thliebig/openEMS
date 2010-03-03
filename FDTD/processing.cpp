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
	if (Op->SnapToMesh(dstart,start)==false) cerr << "Processing::DefineStartStopCoord: Warning: Snapping problem, check start value!!" << endl;
	if (Op->SnapToMesh(dstop,stop)==false) cerr << "Processing::DefineStartStopCoord: Warning: Snapping problem, check stop value!!" << endl;
}

