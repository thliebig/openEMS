#ifndef PROCESSCURRENT_H
#define PROCESSCURRENT_H

#include "processing.h"

class ProcessCurrent : public Processing
{
public:
	ProcessCurrent(Operator* op, Engine* eng);
	virtual ~ProcessCurrent();

	virtual void OpenFile(string outfile);

	virtual void DefineStartStopCoord(double* dstart, double* dstop);

	virtual int Process();

protected:
	ofstream file;

	vector<FDTD_FLOAT> v_current;
};

#endif // PROCESSCURRENT_H
