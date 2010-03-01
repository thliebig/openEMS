#ifndef PROCESSVOLTAGE_H
#define PROCESSVOLTAGE_H

#include "processing.h"

class ProcessVoltage : public Processing
{
public:
	ProcessVoltage(Operator* op, Engine* eng);
	virtual ~ProcessVoltage();

	void DefineStartStopCoord(double* dstart, double* dstop);

	virtual void Process();

protected:
	unsigned int start[3];
	unsigned int stop[3];

	vector<FDTD_FLOAT> voltages;
};

#endif // PROCESSVOLTAGE_H
