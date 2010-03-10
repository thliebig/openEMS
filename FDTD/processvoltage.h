#ifndef PROCESSVOLTAGE_H
#define PROCESSVOLTAGE_H

#include "processing.h"

//! Process voltage along a line from start to stop coordinates. ATM integration along the axis e.g.: in x, then y then z direction (Future: diagonal integration)
class ProcessVoltage : public Processing
{
public:
	ProcessVoltage(Operator* op, Engine* eng);
	virtual ~ProcessVoltage();

	virtual void OpenFile(string outfile);

	virtual int Process();

protected:
	ofstream file;

	vector<FDTD_FLOAT> voltages;
};

#endif // PROCESSVOLTAGE_H
