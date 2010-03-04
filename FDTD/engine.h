#ifndef ENGINE_H
#define ENGINE_H

#include "operator.h"

class Engine
{
	friend class ProcessVoltage;
	friend class ProcessCurrent;
	friend class ProcessFieldsTD;
public:
	Engine(Operator* op);
	virtual ~Engine();

	virtual void Init();
	virtual void Reset();

	//!Iterate a number of timesteps
	virtual bool IterateTS(unsigned int iterTS);

	unsigned int GetNumberOfTimesteps() {return numTS;};

protected:
	Operator* Op;

	FDTD_FLOAT**** volt;
	FDTD_FLOAT**** curr;
	unsigned int numTS;
};

#endif // ENGINE_H
