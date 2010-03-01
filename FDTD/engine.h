#ifndef ENGINE_H
#define ENGINE_H

#include "operator.h"

class Engine
{
public:
	Engine(Operator* op);
	virtual ~Engine();

	virtual void Init();
	virtual void Reset();

	//!Iterate a number of timesteps
	virtual bool IterateTS(unsigned int iterTS);

protected:
	Operator* Op;

	FDTD_FLOAT**** volt;
	FDTD_FLOAT**** curr;
	unsigned int numTS;
};

#endif // ENGINE_H
