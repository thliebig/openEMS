#ifndef PROCESSING_H
#define PROCESSING_H

#include <iostream>
#include <fstream>
#include "operator.h"
#include "engine.h"

class Processing
{
public:
	virtual ~Processing();

	virtual void DefineStartStopCoord(double* dstart, double* dstop);

	virtual void Process() {};

	//! If Disabled Process() will do nothing...
	virtual void SetEnable(bool val) {Enabled=val;}
	//! If Disabled Process() will do nothing...
	virtual bool GetEnable() {return Enabled;}

protected:
	Processing(Operator* op, Engine* eng);
	Operator* Op;
	Engine* Eng;

	bool Enabled;

	unsigned int start[3];
	unsigned int stop[3];
};

#endif // PROCESSING_H
