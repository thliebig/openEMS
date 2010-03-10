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

	void SetProcessInterval(unsigned int interval) {ProcessInterval=interval;}
	virtual int Process() {return GetNextInterval();}

	//! If Disabled Process() will do nothing...
	virtual void SetEnable(bool val) {Enabled=val;}
	//! If Disabled Process() will do nothing...
	virtual bool GetEnable() {return Enabled;}

protected:
	Processing(Operator* op, Engine* eng);
	Operator* Op;
	Engine* Eng;

	bool Enabled;

	bool CheckTimestep();
	int GetNextInterval();
	unsigned int ProcessInterval;

	unsigned int start[3];
	unsigned int stop[3];
};

class ProcessingArray
{
public:
	ProcessingArray() {};
	~ProcessingArray() {};

	void AddProcessing(Processing* proc);

	//! Deletes all given processing's, can be helpful, but use carefull!!!
	void DeleteAll();

	int Process();

protected:
	vector<Processing*> ProcessArray;
};

#endif // PROCESSING_H
