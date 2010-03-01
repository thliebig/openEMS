#ifndef PROCESSING_H
#define PROCESSING_H

#include <iostream>
#include <fstream>
#include "operator.h"
#include "engine.h"

class Processing
{
public:
	Processing(Operator* op, Engine* eng);
	virtual ~Processing();

	virtual void OpenFile(string outfile);

	virtual void Process() {};

protected:
	Operator* Op;
	Engine* Eng;
	ofstream file;
};

#endif // PROCESSING_H
