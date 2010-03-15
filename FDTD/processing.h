/*
*	Copyright (C) 2010 Thorsten Liebig (Thorsten.Liebig@gmx.de)
*
*	This program is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
*	(at your option) any later version.
*
*	This program is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*
*	You should have received a copy of the GNU General Public License
*	along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

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

	virtual void DefineStartStopCoord(double* dstart, double* dstop);

	void SetProcessInterval(unsigned int interval) {ProcessInterval=interval;}
	virtual int Process() {return GetNextInterval();}

	//! If Disabled Process() will do nothing...
	virtual void SetEnable(bool val) {Enabled=val;}
	//! If Disabled Process() will do nothing...
	virtual bool GetEnable() {return Enabled;}

protected:
	Operator* Op;
	Engine* Eng;

	bool Enabled;

	bool CheckTimestep();
	int GetNextInterval();
	unsigned int ProcessInterval;

	unsigned int start[3];
	unsigned int stop[3];

	double CalcLineIntegral(unsigned int* start, unsigned int* stop, int field);
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
