#include "processing.h"

Processing::Processing(Operator* op, Engine* eng)
{
	Op=op;
	Eng=eng;
	Enabled = true;
}

Processing::~Processing()
{
}

bool Processing::CheckTimestep()
{
	if (Eng->GetNumberOfTimesteps()%ProcessInterval==0) return true;
	return false;
}

int Processing::GetNextInterval()
{
	if (Enabled==false) return -1;
	return ProcessInterval - Eng->GetNumberOfTimesteps()%ProcessInterval;
}

void Processing::DefineStartStopCoord(double* dstart, double* dstop)
{
	if (Op->SnapToMesh(dstart,start)==false) cerr << "Processing::DefineStartStopCoord: Warning: Snapping problem, check start value!!" << endl;
	if (Op->SnapToMesh(dstop,stop)==false) cerr << "Processing::DefineStartStopCoord: Warning: Snapping problem, check stop value!!" << endl;
}


void ProcessingArray::AddProcessing(Processing* proc)
{
	ProcessArray.push_back(proc);
}

void ProcessingArray::DeleteAll()
{
	for (size_t i=0;i<ProcessArray.size();++i)
	{
		delete ProcessArray.at(i);
	}
	ProcessArray.clear();
}

int ProcessingArray::Process()
{
	int nextProcess=1e100;
	//this could be done nicely in parallel??
	for (size_t i=0;i<ProcessArray.size();++i)
	{
		int step = ProcessArray.at(i)->Process();
		if ((step>0) && (step<nextProcess))
			nextProcess=step;
	}
	return nextProcess;
}
