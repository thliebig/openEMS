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

double Processing::CalcLineIntegral(unsigned int* start, unsigned int* stop, int field)
{
	double result=0;
	FDTD_FLOAT**** array;
	if (field==0)
		array=Eng->volt;
	else if (field==1)
		array=Eng->curr;
	else return 0.0;

	unsigned int pos[3]={start[0],start[1],start[2]};
//	cerr << Eng->volt[1][pos[0]][pos[1]][pos[2]] << endl;
	for (int n=0;n<3;++n)
	{
		if (start[n]<stop[n])
		{
			for (;pos[n]<stop[n];++pos[n])
			{
				result+=array[n][pos[0]][pos[1]][pos[2]];
//				cerr << n << " " << pos[0] << " " << pos[1] << " " <<  pos[2] << " " << Eng->volt[n][pos[0]][pos[1]][pos[2]] << endl;
			}
		}
		else
		{
			for (;pos[n]>stop[n];--pos[n])
				result-=array[n][pos[0]][pos[1]][pos[2]];
		}
	}
	return result;
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
