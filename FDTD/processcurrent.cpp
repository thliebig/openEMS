#include "processcurrent.h"

ProcessCurrent::ProcessCurrent(Operator* op, Engine* eng) : Processing(op, eng)
{
}

ProcessCurrent::~ProcessCurrent()
{
}

void ProcessCurrent::OpenFile(string outfile)
{
	if (file.is_open()) file.close();

	file.open(outfile.c_str());
	if (file.is_open()==false)
	{
		cerr << "Can't open file: " << outfile << endl;
		return;
	}
}

void ProcessCurrent::DefineStartStopCoord(double* dstart, double* dstop)
{
	if (Op->SnapToMesh(dstart,start,true)==false) cerr << "ProcessCurrent::DefineStartStopCoord: Warning: Snapping problem, check start value!!" << endl;
	if (Op->SnapToMesh(dstop,stop,true)==false) cerr << "ProcessCurrent::DefineStartStopCoord: Warning: Snapping problem, check stop value!!" << endl;
}

void ProcessCurrent::Process()
{
	if (Enabled==false) return;
	FDTD_FLOAT current=0;
//	FDTD_FLOAT help=0;
	double sign[3]={1,1,1};
	unsigned int pos[3]={start[0],start[1],start[2]};
	double loc_start[]={start[0],start[1],start[2]};
	double loc_stop[]={stop[0],stop[1],stop[2]};

	for (int n=0;n<3;++n)
	{
		if (start[n]>stop[n])
		{
			unsigned int help=start[n];
			start[n]=stop[n];
			stop[n]=help;
		}
	}

	//x-current
	for (int i=start[0];i<stop[0];++i)
		current+=Eng->curr[0][i][start[1]][start[2]];
	//y-current
	for (int i=start[1];i<stop[1];++i)
		current+=Eng->curr[1][stop[0]][i][start[2]];
	//z-current
	for (int i=start[2];i<stop[2];++i)
		current+=Eng->curr[2][stop[0]][stop[1]][i];
	//x-current
	for (int i=start[0];i<stop[0];++i)
		current-=Eng->curr[0][i][stop[1]][stop[2]];
	//y-current
	for (int i=start[1];i<stop[1];++i)
		current-=Eng->curr[1][start[0]][i][stop[2]];
	//z-current
	for (int i=start[2];i<stop[2];++i)
		current-=Eng->curr[2][start[0]][start[1]][i];

//	cerr << "ts: " << Eng->numTS << " i: " << current << endl;
	v_current.push_back(current);
	file << (double)Eng->numTS*Op->GetTimestep() << "\t" << current << endl;
}
