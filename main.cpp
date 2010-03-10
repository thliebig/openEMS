#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include "tools/array_ops.h"
#include "FDTD/operator.h"
#include "FDTD/engine.h"
#include "FDTD/processvoltage.h"
#include "FDTD/processcurrent.h"
#include "FDTD/processfields_td.h"

#include "ContinuousStructure.h"

#include "examples/FDTD_examples.h"

using namespace std;

int main(int argc, char *argv[])
{
	time_t startTime=time(NULL);

	//*************** setup/read geometry ************//
	cerr << "Create Geometry..." << endl;
	ContinuousStructure CSX;

	BuildPlaneWave(CSX);
//	BuildMSL(CSX);

	//*************** setup operator ************//
	cerr << "Create Operator..." << endl;
	Operator cop;
	cop.SetGeometryCSX(&CSX);

	cop.CalcECOperator();

	double fmax=1e9;
	cop.CalcGaussianPulsExcitation(fmax/2,fmax/2);

	time_t OpDoneTime=time(NULL);

	cop.ShowSize();
	bool bnd[] = {1,1,0,0,0,0};
	cop.ApplyMagneticBC(bnd);

	cerr << "Nyquist number of timesteps: " << cop.GetNyquistNum(fmax) << endl;
	unsigned int Nyquist = cop.GetNyquistNum(fmax);

	cerr << "Time for operator: " << difftime(OpDoneTime,startTime) << endl;

	//create FDTD engine
	Engine eng(&cop);

	time_t currTime = time(NULL);

	//*************** setup processing ************//
	ProcessingArray PA;

	ProcessVoltage* PV = new ProcessVoltage(&cop,&eng);
	PA.AddProcessing(PV);
	PV->SetProcessInterval(Nyquist/3); //three times as often as nyquist dictates
	PV->OpenFile("tmp/u1");
	double start[]={0,50,0};
	double stop[]={0,0,0};
	PV->DefineStartStopCoord(start,stop);

	ProcessCurrent* PCurr = new ProcessCurrent(&cop,&eng);
	PA.AddProcessing(PCurr);
	PCurr->SetProcessInterval(Nyquist/3); //three times as often as nyquist dictates
	PCurr->OpenFile("tmp/i1");
	start[0]=-50;start[1]=40;start[2]=-0;
	stop[0]=50;stop[1]=60;stop[2]=-0;
	PCurr->DefineStartStopCoord(start,stop);

	unsigned int maxIter = 1800;

	bool EnableDump = true;

	vector<CSProperties*> DumpProps = CSX.GetPropertyByType(CSProperties::DUMPBOX);
	for (size_t i=0;i<DumpProps.size();++i)
	{
		ProcessFieldsTD* ProcTD = new ProcessFieldsTD(&cop,&eng);
		ProcTD->SetEnable(EnableDump);
		ProcTD->SetProcessInterval(Nyquist/2); //twice as often as nyquist dictates

		//only looking for one prim atm
		CSPrimitives* prim = DumpProps.at(i)->GetPrimitive(0);
		if (prim==NULL)
			delete 	ProcTD;
		else
		{
			bool acc;
			double* bnd = prim->GetBoundBox(acc);
			start[0]= bnd[0];start[1]=bnd[2];start[2]=bnd[4];
			stop[0] = bnd[1];stop[1] =bnd[3];stop[2] =bnd[5];
			ProcTD->DefineStartStopCoord(start,stop);
			CSPropDumpBox* db = DumpProps.at(i)->ToDumpBox();
			if (db)
			{
				ProcTD->SetDumpType(db->GetDumpType());
				ProcTD->SetDumpMode(db->GetDumpMode());
				ProcTD->SetFilePattern(db->GetName());
				PA.AddProcessing(ProcTD);
			}
			else
				delete 	ProcTD;
		}
	}

	int step=PA.Process();
	if ((step<0) || (step>maxIter)) step=maxIter;
	//*************** simulate ************//
	while (eng.GetNumberOfTimesteps()<maxIter)
	{
		eng.IterateTS(step);
		step=PA.Process();
//		cerr << " do " << step << " steps; current: " << eng.GetNumberOfTimesteps() << endl;
		if ((step<0) || (step>maxIter - eng.GetNumberOfTimesteps())) step=maxIter - eng.GetNumberOfTimesteps();
	}

	//*************** postproc ************//
	time_t prevTime = currTime;
	currTime = time(NULL);

	double t_diff = difftime(currTime,prevTime);

	cerr << "Time for " << eng.GetNumberOfTimesteps() << " iterations with " << cop.GetNumberCells() << " cells : " << t_diff << " sec" << endl;
	cerr << "Speed: " << (double)cop.GetNumberCells()*(double)eng.GetNumberOfTimesteps()/t_diff/1e6 << " MCells/s " << endl;

	//cleanup
	PA.DeleteAll();

    return 0;
}
