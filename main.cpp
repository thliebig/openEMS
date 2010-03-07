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

	//BuildPlaneWave(CSX);
	BuildMSL(CSX);

	//*************** setup operator ************//
	cerr << "Create Operator..." << endl;
	Operator cop;
	cop.SetGeometryCSX(&CSX);

	cop.CalcECOperator();

	double fmax=1e9;
	cop.CalcGaussianPulsExcitation(fmax/2,fmax/2);

	time_t OpDoneTime=time(NULL);

	cop.ShowSize();
	bool bnd[] = {1,1,0,1,0,1};
	cop.ApplyMagneticBC(bnd);

	cerr << "Nyquist number of timesteps: " << cop.GetNyquistNum(fmax) << endl;
	unsigned int NrIter = cop.GetNyquistNum(fmax)/3;

	cerr << "Time for operator: " << difftime(OpDoneTime,startTime) << endl;

	//create FDTD engine
	Engine eng(&cop);

	time_t currTime = time(NULL);

	//*************** setup processing ************//
	ProcessVoltage PV(&cop,&eng);
	PV.OpenFile("tmp/u1");
	double start[]={0,50,0};
	double stop[]={0,0,0};
	PV.DefineStartStopCoord(start,stop);

	ProcessCurrent PCurr(&cop,&eng);
	PCurr.OpenFile("tmp/i1");
	start[0]=-50;start[1]=40;start[2]=-0;
	stop[0]=50;stop[1]=60;stop[2]=-0;
	PCurr.DefineStartStopCoord(start,stop);

	unsigned int maxIter = 1800;

	vector<CSProperties*> DumpProps = CSX.GetPropertyByType(CSProperties::DUMPBOX);
	vector<ProcessFieldsTD*> TD_Dump;
	ProcessFieldsTD PETD(&cop,&eng);
	if (DumpProps.size()>0)
	{
		CSPrimitives* prim = DumpProps.at(0)->GetPrimitive(0);
		bool acc;
		double* bnd = prim->GetBoundBox(acc);
		start[0]= bnd[0];start[1]=bnd[2];start[2]=bnd[4];
		stop[0] = bnd[1];stop[1] =bnd[3];stop[2] =bnd[5];
		CSPropDumpBox* db = DumpProps.at(0)->ToDumpBox();
		if (db)
		{
			if (db->GetEFieldDump())
			{
				PETD.SetDumpType(0);
				PETD.SetFilePattern("tmp/Et_");
			}
			else
				PETD.SetEnable(false);
		}
		else
			PETD.SetEnable(false);
		PETD.DefineStartStopCoord(start,stop);
	}

	PETD.SetEnable(true);

	PV.Process();
	PCurr.Process();
	PETD.Process();

	//*************** simulate ************//
	for (unsigned int i=0;i<maxIter;i+=NrIter)
	{
		eng.IterateTS(NrIter);
		PV.Process();
		PCurr.Process();
		PETD.Process();
	}

	//*************** postproc ************//
	time_t prevTime = currTime;
	currTime = time(NULL);

	double t_diff = difftime(currTime,prevTime);

	cerr << "Time for " << eng.GetNumberOfTimesteps() << " iterations with " << cop.GetNumberCells() << " cells : " << t_diff << " sec" << endl;
	cerr << "Speed: " << (double)cop.GetNumberCells()*(double)eng.GetNumberOfTimesteps()/t_diff/1e6 << " MCells/s " << endl;
    return 0;
}
