#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include "tools/array_ops.h"
#include "FDTD/cartoperator.h"
#include "FDTD/engine.h"
#include "FDTD/processvoltage.h"
#include "FDTD/processfields_td.h"

#include "ContinuousStructure.h"

using namespace std;

void BuildMSL(ContinuousStructure &CSX);
void BuildDipol(ContinuousStructure &CSX);

int main(int argc, char *argv[])
{
	time_t startTime=time(NULL);

	//*************** setup/read geometry ************//
	ContinuousStructure CSX;
	BuildDipol(CSX);

	//*************** setup operator ************//
	CartOperator cop;
	cop.SetGeometryCSX(&CSX);

	cop.CalcECOperator();

	double fmax=1e9;
	cop.CalcGaussianPulsExcitation(fmax/2,fmax/2);

	time_t OpDoneTime=time(NULL);

	cop.ShowSize();

	cerr << "Nyquist number of timesteps: " << cop.GetNyquistNum(fmax) << endl;
	unsigned int NrIter = cop.GetNyquistNum(fmax)/3;

	cerr << "Time for operator: " << difftime(OpDoneTime,startTime) << endl;

	//create FDTD engine
	Engine eng(&cop);

	time_t currTime = time(NULL);

	//*************** setup processing ************//
	ProcessVoltage PV(&cop,&eng);
	PV.OpenFile("tmp/u1");
	double start[]={-0,-75,0};
	double stop[]={-0,75,0};
	PV.DefineStartStopCoord(start,stop);
	unsigned int maxIter = 5000;

	ProcessFieldsTD PETD(&cop,&eng);
	start[0]=-1000;start[1]=0;start[2]=-1000;
	stop[0]=1000;stop[1]=0;stop[2]=1000;
	PETD.SetFilePattern("tmp/Et_");
	PETD.DefineStartStopCoord(start,stop);

	PV.Process();
	PETD.Process();

	//*************** simulate ************//
	for (unsigned int i=0;i<maxIter;i+=NrIter)
	{
		eng.IterateTS(NrIter);
		PV.Process();
		PETD.Process();
	}

	//*************** postproc ************//
	time_t prevTime = currTime;
	currTime = time(NULL);

	double t_diff = difftime(currTime,prevTime);

	cerr << "Time for " << eng.GetNumberOfTimesteps() << " iterations with " << cop.GetNumberCells() << " cells : " << t_diff << " sec" << endl;
	cerr << "Speed (MCells/s): " << (double)cop.GetNumberCells()*(double)eng.GetNumberOfTimesteps()/t_diff/1e6 << endl;
    return 0;
}

void BuildDipol(ContinuousStructure &CSX)
{
	CSPropMaterial* mat = new CSPropMaterial(CSX.GetParameterSet());
	mat->SetKappa(0.001);
	CSX.AddProperty(mat);

	CSPrimBox* matbox = new CSPrimBox(CSX.GetParameterSet(),mat);
	matbox->SetCoord(0,-1000.0);matbox->SetCoord(1,1000.0);
	matbox->SetCoord(2,-1000.0);matbox->SetCoord(3,1000.0);
	matbox->SetCoord(4,-1000.0);matbox->SetCoord(5,1000.0);
	CSX.AddPrimitive(matbox);

	CSPropElectrode* elec = new CSPropElectrode(CSX.GetParameterSet());
	elec->SetExcitation(1,1);
	elec->SetExcitType(1);
	elec->SetActiveDir(0,0);//disable x
	elec->SetActiveDir(0,2);//disable z
//	elec->SetDelay(2.0e-9);
	CSX.AddProperty(elec);

	CSPrimBox* box = new CSPrimBox(CSX.GetParameterSet(),elec);
	box->SetCoord(0,-1.0);box->SetCoord(1,1.0);
	box->SetCoord(2,-75.0);box->SetCoord(3,75.0);
	box->SetCoord(4,-1.0);box->SetCoord(5,1.0);
	CSX.AddPrimitive(box);

	CSRectGrid* grid = CSX.GetGrid();

	for (int n=-1000;n<=1000;n+=20)
		grid->AddDiscLine(2,(double)n);
	for (int n=-1000;n<=1000;n+=20)
		grid->AddDiscLine(0,(double)n);
	for (int n=-1000;n<=1000;n+=20)
		grid->AddDiscLine(1,(double)n);

	grid->SetDeltaUnit(1e-3);

	CSX.Write2XML("tmp/Dipol.xml");
}

void BuildMSL(ContinuousStructure &CSX)
{
//	CSPropMaterial* mat = new CSPropMaterial(CSX.GetParameterSet());
//	mat->SetEpsilon(3.6);
//	CSX.AddProperty(mat);
//
//	CSPrimBox* box = new CSPrimBox(CSX.GetParameterSet(),mat);
//	box->SetCoord(0,-200.0);box->SetCoord(1,200.0);
//	box->SetCoord(2,0.0);box->SetCoord(3,50.0);
//	box->SetCoord(4,-1000.0);box->SetCoord(5,1000.0);
//	CSX.AddPrimitive(box);
//
//	CSPropMaterial* MSL = new CSPropMaterial(CSX.GetParameterSet());
//	MSL->SetKappa(56e6);
//	CSX.AddProperty(MSL);
//
//	box = new CSPrimBox(CSX.GetParameterSet(),MSL);
//	box->SetCoord(0,-20.0);box->SetCoord(1,20.0);
//	box->SetCoord(2,0.0);box->SetCoord(3,50.0);
//	box->SetCoord(4,-1000.0);box->SetCoord(5,1000.0);
//	CSX.AddPrimitive(box);

	CSPropElectrode* elec = new CSPropElectrode(CSX.GetParameterSet());
	elec->SetExcitation(1,1);
	elec->SetExcitType(0);
//	elec->SetDelay(2.0e-9);
	CSX.AddProperty(elec);

	CSPrimBox* box = new CSPrimBox(CSX.GetParameterSet(),elec);
	box->SetCoord(0,-1.0);box->SetCoord(1,1.0);
	box->SetCoord(2,0.0);box->SetCoord(3,50.0);
	box->SetCoord(4,-1.0);box->SetCoord(5,1.0);
	CSX.AddPrimitive(box);

	CSRectGrid* grid = CSX.GetGrid();

	for (int n=-300;n<=300;n+=20)
		grid->AddDiscLine(0,(double)n);
	for (int n=0;n<=300;n+=20)
		grid->AddDiscLine(1,(double)n);
	for (int n=-1000;n<=1000;n+=20)
		grid->AddDiscLine(2,(double)n);

	grid->SetDeltaUnit(1e-3);

	CSX.Write2XML("tmp/MSL.xml");
}
