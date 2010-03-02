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
void BuildPlaneWave(ContinuousStructure &CSX);

int main(int argc, char *argv[])
{
	time_t startTime=time(NULL);

	//*************** setup/read geometry ************//
	ContinuousStructure CSX;

	BuildPlaneWave(CSX);

	//*************** setup operator ************//
	CartOperator cop;
	cop.SetGeometryCSX(&CSX);

	cop.CalcECOperator();

	double fmax=1e9;
	cop.CalcGaussianPulsExcitation(fmax/2,fmax/8);

	time_t OpDoneTime=time(NULL);

	cop.ShowSize();
	bool bnd[] = {1,1,0,0,0,0};
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
	double start[]={-100,-250,0};
	double stop[]={-100,250,0};
	PV.DefineStartStopCoord(start,stop);
	unsigned int maxIter = 2000;

	ProcessFieldsTD PETD(&cop,&eng);
	start[0]=-480;start[1]=0;start[2]=-4000;
	stop[0]=480;stop[1]=0;stop[2]=4000;
//	start[0]=-250;start[1]=-250;start[2]=0;
//	stop[0]=250;stop[1]=250;stop[2]=0;
	PETD.SetDumpType(0);
	PETD.SetFilePattern("tmp/Et_");
	PETD.DefineStartStopCoord(start,stop);

	ProcessFieldsTD PHTD(&cop,&eng);
	start[0]=-1000;start[2]=0;start[1]=-1000;
	stop[0]=1000;stop[2]=0;stop[1]=1000;
	PHTD.SetDumpType(1);
	PHTD.SetFilePattern("tmp/Ht_");
	PHTD.DefineStartStopCoord(start,stop);

	PETD.SetEnable(true);
	PHTD.SetEnable(false);

	PV.Process();
	PETD.Process();
	PHTD.Process();	

	//*************** simulate ************//
	for (unsigned int i=0;i<maxIter;i+=NrIter)
	{
		eng.IterateTS(NrIter);
		PV.Process();
		PETD.Process();
		PHTD.Process();
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
	box->SetCoord(0,-10.0);box->SetCoord(1,10.0);
	box->SetCoord(2,-75.0);box->SetCoord(3,75.0);
	box->SetCoord(4,-10.0);box->SetCoord(5,10.0);
	CSX.AddPrimitive(box);

	CSRectGrid* grid = CSX.GetGrid();

	for (int n=-990;n<=990;n+=20)
		grid->AddDiscLine(2,(double)n);
	for (int n=-990;n<=990;n+=20)
		grid->AddDiscLine(0,(double)n);
	for (int n=-990;n<=990;n+=20)
		grid->AddDiscLine(1,(double)n);

	grid->SetDeltaUnit(1e-3);

	CSX.Write2XML("tmp/Dipol.xml");
}

void BuildPlaneWave(ContinuousStructure &CSX)
{
//	CSPropMaterial* mat = new CSPropMaterial(CSX.GetParameterSet());
////	mat->SetKappa(0.001);
//	CSX.AddProperty(mat);
//
//	CSPrimBox* matbox = new CSPrimBox(CSX.GetParameterSet(),mat);
//	matbox->SetCoord(0,-1000.0);matbox->SetCoord(1,1000.0);
//	matbox->SetCoord(2,-1000.0);matbox->SetCoord(3,1000.0);
//	matbox->SetCoord(4,-4000.0);matbox->SetCoord(5,4000.0);
//	CSX.AddPrimitive(matbox);

	CSPropElectrode* elec = new CSPropElectrode(CSX.GetParameterSet());
	elec->SetExcitation(1,1);
	elec->SetExcitType(0);
	elec->SetActiveDir(0,0);//disable x
	elec->SetActiveDir(0,2);//disable z
//	elec->SetDelay(2.0e-9);
	CSX.AddProperty(elec);

	CSPrimBox* box = new CSPrimBox(CSX.GetParameterSet(),elec);
	box->SetCoord(0,-500.0);box->SetCoord(1,500.0);
	box->SetCoord(2,-500.0);box->SetCoord(3,500.0);
	box->SetCoord(4,-4000.0);box->SetCoord(5,-4000.0);
	CSX.AddPrimitive(box);

	CSPropMetal* metal = new CSPropMetal(CSX.GetParameterSet());
	CSX.AddProperty(metal);
	CSPrimCylinder* cyl = new CSPrimCylinder(CSX.GetParameterSet(),metal);
	cyl->SetRadius(100);
	cyl->SetCoord(0,0.0);cyl->SetCoord(1,0.0);
	cyl->SetCoord(2,-250.0);cyl->SetCoord(3,250.0);
	cyl->SetCoord(4,-0000.0);cyl->SetCoord(5,-0000.0);
	CSX.AddPrimitive(cyl);

	CSRectGrid* grid = CSX.GetGrid();

	for (int n=-500;n<=500;n+=20)
		grid->AddDiscLine(0,(double)n);
	for (int n=-500;n<=500;n+=20)
		grid->AddDiscLine(1,(double)n);
	for (int n=-4000;n<=4000;n+=20)
		grid->AddDiscLine(2,(double)n);

	grid->SetDeltaUnit(1e-3);

	CSX.Write2XML("tmp/PlaneWave.xml");
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
