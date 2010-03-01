#include <iostream>
#include <time.h>
#include "FDTD/cartoperator.h"
#include "FDTD/engine.h"

#include "ContinuousStructure.h"

using namespace std;

void BuildMSL(ContinuousStructure &CSX);

int main(int argc, char *argv[])
{
	ContinuousStructure CSX;

	time_t startTime=time(NULL);

//	CSX.ReadFromXML("csx-files/MSL.xml");
	BuildMSL(CSX);

	CartOperator cop;
	cop.SetGeometryCSX(&CSX);

	cop.CalcECOperator();

	cop.CalcGaussianPulsExcitation(1e9,1e9);

	time_t OpDoneTime=time(NULL);

	cop.ShowSize();

	cerr << "Time for operator: " << difftime(OpDoneTime,startTime) << endl;

	Engine eng(&cop);

	time_t currTime = time(NULL);

	unsigned int NrIter = 500;

	eng.IterateTS(NrIter);

	time_t prevTime = currTime;
	currTime = time(NULL);

	double t_diff = difftime(currTime,prevTime);

	cerr << "Time for " << NrIter << " iterations with " << cop.GetNumberCells() << " cells : " << t_diff << " sec" << endl;
	cerr << "Speed (MCells/s): " << (double)cop.GetNumberCells()*(double)NrIter/t_diff/1e6 << endl;
    return 0;
}

void BuildMSL(ContinuousStructure &CSX)
{
	CSPropMaterial* mat = new CSPropMaterial(CSX.GetParameterSet());
	mat->SetEpsilon(3.6);
	CSX.AddProperty(mat);

	CSPrimBox* box = new CSPrimBox(CSX.GetParameterSet(),mat);
	box->SetCoord(0,-200.0);box->SetCoord(1,200.0);
	box->SetCoord(2,0.0);box->SetCoord(3,50.0);
	box->SetCoord(4,-1000.0);box->SetCoord(5,1000.0);
	CSX.AddPrimitive(box);

	CSPropMaterial* MSL = new CSPropMaterial(CSX.GetParameterSet());
	MSL->SetKappa(56e6);
	CSX.AddProperty(MSL);

	box = new CSPrimBox(CSX.GetParameterSet(),MSL);
	box->SetCoord(0,-20.0);box->SetCoord(1,20.0);
	box->SetCoord(2,0.0);box->SetCoord(3,50.0);
	box->SetCoord(4,-1000.0);box->SetCoord(5,1000.0);
	CSX.AddPrimitive(box);

	CSPropElectrode* elec = new CSPropElectrode(CSX.GetParameterSet());
	elec->SetExcitation(1,1);
	elec->SetExcitType(0);
//	elec->SetDelay(2.0e-9);
	CSX.AddProperty(elec);

	box = new CSPrimBox(CSX.GetParameterSet(),elec);
	box->SetCoord(0,-20.0);box->SetCoord(1,20.0);
	box->SetCoord(2,0.0);box->SetCoord(3,50.0);
	box->SetCoord(4,0.0);box->SetCoord(5,0.0);
	CSX.AddPrimitive(box);

	CSRectGrid* grid = CSX.GetGrid();

	for (int n=-1000;n<=1000;n+=20)
		grid->AddDiscLine(2,(double)n);
	for (int n=-200;n<=200;n+=10)
		grid->AddDiscLine(0,(double)n);
	for (int n=0;n<=300;n+=10)
		grid->AddDiscLine(1,(double)n);

	grid->SetDeltaUnit(1e-3);

	CSX.Write2XML("csx-files/MSL.xml");
}
