#include "FDTD_examples.h"
#include "../tools/constants.h"

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
	double width = 1000;
	double hight = 1000;
	double length = 4000;
	double abs_l = 200;

	CSPrimBox* box = NULL;
	//fake pml....
	CSPropMaterial* mat = new CSPropMaterial(CSX.GetParameterSet());
//	mat->SetEpsilon(3.6);
	double finalKappa = 0.3/pow(abs_l,4);
	mat->SetKappa(finalKappa);
	std::ostringstream fct;
	fct << "pow(abs(z)-" << length/2.0-abs_l << ",4)";
	mat->SetKappaWeightFunction(fct.str(),0);
	mat->SetKappaWeightFunction(fct.str(),1);
	mat->SetKappaWeightFunction(fct.str(),2);
	mat->SetSigma(finalKappa*__MUE0__/__EPS0__);
	mat->SetSigmaWeightFunction(fct.str(),0);
	mat->SetSigmaWeightFunction(fct.str(),1);
	mat->SetSigmaWeightFunction(fct.str(),2);
	CSX.AddProperty(mat);
	box = new CSPrimBox(CSX.GetParameterSet(),mat);
	box->SetCoord(0,width/-2.0);box->SetCoord(1,width/2.0);
	box->SetCoord(2,hight/-2.0);box->SetCoord(3,hight/2.0);
	box->SetCoord(4,length/2.0-abs_l); box->SetCoord(5,length/2.0);
	CSX.AddPrimitive(box);
	box = new CSPrimBox(CSX.GetParameterSet(),mat);
	box->SetCoord(0,width/-2.0);box->SetCoord(1,width/2.0);
	box->SetCoord(2,hight/-2.0);box->SetCoord(3,hight/2.0);
	box->SetCoord(4,length/-2.0+abs_l); box->SetCoord(5,length/-2.0);
	CSX.AddPrimitive(box);

	CSPropMaterial* mat2 = new CSPropMaterial(CSX.GetParameterSet());
	mat2->SetEpsilon(2);
	CSX.AddProperty(mat2);
	box = new CSPrimBox(CSX.GetParameterSet(),mat2);
	box->SetCoord(0,width/-2.0);box->SetCoord(1,width/2.0);
	box->SetCoord(2,hight/-2.0);box->SetCoord(3,hight/2.0);
	box->SetCoord(4,length/-4.0);box->SetCoord(5,length/4.0);
	CSX.AddPrimitive(box);

	CSPropElectrode* elec = new CSPropElectrode(CSX.GetParameterSet());
	elec->SetExcitation(1,1);
	elec->SetExcitType(0);
	elec->SetActiveDir(0,0);//disable x
	elec->SetActiveDir(0,2);//disable z
//	elec->SetDelay(2.0e-9);
	CSX.AddProperty(elec);

	box = new CSPrimBox(CSX.GetParameterSet(),elec);
	box->SetCoord(0,width/-2.0);box->SetCoord(1,width/2.0);
	box->SetCoord(2,hight/-2.0);box->SetCoord(3,hight/2.0);
	box->SetCoord(4,0.0);box->SetCoord(5,0.0);
	CSX.AddPrimitive(box);

	//E-field dump
	CSPropDumpBox* Edump = new CSPropDumpBox(CSX.GetParameterSet());
	Edump->SetDumpType(0);
	Edump->SetName("tmp/Et_");
	CSX.AddProperty(Edump);
	box = new CSPrimBox(CSX.GetParameterSet(),Edump);
	box->SetCoord(0,width/-3.0);box->SetCoord(1,width/3.0);
	box->SetCoord(2,0.0);box->SetCoord(3,0.0);
	box->SetCoord(4,length/-2.0+abs_l);box->SetCoord(5,length/2.0-abs_l);
	CSX.AddPrimitive(box);

	//H-field dump
	CSPropDumpBox* Hdump = new CSPropDumpBox(CSX.GetParameterSet());
	Hdump->SetDumpType(1);
	Hdump->SetName("tmp/Ht_");
	CSX.AddProperty(Hdump);
	box = new CSPrimBox(CSX.GetParameterSet(),Hdump);
	box->SetCoord(0,width/-3.0);box->SetCoord(1,width/3.0);
	box->SetCoord(2,0.0);box->SetCoord(3,0.0);
	box->SetCoord(4,length/-2.0+abs_l);box->SetCoord(5,length/2.0-abs_l);
	CSX.AddPrimitive(box);

//	CSPropMetal* metal = new CSPropMetal(CSX.GetParameterSet());
//	CSX.AddProperty(metal);
//	CSPrimCylinder* cyl = new CSPrimCylinder(CSX.GetParameterSet(),metal);
//	cyl->SetRadius(100);
//	cyl->SetCoord(0,0.0);cyl->SetCoord(1,0.0);
//	cyl->SetCoord(2,-250.0);cyl->SetCoord(3,250.0);
//	cyl->SetCoord(4,-0000.0);cyl->SetCoord(5,-0000.0);
//	CSX.AddPrimitive(cyl);

	CSRectGrid* grid = CSX.GetGrid();

	for (int n=width/-2.0;n<=width/2;n+=20)
		grid->AddDiscLine(0,(double)n);
	for (int n=hight/-2.0;n<=hight/2.0;n+=20)
		grid->AddDiscLine(1,(double)n);
	for (int n=length/-2.0;n<=length/2.0;n+=20)
		grid->AddDiscLine(2,(double)n);

	grid->SetDeltaUnit(1e-3);

	CSX.Write2XML("tmp/PlaneWave.xml");
}

void BuildMSL(ContinuousStructure &CSX)
{
	double width = 1000;
	double hight = 500;
	double length = 2000;
	double abs_l = 200;

	//substrate....
	CSPropMaterial* mat = new CSPropMaterial(CSX.GetParameterSet());
//	mat->SetEpsilon(3.6);
	double finalKappa = 0.3/pow(abs_l,4);
	mat->SetKappa(finalKappa);
	std::ostringstream fct;
	fct << "pow(abs(z)-" << length/2.0-abs_l << ",4)";
	mat->SetKappaWeightFunction(fct.str(),0);
	mat->SetKappaWeightFunction(fct.str(),1);
	mat->SetKappaWeightFunction(fct.str(),2);
	mat->SetSigma(finalKappa*__MUE0__/__EPS0__);
	mat->SetSigmaWeightFunction(fct.str(),0);
	mat->SetSigmaWeightFunction(fct.str(),1);
	mat->SetSigmaWeightFunction(fct.str(),2);
	CSX.AddProperty(mat);
	CSPrimBox* box = new CSPrimBox(CSX.GetParameterSet(),mat);
	box->SetCoord(0,width/-2.0);box->SetCoord(1,width/2.0);
	box->SetCoord(2,0.0);box->SetCoord(3,hight);
	box->SetCoord(4,length/2.0-abs_l); box->SetCoord(5,length/2.0);
	CSX.AddPrimitive(box);
//	box = new CSPrimBox(CSX.GetParameterSet(),mat);
//	box->SetCoord(0,width/-2.0);box->SetCoord(1,width/2.0);
//	box->SetCoord(2,0.0);box->SetCoord(3,hight);
//	box->SetCoord(4,length/-2.0+abs_l); box->SetCoord(5,length/-2.0);
//	CSX.AddPrimitive(box);

	//MSL
	CSPropMetal* MSL = new CSPropMetal(CSX.GetParameterSet());
	CSX.AddProperty(MSL);
	box = new CSPrimBox(CSX.GetParameterSet(),MSL);
	box->SetCoord(0,-40.0);box->SetCoord(1,40.0);
	box->SetCoord(2,50.0);box->SetCoord(3,50.0);
	box->SetCoord(4,length/-2);box->SetCoord(5,length/2.0);
	box->SetPriority(100);
	CSX.AddPrimitive(box);

	//MSL excite...
	CSPropElectrode* elec = new CSPropElectrode(CSX.GetParameterSet());
	elec->SetExcitation(-1,1);
	elec->SetExcitType(0);
//	elec->SetDelay(2.0e-9);
	CSX.AddProperty(elec);
	box = new CSPrimBox(CSX.GetParameterSet(),elec);
	box->SetCoord(0,-40.0);box->SetCoord(1,40.0);
	box->SetCoord(2,0.0);box->SetCoord(3,50.0);
	box->SetCoord(4,length/-2.0);box->SetCoord(5,length/-2.0);
	CSX.AddPrimitive(box);

	//E-field dump
	CSPropDumpBox* Edump = new CSPropDumpBox(CSX.GetParameterSet());
	Edump->SetDumpType(0);
	Edump->SetName("tmp/Et_");
	CSX.AddProperty(Edump);
	box = new CSPrimBox(CSX.GetParameterSet(),Edump);
	box->SetCoord(0,width/-2.0);box->SetCoord(1,width/2.0);
	box->SetCoord(2,25.0);box->SetCoord(3,25.);
	box->SetCoord(4,length/-2.0);box->SetCoord(5,length/2.0);
	CSX.AddPrimitive(box);

	CSRectGrid* grid = CSX.GetGrid();

	for (double n=width/-2.0;n<=width/2;n+=20)
		grid->AddDiscLine(0,n);
	for (double n=0;n<=500;n+=10)
		grid->AddDiscLine(1,n);
	for (double n=length/-2.0;n<=length/2.0;n+=20)
		grid->AddDiscLine(2,n);

	grid->SetDeltaUnit(1e-3);

	CSX.Write2XML("tmp/MSL.xml");
}


