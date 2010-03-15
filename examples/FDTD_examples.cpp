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

#include "FDTD_examples.h"
#include "../tools/constants.h"

void BuildDipol(const char* filename)
{
	int maxIter = 1000;
	double f0=0.5e9;
	double fc=0.5e9;
	int Excit_Type=0;
	int bounds[] = {1,1,0,0,0,0};

	cerr << "Create Geometry..." << endl;
	ContinuousStructure CSX;

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

	//*************** Create XML file **********************
	TiXmlDocument doc(filename);
	doc.InsertEndChild(TiXmlDeclaration("1.0","ISO-8859-1","yes"));

	TiXmlElement openEMS("openEMS");

	TiXmlElement FDTD_Opts("FDTD");
	FDTD_Opts.SetAttribute("NumberOfTimesteps",maxIter);

	TiXmlElement Excite("Excitation");
	Excite.SetAttribute("Type",Excit_Type);
	Excite.SetAttribute("f0",f0);
	Excite.SetAttribute("fc",fc);
	FDTD_Opts.InsertEndChild(Excite);

	TiXmlElement BC("BoundaryCond");
	BC.SetAttribute("xmin",bounds[0]);
	BC.SetAttribute("xmax",bounds[1]);
	BC.SetAttribute("ymin",bounds[2]);
	BC.SetAttribute("ymax",bounds[3]);
	BC.SetAttribute("zmin",bounds[4]);
	BC.SetAttribute("zmax",bounds[5]);
	FDTD_Opts.InsertEndChild(BC);

	openEMS.InsertEndChild(FDTD_Opts);

	if (CSX.Write2XML(&openEMS,true)==false)
	{
		cerr << "writing failed" << endl;
		exit(-1);
	}

	doc.InsertEndChild(openEMS);

	doc.SaveFile();
}

void BuildPlaneWave(const char* filename)
{
	int maxIter = 1000;
	double f0=0.5e9;
	double fc=0.5e9;
	int Excit_Type=0;
	int bounds[] = {1,1,0,0,0,0};

	cerr << "Create Geometry..." << endl;
	ContinuousStructure CSX;

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
	Edump->SetName("Et_");
	CSX.AddProperty(Edump);
	box = new CSPrimBox(CSX.GetParameterSet(),Edump);
	box->SetCoord(0,width/-3.0);box->SetCoord(1,width/3.0);
	box->SetCoord(2,0.0);box->SetCoord(3,0.0);
	box->SetCoord(4,length/-2.0+abs_l);box->SetCoord(5,length/2.0-abs_l);
	CSX.AddPrimitive(box);

	//H-field dump
	CSPropDumpBox* Hdump = new CSPropDumpBox(CSX.GetParameterSet());
	Hdump->SetDumpType(1);
	Hdump->SetName("Ht_");
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

	//*************** Create XML file **********************
	TiXmlDocument doc(filename);
	doc.InsertEndChild(TiXmlDeclaration("1.0","ISO-8859-1","yes"));

	TiXmlElement openEMS("openEMS");

	TiXmlElement FDTD_Opts("FDTD");
	FDTD_Opts.SetAttribute("NumberOfTimesteps",maxIter);

	TiXmlElement Excite("Excitation");
	Excite.SetAttribute("Type",Excit_Type);
	Excite.SetAttribute("f0",f0);
	Excite.SetAttribute("fc",fc);
	FDTD_Opts.InsertEndChild(Excite);

	TiXmlElement BC("BoundaryCond");
	BC.SetAttribute("xmin",bounds[0]);
	BC.SetAttribute("xmax",bounds[1]);
	BC.SetAttribute("ymin",bounds[2]);
	BC.SetAttribute("ymax",bounds[3]);
	BC.SetAttribute("zmin",bounds[4]);
	BC.SetAttribute("zmax",bounds[5]);
	FDTD_Opts.InsertEndChild(BC);

	openEMS.InsertEndChild(FDTD_Opts);

	if (CSX.Write2XML(&openEMS,true)==false)
	{
		cerr << "writing failed" << endl;
		exit(-1);
	}

	doc.InsertEndChild(openEMS);

	doc.SaveFile();
}

void BuildMSL(const char* filename)
{
	int maxIter = 1000;
	double f0=0.5e9;
	double fc=0.5e9;
	int Excit_Type=0;
	int bounds[] = {1,1,0,0,0,0};

	cerr << "Create Geometry..." << endl;
	ContinuousStructure CSX;

	double width = 1000;
	double height = 500;
	double length = 2000;
	double abs_l = 200;
	double MSL_height=50;
	double MSL_width=80;
	double delta[] = {20,10,20};

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
	box->SetCoord(2,0.0);box->SetCoord(3,height);
	box->SetCoord(4,length/2.0-abs_l); box->SetCoord(5,length/2.0);
	CSX.AddPrimitive(box);
//	box = new CSPrimBox(CSX.GetParameterSet(),mat);
//	box->SetCoord(0,width/-2.0);box->SetCoord(1,width/2.0);
//	box->SetCoord(2,0.0);box->SetCoord(3,hight);
//	box->SetCoord(4,length/-2.0+abs_l); box->SetCoord(5,length/-2.0);
//	CSX.AddPrimitive(box);

	//MSL
	CSProperties* MSL = NULL;
	CSPropMaterial* MSL_mat = new CSPropMaterial(CSX.GetParameterSet());
	MSL_mat->SetKappa(56e6);
	MSL = MSL_mat;
//	MSL = new CSPropMetal(CSX.GetParameterSet());

	CSX.AddProperty(MSL);
	box = new CSPrimBox(CSX.GetParameterSet(),MSL);
	box->SetCoord(0,MSL_width/-2.0);box->SetCoord(1,MSL_width/2.0);
	box->SetCoord(2,MSL_height);box->SetCoord(3,MSL_height+delta[1]);
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

	CSPropDumpBox* Edump = NULL;
//	//E-field dump xz
//	Edump = new CSPropDumpBox(CSX.GetParameterSet());
//	Edump->SetDumpType(0);
//	Edump->SetName("Et_xz_");
//	CSX.AddProperty(Edump);
//	box = new CSPrimBox(CSX.GetParameterSet(),Edump);
//	box->SetCoord(0,width/-2.0);box->SetCoord(1,width/2.0);
//	box->SetCoord(2,25.0);box->SetCoord(3,25.);
//	box->SetCoord(4,length/-2.0);box->SetCoord(5,length/2.0);
//	CSX.AddPrimitive(box);
//
//	//E-field dump xy
//	Edump = new CSPropDumpBox(CSX.GetParameterSet());
//	Edump->SetDumpType(0);
//	Edump->SetName("Et_xy_");
//	CSX.AddProperty(Edump);
//	box = new CSPrimBox(CSX.GetParameterSet(),Edump);
//	box->SetCoord(0,width/-2.0);box->SetCoord(1,width/2.0);
//	box->SetCoord(2,0.0);box->SetCoord(3,height);
//	box->SetCoord(4,0.0);box->SetCoord(5,0.0);
//	CSX.AddPrimitive(box);

	//E-field dump 3D
	Edump = new CSPropDumpBox(CSX.GetParameterSet());
	Edump->SetDumpType(0);
	Edump->SetDumpMode(2); //cell interpolated dump
	Edump->SetName("Et_");
	CSX.AddProperty(Edump);
	box = new CSPrimBox(CSX.GetParameterSet(),Edump);
	box->SetCoord(0,MSL_width*-1.5);box->SetCoord(1,MSL_width*1.5);
	box->SetCoord(2,0.0);box->SetCoord(3,MSL_height*1.5);
	box->SetCoord(4,length/-2.0);box->SetCoord(5,length/2.0);
	CSX.AddPrimitive(box);

	//voltage calc
	CSPropProbeBox* volt = new CSPropProbeBox(CSX.GetParameterSet());
	volt->SetProbeType(0);
	volt->SetName("u1");
	CSX.AddProperty(volt);
	box = new CSPrimBox(CSX.GetParameterSet(),volt);
	box->SetCoord(0,0.0);box->SetCoord(1,0.0);
	box->SetCoord(2,MSL_height);box->SetCoord(3,0.0);
	box->SetCoord(4,0.0);box->SetCoord(5,0.0);
	CSX.AddPrimitive(box);

	//current calc
	CSPropProbeBox* curr = new CSPropProbeBox(CSX.GetParameterSet());
	curr->SetProbeType(1);
	curr->SetName("i1");
	CSX.AddProperty(curr);
	box = new CSPrimBox(CSX.GetParameterSet(),curr);
	box->SetCoord(0,MSL_width*-1.5);box->SetCoord(1,MSL_width*1.5);
	box->SetCoord(2,MSL_height/2.0);box->SetCoord(3,MSL_height*1.5);
	box->SetCoord(4,0.0);box->SetCoord(5,0.0);
	CSX.AddPrimitive(box);

	CSRectGrid* grid = CSX.GetGrid();

	for (double n=width/-2.0;n<=width/2;n+=delta[0])
		grid->AddDiscLine(0,n);
	for (double n=0;n<=height;n+=delta[1])
		grid->AddDiscLine(1,n);
	for (double n=length/-2.0;n<=length/2.0;n+=delta[2])
		grid->AddDiscLine(2,n);

	grid->SetDeltaUnit(1e-3);

	//*************** Create XML file **********************
	TiXmlDocument doc(filename);
	doc.InsertEndChild(TiXmlDeclaration("1.0","ISO-8859-1","yes"));

	TiXmlElement openEMS("openEMS");

	TiXmlElement FDTD_Opts("FDTD");
	FDTD_Opts.SetAttribute("NumberOfTimesteps",maxIter);

	TiXmlElement Excite("Excitation");
	Excite.SetAttribute("Type",Excit_Type);
	Excite.SetAttribute("f0",f0);
	Excite.SetAttribute("fc",fc);
	FDTD_Opts.InsertEndChild(Excite);

	TiXmlElement BC("BoundaryCond");
	BC.SetAttribute("xmin",bounds[0]);
	BC.SetAttribute("xmax",bounds[1]);
	BC.SetAttribute("ymin",bounds[2]);
	BC.SetAttribute("ymax",bounds[3]);
	BC.SetAttribute("zmin",bounds[4]);
	BC.SetAttribute("zmax",bounds[5]);
	FDTD_Opts.InsertEndChild(BC);

	openEMS.InsertEndChild(FDTD_Opts);

	if (CSX.Write2XML(&openEMS,true)==false)
	{
		cerr << "writing failed" << endl;
		exit(-1);
	}

	doc.InsertEndChild(openEMS);

	doc.SaveFile();
}

void BuildCoaxial_Cartesian(const char* filename)
{
	int maxIter = 1000;
	double f0=1e9;
	double fc=1e9;
	int Excit_Type=0;
	int bounds[] = {0,0,0,0,0,0};

	cerr << "Create Geometry..." << endl;
	ContinuousStructure CSX;

	double rad[] = {100, 230};
	double length[] = {-500,1500};
	double abs_l = 200;

	double delta[] = {5,5,10};

	CSPrimBox* box = NULL;
	//fake pml....
	CSPropMaterial* mat = new CSPropMaterial(CSX.GetParameterSet());
//	mat->SetEpsilon(3.6);
	double finalKappa = 0.3/pow(abs_l,4);
	mat->SetKappa(finalKappa);
	std::ostringstream fct;
	fct << "pow(abs(z)-" << length[1]-abs_l << ",4)";
	mat->SetKappaWeightFunction(fct.str(),0);
	mat->SetKappaWeightFunction(fct.str(),1);
	mat->SetKappaWeightFunction(fct.str(),2);
	mat->SetSigma(finalKappa*__MUE0__/__EPS0__);
	mat->SetSigmaWeightFunction(fct.str(),0);
	mat->SetSigmaWeightFunction(fct.str(),1);
	mat->SetSigmaWeightFunction(fct.str(),2);
	CSX.AddProperty(mat);
	
	box = new CSPrimBox(CSX.GetParameterSet(),mat);
	box->SetCoord(0,-1.0*rad[1]);box->SetCoord(1,1.0*rad[1]);
	box->SetCoord(2,-1.0*rad[1]);box->SetCoord(3,1.0*rad[1]);
	box->SetCoord(4,length[1]-abs_l); box->SetCoord(5,length[1]);
	box->SetPriority(10);
	CSX.AddPrimitive(box);

	CSPropElectrode* elec = new CSPropElectrode(CSX.GetParameterSet());
	elec->SetExcitation(1.0,0);
	elec->SetExcitation(1.0,1);
	elec->SetWeightFunction("x/pow(rho,2)",0);
	elec->SetWeightFunction("y/pow(rho,2)",1);
	elec->SetExcitType(0);
//	elec->SetActiveDir(0,0);//disable x
//	elec->SetActiveDir(0,2);//disable z
//	elec->SetDelay(2.0e-9);
	CSX.AddProperty(elec);

//	double coords[] = {-100,-100,0};
//	cerr << elec->GetWeightedExcitation(0,coords) << endl;
//	cerr << elec->GetWeightedExcitation(1,coords) << endl;
//	exit(0);

	box = new CSPrimBox(CSX.GetParameterSet(),elec);
	box->SetCoord(0,-1.0*rad[1]);box->SetCoord(1,1.0*rad[1]);
	box->SetCoord(2,-1.0*rad[1]);box->SetCoord(3,1.0*rad[1]);
	box->SetCoord(4,length[0]);box->SetCoord(5,length[0]);
	box->SetPriority(5);
	CSX.AddPrimitive(box);

	//E-field dump
	CSPropDumpBox* Edump = new CSPropDumpBox(CSX.GetParameterSet());
	Edump->SetDumpType(0);
	Edump->SetDumpMode(2);
	Edump->SetName("Et_");
	CSX.AddProperty(Edump);
	box = new CSPrimBox(CSX.GetParameterSet(),Edump);
	box->SetCoord(0,-1*rad[1]);box->SetCoord(1,rad[1]);
	box->SetCoord(2,-0*rad[1]);box->SetCoord(3,0*rad[1]);
	box->SetCoord(4,length[0]);box->SetCoord(5,length[1]);
	CSX.AddPrimitive(box);

	//voltage calc
	CSPropProbeBox* volt = new CSPropProbeBox(CSX.GetParameterSet());
	volt->SetProbeType(0);
	volt->SetName("u1");
	CSX.AddProperty(volt);
	box = new CSPrimBox(CSX.GetParameterSet(),volt);
	box->SetCoord(0,rad[0]);box->SetCoord(1,rad[1]);
	box->SetCoord(2,0.0);box->SetCoord(3,0.0);
	box->SetCoord(4,0.0);box->SetCoord(5,0.0);
	CSX.AddPrimitive(box);

	//current calc
	CSPropProbeBox* curr = new CSPropProbeBox(CSX.GetParameterSet());
	curr->SetProbeType(1);
	curr->SetName("i1");
	CSX.AddProperty(curr);
	box = new CSPrimBox(CSX.GetParameterSet(),curr);
	box->SetCoord(0,-1.5*rad[0]);box->SetCoord(1,1.5*rad[0]);
	box->SetCoord(2,-1.5*rad[0]);box->SetCoord(3,1.5*rad[0]);
	box->SetCoord(4,0.0);box->SetCoord(5,0.0);
	CSX.AddPrimitive(box);

	CSPropMaterial* metal = new CSPropMaterial(CSX.GetParameterSet());
	metal->SetKappa(56e9);
//	CSPropMetal* metal = new CSPropMetal(CSX.GetParameterSet());
	CSX.AddProperty(metal);
	CSPrimCylinder* cyl = new CSPrimCylinder(CSX.GetParameterSet(),metal);
	cyl->SetRadius(rad[0]);
	cyl->SetCoord(0,0.0);cyl->SetCoord(1,0.0);
	cyl->SetCoord(2,0.0);cyl->SetCoord(3,0.0);
	cyl->SetCoord(4,length[0]);cyl->SetCoord(5,length[1]);
	cyl->SetPriority(100);
	CSX.AddPrimitive(cyl);
	box = new CSPrimBox(CSX.GetParameterSet(),metal);
	box->SetCoord(0,-1.0*rad[1]);box->SetCoord(1,1.0*rad[1]);
	box->SetCoord(2,-1.0*rad[1]);box->SetCoord(3,1.0*rad[1]);
	box->SetCoord(4,length[0]);box->SetCoord(5,length[1]);
	box->SetPriority(1);
	CSX.AddPrimitive(box);

	CSPropMaterial* air = new CSPropMaterial(CSX.GetParameterSet());
	CSX.AddProperty(air);
	cyl = new CSPrimCylinder(CSX.GetParameterSet(),air);
	cyl->SetRadius(rad[1]);
	cyl->SetCoord(0,0.0);cyl->SetCoord(1,0.0);
	cyl->SetCoord(2,0.0);cyl->SetCoord(3,0.0);
	cyl->SetCoord(4,length[0]);cyl->SetCoord(5,length[1]);
	cyl->SetPriority(9);
	CSX.AddPrimitive(cyl);

	CSRectGrid* grid = CSX.GetGrid();

	for (int n=-1.0*rad[1];n<=rad[1];n+=delta[0])
		grid->AddDiscLine(0,(double)n);
	for (int n=-1.0*rad[1];n<=rad[1];n+=delta[1])
		grid->AddDiscLine(1,(double)n);
	for (int n=length[0];n<=length[1];n+=delta[2])
		grid->AddDiscLine(2,(double)n);

	grid->SetDeltaUnit(1e-3);

	//*************** Create XML file **********************
	TiXmlDocument doc(filename);
	doc.InsertEndChild(TiXmlDeclaration("1.0","ISO-8859-1","yes"));

	TiXmlElement openEMS("openEMS");

	TiXmlElement FDTD_Opts("FDTD");
	FDTD_Opts.SetAttribute("NumberOfTimesteps",maxIter);

	TiXmlElement Excite("Excitation");
	Excite.SetAttribute("Type",Excit_Type);
	Excite.SetAttribute("f0",f0);
	Excite.SetAttribute("fc",fc);
	FDTD_Opts.InsertEndChild(Excite);

	TiXmlElement BC("BoundaryCond");
	BC.SetAttribute("xmin",bounds[0]);
	BC.SetAttribute("xmax",bounds[1]);
	BC.SetAttribute("ymin",bounds[2]);
	BC.SetAttribute("ymax",bounds[3]);
	BC.SetAttribute("zmin",bounds[4]);
	BC.SetAttribute("zmax",bounds[5]);
	FDTD_Opts.InsertEndChild(BC);

	openEMS.InsertEndChild(FDTD_Opts);

	if (CSX.Write2XML(&openEMS,true)==false)
	{
		cerr << "writing failed" << endl;
		exit(-1);
	}

	doc.InsertEndChild(openEMS);

	doc.SaveFile();
}

void BuildHelix(const char* filename)
{
	int maxIter = 1e5;
	double f0=0.5e9;
	double fc=0.5e9;
	int Excit_Type=0;
	int bounds[] = {1,1,1,1,1,1};

	cerr << "Create Helix Geometry..." << endl;
	ContinuousStructure CSX;

	double feed_length=10;
	double wire_rad = 0.7;
	double coil_rad = 10;
	double coil_length = 50;
	double delta[] = {0.5,0.5,0.5};

	//MSL
	CSPropMaterial* copper = new CSPropMaterial(CSX.GetParameterSet());
	copper->SetKappa(56e6);
	copper->SetName("copper");
	CSX.AddProperty(copper);

	CSPrimUserDefined* helix = new CSPrimUserDefined(CSX.GetParameterSet(),copper);
	helix->SetCoordSystem(CSPrimUserDefined::CYLINDER_SYSTEM);
	helix->SetFunction("(r>9)&(r<11)&(sqrt(pow(x-r*cos(2*pi*z/6.25),2)+pow(y-r*sin(2*pi*z/6.25),2))<2)&(z>0)&(z<50)");
	CSX.AddPrimitive(helix);
	CSPrimCylinder* cyl = new CSPrimCylinder(CSX.GetParameterSet(),copper);
	cyl->SetRadius(wire_rad);
	cyl->SetCoord(0,coil_rad);cyl->SetCoord(1,coil_rad+feed_length);
	cyl->SetCoord(2,0.0);cyl->SetCoord(3,0.0);
	cyl->SetCoord(4,0.0);cyl->SetCoord(5,0.0);
	CSX.AddPrimitive(cyl);
	cyl = new CSPrimCylinder(CSX.GetParameterSet(),copper);
	cyl->SetRadius(wire_rad);
	cyl->SetCoord(0,coil_rad);cyl->SetCoord(1,coil_rad+feed_length);
	cyl->SetCoord(2,0.0);cyl->SetCoord(3,0.0);
	cyl->SetCoord(4,coil_length);cyl->SetCoord(5,coil_length);
	CSX.AddPrimitive(cyl);

	double kappa_resist = (coil_length-2.0*delta[2])/(PI*wire_rad*wire_rad)/50/1e-3;
	CSPropMaterial* Src_Resist = new CSPropMaterial(CSX.GetParameterSet());
	Src_Resist->SetKappa(kappa_resist);
	Src_Resist->SetName("resist");
	CSX.AddProperty(Src_Resist);
	cyl = new CSPrimCylinder(CSX.GetParameterSet(),Src_Resist);
	cyl->SetRadius(wire_rad);
	cyl->SetCoord(0,coil_rad+feed_length);cyl->SetCoord(1,coil_rad+feed_length);
	cyl->SetCoord(2,0.0);cyl->SetCoord(3,0.0);
	cyl->SetCoord(4,0.0);cyl->SetCoord(5,coil_length/2.0-delta[2]);
	CSX.AddPrimitive(cyl);
	cyl = new CSPrimCylinder(CSX.GetParameterSet(),Src_Resist);
	cyl->SetRadius(wire_rad);
	cyl->SetCoord(0,coil_rad+feed_length);cyl->SetCoord(1,coil_rad+feed_length);
	cyl->SetCoord(2,0.0);cyl->SetCoord(3,0.0);
	cyl->SetCoord(4,coil_length/2.0+delta[2]);cyl->SetCoord(5,coil_length);
	CSX.AddPrimitive(cyl);

	CSPropElectrode* elec = new CSPropElectrode(CSX.GetParameterSet());
	elec->SetExcitation(1.0,2);
	elec->SetExcitType(1);
//	elec->SetDelay(2.0e-9);
	CSX.AddProperty(elec);
	cyl = new CSPrimCylinder(CSX.GetParameterSet(),elec);
	cyl->SetRadius(wire_rad);
	cyl->SetCoord(0,coil_rad+feed_length);cyl->SetCoord(1,coil_rad+feed_length);
	cyl->SetCoord(2,0.0);cyl->SetCoord(3,0.0);
	cyl->SetCoord(4,coil_length/2.0-delta[2]);cyl->SetCoord(5,coil_length/2.0+delta[2]);
	CSX.AddPrimitive(cyl);

	CSPropDumpBox* Edump = NULL;
	CSPrimBox* box = NULL;
	//E-field dump xz
//	Edump = new CSPropDumpBox(CSX.GetParameterSet());
//	Edump->SetDumpType(0);
//	Edump->SetName("Et_xz_");
//	CSX.AddProperty(Edump);
//	box = new CSPrimBox(CSX.GetParameterSet(),Edump);
//	box->SetCoord(0,coil_rad/-2.0-25.0);box->SetCoord(1,coil_rad/2.0+25.0+feed_length);
//	box->SetCoord(2,0.0);box->SetCoord(3,0.0);
//	box->SetCoord(4,-25.0);box->SetCoord(5,coil_length+25.0);
//	CSX.AddPrimitive(box);
//
//	//E-field dump xy
//	Edump = new CSPropDumpBox(CSX.GetParameterSet());
//	Edump->SetDumpType(0);
//	Edump->SetName("Et_xy_");
//	CSX.AddProperty(Edump);
//	box = new CSPrimBox(CSX.GetParameterSet(),Edump);
//	box->SetCoord(0,width/-2.0);box->SetCoord(1,width/2.0);
//	box->SetCoord(2,0.0);box->SetCoord(3,height);
//	box->SetCoord(4,0.0);box->SetCoord(5,0.0);
//	CSX.AddPrimitive(box);

//	//E-field dump 3D
//	Edump = new CSPropDumpBox(CSX.GetParameterSet());
//	Edump->SetDumpType(0);
//	Edump->SetDumpMode(2); //cell interpolated dump
//	Edump->SetName("Et_");
//	CSX.AddProperty(Edump);
//	box = new CSPrimBox(CSX.GetParameterSet(),Edump);
//	box->SetCoord(0,MSL_width*-1.5);box->SetCoord(1,MSL_width*1.5);
//	box->SetCoord(2,0.0);box->SetCoord(3,MSL_height*1.5);
//	box->SetCoord(4,length/-2.0);box->SetCoord(5,length/2.0);
//	CSX.AddPrimitive(box);
//
	//voltage calc
	CSPropProbeBox* volt = new CSPropProbeBox(CSX.GetParameterSet());
	volt->SetProbeType(0);
	volt->SetName("u1");
	CSX.AddProperty(volt);
	box = new CSPrimBox(CSX.GetParameterSet(),volt);
	box->SetCoord(0,coil_rad+feed_length);box->SetCoord(1,coil_rad+feed_length);
	box->SetCoord(2,0.0);box->SetCoord(3,0.0);
	box->SetCoord(4,coil_length/2.0+delta[2]);box->SetCoord(5,coil_length/2.0-delta[2]);
	CSX.AddPrimitive(box);

	//current calc
	double curr_dist = 2;
	CSPropProbeBox* curr = new CSPropProbeBox(CSX.GetParameterSet());
	curr->SetProbeType(1);
	curr->SetName("i1");
	CSX.AddProperty(curr);
	box = new CSPrimBox(CSX.GetParameterSet(),curr);
	box->SetCoord(0,coil_rad+feed_length-curr_dist);box->SetCoord(1,coil_rad+feed_length+curr_dist);
	box->SetCoord(2,-curr_dist);box->SetCoord(3,curr_dist);
	box->SetCoord(4,coil_length/2.0);box->SetCoord(5,coil_length/2.0);
	CSX.AddPrimitive(box);

	CSRectGrid* grid = CSX.GetGrid();

	for (double n=coil_rad/-2.0-25.0;n<=coil_rad/2.0+25.0+feed_length;n+=delta[0])
		grid->AddDiscLine(0,n);
	for (double n=coil_rad/-2.0-25.0;n<=coil_rad/2.0+25.0;n+=delta[1])
		grid->AddDiscLine(1,n);
	for (double n=-25.0;n<=coil_length+25.0;n+=delta[2])
		grid->AddDiscLine(2,n);

	grid->SetDeltaUnit(1e-3);

	//*************** Create XML file **********************
	TiXmlDocument doc(filename);
	doc.InsertEndChild(TiXmlDeclaration("1.0","ISO-8859-1","yes"));

	TiXmlElement openEMS("openEMS");

	TiXmlElement FDTD_Opts("FDTD");
	FDTD_Opts.SetAttribute("NumberOfTimesteps",maxIter);

	TiXmlElement Excite("Excitation");
	Excite.SetAttribute("Type",Excit_Type);
	Excite.SetAttribute("f0",f0);
	Excite.SetAttribute("fc",fc);
	FDTD_Opts.InsertEndChild(Excite);

	TiXmlElement BC("BoundaryCond");
	BC.SetAttribute("xmin",bounds[0]);
	BC.SetAttribute("xmax",bounds[1]);
	BC.SetAttribute("ymin",bounds[2]);
	BC.SetAttribute("ymax",bounds[3]);
	BC.SetAttribute("zmin",bounds[4]);
	BC.SetAttribute("zmax",bounds[5]);
	FDTD_Opts.InsertEndChild(BC);

	openEMS.InsertEndChild(FDTD_Opts);

	if (CSX.Write2XML(&openEMS,true)==false)
	{
		cerr << "writing failed" << endl;
		exit(-1);
	}

	doc.InsertEndChild(openEMS);

	doc.SaveFile();
}
