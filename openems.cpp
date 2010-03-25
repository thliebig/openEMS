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

#include "openems.h"
#include <iomanip>
#include "tools/array_ops.h"
#include "FDTD/operator.h"
#include "FDTD/engine.h"
#include "FDTD/processvoltage.h"
#include "FDTD/processcurrent.h"
#include "FDTD/processfields_td.h"
#include <sys/time.h>
#include <time.h>

//external libs
#include "tinyxml.h"
#include "ContinuousStructure.h"

double CalcDiffTime(timeval t1, timeval t2)
{
	double s_diff = difftime(t1.tv_sec,t2.tv_sec);
	s_diff += ((double)t1.tv_usec-(double)t2.tv_usec)*1e-6;
	return s_diff;
}

openEMS::openEMS()
{
	FDTD_Op=NULL;
	FDTD_Eng=NULL;
	PA=NULL;
	Enable_Dumps = true;
	DebugMat = false;
	DebugOp = false;
	endCrit = 1e-6;
}

openEMS::~openEMS()
{
	delete FDTD_Op;
	FDTD_Op=NULL;
	delete FDTD_Eng;
	FDTD_Eng=NULL;
	delete PA;
	PA=NULL;
}

void openEMS::Reset()
{
	delete FDTD_Op;
	FDTD_Op=NULL;
	delete FDTD_Eng;
	FDTD_Eng=NULL;
	if (PA) PA->DeleteAll();
	delete PA;
	PA=NULL;
}

int openEMS::SetupFDTD(const char* file)
{
	if (file==NULL) return -1;
	Reset();
	double f0=0;
	double fc=0;
	int Excit_Type=0;
	int bounds[6];

	time_t startTime=time(NULL);

	TiXmlDocument doc(file);
	if (!doc.LoadFile())
	{
		cerr << "openEMS: Error File-Loading failed!!! File: " << file << endl;
		exit(-1);
	}

	cout << "Read openEMS Settings..." << endl;
	TiXmlElement* openEMSxml = doc.FirstChildElement("openEMS");
	if (openEMSxml==NULL)
	{
		cerr << "Can't read openEMS ... " << endl;
		exit(-1);
	}

	TiXmlElement* FDTD_Opts = openEMSxml->FirstChildElement("FDTD");
	if (FDTD_Opts==NULL)
	{
		cerr << "Can't read openEMS FDTD Settings... " << endl;
		exit(-1);
	}
	FDTD_Opts->QueryIntAttribute("NumberOfTimesteps",&NrTS);
	FDTD_Opts->QueryDoubleAttribute("endCriteria",&endCrit);
	if (endCrit==0)
		endCrit=1e-6;

	TiXmlElement* Excite = FDTD_Opts->FirstChildElement("Excitation");
	if (Excite==NULL)
	{
		cerr << "Can't read openEMS Excitation Settings... " << endl;
		exit(-2);
	}
	Excite->QueryIntAttribute("Type",&Excit_Type);
	if (Excit_Type==0)
	{
		Excite->QueryDoubleAttribute("f0",&f0);
		Excite->QueryDoubleAttribute("fc",&fc);
	}
	else if (Excit_Type==1)
	{
		Excite->QueryDoubleAttribute("f0",&f0);
	}

	TiXmlElement* BC = FDTD_Opts->FirstChildElement("BoundaryCond");
	if (BC==NULL)
	{
		cerr << "Can't read openEMS boundary cond Settings... " << endl;
		exit(-3);
	}
	BC->QueryIntAttribute("xmin",&bounds[0]);
	BC->QueryIntAttribute("xmax",&bounds[1]);
	BC->QueryIntAttribute("ymin",&bounds[2]);
	BC->QueryIntAttribute("ymax",&bounds[3]);
	BC->QueryIntAttribute("zmin",&bounds[4]);
	BC->QueryIntAttribute("zmax",&bounds[5]);

	cout << "Read Geometry..." << endl;
	ContinuousStructure CSX;
	string EC(CSX.ReadFromXML(openEMSxml));
	if (EC.empty()==false)
	{
		cerr << EC << endl;
		return(-2);
	}

	bool PMC[6];
	for (int n=0;n<6;++n)
		PMC[n]=(bounds[n]==1);

	//*************** setup operator ************//
	cout << "Create Operator..." << endl;
	FDTD_Op = new Operator();
	if (FDTD_Op->SetGeometryCSX(&CSX)==false) return(-1);

	if (DebugMat)
	{
		FDTD_Op->DumpMaterial2File("material_dump.vtk");
	}
	FDTD_Op->CalcECOperator();
	if (DebugOp)
	{
		FDTD_Op->DumpOperator2File("operator_dump.vtk");
	}

	fc = 0;
	if (Excit_Type==0)
		FDTD_Op->CalcGaussianPulsExcitation(f0,fc);
	else if (Excit_Type==1)
		FDTD_Op->CalcSinusExcitation(f0,NrTS);
	else
	{
		cerr << "openEMS: Excitation type is unknown" << endl;
		exit(-1);
	}

	time_t OpDoneTime=time(NULL);

	FDTD_Op->ShowSize();

	FDTD_Op->ApplyMagneticBC(PMC);

	cout << "Nyquist number of timesteps: " << FDTD_Op->GetNyquistNum(f0+fc) << endl;
	unsigned int Nyquist = FDTD_Op->GetNyquistNum(f0+fc);

	cout << "Creation time for operator: " << difftime(OpDoneTime,startTime) << " s" << endl;

	//create FDTD engine
	FDTD_Eng = new Engine(FDTD_Op);

	time_t currTime = time(NULL);

	//*************** setup processing ************//
	cout << "Setting up processing..." << endl;
	PA = new ProcessingArray(Nyquist);

	double start[3];
	double stop[3];
	vector<CSProperties*> Probes = CSX.GetPropertyByType(CSProperties::PROBEBOX);
	for (size_t i=0;i<Probes.size();++i)
	{
		//only looking for one prim atm
		CSPrimitives* prim = Probes.at(i)->GetPrimitive(0);
		if (prim!=NULL)
		{
			bool acc;
			double* bnd = prim->GetBoundBox(acc,true);
			start[0]= bnd[0];start[1]=bnd[2];start[2]=bnd[4];
			stop[0] = bnd[1];stop[1] =bnd[3];stop[2] =bnd[5];
			CSPropProbeBox* pb = Probes.at(i)->ToProbeBox();
			Processing* proc = NULL;
			if (pb)
			{
				if (pb->GetProbeType()==0)
				{
					ProcessVoltage* procVolt = new ProcessVoltage(FDTD_Op,FDTD_Eng);
					procVolt->OpenFile(pb->GetName());
					proc=procVolt;
				}
				if (pb->GetProbeType()==1)
				{
					ProcessCurrent* procCurr = new ProcessCurrent(FDTD_Op,FDTD_Eng);
					procCurr->OpenFile(pb->GetName());
					proc=procCurr;
				}
				proc->SetProcessInterval(Nyquist/3); //three times as often as nyquist dictates
				proc->DefineStartStopCoord(start,stop);
				PA->AddProcessing(proc);
			}
			else
				delete 	proc;
		}
	}

	vector<CSProperties*> DumpProps = CSX.GetPropertyByType(CSProperties::DUMPBOX);
	for (size_t i=0;i<DumpProps.size();++i)
	{
		ProcessFieldsTD* ProcTD = new ProcessFieldsTD(FDTD_Op,FDTD_Eng);
		ProcTD->SetEnable(Enable_Dumps);
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
			CSPropDumpBox* db = DumpProps.at(i)->ToDumpBox();
			if (db)
			{
				ProcTD->SetDumpType(db->GetDumpType());
				ProcTD->SetDumpMode(db->GetDumpMode());
				ProcTD->SetFilePattern(db->GetName());
				ProcTD->DefineStartStopCoord(start,stop);
				PA->AddProcessing(ProcTD);
			}
			else
				delete 	ProcTD;
		}
	}
	return 0;
}

void openEMS::RunFDTD()
{
	cout << "Running FDTD engine... this may take a while... grab a coup of coffee?!?" << endl;

	timeval currTime;
	gettimeofday(&currTime,NULL);
	timeval startTime = currTime;
	timeval prevTime= currTime;
	ProcessFields ProcField(FDTD_Op,FDTD_Eng);
	double maxE=0,currE=0;
	double change=1;
	int prevTS=0,currTS=0;
	double speed = (double)FDTD_Op->GetNumberCells()/1e6;
	double t_diff;
	//*************** simulate ************//
	int step=PA->Process();
	if ((step<0) || (step>NrTS)) step=NrTS;
	while ((FDTD_Eng->GetNumberOfTimesteps()<NrTS) && (change>endCrit))
	{
		FDTD_Eng->IterateTS(step);
		step=PA->Process();
//		cout << " do " << step << " steps; current: " << eng.GetNumberOfTimesteps() << endl;
		currTS = FDTD_Eng->GetNumberOfTimesteps();
		if ((step<0) || (step>NrTS - currTS)) step=NrTS - currTS;

		gettimeofday(&currTime,NULL);

		t_diff = CalcDiffTime(currTime,prevTime);
		if (t_diff>4)
		{
			currE = ProcField.CalcTotalEnergy();
			if (currE>maxE)
				maxE=currE;
			cout << "[@" << setw(8) << (int)CalcDiffTime(currTime,startTime)  <<  "s] Timestep: " << setw(12)  << currTS << " (" << setw(6) << setprecision(2) << std::fixed << (double)currTS/(double)NrTS*100.0  << "%)" ;
			cout << " with currently " << setw(6) << setprecision(1) << std::fixed << speed*(double)(currTS-prevTS)/t_diff << " MCells/s" ;
			if (maxE)
				change = currE/maxE;
			cout << " --- Energy: ~" << setw(6) << setprecision(2) << std::scientific << currE << " (decrement: " << setw(6)  << setprecision(2) << std::fixed << fabs(10.0*log10(change)) << "dB)" << endl;
			prevTime=currTime;
			prevTS=currTS;
		}
	}

	//*************** postproc ************//
	prevTime = currTime;
	gettimeofday(&currTime,NULL);

	t_diff = CalcDiffTime(currTime,startTime);

	cout << "Time for " << FDTD_Eng->GetNumberOfTimesteps() << " iterations with " << FDTD_Op->GetNumberCells() << " cells : " << t_diff << " sec" << endl;
	cout << "Speed: " << speed*(double)FDTD_Eng->GetNumberOfTimesteps()/t_diff << " MCells/s " << endl;
}
