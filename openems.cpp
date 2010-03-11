#include "openems.h"
#include "tools/array_ops.h"
#include "FDTD/operator.h"
#include "FDTD/engine.h"
#include "FDTD/processvoltage.h"
#include "FDTD/processcurrent.h"
#include "FDTD/processfields_td.h"

//external libs
#include "tinyxml.h"
#include "ContinuousStructure.h"

openEMS::openEMS()
{
	FDTD_Op=NULL;
	FDTD_Eng=NULL;
	PA=NULL;
	Enable_Dumps = true;
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

	cerr << "Read openEMS Settings..." << endl;
	TiXmlElement* FDTD_Opts = doc.FirstChildElement("openEMS-Parameter");
	if (FDTD_Opts==NULL)
	{
		cerr << "Can't read openEMS Settings... " << endl;
		exit(-1);
	}
	FDTD_Opts->QueryIntAttribute("NumberOfTimesteps",&NrTS);

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

	cerr << "Read Geometry..." << endl;
	ContinuousStructure CSX;
	string EC(CSX.ReadFromXML(&doc));
	if (EC.empty()==false)
	{
		cerr << EC << endl;
		return(-2);
	}

	bool PMC[6];
	for (int n=0;n<6;++n)
		PMC[n]=(bounds[n]==1);

	//*************** setup operator ************//
	cerr << "Create Operator..." << endl;
	FDTD_Op = new Operator();
	if (FDTD_Op->SetGeometryCSX(&CSX)==false) return(-1);

	FDTD_Op->CalcECOperator();

	if (Excit_Type==0)
		FDTD_Op->CalcGaussianPulsExcitation(f0,fc);
	else
	{
		cerr << "openEMS: Excitation type is unknown" << endl;
		exit(-1);
	}

	time_t OpDoneTime=time(NULL);

	FDTD_Op->ShowSize();

	FDTD_Op->ApplyMagneticBC(PMC);

	cerr << "Nyquist number of timesteps: " << FDTD_Op->GetNyquistNum(f0+fc) << endl;
	unsigned int Nyquist = FDTD_Op->GetNyquistNum(f0+fc);

	cerr << "Time for operator: " << difftime(OpDoneTime,startTime) << endl;

	//create FDTD engine
	FDTD_Eng = new Engine(FDTD_Op);

	time_t currTime = time(NULL);

	//*************** setup processing ************//
	PA = new ProcessingArray();

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
	time_t currTime = time(NULL);
	//*************** simulate ************//
	int step=PA->Process();
	if ((step<0) || (step>NrTS)) step=NrTS;
	while (FDTD_Eng->GetNumberOfTimesteps()<NrTS)
	{
		FDTD_Eng->IterateTS(step);
		step=PA->Process();
//		cerr << " do " << step << " steps; current: " << eng.GetNumberOfTimesteps() << endl;
		if ((step<0) || (step>NrTS - FDTD_Eng->GetNumberOfTimesteps())) step=NrTS - FDTD_Eng->GetNumberOfTimesteps();
	}

	//*************** postproc ************//
	time_t prevTime = currTime;
	currTime = time(NULL);

	double t_diff = difftime(currTime,prevTime);

	cerr << "Time for " << FDTD_Eng->GetNumberOfTimesteps() << " iterations with " << FDTD_Op->GetNumberCells() << " cells : " << t_diff << " sec" << endl;
	cerr << "Speed: " << (double)FDTD_Op->GetNumberCells()*(double)FDTD_Eng->GetNumberOfTimesteps()/t_diff/1e6 << " MCells/s " << endl;
}
