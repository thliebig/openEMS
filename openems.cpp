/*
*	Copyright (C) 2010-2015 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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
#include <iostream>
#include <fstream>
#include "tools/array_ops.h"
#include "tools/useful.h"
#include "FDTD/operator_cylinder.h"
#include "FDTD/operator_cylindermultigrid.h"
#include "FDTD/engine_multithread.h"
#include "FDTD/operator_multithread.h"
#include "FDTD/extensions/operator_ext_excitation.h"
#include "FDTD/extensions/operator_ext_tfsf.h"
#include "FDTD/extensions/operator_ext_mur_abc.h"
#include "FDTD/extensions/operator_ext_upml.h"
#include "FDTD/extensions/operator_ext_lorentzmaterial.h"
#include "FDTD/extensions/operator_ext_conductingsheet.h"
#include "FDTD/extensions/operator_ext_steadystate.h"
#include "FDTD/extensions/engine_ext_steadystate.h"
#include "FDTD/engine_interface_fdtd.h"
#include "FDTD/engine_interface_cylindrical_fdtd.h"
#include "Common/processvoltage.h"
#include "Common/processcurrent.h"
#include "Common/processfieldprobe.h"
#include "Common/processmodematch.h"
#include "Common/processfields_td.h"
#include "Common/processfields_fd.h"
#include "Common/processfields_sar.h"
#include <hdf5.h>            // only for H5get_libversion()
#include <boost/version.hpp> // only for BOOST_LIB_VERSION
#include <vtkVersion.h>

//external libs
#include "tinyxml.h"
#include "ContinuousStructure.h"
#include "CSPropProbeBox.h"
#include "CSPropDumpBox.h"

using namespace std;

double CalcDiffTime(timeval t1, timeval t2)
{
	double s_diff = t1.tv_sec - t2.tv_sec;
	s_diff += (t1.tv_usec-t2.tv_usec)*1e-6;
	return s_diff;
}

openEMS::openEMS()
{
	FDTD_Op=NULL;
	FDTD_Eng=NULL;
	Eng_Ext_SSD=NULL;
	m_CSX=NULL;
	PA=NULL;
	CylinderCoords = false;
	Enable_Dumps = true;
	DebugMat = false;
	DebugOp = false;
	m_debugCSX = false;
	m_debugBox = m_debugPEC = m_no_simulation = false;
	m_DumpStats = false;
	endCrit = 1e-6;
	m_OverSampling = 4;
	m_CellConstantMaterial=false;

	m_engine = EngineType_Multithreaded; //default engine type
	m_engine_numThreads = 0;

	m_Abort = false;
	m_Exc = 0;
}

openEMS::~openEMS()
{
	Reset();
}

void openEMS::Reset()
{
	if (PA) PA->DeleteAll();
	delete PA;
	PA=0;
	delete FDTD_Eng;
	FDTD_Eng=0;
	delete FDTD_Op;
	FDTD_Op=0;
	delete m_CSX;
	m_CSX=0;
	delete m_Exc;
	m_Exc=0;
}

//! \brief processes a command line argument
//! \return true if argument is known
//! \return false if argument is unknown
bool openEMS::parseCommandLineArgument( const char *argv )
{
	if (!argv)
		return false;

	if (strcmp(argv,"--disable-dumps")==0)
	{
		cout << "openEMS - disabling all field dumps" << endl;
		SetEnableDumps(false);
		return true;
	}
	else if (strcmp(argv,"--debug-material")==0)
	{
		cout << "openEMS - dumping material to 'material_dump.vtk'" << endl;
		DebugMaterial();
		return true;
	}
	else if (strcmp(argv,"--debug-operator")==0)
	{
		cout << "openEMS - dumping operator to 'operator_dump.vtk'" << endl;
		DebugOperator();
		return true;
	}
	else if (strcmp(argv,"--debug-boxes")==0)
	{
		cout << "openEMS - dumping boxes to 'box_dump*.vtk'" << endl;
		DebugBox();
		return true;
	}
	else if (strcmp(argv,"--debug-PEC")==0)
	{
		cout << "openEMS - dumping PEC info to 'PEC_dump.vtk'" << endl;
		m_debugPEC = true;
		return true;
	}
	else if (strcmp(argv,"--debug-CSX")==0)
	{
		cout << "openEMS - dumping CSX geometry to 'debugCSX.xml'" << endl;
		m_debugCSX = true;
		return true;
	}
	else if (strcmp(argv,"--engine=basic")==0)
	{
		cout << "openEMS - enabled basic engine" << endl;
		m_engine = EngineType_Basic;
		return true;
	}
	else if (strcmp(argv,"--engine=sse")==0)
	{
		cout << "openEMS - enabled sse engine" << endl;
		m_engine = EngineType_SSE;
		return true;
	}
	else if (strcmp(argv,"--engine=sse-compressed")==0)
	{
		cout << "openEMS - enabled compressed sse engine" << endl;
		m_engine = EngineType_SSE_Compressed;
		return true;
	}
	else if (strcmp(argv,"--engine=multithreaded")==0)
	{
		cout << "openEMS - enabled multithreading" << endl;
		m_engine = EngineType_Multithreaded;
		return true;
	}
	else if (strncmp(argv,"--numThreads=",13)==0)
	{
		m_engine_numThreads = atoi(argv+13);
		cout << "openEMS - fixed number of threads: " << m_engine_numThreads << endl;
		return true;
	}
	else if (strcmp(argv,"--engine=fastest")==0)
	{
		cout << "openEMS - enabled multithreading engine" << endl;
		m_engine = EngineType_Multithreaded;
		return true;
	}
	else if (strcmp(argv,"--no-simulation")==0)
	{
		cout << "openEMS - disabling simulation => preprocessing only" << endl;
		m_no_simulation = true;
		return true;
	}
	else if (strcmp(argv,"--dump-statistics")==0)
	{
		cout << "openEMS - dump simulation statistics to '" << __OPENEMS_RUN_STAT_FILE__ << "' and '" << __OPENEMS_STAT_FILE__ << "'" << endl;
		m_DumpStats = true;
		return true;
	}

	return false;
}

string openEMS::GetExtLibsInfo()
{
	stringstream str;

	str << "\tUsed external libraries:" << endl;
	str << "\t\t" << ContinuousStructure::GetInfoLine(true) << endl;

	// libhdf5
	unsigned int major, minor, release;
	if (H5get_libversion( &major, &minor, &release ) >= 0)
	{
		str << "\t\t" << "hdf5   -- Version: " << major << '.' << minor << '.' << release << endl;
		str << "\t\t" << "          compiled against: " H5_VERS_INFO << endl;
	}

	// tinyxml
	str << "\t\t" << "tinyxml -- compiled against: " << TIXML_MAJOR_VERSION << '.' << TIXML_MINOR_VERSION << '.' << TIXML_PATCH_VERSION << endl;

	// fparser
	str << "\t\t" << "fparser" << endl;

	// boost
	str << "\t\t" << "boost  -- compiled against: " << BOOST_LIB_VERSION << endl;

	//vtk
	str << "\t\t" << "vtk -- Version: " << vtkVersion::GetVTKMajorVersion() << "." << vtkVersion::GetVTKMinorVersion() << "." << vtkVersion::GetVTKBuildVersion() << endl;
	str << "\t\t" << "       compiled against: " << VTK_VERSION << endl;

	return str.str();
}

bool openEMS::SetupBoundaryConditions(TiXmlElement* BC)
{
	int EC; //error code of tinyxml
	int bounds[6] = {0,0,0,0,0,0};   //default boundary cond. (PEC)
	unsigned int pml_size[6] = {8,8,8,8,8,8}; //default pml size
	string s_bc;
	const char* tmp = BC->Attribute("PML_Grading");
	string pml_gradFunc;
	if (tmp)
		pml_gradFunc = string(tmp);

	string bound_names[] = {"xmin","xmax","ymin","ymax","zmin","zmax"};

	for (int n=0; n<6; ++n)
	{
		EC = BC->QueryIntAttribute(bound_names[n].c_str(),&bounds[n]);
		if (EC==TIXML_SUCCESS)
			continue;
		if (EC==TIXML_WRONG_TYPE)
		{
			tmp = BC->Attribute(bound_names[n].c_str());
			if (tmp)
				s_bc = string(tmp);
			if (s_bc=="PEC")
				bounds[n] = 0;
			else if (s_bc=="PMC")
				bounds[n] = 1;
			else if (s_bc=="MUR")
				bounds[n] = 2;
			else if (strncmp(s_bc.c_str(),"PML_=",4)==0)
			{
				bounds[n] = 3;
				pml_size[n] = atoi(s_bc.c_str()+4);
			}
			else
				cerr << "openEMS::SetupBoundaryConditions: Warning,  boundary condition for \"" << bound_names[n] << "\" unknown... set to PEC " << endl;
		}
		else
			cerr << "openEMS::SetupBoundaryConditions: Warning, boundary condition for \"" << bound_names[n] << "\" not found... set to PEC " << endl;
	}

	FDTD_Op->SetBoundaryCondition(bounds); //operator only knows about PEC and PMC, everything else is defined by extensions (see below)

	/**************************** create all operator/engine extensions here !!!! **********************************/
	//Mur-ABC, defined as extension to the operator
	double mur_v_ph = 0;
	//read general mur phase velocity
	if (BC->QueryDoubleAttribute("MUR_PhaseVelocity",&mur_v_ph) != TIXML_SUCCESS)
		mur_v_ph = -1;
	string mur_v_ph_names[6] = {"MUR_PhaseVelocity_xmin", "MUR_PhaseVelocity_xmax", "MUR_PhaseVelocity_ymin", "MUR_PhaseVelocity_ymax", "MUR_PhaseVelocity_zmin", "MUR_PhaseVelocity_zmax"};
	for (int n=0; n<6; ++n)
	{
		FDTD_Op->SetBCSize(n, 0);
		if (bounds[n]==2) //Mur-ABC
		{
			FDTD_Op->SetBCSize(n, 1);
			Operator_Ext_Mur_ABC* op_ext_mur = new Operator_Ext_Mur_ABC(FDTD_Op);
			op_ext_mur->SetDirection(n/2,n%2);
			double v_ph = 0;
			//read special mur phase velocity or assign general phase velocity
			if (BC->QueryDoubleAttribute(mur_v_ph_names[n].c_str(),&v_ph) == TIXML_SUCCESS)
				op_ext_mur->SetPhaseVelocity(v_ph);
			else if (mur_v_ph>0)
				op_ext_mur->SetPhaseVelocity(mur_v_ph);
			FDTD_Op->AddExtension(op_ext_mur);
		}
		if (bounds[n]==3)
			FDTD_Op->SetBCSize(n, pml_size[n]);
	}

	//create the upml
	Operator_Ext_UPML::Create_UPML(FDTD_Op,bounds,pml_size,pml_gradFunc);

	return true;
}

Engine_Interface_FDTD* openEMS::NewEngineInterface(int multigridlevel)
{
	Operator_CylinderMultiGrid* op_cyl_mg = dynamic_cast<Operator_CylinderMultiGrid*>(FDTD_Op);
	while (op_cyl_mg && multigridlevel>0)
	{
		int mgl = op_cyl_mg->GetMultiGridLevel();
		if (mgl==multigridlevel)
		{
			if (g_settings.GetVerboseLevel()>0)
				cout << __func__ << ": Operator with requested multi-grid level found." << endl;
			return new Engine_Interface_Cylindrical_FDTD(op_cyl_mg);
		}
		Operator_Cylinder* op_cyl_inner = op_cyl_mg->GetInnerOperator();
		op_cyl_mg = dynamic_cast<Operator_CylinderMultiGrid*>(op_cyl_inner);
		if (op_cyl_mg==NULL) //inner most operator reached
		{
			if (g_settings.GetVerboseLevel()>0)
				cout << __func__ << ": Operator with highest multi-grid level chosen." << endl;
			return new Engine_Interface_Cylindrical_FDTD(op_cyl_inner);
		}
		// try next level
	}
	Operator_Cylinder* op_cyl = dynamic_cast<Operator_Cylinder*>(FDTD_Op);
	if (op_cyl)
		return new Engine_Interface_Cylindrical_FDTD(op_cyl);
	Operator_sse* op_sse = dynamic_cast<Operator_sse*>(FDTD_Op);
	if (op_sse)
		return new Engine_Interface_SSE_FDTD(op_sse);
	return new Engine_Interface_FDTD(FDTD_Op);
}

bool openEMS::SetupProcessing()
{
	//*************** setup processing ************//
	if (g_settings.GetVerboseLevel()>0)
		cout << "Setting up processing..." << endl;

	unsigned int Nyquist = FDTD_Op->GetExcitationSignal()->GetNyquistNum();
	PA = new ProcessingArray(Nyquist);

	double start[3];
	double stop[3];
	bool l_MultiBox = false;
	vector<CSProperties*> Probes = m_CSX->GetPropertyByType(CSProperties::PROBEBOX);
	for (size_t i=0; i<Probes.size(); ++i)
	{
		//check whether one or more probe boxes are defined
		l_MultiBox =  (Probes.at(i)->GetQtyPrimitives()>1);

		for (size_t nb=0; nb<Probes.at(i)->GetQtyPrimitives(); ++nb)
		{
			CSPrimitives* prim = Probes.at(i)->GetPrimitive(nb);
			if (prim!=NULL)
			{
				double bnd[6] = {0,0,0,0,0,0};
				prim->GetBoundBox(bnd,true);
				start[0]= bnd[0];
				start[1]=bnd[2];
				start[2]=bnd[4];
				stop[0] = bnd[1];
				stop[1] =bnd[3];
				stop[2] =bnd[5];
				CSPropProbeBox* pb = Probes.at(i)->ToProbeBox();
				ProcessIntegral* proc = NULL;
				if (pb)
				{
					if (pb->GetProbeType()==0)
					{
						ProcessVoltage* procVolt = new ProcessVoltage(NewEngineInterface());
						proc=procVolt;
					}
					else if (pb->GetProbeType()==1)
					{
						ProcessCurrent* procCurr = new ProcessCurrent(NewEngineInterface());
						proc=procCurr;
					}
					else if (pb->GetProbeType()==2)
						proc = new ProcessFieldProbe(NewEngineInterface(),0);
					else if (pb->GetProbeType()==3)
						proc = new ProcessFieldProbe(NewEngineInterface(),1);
					else if ((pb->GetProbeType()==10) || (pb->GetProbeType()==11))
					{
						ProcessModeMatch* pmm = new ProcessModeMatch(NewEngineInterface());
						pmm->SetFieldType(pb->GetProbeType()-10);
						pmm->SetModeFunction(0,pb->GetAttributeValue("ModeFunctionX"));
						pmm->SetModeFunction(1,pb->GetAttributeValue("ModeFunctionY"));
						pmm->SetModeFunction(2,pb->GetAttributeValue("ModeFunctionZ"));
						proc = pmm;
					}
					else
					{
						cerr << "openEMS::SetupFDTD: Warning: Probe type " << pb->GetProbeType() << " of property '" << pb->GetName() << "' is unknown..." << endl;
						continue;
					}
					if (CylinderCoords)
						proc->SetMeshType(Processing::CYLINDRICAL_MESH);
					if ((pb->GetProbeType()==1) || (pb->GetProbeType()==3) || (pb->GetProbeType()==11))
					{
						proc->SetDualTime(true);
						proc->SetDualMesh(true);
					}
					proc->SetProcessInterval(Nyquist/m_OverSampling);
					if (pb->GetStartTime()>0 || pb->GetStopTime()>0)
						proc->SetProcessStartStopTime(pb->GetStartTime(), pb->GetStopTime());
					proc->AddFrequency(pb->GetFDSamples());
					proc->GetNormalDir(pb->GetNormalDir());
					if (l_MultiBox==false)
						proc->SetName(pb->GetName());
					else
						proc->SetName(pb->GetName(),nb);
					proc->DefineStartStopCoord(start,stop);
					if (g_settings.showProbeDiscretization())
						proc->ShowSnappedCoords();
					proc->SetWeight(pb->GetWeighting());
					PA->AddProcessing(proc);
					prim->SetPrimitiveUsed(true);
				}
				else
					delete 	proc;
			}
		}
	}

	vector<CSProperties*> DumpProps = m_CSX->GetPropertyByType(CSProperties::DUMPBOX);
	for (size_t i=0; i<DumpProps.size(); ++i)
	{
		ProcessFields* ProcField=NULL;

		//check whether one or more probe boxes are defined
		l_MultiBox =  (DumpProps.at(i)->GetQtyPrimitives()>1);

		for (size_t nb=0; nb<DumpProps.at(i)->GetQtyPrimitives(); ++nb)
		{

			CSPrimitives* prim = DumpProps.at(i)->GetPrimitive(nb);
			if (prim!=NULL)
			{
				double bnd[6] = {0,0,0,0,0,0};
				prim->GetBoundBox(bnd,true);
				start[0]= bnd[0];
				start[1]=bnd[2];
				start[2]=bnd[4];
				stop[0] = bnd[1];
				stop[1] =bnd[3];
				stop[2] =bnd[5];
				CSPropDumpBox* db = DumpProps.at(i)->ToDumpBox();
				if (db)
				{
					if ((db->GetDumpType()>=0) && (db->GetDumpType()<=3))
						ProcField = new ProcessFieldsTD(NewEngineInterface(db->GetMultiGridLevel()));
					else if ((db->GetDumpType()>=10) && (db->GetDumpType()<=13))
						ProcField = new ProcessFieldsFD(NewEngineInterface(db->GetMultiGridLevel()));
					else if ( ((db->GetDumpType()>=20) && (db->GetDumpType()<=22)) || (db->GetDumpType()==29) )
					{
						ProcessFieldsSAR* procSAR = new ProcessFieldsSAR(NewEngineInterface(db->GetMultiGridLevel()));
						ProcField = procSAR;
						string method = db->GetAttributeValue("SAR_Method");
						if (!method.empty())
							procSAR->SetSARAveragingMethod(method);
						// use (center)-cell based conductivity only
						procSAR->SetUseCellConductivity(true);
					}
					else
						cerr << "openEMS::SetupFDTD: unknown dump box type... skipping!" << endl;
					if (ProcField)
					{
						ProcField->SetEnable(Enable_Dumps);
						ProcField->SetProcessInterval(Nyquist/m_OverSampling);
						if (db->GetStopTime()>0 || db->GetStartTime()>0)
							ProcField->SetProcessStartStopTime(db->GetStartTime(), db->GetStopTime());
						if ((db->GetDumpType()==1) || (db->GetDumpType()==11))
						{
							ProcField->SetDualTime(true);
							//make dualMesh the default mesh for h-field dumps, maybe overwritten by interpolation type (node-interpolation)
							ProcField->SetDualMesh(true);
						}
						if (db->GetDumpType()>=10)
						{
							ProcField->AddFrequency(db->GetFDSamples());
							ProcField->SetDumpType((ProcessFields::DumpType)(db->GetDumpType()-10));
						}
						else
							ProcField->SetDumpType((ProcessFields::DumpType)db->GetDumpType());

						if (db->GetDumpType()==20)
							ProcField->SetDumpType(ProcessFields::SAR_LOCAL_DUMP);
						if (db->GetDumpType()==21)
							ProcField->SetDumpType(ProcessFields::SAR_1G_DUMP);
						if (db->GetDumpType()==22)
							ProcField->SetDumpType(ProcessFields::SAR_10G_DUMP);
						if (db->GetDumpType()==29)
							ProcField->SetDumpType(ProcessFields::SAR_RAW_DATA);

						//SetupMaterialStorages() has previewed storage needs... refresh here to prevent cleanup!!!
						if ( ProcField->NeedConductivity() && Enable_Dumps )
							FDTD_Op->SetMaterialStoreFlags(1,true);

						ProcField->SetDumpMode((Engine_Interface_Base::InterpolationType)db->GetDumpMode());
						ProcField->SetFileType((ProcessFields::FileType)db->GetFileType());
						if (CylinderCoords)
							ProcField->SetMeshType(Processing::CYLINDRICAL_MESH);
						if (db->GetSubSampling())
							for (int n=0; n<3; ++n)
								ProcField->SetSubSampling(db->GetSubSampling(n),n);
						if (db->GetOptResolution())
							for (int n=0; n<3; ++n)
								ProcField->SetOptResolution(db->GetOptResolution(n),n);

						if (l_MultiBox==false)
							ProcField->SetName(db->GetName());
						else
							ProcField->SetName(db->GetName(),nb);

						ProcField->SetFileName(ProcField->GetName());
						ProcField->DefineStartStopCoord(start,stop);
						if (g_settings.showProbeDiscretization())
							ProcField->ShowSnappedCoords();
						PA->AddProcessing(ProcField);
						prim->SetPrimitiveUsed(true);
					}
				}
			}
		}
	}

	return true;
}

bool openEMS::SetupMaterialStorages()
{
	vector<CSProperties*> DumpProps = m_CSX->GetPropertyByType(CSProperties::DUMPBOX);
	for (size_t i=0; i<DumpProps.size(); ++i)
	{
		CSPropDumpBox* db = DumpProps.at(i)->ToDumpBox();
		if (!db)
			continue;
		if (db->GetQtyPrimitives()==0)
			continue;
		//check for current density dump types
		if ( ((db->GetDumpType()==2) || (db->GetDumpType()==12) || // current density storage
			  (db->GetDumpType()==20) || (db->GetDumpType()==21) || (db->GetDumpType()==22)) && // SAR dump types
			  Enable_Dumps )
			FDTD_Op->SetMaterialStoreFlags(1,true); //tell operator to store kappa material data
	}
	return true;
}

bool openEMS::SetupOperator(TiXmlElement* FDTD_Opts)
{
	if (CylinderCoords)
	{
		const char* radii = FDTD_Opts->Attribute("MultiGrid");
		if (radii)
		{
			string rad(radii);
			FDTD_Op = Operator_CylinderMultiGrid::New(SplitString2Double(rad,','),m_engine_numThreads);
			if (FDTD_Op==NULL)
				FDTD_Op = Operator_Cylinder::New(m_engine_numThreads);
		}
		else
			FDTD_Op = Operator_Cylinder::New(m_engine_numThreads);
	}
	else if (m_engine == EngineType_SSE)
	{
		FDTD_Op = Operator_sse::New();
	}
	else if (m_engine == EngineType_SSE_Compressed)
	{
		FDTD_Op = Operator_SSE_Compressed::New();
	}
	else if (m_engine == EngineType_Multithreaded)
	{
		FDTD_Op = Operator_Multithread::New(m_engine_numThreads);
	}
	else
	{
		FDTD_Op = Operator::New();
	}
	return true;
}


int openEMS::SetupFDTD(const char* file)
{
	if (file==NULL) return -1;
	Reset();

	if (g_settings.GetVerboseLevel()>0)
		cout << "Read openEMS xml file: " << file << " ..." << endl;

	timeval startTime;
	gettimeofday(&startTime,NULL);


	TiXmlDocument doc(file);
	if (!doc.LoadFile())
	{
		cerr << "openEMS: Error File-Loading failed!!! File: " << file << endl;
		exit(-1);
	}

	if (g_settings.GetVerboseLevel()>0)
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
	double dhelp=0;
	FDTD_Opts->QueryDoubleAttribute("NumberOfTimesteps",&dhelp);
	if (dhelp<0)
		NrTS=0;
	else
		NrTS = (unsigned int)dhelp;

	int ihelp = 0;
	FDTD_Opts->QueryIntAttribute("CylinderCoords",&ihelp);
	if (ihelp==1)
		CylinderCoords = true;

	FDTD_Opts->QueryDoubleAttribute("endCriteria",&endCrit);
	if (endCrit==0)
		endCrit=1e-6;

	FDTD_Opts->QueryIntAttribute("OverSampling",&m_OverSampling);
	if (m_OverSampling<2)
		m_OverSampling=2;

	TiXmlElement* BC = FDTD_Opts->FirstChildElement("BoundaryCond");
	if (BC==NULL)
	{
		cerr << "Can't read openEMS boundary cond Settings... " << endl;
		exit(-3);
	}

	if (g_settings.GetVerboseLevel()>0)
		cout << "Read Geometry..." << endl;
	m_CSX = new ContinuousStructure();
	string EC(m_CSX->ReadFromXML(openEMSxml));
	if (EC.empty()==false)
	{
		cerr << EC << endl;
//		return(-2);
	}

	if (g_settings.GetVerboseLevel()>2)
		m_CSX->ShowPropertyStatus(cerr);

	if (CylinderCoords)
		if (m_CSX->GetCoordInputType()!=CYLINDRICAL)
		{
			cerr << "openEMS::SetupFDTD: Warning: Coordinate system found in the CSX file is not a cylindrical. Forcing to cylindrical coordinate system!" << endl;
			m_CSX->SetCoordInputType(CYLINDRICAL); //tell CSX to use cylinder-coords
		}

	if (m_debugCSX)
		m_CSX->Write2XML("debugCSX.xml");

	//*************** setup operator ************//
	if (SetupOperator(FDTD_Opts)==false)
		return 2;

	// default material averaging is quarter cell averaging
	FDTD_Op->SetQuarterCellMaterialAvg();

	// check for cell constant material averaging
	if (FDTD_Opts->QueryIntAttribute("CellConstantMaterial",&ihelp)==TIXML_SUCCESS)
		m_CellConstantMaterial=(ihelp==1);
	if (m_CellConstantMaterial)
	{
		FDTD_Op->SetCellConstantMaterial();
		if (g_settings.GetVerboseLevel()>0)
			cout << "Enabling constant cell material assumption." << endl;
	}

	m_Exc = new Excitation();
	if (!m_Exc->setupExcitation(FDTD_Opts->FirstChildElement("Excitation")))
		exit(2);
	FDTD_Op->SetExcitationSignal(m_Exc);
	FDTD_Op->AddExtension(new Operator_Ext_Excitation(FDTD_Op));
	if (!CylinderCoords)
		FDTD_Op->AddExtension(new Operator_Ext_TFSF(FDTD_Op));

	if (FDTD_Op->SetGeometryCSX(m_CSX)==false) return(2);

	SetupBoundaryConditions(BC);

	int TS_method=0;
	if (FDTD_Opts->QueryIntAttribute("TimeStepMethod",&TS_method)==TIXML_SUCCESS)
		FDTD_Op->SetTimeStepMethod(TS_method);

	double timestep=0;
	FDTD_Opts->QueryDoubleAttribute("TimeStep",&timestep);
	if (timestep)
		FDTD_Op->SetTimestep(timestep);
	double timestepfactor=0;
	if (FDTD_Opts->QueryDoubleAttribute("TimeStepFactor",&timestepfactor)==TIXML_SUCCESS)
		FDTD_Op->SetTimestepFactor(timestepfactor);

	// Is a steady state detection requested
	Operator_Ext_SteadyState* Op_Ext_SSD = NULL;
	if (m_Exc->GetSignalPeriod()>0)
	{
		cout << "Create a steady state detection using a period of " << m_Exc->GetSignalPeriod() << " s" << endl;
		Op_Ext_SSD = new Operator_Ext_SteadyState(FDTD_Op, m_Exc->GetSignalPeriod());
		unsigned int pos[3];
		for (int p=0;p<3;++p)
			pos[p] = FDTD_Op->GetNumberOfLines(p)/2;
		Op_Ext_SSD->Add_E_Probe(pos, 0);
		Op_Ext_SSD->Add_E_Probe(pos, 1);
		Op_Ext_SSD->Add_E_Probe(pos, 2);

		for (int n=0;n<3;++n)
		{
			for (int p=0;p<3;++p)
				pos[p] = FDTD_Op->GetNumberOfLines(p)/2;
			pos[n] = 10;
			Op_Ext_SSD->Add_E_Probe(pos, 0);
			Op_Ext_SSD->Add_E_Probe(pos, 1);
			Op_Ext_SSD->Add_E_Probe(pos, 2);

			pos[n] = FDTD_Op->GetNumberOfLines(n)-10;
			Op_Ext_SSD->Add_E_Probe(pos, 0);
			Op_Ext_SSD->Add_E_Probe(pos, 1);
			Op_Ext_SSD->Add_E_Probe(pos, 2);
		}
		FDTD_Op->AddExtension(Op_Ext_SSD);
	}

	if ((m_CSX->GetQtyPropertyType(CSProperties::LORENTZMATERIAL)>0) || (m_CSX->GetQtyPropertyType(CSProperties::DEBYEMATERIAL)>0))
		FDTD_Op->AddExtension(new Operator_Ext_LorentzMaterial(FDTD_Op));
	if (m_CSX->GetQtyPropertyType(CSProperties::CONDUCTINGSHEET)>0)
	{
		double f_max=0;
		FDTD_Opts->QueryDoubleAttribute("f_max",&f_max);
		if (f_max>0)
			FDTD_Op->AddExtension(new Operator_Ext_ConductingSheet(FDTD_Op,f_max));
		else
			cerr << __func__ << ": Error, max. frequency is <=0, disabling conducting sheet material..." << endl;
	}

	//check all properties to request material storage during operator creation...
	SetupMaterialStorages();

	/*******************   create the EC-FDTD operator *****************************/
	Operator::DebugFlags debugFlags = Operator::None;
	if (DebugMat)
		debugFlags |= Operator::debugMaterial;
	if (DebugOp)
		debugFlags |= Operator::debugOperator;
	if (m_debugPEC)
		debugFlags |= Operator::debugPEC;

	FDTD_Op->CalcECOperator( debugFlags );
	/*******************************************************************************/

	//reset flags for material storage, if no dump-box resets it to true, it will be cleaned up...
	FDTD_Op->SetMaterialStoreFlags(0,false);
	FDTD_Op->SetMaterialStoreFlags(1,false);
	FDTD_Op->SetMaterialStoreFlags(2,false);
	FDTD_Op->SetMaterialStoreFlags(3,false);

	double maxTime=0;
	FDTD_Opts->QueryDoubleAttribute("MaxTime",&maxTime);
	unsigned int maxTime_TS = (unsigned int)(maxTime/FDTD_Op->GetTimestep());
	if ((maxTime_TS>0) && (maxTime_TS<NrTS))
		NrTS = maxTime_TS;

	if (!m_Exc->buildExcitationSignal(NrTS))
		exit(2);
	m_Exc->DumpVoltageExcite("et");
	m_Exc->DumpCurrentExcite("ht");

	timeval OpDoneTime;
	gettimeofday(&OpDoneTime,NULL);

	if (g_settings.GetVerboseLevel()>0)
	{
		FDTD_Op->ShowStat();
		FDTD_Op->ShowExtStat();
		cout << "Creation time for operator: " << CalcDiffTime(OpDoneTime,startTime) << " s" << endl;
	}
	cout << "FDTD simulation size: " << FDTD_Op->GetNumberOfLines(0) << "x" << FDTD_Op->GetNumberOfLines(1) << "x" << FDTD_Op->GetNumberOfLines(2) << " --> "  << FDTD_Op->GetNumberCells() << " FDTD cells " << endl;
	cout << "FDTD timestep is: " <<FDTD_Op->GetTimestep()  << " s; Nyquist rate: " <<  m_Exc->GetNyquistNum() << " timesteps @" << CalcNyquistFrequency(m_Exc->GetNyquistNum(),FDTD_Op->GetTimestep()) << " Hz" << endl;
	if (m_Exc->GetNyquistNum()>1000)
		cerr << "openEMS::SetupFDTD: Warning, the timestep seems to be very small --> long simulation. Check your mesh!?" << endl;

	cout << "Excitation signal length is: " <<  m_Exc->GetLength() << " timesteps (" <<  m_Exc->GetLength()*FDTD_Op->GetTimestep() << "s)" << endl;
	cout << "Max. number of timesteps: " << NrTS << " ( --> " << (double)NrTS/(double)(m_Exc->GetLength()) << " * Excitation signal length)" << endl;
	if ( ((double)NrTS/(double)m_Exc->GetLength() < 3) && (m_Exc->GetExciteType()==0))
		cerr << "openEMS::SetupFDTD: Warning, max. number of timesteps is smaller than three times the excitation. " << endl << \
				"\tYou may want to choose a higher number of max. timesteps... " << endl;

	if (m_no_simulation)
	{
		// simulation was disabled (to generate debug output only)
		return 1;
	}

	//create FDTD engine
	FDTD_Eng = FDTD_Op->CreateEngine();

	if (Op_Ext_SSD)
	{
		Eng_Ext_SSD = dynamic_cast<Engine_Ext_SteadyState*>(Op_Ext_SSD->GetEngineExtention());
		Eng_Ext_SSD->SetEngineInterface(this->NewEngineInterface());
	}

	//setup all processing classes
	if (SetupProcessing()==false)
		return 2;

	// Cleanup all unused material storages...
	FDTD_Op->CleanupMaterialStorage();

	//check and warn for unused properties and primitives
	m_CSX->WarnUnusedPrimitves(cerr);

	// dump all boxes (voltage, current, fields, ...)
	if (m_debugBox)
	{
		PA->DumpBoxes2File("box_dump_");
	}

	return 0;
}

string FormatTime(int sec)
{
	stringstream ss;
	if (sec<60)
	{
		ss << setw(9) << sec << "s";
		return ss.str();
	}
	if (sec<3600)
	{
		ss << setw(6) << sec/60 << "m" << setw(2) << setfill('0') << sec%60 << "s";
		return ss.str();
	}
	ss << setw(3) << sec/3600 << "h" << setw(2) << setfill('0')  << (sec%3600)/60 << "m" << setw(2) << setfill('0')  << sec%60 << "s";
	return ss.str();
}

bool openEMS::CheckAbortCond()
{
	if (m_Abort) //abort was set externally
		return true;

	//check whether the file "ABORT" exist in current working directory
	ifstream ifile("ABORT");
	if (ifile)
	{
		ifile.close();
		cerr << "openEMS::CheckAbortCond(): Found file \"ABORT\", aborting simulation..." << endl;
		return true;
	}

	return false;
}

void openEMS::RunFDTD()
{
	cout << "Running FDTD engine... this may take a while... grab a cup of coffee?!?" << endl;

	//special handling of a field processing, needed to realize the end criteria...
	ProcessFields* ProcField = new ProcessFields(NewEngineInterface());
	PA->AddProcessing(ProcField);
	double maxE=0,currE=0;

	//init processings
	PA->InitAll();

	//add all timesteps to end-crit field processing with max excite amplitude
	unsigned int maxExcite = FDTD_Op->GetExcitationSignal()->GetMaxExcitationTimestep();
//	for (unsigned int n=0; n<FDTD_Op->Exc->Volt_Count; ++n)
//		ProcField->AddStep(FDTD_Op->Exc->Volt_delay[n]+maxExcite);
	ProcField->AddStep(maxExcite);

	double change=1;
	int prevTS=0,currTS=0;
	double numCells = FDTD_Op->GetNumberCells();
	double speed = 0;
	double t_diff;
	double t_run;

	timeval currTime;
	gettimeofday(&currTime,NULL);
	timeval startTime = currTime;
	timeval prevTime= currTime;

	if (m_DumpStats)
		InitRunStatistics(__OPENEMS_RUN_STAT_FILE__);
	//*************** simulate ************//

	PA->PreProcess();
	int step=PA->Process();
	if ((step<0) || (step>(int)NrTS)) step=NrTS;
	while ((FDTD_Eng->GetNumberOfTimesteps()<NrTS) && (change>endCrit) && !CheckAbortCond())
	{
		FDTD_Eng->IterateTS(step);
		step=PA->Process();

		if ((Eng_Ext_SSD==NULL) && ProcField->CheckTimestep())
		{
			currE = ProcField->CalcTotalEnergyEstimate();
			if (currE>maxE)
				maxE=currE;
		}

//		cout << " do " << step << " steps; current: " << eng.GetNumberOfTimesteps() << endl;
		currTS = FDTD_Eng->GetNumberOfTimesteps();
		if ((step<0) || (step>(int)(NrTS - currTS))) step=NrTS - currTS;

		gettimeofday(&currTime,NULL);

		t_diff = CalcDiffTime(currTime,prevTime);

		if (t_diff>4)
		{
			t_run = CalcDiffTime(currTime,startTime);
			speed = numCells*(currTS-prevTS)/t_diff;
			cout << "[@" <<  FormatTime(t_run) <<  "] Timestep: " << setw(12)  << currTS ;
			cout << " || Speed: " << setw(6) << setprecision(1) << std::fixed << speed*1e-6 << " MC/s (" <<  setw(4) << setprecision(3) << std::scientific << t_diff/(currTS-prevTS) << " s/TS)" ;
			if (Eng_Ext_SSD==NULL)
			{
				currE = ProcField->CalcTotalEnergyEstimate();
				if (currE>maxE)
					maxE=currE;
				if (maxE)
					change = currE/maxE;
				cout << " || Energy: ~" << setw(6) << setprecision(2) << std::scientific << currE << " (-" << setw(5)  << setprecision(2) << std::fixed << fabs(10.0*log10(change)) << "dB)" << endl;
			}
			else
			{
				change = Eng_Ext_SSD->GetLastDiff();
				cout << " || SteadyState: " << setw(6) << setprecision(2) << std::fixed << 10.0*log10(change) << " dB" << endl;
			}
			prevTime=currTime;
			prevTS=currTS;

			PA->FlushNext();

			if (m_DumpStats)
				DumpRunStatistics(__OPENEMS_RUN_STAT_FILE__, t_run, currTS, speed, currE);
		}
	}
	if ((change>endCrit) && (FDTD_Op->GetExcitationSignal()->GetExciteType()==0))
		cerr << "RunFDTD: Warning: Max. number of timesteps was reached before the end-criteria of -" << fabs(10.0*log10(endCrit)) << "dB was reached... " << endl << \
				"\tYou may want to choose a higher number of max. timesteps... " << endl;

	gettimeofday(&currTime,NULL);
	t_diff = CalcDiffTime(currTime,startTime);

	cout << "Time for " << FDTD_Eng->GetNumberOfTimesteps() << " iterations with " << FDTD_Op->GetNumberCells() << " cells : " << t_diff << " sec" << endl;
	cout << "Speed: " << numCells*(double)FDTD_Eng->GetNumberOfTimesteps()/t_diff*1e-6 << " MCells/s " << endl;

	if (m_DumpStats)
		DumpStatistics(__OPENEMS_STAT_FILE__, t_diff);

	//*************** postproc ************//
	PA->PostProcess();
}

bool openEMS::DumpStatistics(const string& filename, double time)
{
	ofstream stat_file;
	stat_file.open(filename.c_str());

	if (!stat_file.is_open())
	{
		cerr << "openEMS::DumpStatistics: Error, opening file failed..." << endl;
		return false;
	}
	stat_file << std::setprecision( 16 );
	stat_file << FDTD_Op->GetNumberCells() << "\t% number of cells" << endl;
	stat_file << FDTD_Op->GetTimestep() << "\t% timestep (s)" << endl;
	stat_file << FDTD_Eng->GetNumberOfTimesteps() << "\t% number of iterations" << endl;
	stat_file << FDTD_Eng->GetNumberOfTimesteps()*FDTD_Op->GetTimestep() << "\t% total numercial time (s)" << endl;
	stat_file << time << "\t% simulation time (s)" << endl;
	stat_file << (double)FDTD_Op->GetNumberCells()*(double)FDTD_Eng->GetNumberOfTimesteps()/time << "\t% speed (cells/s)" << endl;

	stat_file.close();
	return true;
}

bool openEMS::InitRunStatistics(const string& filename)
{
	ofstream stat_file;
	stat_file.open(filename.c_str(), ios_base::out);

	if (!stat_file.is_open())
	{
		cerr << "openEMS::InitRunStatistics: Error, opening file failed..." << endl;
		return false;
	}
	stat_file << "%time\ttimestep\tspeed\tenergy" << endl;
	stat_file.close();
	return true;
}

bool openEMS::DumpRunStatistics(const string& filename, double time, unsigned int ts, double speed, double energy)
{
	ofstream stat_file;
	stat_file.open(filename.c_str(), ios_base::app);

	if (!stat_file.is_open())
	{
		cerr << "openEMS::DumpRunStatistics: Error, opening file failed..." << endl;
		return false;
	}
	stat_file << time << "\t" << ts << "\t" << speed << "\t" << energy << endl;
	stat_file.close();
	return true;
}
