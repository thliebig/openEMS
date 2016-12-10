/*
*	Copyright (C) 2011 Thorsten Liebig (Thorsten.Liebig@gmx.de)
*
*	This program is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
*	(at your option) any later version.
*
*	This program is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY{} without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*
*	You should have received a copy of the GNU General Public License
*	along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "openems_fdtd_mpi.h"
#include "FDTD/engine_interface_fdtd.h"
#include "FDTD/operator_mpi.h"
#include "FDTD/operator_cylinder.h"
#include "FDTD/engine_mpi.h"
#include "Common/processfields.h"
#include "Common/processintegral.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include "mpi.h"
#include "tools/useful.h"
#include "tinyxml.h"

openEMS_FDTD_MPI::openEMS_FDTD_MPI(bool m_MPI_Debug) : openEMS()
{
	m_MyID = MPI::COMM_WORLD.Get_rank();
	m_NumProc = MPI::COMM_WORLD.Get_size();

	m_MaxEnergy = 0;
	m_EnergyDecrement = 1;
	m_MPI_Op = NULL;

	if (m_NumProc>1)
		m_MPI_Enabled=true;
	else
		m_MPI_Enabled=false;

	if (m_MyID==0)
	{
		m_Gather_Buffer = new int[m_NumProc];
		m_Energy_Buffer = new double[m_NumProc];
	}
	else
	{
		m_Gather_Buffer = NULL;
		m_Energy_Buffer = NULL;
	}

	m_MPI_Elem = NULL;
	m_Original_Grid = NULL;

	//redirect output to file for all ranks > 0
	if ((m_MyID>0) && (m_MPI_Debug==false))
	{
		stringstream out_name;
		out_name << "ID" << m_MyID << "_" << "output.txt";
		m_Output = new ofstream();
		m_Output->open(out_name.str().c_str());
		cout.rdbuf(m_Output->rdbuf());
		cerr.rdbuf(m_Output->rdbuf());
	}
	else
		m_Output = NULL;
}

openEMS_FDTD_MPI::~openEMS_FDTD_MPI()
{
	delete[] m_Gather_Buffer;
	m_Gather_Buffer = NULL;
	delete[] m_Energy_Buffer;
	m_Energy_Buffer = NULL;
	delete m_Original_Grid;
	m_Original_Grid = NULL;
	delete m_Output;
	m_Output=NULL;
}

string openEMS_FDTD_MPI::GetExtLibsInfo()
{
	stringstream str;

	str << openEMS::GetExtLibsInfo();

	// MPI
	str << "\t\t" << "MPI -- Version: " << MPI_VERSION << "." << MPI_SUBVERSION << endl;
	str << "\t\t" << "       compiled against: ";
#ifdef MPICH2_VERSION
	str << "MPICH2 " << MPICH2_VERSION << endl;
#endif

#ifdef OMPI_MAJOR_VERSION
	str << "openMPI" << OMPI_MAJOR_VERSION << "." << OMPI_MINOR_VERSION << "." << OMPI_RELEASE_VERSION << endl;
#endif

	return str.str();
}


bool openEMS_FDTD_MPI::parseCommandLineArgument( const char *argv )
{
	if (!argv)
		return false;

	bool ret = openEMS::parseCommandLineArgument( argv );

	if (ret)
		return ret;

	if (strcmp(argv,"--engine=MPI")==0)
	{
		cout << "openEMS_FDTD_MPI - enabled MPI parallel processing" << endl;
		m_engine = EngineType_MPI;
		return true;
	}

	return false;
}

bool openEMS_FDTD_MPI::Parse_XML_FDTDSetup(TiXmlElement* FDTD_Opts)
{
	m_MPI_Elem = FDTD_Opts->FirstChildElement("MPI");
	if (!m_MPI_Enabled)
	{
		if ((m_MPI_Elem!=NULL))
			cerr << "openEMS_FDTD_MPI::SetupMPI: Warning: Number of MPI processes is 1, skipping MPI engine... " << endl;
		return openEMS::Parse_XML_FDTDSetup(FDTD_Opts);
	}

	if (m_MPI_Elem==NULL)
	{
		MPI_Barrier(MPI_COMM_WORLD);
		if (m_MyID==0)
			cerr << "openEMS_FDTD_MPI::SetupMPI: Error: no MPI settings found, exiting MPI engine... " << endl;
		exit(-1);
	}

	CSRectGrid* grid = m_CSX->GetGrid();
	delete m_Original_Grid;
	m_Original_Grid = CSRectGrid::Clone(grid);

	string arg_Pos_Names[] = {"SplitPos_X", "SplitPos_Y", "SplitPos_Z"};
	string arg_N_Names[] = {"SplitN_X", "SplitN_Y", "SplitN_Z"};
	const char* tmp = NULL;
	for (int n=0;n<3;++n)
	{
		m_SplitNumber[n].clear();
		m_SplitNumber[n].push_back(0);
		tmp = m_MPI_Elem->Attribute(arg_Pos_Names[n].c_str());
		if (tmp) //check if a split position is requested
		{
			vector<double> SplitLines = SplitString2Double(tmp, ',');
			bool inside;
			unsigned int line;
			for (size_t lineN = 0; lineN<SplitLines.size();++lineN)
			{
				line = m_Original_Grid->Snap2LineNumber(n, SplitLines.at(lineN), inside);
				if (inside)
					m_SplitNumber[n].push_back(line);
			}
		}
		else //check if a number of splits is requested
		{
			int SplitN=0;
			if (m_MPI_Elem->QueryIntAttribute( arg_N_Names[n].c_str(), &SplitN) == TIXML_SUCCESS)
			{
				if (SplitN>1)
				{

					vector<unsigned int> jobs = AssignJobs2Threads(m_Original_Grid->GetQtyLines(n)-1, SplitN, true);
					unsigned int line=0;
					for (size_t i = 0; i<jobs.size()-1;++i)
					{
						line += jobs.at(i);
						m_SplitNumber[n].push_back(line);
					}
				}
			}
		}

		m_SplitNumber[n].push_back(m_Original_Grid->GetQtyLines(n)-1);
		unique(m_SplitNumber[n].begin(), m_SplitNumber[n].end());
	}

	return openEMS::Parse_XML_FDTDSetup(FDTD_Opts);
}

bool openEMS_FDTD_MPI::SetupMPI()
{
	if (!m_MPI_Enabled)
		return true;

	MPI_Barrier(MPI_COMM_WORLD);

	//validate number of processes
	unsigned int numProcs = (m_SplitNumber[0].size()-1)*(m_SplitNumber[1].size()-1)*(m_SplitNumber[2].size()-1);
	if (numProcs!=m_NumProc)
	{
		if (m_MyID==0)
			cerr << "openEMS_FDTD_MPI::SetupMPI: Error: Requested splits require " << numProcs << " processes, but " << m_NumProc << " were found! Exit! " << endl;
		exit(10);
	}

	//create process table
	unsigned int procN = 0;
	unsigned int splits[] = {(unsigned int)m_SplitNumber[0].size()-1, (unsigned int)m_SplitNumber[1].size()-1, (unsigned int)m_SplitNumber[2].size()-1};
	m_MPI_Op->SetSplitNumbers(0,splits[0]);
	m_MPI_Op->SetSplitNumbers(1,splits[1]);
	m_MPI_Op->SetSplitNumbers(2,splits[2]);
	unsigned int*** procTable=Create3DArray<unsigned int>(splits);
	for (size_t i=0;i<m_SplitNumber[0].size()-1;++i)
		for (size_t j=0;j<m_SplitNumber[1].size()-1;++j)
			for (size_t k=0;k<m_SplitNumber[2].size()-1;++k)
			{
				procTable[i][j][k] = procN;
				++procN;
			}
	m_MPI_Op->SetProcessTable(procTable);

	CSRectGrid* grid = m_CSX->GetGrid();
	//assign mesh and neighbors to this process
	for (size_t i=0;i<m_SplitNumber[0].size()-1;++i)
	{
		for (size_t j=0;j<m_SplitNumber[1].size()-1;++j)
		{
			for (size_t k=0;k<m_SplitNumber[2].size()-1;++k)
			{
				if (procTable[i][j][k] == m_MyID)
				{
					m_MPI_Op->SetProcessTablePosition(0,i);
					m_MPI_Op->SetProcessTablePosition(1,j);
					m_MPI_Op->SetProcessTablePosition(2,k);

					grid->ClearLines(0);
					grid->ClearLines(1);
					grid->ClearLines(2);

					for (unsigned int n=m_SplitNumber[0].at(i);n<=m_SplitNumber[0].at(i+1);++n)
						grid->AddDiscLine(0, m_Original_Grid->GetLine(0,n) );
					for (unsigned int n=m_SplitNumber[1].at(j);n<=m_SplitNumber[1].at(j+1);++n)
						grid->AddDiscLine(1, m_Original_Grid->GetLine(1,n) );
					for (unsigned int n=m_SplitNumber[2].at(k);n<=m_SplitNumber[2].at(k+1);++n)
						grid->AddDiscLine(2, m_Original_Grid->GetLine(2,n) );

					m_MPI_Op->SetSplitPos(0,m_SplitNumber[0].at(i));
					m_MPI_Op->SetSplitPos(1,m_SplitNumber[1].at(j));
					m_MPI_Op->SetSplitPos(2,m_SplitNumber[2].at(k));

					if (i>0)
						m_MPI_Op->SetNeighborDown(0,procTable[i-1][j][k]);
					if (i<m_SplitNumber[0].size()-2)
					{
						//add one additional line
						grid->AddDiscLine(0, m_Original_Grid->GetLine(0,m_SplitNumber[0].at(i+1)+1 ));
						m_MPI_Op->SetNeighborUp(0,procTable[i+1][j][k]);
					}

					if (j>0)
						m_MPI_Op->SetNeighborDown(1,procTable[i][j-1][k]);
					if (j<m_SplitNumber[1].size()-2)
					{
						//add one additional line
						grid->AddDiscLine(1, m_Original_Grid->GetLine(1,m_SplitNumber[1].at(j+1)+1 ));
						m_MPI_Op->SetNeighborUp(1,procTable[i][j+1][k]);
					}

					if (k>0)
						m_MPI_Op->SetNeighborDown(2,procTable[i][j][k-1]);
					if (k<m_SplitNumber[2].size()-2)
					{
						//add one additional line
						grid->AddDiscLine(2, m_Original_Grid->GetLine(2,m_SplitNumber[2].at(k+1)+1 ));
						m_MPI_Op->SetNeighborUp(2,procTable[i][j][k+1]);
					}

				}
			}
		}
	}

	m_MPI_Op->SetOriginalMesh(m_Original_Grid);

	m_MPI_Op->SetTag(0);

	return true;
}


bool openEMS_FDTD_MPI::SetupOperator()
{
	bool ret = true;
	if (m_engine == EngineType_MPI)
	{
		FDTD_Op = Operator_MPI::New();
	}
	else
	{
		ret = openEMS::SetupOperator();
	}

	m_MPI_Op = dynamic_cast<Operator_MPI*>(FDTD_Op);

	if ((m_MPI_Enabled) && (m_MPI_Op==NULL))
	{
		cerr << "openEMS_FDTD_MPI::SetupOperator: Error: MPI is enabled but requested engine does not support MPI... EXIT!!!" << endl;
		MPI_Barrier(MPI_COMM_WORLD);
		exit(0);
	}

	ret &= SetupMPI();

	return ret;
}

unsigned int openEMS_FDTD_MPI::GetNextStep()
{
	//start processing and get local next step
	int step=PA->Process();
	double currTS = FDTD_Eng->GetNumberOfTimesteps();
	if ((step<0) || (step>(int)(NrTS - currTS))) step=NrTS - currTS;

	int local_step=step;

	//find the smallest next step requestes by all processings
	MPI_Reduce(&local_step, &step, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
	//send the smallest next step to all
	MPI_Bcast(&step, 1, MPI_INT, 0, MPI_COMM_WORLD);

	return step;
}

bool openEMS_FDTD_MPI::CheckEnergyCalc()
{
	int local_Check = (int)m_ProcField->CheckTimestep();
	int result;

	//check if some process request an energy calculation --> the sum is larger than 0
	MPI_Reduce(&local_Check, &result, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	//send result to all
	MPI_Bcast(&result, 1, MPI_INT, 0, MPI_COMM_WORLD);

	//calc energy if result is non-zero
	return result>0;
}

double openEMS_FDTD_MPI::CalcEnergy()
{
	double energy = 0;
	double loc_energy= m_ProcField->CalcTotalEnergyEstimate();

	//calc the sum of all local energies
	MPI_Reduce(&loc_energy, &energy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	//send sum-energy to all processes
	MPI_Bcast(&energy, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (energy>m_MaxEnergy)
		m_MaxEnergy = energy;
	if (m_MaxEnergy)
		m_EnergyDecrement = energy/m_MaxEnergy;

	return energy;
}

bool openEMS_FDTD_MPI::SetupProcessing()
{
	bool ret = openEMS::SetupProcessing();

	//search for active processings in different processes
	size_t numProc = PA->GetNumberOfProcessings();
	int active=0;
	bool deactivate = false;
	bool rename = false;
	for (size_t n=0;n<numProc;++n)
	{
		Processing* proc = PA->GetProcessing(n);
		int isActive = (int)proc->GetEnable();
		//sum of all active processings
		MPI_Reduce(&isActive, &active, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		deactivate = false;
		rename = false;
		if ((m_MyID==0) && (active>1)) //more than one active processing...
		{
			deactivate = true; //default
			if (dynamic_cast<ProcessIntegral*>(proc)!=NULL)
			{
				//type is integral processing --> disable! Needs to be fixed!
				cerr << "openEMS_FDTD_MPI::SetupProcessing(): Warning: Processing: " << proc->GetName() << " occures multiple times and is being deactivated..." << endl;
				cerr << "openEMS_FDTD_MPI::SetupProcessing(): Note: Processing: Make sure that there are no splits inside probes or sources." << endl;
				deactivate = true;
				rename = false;
			}
			if (dynamic_cast<ProcessFields*>(proc)!=NULL)
			{
				//type is field processing --> renameing! Needs to be fixed!
				cerr << "openEMS_FDTD_MPI::SetupProcessing(): Warning: Processing: " << proc->GetName() << " occures multiple times and is being renamed..." << endl;
				deactivate = false;
				rename = true;
			}
		}
		//broadcast informations to all
		MPI_Bcast(&deactivate, 1, MPI::BOOL, 0, MPI_COMM_WORLD);
		MPI_Bcast(&rename, 1, MPI::BOOL, 0, MPI_COMM_WORLD);
		if (deactivate)
			proc->SetEnable(false);
		if (rename)
		{
			ProcessFields* ProcField = dynamic_cast<ProcessFields*>(proc);
			if (ProcField)
			{
				stringstream name_ss;
				name_ss << "ID" << m_MyID << "_" << ProcField->GetName();
				ProcField->SetName(name_ss.str());
				ProcField->SetFileName(name_ss.str());
			}
		}
	}
	return ret;
}

int openEMS_FDTD_MPI::SetupFDTD()
{
	return openEMS::SetupFDTD();
}

void openEMS_FDTD_MPI::RunFDTD()
{
	if (!m_MPI_Enabled)
		return openEMS::RunFDTD();

	cout << "Running MPI-FDTD engine... this may take a while... grab a cup of coffee?!?" << endl;

	//get the sum of all cells
	unsigned int local_NrCells=FDTD_Op->GetNumberCells();
	MPI_Reduce(&local_NrCells, &m_NumberCells, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(&m_NumberCells, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

	//special handling of a field processing, needed to realize the end criteria...
	m_ProcField = new ProcessFields(NewEngineInterface());
	PA->AddProcessing(m_ProcField);

	//init processings
	PA->InitAll();

	double currE=0;

	//add all timesteps to end-crit field processing with max excite amplitude
	unsigned int maxExcite = FDTD_Op->GetExcitationSignal()->GetMaxExcitationTimestep();
//	for (unsigned int n=0; n<FDTD_Op->Exc->Volt_Count; ++n)
//		m_ProcField->AddStep(FDTD_Op->Exc->Volt_delay[n]+maxExcite);
	m_ProcField->AddStep(maxExcite);

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
	int step = GetNextStep();

	while ((step>0) && !CheckAbortCond())
	{
		FDTD_Eng->IterateTS(step);
		step = GetNextStep();

		currTS = FDTD_Eng->GetNumberOfTimesteps();

		currE = 0;
		gettimeofday(&currTime,NULL);
		t_diff = CalcDiffTime(currTime,prevTime);

		if (CheckEnergyCalc())
			currE = CalcEnergy();

		//make sure all processes are at the same simulation time
		MPI_Bcast(&t_diff, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		if (t_diff>4)
		{
			if (currE==0)
				currE = CalcEnergy();
			if (m_MyID==0)
			{
				t_run = CalcDiffTime(currTime,startTime);
				speed = numCells*(currTS-prevTS)/t_diff;
				cout << "[@" << FormatTime(t_run)  <<  "] Timestep: " << setw(12)  << currTS ;
				cout << " || Speed: " << setw(6) << setprecision(1) << std::fixed << speed*1e-6 << " MC/s (" <<  setw(4) << setprecision(3) << std::scientific << t_diff/(currTS-prevTS) << " s/TS)" ;
				cout << " || Energy: ~" << setw(6) << setprecision(2) << std::scientific << currE << " (-" << setw(5)  << setprecision(2) << std::fixed << fabs(10.0*log10(m_EnergyDecrement)) << "dB)" << endl;

				//set step to zero to abort simulation and send to all
				if (m_EnergyDecrement<endCrit)
					step=0;

				if (m_DumpStats)
					DumpRunStatistics(__OPENEMS_RUN_STAT_FILE__, t_run, currTS, speed, currE);
			}

			MPI_Bcast(&step, 1, MPI_INT, 0, MPI_COMM_WORLD);
			
			prevTime=currTime;
			prevTS=currTS;

			PA->FlushNext();
		}
	}
	if ((m_MyID==0) && (m_EnergyDecrement>endCrit) && (FDTD_Op->GetExcitationSignal()->GetExciteType()==0))
		cerr << "RunFDTD: max. number of timesteps was reached before the end-criteria of -" << fabs(10.0*log10(endCrit)) << "dB was reached... " << endl << \
				"\tYou may want to choose a higher number of max. timesteps... " << endl;

	gettimeofday(&currTime,NULL);

	t_diff = CalcDiffTime(currTime,startTime);

	if (m_MyID==0)
	{
		cout << "Time for " << FDTD_Eng->GetNumberOfTimesteps() << " iterations with " << FDTD_Op->GetNumberCells() << " cells : " << t_diff << " sec" << endl;
		cout << "Speed: " << numCells*(double)FDTD_Eng->GetNumberOfTimesteps()/t_diff*1e-6 << " MCells/s " << endl;

		if (m_DumpStats)
			DumpStatistics(__OPENEMS_STAT_FILE__, t_diff);
	}

	//*************** postproc ************//
	PA->PostProcess();
}
