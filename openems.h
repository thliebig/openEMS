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

#ifndef OPENEMS_H
#define OPENEMS_H

#include <sstream>
#ifndef __GNUC__
#include <Winsock2.h> // for struct timeval
#else
#include <sys/time.h>
#endif
#include <time.h>
#include <vector>

#include "openems_global.h"

#define __OPENEMS_STAT_FILE__ "openEMS_stats.txt"
#define __OPENEMS_RUN_STAT_FILE__ "openEMS_run_stats.txt"

class Operator;
class Engine;
class Engine_Interface_FDTD;
class ProcessingArray;
class TiXmlElement;
class ContinuousStructure;
class Engine_Interface_FDTD;
class Excitation;
class Engine_Ext_SteadyState;

double CalcDiffTime(timeval t1, timeval t2);
std::string FormatTime(int sec);

class OPENEMS_EXPORT openEMS
{
public:
	openEMS();
	virtual ~openEMS();

	virtual bool parseCommandLineArgument( const char *argv );
	static void showUsage();

	bool ParseFDTDSetup(std::string file);
	virtual bool Parse_XML_FDTDSetup(TiXmlElement* openEMSxml);
	virtual int SetupFDTD();
	virtual void RunFDTD();

	void Reset();

	void SetNumberOfTimeSteps(unsigned int val) {NrTS=val;}
	void SetEnableDumps(bool val) {Enable_Dumps=val;}
	void SetEndCriteria(double val) {endCrit=val;}
	void SetOverSampling(int val) {m_OverSampling=val;}
	void SetCellConstantMaterial(bool val) {m_CellConstantMaterial=val;}

	void SetCylinderCoords(bool val) {CylinderCoords=val;}
	void SetupCylinderMultiGrid(std::vector<double> val) {m_CC_MultiGrid=val;}
	void SetupCylinderMultiGrid(std::string val);

	void SetTimeStepMethod(int val) {m_TS_method=val;}
	void SetTimeStep(double val) {m_TS=val;}
	void SetTimeStepFactor(double val) {m_TS_fac=val;}
	void SetMaxTime(double val) {m_maxTime=val;}

	void SetNumberOfThreads(int val);

	void DebugMaterial() {DebugMat=true;}
	void DebugOperator() {DebugOp=true;}
	void DebugBox() {m_debugBox=true;}
	void DebugPEC() {m_debugPEC=true;}
	void DebugCSX() {m_debugCSX=true;}

	void Set_BC_Type(int idx, int type);
	int Get_BC_Type(int idx);
	void Set_BC_PML(int idx, unsigned int size);
	int Get_PML_Size(int idx);
	void Set_Mur_PhaseVel(int idx, double val);

	//! Get information about external libs used by openEMS
	static std::string GetExtLibsInfo(std::string prefix="\t");

	//! Get welcome screen for openEMS
	static void WelcomeScreen();

	//! Set this to about FDTD iteration process
	void SetAbort(bool val) {m_Abort=val;}
	//! Check for abort conditions
	bool CheckAbortCond();

	void SetGaussExcite(double f0, double fc);
	void SetSinusExcite(double f0);
	void SetDiracExcite(double f_max);
	void SetStepExcite(double f_max);

	Excitation* InitExcitation();

	void SetCSX(ContinuousStructure* csx);

	Engine_Interface_FDTD* NewEngineInterface(int multigridlevel = 0);

	void SetVerboseLevel(int level);

protected:
	bool CylinderCoords;
	std::vector<double> m_CC_MultiGrid;

	ContinuousStructure* m_CSX;

	//! Number of Timesteps
	unsigned int NrTS;
	int m_TS_method;
	double m_TS;
	double m_TS_fac;
	double m_maxTime;

	// some command line flags
	bool Enable_Dumps;
	bool DebugMat;
	bool DebugOp;
	bool m_debugCSX;
	bool m_DumpStats;
	bool m_debugBox, m_debugPEC, m_no_simulation;

	double endCrit;
	int m_OverSampling;
	bool m_CellConstantMaterial;
	Operator* FDTD_Op;
	Engine* FDTD_Eng;
	Engine_Ext_SteadyState* Eng_Ext_SSD;
	ProcessingArray* PA;

	Excitation* m_Exc;

	bool m_Abort;

#ifdef MPI_SUPPORT
	enum EngineType {EngineType_Basic, EngineType_SSE, EngineType_SSE_Compressed, EngineType_Multithreaded, EngineType_MPI};
#else
	enum EngineType {EngineType_Basic, EngineType_SSE, EngineType_SSE_Compressed, EngineType_Multithreaded};
#endif
	EngineType m_engine;
	unsigned int m_engine_numThreads;

	//! Setup an operator matching the requested engine
	virtual bool SetupOperator();

	//! Read boundary conditions from xml element and apply to FDTD operator
	bool SetupBoundaryConditions();
	int m_BC_type[6];
	unsigned int m_PML_size[6];
	double m_Mur_v_ph[6];

	//! Check whether or not the FDTD-Operator has to store material data.
	bool SetupMaterialStorages();

	//! Setup all processings.
	virtual bool SetupProcessing();

	//! Dump statistics to file
	virtual bool DumpStatistics(const std::string& filename, double time);

	//! Dump run statistivs to file
	virtual bool InitRunStatistics(const std::string& filename);
	virtual bool DumpRunStatistics(const std::string& filename, double time, unsigned int ts, double speed, double energy);
};

#endif // OPENEMS_H
