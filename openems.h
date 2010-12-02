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

using namespace std;

class Operator;
class Engine;
class Engine_Interface_FDTD;
class ProcessingArray;
class TiXmlElement;

class openEMS
{
public:
    openEMS();
	~openEMS();

	bool parseCommandLineArgument( const char *argv );

	int SetupFDTD(const char* file);

	void RunFDTD();

	void Reset();

	void SetEnableDumps(bool val) {Enable_Dumps=val;}
	void SetEndCriteria(double val) {endCrit=val;}

	void DebugMaterial() {DebugMat=true;}
	void DebugOperator() {DebugOp=true;}
	void DebugBox() {m_debugBox=true;}

	//! Get informations about external libs used by openEMS
	static string GetExtLibsInfo();

	//! Set this to about FDTD iteration process
	void SetAbort(bool val) {m_Abort=val;}
	//! Check for abort conditions
	bool CheckAbortCond();

protected:
	bool CylinderCoords;

	//! Number of Timesteps
	unsigned int NrTS;
	bool Enable_Dumps;
	bool DebugMat;
	bool DebugOp;
	bool m_debugCSX;
	bool m_debugBox, m_debugPEC, m_no_simulation;
	double endCrit;
	int m_OverSampling;
	Operator* FDTD_Op;
	Engine* FDTD_Eng;
	ProcessingArray* PA;

	bool m_Abort;

	enum EngineType {EngineType_Standard, EngineType_SSE, EngineType_SSE_Compressed, EngineType_Multithreaded};
	EngineType m_engine;
	unsigned int m_engine_numThreads;

	//! Read boundary conditions from xml element and apply to FDTD operator
	bool SetupBoundaryConditions(TiXmlElement* BC);
};

#endif // OPENEMS_H
