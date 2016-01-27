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

#ifndef OPENEMS_FDTD_MPI_H
#define OPENEMS_FDTD_MPI_H

#include "openems.h"

class ProcessFields;
class Operator_MPI;
class CSRectGrid;

class openEMS_FDTD_MPI : public openEMS
{
public:
	openEMS_FDTD_MPI(bool m_MPI_Debug=false);
	virtual ~openEMS_FDTD_MPI();

	virtual bool Parse_XML_FDTDSetup(TiXmlElement* FDTD_Opts);
	virtual int SetupFDTD();
	virtual void RunFDTD();

	virtual bool parseCommandLineArgument( const char *argv );

	static std::string GetExtLibsInfo();

protected:
	Operator_MPI* m_MPI_Op;
	CSRectGrid* m_Original_Grid;
	unsigned int m_MyID;
	unsigned int m_NumProc;
	bool m_MPI_Enabled;
	unsigned int m_NumberCells;

	std::vector<unsigned int> m_SplitNumber[3];
	TiXmlElement* m_MPI_Elem;
	virtual bool SetupMPI();
	virtual bool SetupOperator();

	int* m_Gather_Buffer;
	unsigned int GetNextStep();

	ProcessFields* m_ProcField;
	double m_MaxEnergy;
	double m_EnergyDecrement;
	double* m_Energy_Buffer;
	//! Check if energy calc is requested...
	bool CheckEnergyCalc();
	//! Calc energy in all processes and add up
	double CalcEnergy();

	virtual bool SetupProcessing();

	//output redirection to file for ranks > 0
	std::ofstream* m_Output;
};

#endif // OPENEMS_FDTD_MPI_H
