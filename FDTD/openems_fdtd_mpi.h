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
    openEMS_FDTD_MPI();
	virtual ~openEMS_FDTD_MPI();

	virtual void RunFDTD();

	virtual bool parseCommandLineArgument( const char *argv );

protected:
	Operator_MPI* m_MPI_Op;
	CSRectGrid* m_Original_Grid;
	int m_MyID;
	int m_NumProc;
	bool m_MPI_Enabled;
	unsigned int m_NumberCells;

	virtual bool SetupMPI(TiXmlElement* FDTD_Opts);
	virtual bool SetupOperator(TiXmlElement* FDTD_Opts);

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
};

#endif // OPENEMS_FDTD_MPI_H
