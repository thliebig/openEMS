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

#ifndef OPERATOR_MPI_H
#define OPERATOR_MPI_H

#include "operator_sse_compressed.h"

class Operator_MPI : public Operator_SSE_Compressed
{
	friend class Engine_MPI;
public:
	//! Create a new operator
	static Operator_MPI* New();
	virtual ~Operator_MPI();

	bool GetMPIEnabled() const {return m_MPI_Enabled;}

	virtual bool SetGeometryCSX(ContinuousStructure* geo);

	virtual void SetBoundaryCondition(int* BCs);
	virtual void ApplyElectricBC(bool* dirs);

	virtual Engine* CreateEngine() const;

protected:
	Operator_MPI();
	bool m_MPI_Enabled;
	virtual void Init();
	void Delete();
	virtual void Reset();

	virtual double CalcTimestep();

	int m_MyID;
	int m_NumProc;
	int m_MyTag;
	char* m_Processor_Name;

	//the up and down neighbors, -1 if non for the given direction
	int m_NeighborUp[3];
	int m_NeighborDown[3];

	string PrependRank(string name);

	virtual void DumpOperator2File(string filename);
	virtual void DumpMaterial2File(string filename);
	virtual void DumpPEC2File( string filename );
};

#endif // OPERATOR_MPI_H
