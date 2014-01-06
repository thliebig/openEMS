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

#include "operator_mpi.h"
#include "operator_sse_compressed.h"
#include "engine_sse_compressed.h"
#include "engine_mpi.h"
#include "extensions/operator_extension.h"
#include "tools/array_ops.h"
#include "tools/useful.h"
#include "mpi.h"

Operator_MPI* Operator_MPI::New()
{
	cout << "Create FDTD operator (compressed SSE + MPI)" << endl;
	Operator_MPI* op = new Operator_MPI();
	op->Init();
	return op;
}

Operator_MPI::Operator_MPI() : Operator_SSE_Compressed()
{
	m_NumProc = MPI::COMM_WORLD.Get_size();

	//enabled only if more than one process is active
	m_MPI_Enabled = m_NumProc>1;
}

Operator_MPI::~Operator_MPI()
{
	Delete();
}

double Operator_MPI::CalcTimestep()
{
	double ret = Operator::CalcTimestep();

	if (!m_MPI_Enabled)
		return ret;

	double local_dT = dT;
	//find the smallest time-step requestes by all processings
	MPI_Reduce(&local_dT, &dT, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
	//send the smallest time-step to all
	MPI_Bcast(&dT, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	return ret;
}

void Operator_MPI::SetBoundaryCondition(int* BCs)
{
	if (!m_MPI_Enabled)
		return Operator_SSE_Compressed::SetBoundaryCondition(BCs);

	//set boundary conditions on MPI interfaces to PEC, ApplyElectricBC will handle proper interface handling...
	for (int n=0;n<3;++n)
	{
		if (m_NeighborUp[n]>=0)
			BCs[2*n+1] = 0;
		if (m_NeighborDown[n]>=0)
			BCs[2*n] = 0;
	}
	Operator_SSE_Compressed::SetBoundaryCondition(BCs);
}

Engine* Operator_MPI::CreateEngine()
{
	if (m_MPI_Enabled)
		m_Engine = Engine_MPI::New(this);
	else
		m_Engine = Engine_SSE_Compressed::New(this);
	return m_Engine;
}

void Operator_MPI::SetNeighborUp(int ny, int id)
{
	if ((ny<0) || (ny>2))
		return;
	m_NeighborUp[ny]=id;
}

void Operator_MPI::SetNeighborDown(int ny, int id)
{
	if ((ny<0) || (ny>2))
		return;
	m_NeighborDown[ny]=id;
}

void Operator_MPI::Init()
{
	Operator_SSE_Compressed::Init();

	m_MyTag = 0;

	for (int i=0;i<3;++i)
	{
		m_NeighborUp[i]=-1;
		m_NeighborDown[i]=-1;
		m_OrigDiscLines[i]=NULL;
	}

	m_ProcTable = NULL;
	m_SplitNumber[0]=0;
	m_SplitNumber[1]=0;
	m_SplitNumber[2]=0;

	int  namelen;
	m_NumProc = MPI::COMM_WORLD.Get_size();
	m_MyID = MPI::COMM_WORLD.Get_rank();

	m_Processor_Name = new char[MPI_MAX_PROCESSOR_NAME];
	MPI::Get_processor_name(m_Processor_Name,namelen);

	if ((m_MPI_Enabled) && (g_settings.GetVerboseLevel()>0))
		cerr << "Operator_MPI::Init(): Running on " << m_Processor_Name << endl;
}

void Operator_MPI::Delete()
{
	delete[] m_Processor_Name;
	m_Processor_Name = NULL;
	for (int i=0;i<3;++i)
	{
		delete[] m_OrigDiscLines[i];
		m_OrigDiscLines[i] = NULL;
	}
	Delete3DArray(m_ProcTable,m_SplitNumber);
	m_ProcTable=NULL;
}

void Operator_MPI::Reset()
{
	Delete();
	Operator_SSE_Compressed::Reset();
}

void Operator_MPI::SetOriginalMesh(CSRectGrid* orig_Mesh)
{
	for (int n=0;n<3;++n)
	{
		delete[] m_OrigDiscLines[n];
		m_OrigDiscLines[n] = orig_Mesh->GetLines(n,NULL,m_OrigNumLines[n]);
	}
}

unsigned int Operator_MPI::GetNumberOfLines(int ny, bool fullMesh) const
{
	if (fullMesh)
		return Operator_SSE_Compressed::GetNumberOfLines(ny,fullMesh);

	if ((!m_MPI_Enabled) || (m_NeighborUp[ny]<0))
		return Operator_SSE_Compressed::GetNumberOfLines(ny,fullMesh);

	return Operator_SSE_Compressed::GetNumberOfLines(ny)-1;
}

void Operator_MPI::AddExtension(Operator_Extension* op_ext)
{
	if (m_MPI_Enabled==false)
		return Operator_SSE_Compressed::AddExtension(op_ext);

	if (op_ext->IsMPISave())
		Operator_SSE_Compressed::AddExtension(op_ext);
	else
	{
		cerr << "Operator_MPI::AddExtension: Warning: Operator extension \"" << op_ext->GetExtensionName() << "\" is not compatible with MPI!! skipping...!" << endl;
		delete op_ext;
	}
}


string Operator_MPI::PrependRank(string name)
{
	stringstream out_name;
	if (m_MPI_Enabled)
		out_name << "ID" << m_MyID << "_" << name;
	else
		out_name << name;
	return out_name.str();
}

void Operator_MPI::DumpOperator2File(string filename)
{
	Operator_SSE_Compressed::DumpOperator2File(PrependRank(filename));
}

void Operator_MPI::DumpMaterial2File(string filename)
{
	Operator_SSE_Compressed::DumpMaterial2File(PrependRank(filename));
}

void Operator_MPI::DumpPEC2File(string filename , unsigned int *range)
{
	Operator_SSE_Compressed::DumpPEC2File(PrependRank(filename), range);
}
