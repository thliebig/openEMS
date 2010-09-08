/*
*	Copyright (C) 2010 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#include "operator_cylindermultigrid.h"
#include "engine_cylindermultigrid.h"
#include "operator_ext_cylinder.h"

Operator_CylinderMultiGrid::Operator_CylinderMultiGrid(vector<double> Split_Radii) : Operator_Cylinder()
{
	m_Split_Radii = Split_Radii;
	m_Split_Rad = m_Split_Radii.back();
	m_Split_Radii.pop_back();
}

Operator_CylinderMultiGrid::~Operator_CylinderMultiGrid()
{
}

Operator_CylinderMultiGrid* Operator_CylinderMultiGrid::New(vector<double> Split_Radii, unsigned int numThreads)
{
	cout << "Create cylindrical multi grid FDTD operator " << endl;
	Operator_CylinderMultiGrid* op = new Operator_CylinderMultiGrid(Split_Radii);
	op->setNumThreads(numThreads);
	op->Init();

	return op;
}

Engine* Operator_CylinderMultiGrid::CreateEngine() const
{
	Engine_CylinderMultiGrid* eng = Engine_CylinderMultiGrid::New(this,m_numThreads);
	return eng;
}

double Operator_CylinderMultiGrid::GetNumberCells() const
{
	if (numLines)
		return (numLines[0]-m_Split_Pos)*(numLines[1])*(numLines[2]) + m_InnerOp->GetNumberCells();
	return 0;
}

bool Operator_CylinderMultiGrid::SetGeometryCSX(ContinuousStructure* geo)
{
	if (Operator_Cylinder::SetGeometryCSX(geo)==false)
		return false;

	if (numLines[1]%2 != 1)
	{
		cerr << "Operator_CylinderMultiGrid::SetGeometryCSX: Error, number of line in alpha direction must be odd... found: " << numLines[1] << endl;
		exit(0);
	}

	//check if mesh is homogenous in alpha-direction
	double diff=discLines[1][1]-discLines[1][0];
	for (unsigned int n=2;n<numLines[1];++n)
	{
		if ( fabs((discLines[1][n]-discLines[1][n-1]) - diff)/diff > 1e-10)
		{
			cerr << "Operator_CylinderMultiGrid::SetGeometryCSX: Error, mesh has to be homogenous in alpha direction for multi grid engine, violation found at: " << n << endl;
			exit(0);
		}
	}

	m_Split_Pos = 0;
	for (unsigned int n=0;n<numLines[0];++n)
	{
		if (m_Split_Rad < discLines[0][n])
		{
			m_Split_Pos = n;
			cout << "Operator_CylinderMultiGrid::SetGeometryCSX: Found mesh split position @" << m_Split_Pos << endl;
			m_Split_Rad = discLines[0][n];
			break;
		}
	}
	if ((m_Split_Pos<4) || (m_Split_Pos>numLines[0]-4))
	{
		cerr << "Operator_CylinderMultiGrid::SetGeometryCSX: Error, split invalid..." << endl;
		return false;
	}

	CSRectGrid* grid = geo->GetGrid();

	grid->ClearLines(0);
	grid->ClearLines(1);
	for (unsigned int n=0; n<m_Split_Pos ;++n)
		grid->AddDiscLine(0,discLines[0][n]);
	for (unsigned int n=0; n<numLines[1]; n+=2)
		grid->AddDiscLine(1,discLines[1][n]);

	if (m_InnerOp->SetGeometryCSX(CSX)==false)
		return false;

	//restore grid to original mesh
	grid->ClearLines(0);
	grid->ClearLines(1);
	for (unsigned int n=0; n<numLines[0]; ++n)
		grid->AddDiscLine(0,discLines[0][n]);
	for (unsigned int n=0; n<numLines[1]; ++n)
		grid->AddDiscLine(1,discLines[1][n]);

	return true;
}

void Operator_CylinderMultiGrid::Init()
{
	Operator_Cylinder::Init();

	if (m_Split_Radii.empty())
		m_InnerOp = Operator_Cylinder::New(m_numThreads);
	else
		m_InnerOp = Operator_CylinderMultiGrid::New(m_Split_Radii,m_numThreads);
}

void Operator_CylinderMultiGrid::CalcStartStopLines(unsigned int &numThreads, vector<unsigned int> &start, vector<unsigned int> &stop) const
{
	if (numLines[0]<numThreads) //in case more threads requested as lines in r-direction, reduce number of worker threads
		numThreads = numLines[0];

	unsigned int linesPerThread = round((float)(numLines[0] - m_Split_Pos + 1) / (float)numThreads);
	if ((numThreads-1) * linesPerThread >= (numLines[0] - m_Split_Pos + 1))
		--numThreads;

	start.resize(numThreads);
	stop.resize(numThreads);

	for (unsigned int n=0; n<numThreads; n++)
	{
		start.at(n) = n * linesPerThread + m_Split_Pos - 1;
		stop.at(n) = (n+1) * linesPerThread - 1 + m_Split_Pos - 1;
		if (n == numThreads-1) // last thread
			stop.at(n) = numLines[0]-1;
	}
}


int Operator_CylinderMultiGrid::CalcECOperator()
{
	int retCode=0;
	if (dT)
		m_InnerOp->SetTimestep(dT);

	//calc inner child first
	m_InnerOp->CalcECOperator();

	dT = m_InnerOp->GetTimestep();

	retCode = Operator_Cylinder::CalcECOperator();
	if (GetTimestepValid()==false)
	{
		cerr << "Operator_CylinderMultiGrid::CalcECOperator(): Warning, timestep invalid... resetting..." << endl;
		dT = opt_dT;
		m_InnerOp->SetTimestep(dT);
		m_InnerOp->CalcECOperator();
		return Operator_Cylinder::CalcECOperator();
	}

	return retCode;
}

bool Operator_CylinderMultiGrid::SetupExcitation(TiXmlElement* Excite, unsigned int maxTS)
{
	if (!m_InnerOp->SetupExcitation(Excite,maxTS))
		return false;
	return Exc->setupExcitation(Excite,maxTS);
}

void Operator_CylinderMultiGrid::Reset()
{
	Operator_Cylinder::Reset();
	m_InnerOp->Reset();
}

void Operator_CylinderMultiGrid::SetBoundaryCondition(int* BCs)
{
	Operator_Cylinder::SetBoundaryCondition(BCs);
	BCs[1] = 0; //always PEC in +r-direction
	m_InnerOp->SetBoundaryCondition(BCs);
}

void Operator_CylinderMultiGrid::AddExtension(Operator_Extension* op_ext)
{
	if (dynamic_cast<Operator_Ext_Cylinder*>(op_ext))
	{
		Operator_Cylinder::AddExtension(op_ext);
		return;
	}
	Operator_Extension* child_Ext = op_ext->Clone(m_InnerOp);
	if (child_Ext==NULL)
	{
		cerr << "Operator_CylinderMultiGrid::AddExtension: Warning, extension: " << op_ext->GetExtensionName() << " can not be cloned for the child operator. Skipping Extension... " << endl;
		return;
	}
	Operator_Cylinder::AddExtension(op_ext);
	//give the copy to child
	m_InnerOp->AddExtension(child_Ext);
}

void Operator_CylinderMultiGrid::ShowStat() const
{
	m_InnerOp->ShowStat();
	m_InnerOp->ShowExtStat();
	Operator_Cylinder::ShowStat();
}
