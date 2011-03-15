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

#include "engine.h"
#include "Common/processfields.h"
#include "operator_cylinder.h"
#include "extensions/operator_extension.h"
#include "extensions/operator_ext_cylinder.h"

Operator_Cylinder* Operator_Cylinder::New(unsigned int numThreads)
{
	cout << "Create cylindrical FDTD operator" << endl;
	Operator_Cylinder* op = new Operator_Cylinder();
	op->setNumThreads(numThreads);
	op->Init();
	return op;
}

Operator_Cylinder::Operator_Cylinder() : Operator_Multithread()
{
	m_MeshType = ProcessFields::CYLINDRICAL_MESH;
}

Operator_Cylinder::~Operator_Cylinder()
{
}

void Operator_Cylinder::Init()
{
	CC_closedAlpha = false;
	CC_R0_included = false;
	Operator_Multithread::Init();
}

int Operator_Cylinder::CalcECOperator( DebugFlags debugFlags )
{
	// debugs only work with the native vector dumps
	bool natDump = g_settings.NativeFieldDumps();
	g_settings.SetNativeFieldDumps(true);
	int rc = Operator_Multithread::CalcECOperator(debugFlags);
	// reset original settings
	g_settings.SetNativeFieldDumps(natDump);
	return rc;
}

inline unsigned int Operator_Cylinder::GetNumberOfLines(int ny) const
{
	//this is necessary for a correct field processing... cylindrical engine has to reset this by adding +1
	if (CC_closedAlpha && ny==1)
		return Operator_Multithread::GetNumberOfLines(ny)-1;

	return Operator_Multithread::GetNumberOfLines(ny);
}

string Operator_Cylinder::GetDirName(int ny) const
{
	if (ny==0) return "rho";
	if (ny==1) return "alpha";
	if (ny==2) return "z";
	return "";
}

double Operator_Cylinder::GetNodeWidth(int ny, const unsigned int pos[3], bool dualMesh) const
{
	if ((ny<0) || (ny>2)) return 0.0;
	if (pos[ny]>=numLines[ny]) return 0.0;
	double width = Operator_Multithread::GetEdgeLength(ny,pos,!dualMesh);
	if (ny==1)
		width *= GetDiscLine(0,pos[0],dualMesh);
	return width;
}

double Operator_Cylinder::GetNodeArea(int ny, const unsigned int pos[3], bool dualMesh) const
{
	if (pos[ny]>=numLines[ny]) return 0.0;
	if (pos[0]>=numLines[0]) return 0.0;
	if (ny==2)
	{
		double da = Operator_Multithread::GetEdgeLength(1,pos,dualMesh)/gridDelta;
		double r1,r2;

		if (dualMesh)
		{
			r1 = GetDiscLine(0,pos[0],false)*gridDelta;
			r2 = r1 + GetEdgeLength(0,pos,false);
		}
		else
		{
			r2 = GetDiscLine(0,pos[0],!dualMesh)*gridDelta;
			r1 = r2 - GetEdgeLength(0,pos,true);
		}

		if (r1<=0)
			return da/2 * pow(r2,2);
		else
			return da/2* (pow(r2,2) - pow(r1,2));
	}

	return Operator_Multithread::GetNodeArea(ny,pos,dualMesh);
}

double Operator_Cylinder::GetEdgeLength(int ny, const unsigned int pos[3], bool dualMesh) const
{
	double length = Operator_Multithread::GetEdgeLength(ny,pos,dualMesh);
	if (ny!=1)
		return length;
	return length * GetDiscLine(0,pos[0],dualMesh);
}

double Operator_Cylinder::GetEdgeArea(int ny, const unsigned int pos[3], bool dualMesh) const
{
	if (ny!=0)
		return GetNodeArea(ny,pos,dualMesh);

	return GetEdgeLength(1,pos,!dualMesh) * GetEdgeLength(2,pos,!dualMesh);
}

bool Operator_Cylinder::SetGeometryCSX(ContinuousStructure* geo)
{
	if (Operator_Multithread::SetGeometryCSX(geo)==false) return false;

	double minmaxA = fabs(discLines[1][numLines[1]-1]-discLines[1][0]);
	if (fabs(minmaxA-2*PI) < (2*PI)/10/numLines[1]) //check minmaxA smaller then a tenth of average alpha-width
	{
		cout << "Operator_Cylinder::SetGeometryCSX: Alpha is a full 2*PI => closed Cylinder..." << endl;
		CC_closedAlpha = true;
		discLines[1][numLines[1]-1] = discLines[1][0] + 2*PI;
		cerr << "Operator_Cylinder::SetGeometryCSX: Warning, not handling the disc-line width and material averaging correctly yet for a closed cylinder..." << endl;
		if (MainOp->GetIndexDelta(1,0)-MainOp->GetIndexDelta(1,numLines[1]-2) > (2*PI)/10/numLines[1])
		{
			cerr << "Operator_Cylinder::SetGeometryCSX: first and last angle delta must be the same... deviation to large..." << MainOp->GetIndexDelta(1,0) - MainOp->GetIndexDelta(1,numLines[1]-2) << endl;
			exit(1);
		}
		if (MainOp->GetIndexDelta(1,0)-MainOp->GetIndexDelta(1,numLines[1]-2) > 0)
		{
			cerr << "Operator_Cylinder::SetGeometryCSX: first and last angle delta must be the same... auto correction of deviation: " << MainOp->GetIndexDelta(1,0) - MainOp->GetIndexDelta(1,numLines[1]-2) << endl;
			discLines[1][numLines[1]-2] = discLines[1][numLines[1]-1]-MainOp->GetIndexDelta(1,0);
		}
	}
	else if (minmaxA>2*PI)
	{
		cerr << "Operator_Cylinder::SetGeometryCSX: Alpha Max-Min must not be larger than 2*PI!!!" << endl;
		Reset();
		return false;
	}
	else
	{
		CC_closedAlpha=false;
	}

	if (discLines[0][0]<0)
	{
		cerr << "Operator_Cylinder::SetGeometryCSX: r<0 not allowed in Cylinder Coordinates!!!" << endl;
		Reset();
		return false;
	}
	else if (discLines[0][0]==0.0)
	{
		cout << "Operator_Cylinder::SetGeometryCSX: r=0 included..." << endl;
		CC_R0_included= true;  //also needed for correct ec-calculation
	}

	if (CC_closedAlpha || CC_R0_included)
		this->AddExtension(new Operator_Ext_Cylinder(this));

	return true;
}

void Operator_Cylinder::ApplyElectricBC(bool* dirs)
{
	if (dirs==NULL) return;
	if (CC_closedAlpha)
	{
		dirs[2]=0;
		dirs[3]=0; //no PEC in alpha directions...
	}
	if (CC_R0_included)
	{
		// E in alpha direction ( aka volt[1][x][y][z] ) is not defined for r==0 --> always zero...
		unsigned int pos[3] = {0,0,0};
		for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
		{
			for (pos[2]=0; pos[2]<numLines[2]; ++pos[2])
			{
				SetVV(1,pos[0],pos[1],pos[2], 0 );
				SetVI(1,pos[0],pos[1],pos[2], 0 );
			}
		}
	}
	Operator_Multithread::ApplyElectricBC(dirs);
}

void Operator_Cylinder::ApplyMagneticBC(bool* dirs)
{
	if (dirs==NULL) return;
	if (CC_closedAlpha)
	{
		dirs[2]=0;
		dirs[3]=0; //no PMC in alpha directions...
	}
	if (CC_R0_included)
	{
		dirs[0]=0;  //no PMC in r_min directions...
	}
	Operator_Multithread::ApplyMagneticBC(dirs);
}

void Operator_Cylinder::AddExtension(Operator_Extension* op_ext)
{
	if (op_ext->IsCylinderCoordsSave())
		m_Op_exts.push_back(op_ext);
	else
		cerr << "Operator_Cylinder::AddExtension: Warning: Operator extension \"" << op_ext->GetExtensionName() << "\" is not compatible with cylinder-coords!! skipping...!" << endl;
}

double Operator_Cylinder::CalcTimestep()
{
	m_TimeStepVar = 1;
	return Operator_Multithread::CalcTimestep();
}
