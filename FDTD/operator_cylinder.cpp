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
#include "engine_cylinder.h"
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
	m_MeshType = CYLINDRICAL;
	m_Cyl_Ext = NULL;
}

Operator_Cylinder::~Operator_Cylinder()
{
}


Engine* Operator_Cylinder::CreateEngine() const
{
	//! create a special cylindrical-engine
	Engine_Cylinder* eng = Engine_Cylinder::New(this, m_numThreads);
	return eng;
}

void Operator_Cylinder::Init()
{
	CC_closedAlpha = false;
	CC_R0_included = false;
	Operator_Multithread::Init();
}

double Operator_Cylinder::GetRawDiscDelta(int ny, const int pos) const
{
	if (CC_closedAlpha && ny==1 && pos==-1)
	{
//		cerr << (discLines[1][numLines[1]-2] - discLines[1][numLines[1]-3]) << " vs "  << Operator_Multithread::GetRawDiscDelta(ny,pos) << endl;
		return (discLines[1][numLines[1]-2] - discLines[1][numLines[1]-3]);
	}

	return Operator_Multithread::GetRawDiscDelta(ny,pos);
}

double Operator_Cylinder::GetMaterial(int ny, const double* coords, int MatType, bool markAsUsed) const
{
	double l_coords[] = {coords[0],coords[1],coords[2]};
	if (CC_closedAlpha && (coords[1]>GetDiscLine(1,0,false)+2*PI))
		l_coords[1]-=2*PI;
	if (CC_closedAlpha && (coords[1]<GetDiscLine(1,0,false)))
		l_coords[1] += 2*PI;
	return Operator_Multithread::GetMaterial(ny,l_coords,MatType,markAsUsed);
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
		return Operator_Multithread::GetNumberOfLines(ny)-2;

	return Operator_Multithread::GetNumberOfLines(ny);
}

string Operator_Cylinder::GetDirName(int ny) const
{
	if (ny==0) return "rho";
	if (ny==1) return "alpha";
	if (ny==2) return "z";
	return "";
}

bool Operator_Cylinder::GetYeeCoords(int ny, unsigned int pos[3], double* coords, bool dualMesh) const
{
	bool ret = Operator_Multithread::GetYeeCoords(ny,pos,coords,dualMesh);

	if (CC_closedAlpha && (coords[1]>=GetDiscLine(1,0,false)+2*PI))
		coords[1]-=2*PI;
	if (CC_closedAlpha && (coords[1]<GetDiscLine(1,0,false)))
		coords[1]+=2*PI;

	return ret;
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

double Operator_Cylinder::GetNodeWidth(int ny, const int pos[3], bool dualMesh) const
{
	if ( (pos[0]<0) || (pos[1]<0 && CC_closedAlpha==false) || (pos[2]<0) )
		return 0.0;

	unsigned int uiPos[]={(unsigned int)pos[0],(unsigned int)pos[1],(unsigned int)pos[2]};
	if (pos[1]<0 && CC_closedAlpha==true)
		uiPos[1]+=numLines[1]-2;

	return GetNodeWidth(ny, uiPos, dualMesh);
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

double Operator_Cylinder::GetNodeArea(int ny, const int pos[3], bool dualMesh) const
{
	if ( (pos[0]<0) || (pos[1]<0 && CC_closedAlpha==false) || (pos[2]<0) )
		return 0.0;

	unsigned int uiPos[]={(unsigned int)pos[0],(unsigned int)pos[1],(unsigned int)pos[2]};
	if (pos[1]<0 && CC_closedAlpha==true)
		uiPos[1]+=numLines[1]-2;

	return GetNodeArea(ny, uiPos, dualMesh);
}

double Operator_Cylinder::GetEdgeLength(int ny, const unsigned int pos[3], bool dualMesh) const
{
	double length = Operator_Multithread::GetEdgeLength(ny,pos,dualMesh);
	if (ny!=1)
		return length;
	return length * GetDiscLine(0,pos[0],dualMesh);
}

double Operator_Cylinder::GetCellVolume(const unsigned int pos[3], bool dualMesh) const
{
	return GetEdgeArea(2,pos,dualMesh)*GetEdgeLength(2,pos,dualMesh);
}

double Operator_Cylinder::GetEdgeArea(int ny, const unsigned int pos[3], bool dualMesh) const
{
	if (ny!=0)
		return GetNodeArea(ny,pos,dualMesh);

	return GetEdgeLength(1,pos,!dualMesh) * GetEdgeLength(2,pos,!dualMesh);
}

bool Operator_Cylinder::SetupCSXGrid(CSRectGrid* grid)
{
	unsigned int alphaNum;
	double* alphaLines = NULL;
	alphaLines = grid->GetLines(1,alphaLines,alphaNum,true);

	double minmaxA = fabs(alphaLines[alphaNum-1]-alphaLines[0]);
	if (fabs(minmaxA-2*PI) < OPERATOR_CYLINDER_CLOSED_ALPHA_THRESHOLD)
	{
		if (g_settings.GetVerboseLevel()>0)
			cout << "Operator_Cylinder::SetupCSXGrid: Alpha is a full 2*PI => closed Cylinder..." << endl;
		CC_closedAlpha = true;
		grid->SetLine(1,alphaNum-1,2*PI+alphaLines[0]);
		grid->AddDiscLine(1,2*PI+alphaLines[1]);
	}
	else if (minmaxA>2*PI)
	{
		cerr << "Operator_Cylinder::SetupCSXGrid: Alpha Max-Min must not be larger than 2*PI!!!" << endl;
		Reset();
		return false;
	}
	else
	{
		CC_closedAlpha=false;
	}

	CC_R0_included = false;
	if (grid->GetLine(0,0)<0)
	{
		cerr << "Operator_Cylinder::SetupCSXGrid: r<0 not allowed in Cylinder Coordinates!!!" << endl;
		Reset();
		return false;
	}
	else if (grid->GetLine(0,0)==0.0)
	{
		if (g_settings.GetVerboseLevel()>0)
			cout << "Operator_Cylinder::SetupCSXGrid: r=0 included..." << endl;
		CC_R0_included = CC_closedAlpha;  //needed for correct ec-calculation, deactivate if closed cylinder is false... --> E_r = 0 anyways

		// use conservative timestep for a mesh including the r==0 singularity
		m_TimeStepVar = 1;
	}

#ifdef MPI_SUPPORT
		// Setup an MPI split in alpha direction for a closed cylinder
		CC_MPI_Alpha = false;
	if ((m_NeighborUp[1]>=0) || (m_NeighborDown[1]>=0)) //check for MPI split in alpha direction
	{
		double minmaxA = 2*PI;// fabs(m_OrigDiscLines[1][m_OrigNumLines[1]-1]-m_OrigDiscLines[1][0]);
		if (fabs(minmaxA-2*PI) < OPERATOR_CYLINDER_CLOSED_ALPHA_THRESHOLD) //check for closed alpha MPI split
		{
			CC_MPI_Alpha = true;
			if (m_OrigDiscLines[0][0]==0)
			{
				cerr << "Operator_Cylinder::SetupCSXGrid: Error: MPI split in alpha direction for closed cylinder including r==0 is currently not supported! Exit!" << endl;
				exit(-2);
			}

			if (m_NeighborUp[1]<0) //check if this process is at the alpha-end
			{
				grid->SetLine(1,alphaNum-1,2*PI+m_OrigDiscLines[1][0]);
				grid->AddDiscLine(1,2*PI+m_OrigDiscLines[1][1]);

				SetNeighborUp(1,m_ProcTable[m_ProcTablePos[0]][0][m_ProcTablePos[2]]);
			}

			if (m_NeighborDown[1]<0) //check if this process is at the alpha-start
			{
				SetNeighborDown(1,m_ProcTable[m_ProcTablePos[0]][m_SplitNumber[1]-1][m_ProcTablePos[2]]);
			}

			//Note: the process table will not reflect this up/down neighbors necessary for a closed cylinder
		}
	}
#endif

	if (Operator_Multithread::SetupCSXGrid(grid)==false)
		return false;

	if (CC_closedAlpha || CC_R0_included)
	{
		m_Cyl_Ext = new Operator_Ext_Cylinder(this);
		this->AddExtension(m_Cyl_Ext);
	}

	return true;
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
	if (op_ext->IsCylinderCoordsSave(CC_closedAlpha, CC_R0_included))
		Operator_Multithread::AddExtension(op_ext);
	else
	{
		cerr << "Operator_Cylinder::AddExtension: Warning: Operator extension \"" << op_ext->GetExtensionName() << "\" is not compatible with cylinder-coords!! skipping...!" << endl;
		delete op_ext;
	}
}
