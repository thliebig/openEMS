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
#include "tools/useful.h"

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


Engine* Operator_Cylinder::CreateEngine()
{
	//! create a special cylindrical-engine
	m_Engine = Engine_Cylinder::New(this, m_numThreads);
	return m_Engine;
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
		return (discLines[1][numLines[1]-2] - discLines[1][numLines[1]-3]);
	if (CC_closedAlpha && ny==1 && (pos==(int)numLines[ny]-1))
		return (discLines[1][2] - discLines[1][1]);

	return Operator_Multithread::GetRawDiscDelta(ny,pos);
}

double Operator_Cylinder::GetMaterial(int ny, const double* coords, int MatType, vector<CSPrimitives*> vPrims, bool markAsUsed) const
{
	double l_coords[] = {coords[0],coords[1],coords[2]};
	if (CC_closedAlpha && (coords[1]>GetDiscLine(1,0,false)+2*PI))
		l_coords[1]-=2*PI;
	if (CC_closedAlpha && (coords[1]<GetDiscLine(1,0,false)))
		l_coords[1] += 2*PI;
	return Operator_Multithread::GetMaterial(ny,l_coords,MatType,vPrims,markAsUsed);
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

double Operator_Cylinder::CalcTimestep()
{
	if (discLines[0][0]==0.0)
		// use conservative timestep for a mesh including the r==0 singularity
		m_TimeStepVar = 1;

	return Operator::CalcTimestep();
}

inline unsigned int Operator_Cylinder::GetNumberOfLines(int ny, bool full) const
{
	if (full)
		return Operator_Multithread::GetNumberOfLines(ny, full);

	//this is necessary for a correct field processing... cylindrical engine has to reset this by adding +1
	if (CC_closedAlpha && ny==1)
		return Operator_Multithread::GetNumberOfLines(ny, true)-2;

	return Operator_Multithread::GetNumberOfLines(ny, full);
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

double Operator_Cylinder::FitToAlphaRange(double a_coord, bool fullMesh) const
{
	double min = GetDiscLine(1,0);
	double max = GetDiscLine(1,GetNumberOfLines(1, fullMesh)-1);
	if ((a_coord>=min) && (a_coord<=max))
		return a_coord;
	while (a_coord<min)
	{
		a_coord+=2*PI;
		if (a_coord>max)
			return a_coord-2*PI;
		if (a_coord>min)
			return a_coord;
	}
	while (a_coord>max)
	{
		a_coord-=2*PI;
		if (a_coord<min)
			return a_coord+2*PI;
		if (a_coord<max)
			return a_coord;
	}
	// this cannot happen
	return a_coord;
}

int Operator_Cylinder::MapAlphaIndex2Range(int pos) const
{
	if (!CC_closedAlpha)
		return pos;
	if (pos<0)
		return (int)numLines[1]+pos-2;
	else if (pos>=(int)numLines[1]-2)
		return pos-(int)numLines[1]+2;
	else
		return pos;
}

bool Operator_Cylinder::GetCellCenterMaterialAvgCoord(const int pos[3], double coord[3]) const
{
	if (!CC_closedAlpha || ((pos[1]>=0) && (pos[1]<(int)numLines[1]-2)))
	{
		return Operator_Multithread::GetCellCenterMaterialAvgCoord(pos, coord);
	}
	if ((pos[0]<0) || (pos[2]<0))
		return false;
	int l_pos[3] = {pos[0], 0, pos[2]};
	l_pos[1] = MapAlphaIndex2Range(pos[1]);
	return Operator_Multithread::GetCellCenterMaterialAvgCoord(l_pos, coord);
}

unsigned int Operator_Cylinder::SnapToMeshLine(int ny, double coord, bool &inside, bool dualMesh, bool fullMesh) const
{
	if (ny==1)
		coord=FitToAlphaRange(coord);
	return Operator_Multithread::SnapToMeshLine(ny, coord, inside, dualMesh, fullMesh);
}

int Operator_Cylinder::SnapBox2Mesh(const double* start, const double* stop, unsigned int* uiStart, unsigned int* uiStop, bool dualMesh, bool fullMesh, int SnapMethod, bool* bStartIn, bool* bStopIn) const
{
	double a_min = GetDiscLine(1,0);
	double a_max = GetDiscLine(1,GetNumberOfLines(1,fullMesh)-1);

	double a_size = stop[1] - start[1];
	double a_center = FitToAlphaRange(0.5*(stop[1]+start[1]));
	double a_start = a_center-a_size/2;
	double a_stop = a_start + a_size;
	if (a_stop>a_max)
		a_stop=a_max;
	if (a_stop<a_min)
		a_stop=a_min;
	if (a_start>a_max)
		a_start=a_max;
	if (a_start<a_min)
		a_start=a_min;

	double l_start[3] = {start[0], a_start, start[2]};
	double l_stop[3]  = {stop[0] , a_stop , stop[2] };
	return Operator_Multithread::SnapBox2Mesh(l_start, l_stop, uiStart, uiStop, dualMesh, fullMesh, SnapMethod, bStartIn, bStopIn);
}

int Operator_Cylinder::SnapLine2Mesh(const double* start, const double* stop, unsigned int* uiStart, unsigned int* uiStop, bool dualMesh, bool fullMesh) const
{
	int ret = Operator_Multithread::SnapLine2Mesh(start, stop, uiStart, uiStop, dualMesh, fullMesh);

	if ((stop[1]>start[1]) && (uiStop[1]<uiStart[1]) && (uiStop[1]==0))
		uiStop[1] = GetNumberOfLines(1, fullMesh)-1-(int)CC_closedAlpha;
	if ((stop[1]<start[1]) && (uiStop[1]>uiStart[1]) && (uiStop[1]==GetNumberOfLines(1, fullMesh)-1-(int)CC_closedAlpha))
		uiStop[1] = 0;

	return ret;
}


Grid_Path Operator_Cylinder::FindPath(double start[], double stop[])
{
	double l_start[3];
	double l_stop[3];

	for (int n=0;n<3;++n)
	{
		l_start[n] = start[n];
		l_stop[n] = stop[n];
	}

	while (fabs(l_stop[1]-l_start[1])>PI)
	{
		if (l_stop[1]>l_start[1])
			l_stop[1]-=2*PI;
		else
			l_stop[1]+=2*PI;
	}

	double help=0;
	if (l_start[1]>l_stop[1])
	{
		for (int n=0;n<3;++n)
		{
			help = l_start[n];
			l_start[n] = l_stop[n];
			l_stop[n] = help;
		}
	}

	double a_start = FitToAlphaRange(l_start[1]);
	double a_stop = FitToAlphaRange(l_stop[1]);

	if (a_stop >= a_start)
	{
		l_start[1] = a_start;
		l_stop[1] = a_stop;
		return Operator_Multithread::FindPath(l_start, l_stop);
	}

	// if a-stop fitted to disc range is now smaller than a-start, it must step over the a-bounds...

	Grid_Path path;
	for (int n=0;n<3;++n)
	{
		if ((l_start[n]<GetDiscLine(n,0)) && (l_stop[n]<GetDiscLine(n,0)))
			return path; //lower bound violation
		if ((l_start[n]>GetDiscLine(n,GetNumberOfLines(n,true)-1)) && (l_stop[n]>GetDiscLine(n,GetNumberOfLines(n,true)-1)))
			return path; //upper bound violation
	}

	if (g_settings.GetVerboseLevel()>2)
		cerr << __func__ << ": A path was leaving the alpha-direction mesh..." << endl;

	// this section comes into play, if the line moves over the angulare mesh-end/start
	// we try to have one part of the path on both "ends" of the mesh and stitch them together

	Grid_Path path1;
	Grid_Path path2;

	// calculate the intersection of the line with the a-max boundary
	double p0[3],p1[3],p2[3];
	for (int n=0;n<3;++n)
	{
		p0[n] = GetDiscLine(n,0);
		p1[n] = p0[n];
		p2[n] = p0[n];
	}
	p0[1] = GetDiscLine(1,GetNumberOfLines(1,true)-1-(int)CC_closedAlpha);
	p1[1] = p0[1];
	p2[1] = p0[1];
	p1[0] = discLines[0][numLines[0]-1];
	p2[2] = discLines[2][numLines[2]-1];

	TransformCoordSystem(p0,p0,m_MeshType,CARTESIAN);
	TransformCoordSystem(p1,p1,m_MeshType,CARTESIAN);
	TransformCoordSystem(p2,p2,m_MeshType,CARTESIAN);

	double c_start[3],c_stop[3];
	TransformCoordSystem(l_start,c_start,m_MeshType,CARTESIAN);
	TransformCoordSystem(l_stop,c_stop,m_MeshType,CARTESIAN);
	double intersect[3];
	double dist;
	int ret = LinePlaneIntersection(p0,p1,p2,c_start,c_stop,intersect,dist);
	if (ret<0)
	{
		cerr << __func__ << ": Error, unable to calculate intersection, this should not happen!" << endl;
		return path; // return empty path;
	}

	if (ret==0)
	{
		TransformCoordSystem(intersect,intersect,CARTESIAN,m_MeshType);
		intersect[1] = GetDiscLine(1,GetNumberOfLines(1,true)-1-(int)CC_closedAlpha);
		l_start[1] = FitToAlphaRange(l_start[1]);
		path1 = Operator::FindPath(l_start, intersect);
		if (g_settings.GetVerboseLevel()>2)
			cerr << __func__ << ": Intersection top: " << intersect[0] << "," << intersect[1] << "," << intersect[2] << endl;
	} //otherwise the path was not intersecting the upper a-bound...

	if (CC_closedAlpha==false)
	{
		for (int n=0;n<3;++n)
		{
			p0[n] = GetDiscLine(n,0);
			p1[n] = p0[n];
			p2[n] = p0[n];
		}
		p1[0] = discLines[0][numLines[0]-1];
		p2[2] = discLines[2][numLines[2]-1];

		TransformCoordSystem(p0,p0,m_MeshType,CARTESIAN);
		TransformCoordSystem(p1,p1,m_MeshType,CARTESIAN);
		TransformCoordSystem(p2,p2,m_MeshType,CARTESIAN);

		TransformCoordSystem(l_start,c_start,m_MeshType,CARTESIAN);
		TransformCoordSystem(l_stop,c_stop,m_MeshType,CARTESIAN);

		ret = LinePlaneIntersection(p0,p1,p2,c_start,c_stop,intersect,dist);
		TransformCoordSystem(intersect,intersect,CARTESIAN,m_MeshType);
	}

	if (ret==0)
	{
		intersect[1] = GetDiscLine(1,0);
		l_stop[1] = FitToAlphaRange(l_stop[1]);
		path2 = Operator::FindPath(intersect, l_stop);
		if (g_settings.GetVerboseLevel()>2)
			cerr << __func__ << ": Intersection bottom: " << intersect[0] << "," << intersect[1] << "," << intersect[2] << endl;
	}

	//combine path
	for (size_t t=0; t<path1.dir.size(); ++t)
	{
		path.posPath[0].push_back(path1.posPath[0].at(t));
		path.posPath[1].push_back(path1.posPath[1].at(t));
		path.posPath[2].push_back(path1.posPath[2].at(t));
		path.dir.push_back(path1.dir.at(t));
	}
	for (size_t t=0; t<path2.dir.size(); ++t)
	{
		path.posPath[0].push_back(path2.posPath[0].at(t));
		path.posPath[1].push_back(path2.posPath[1].at(t));
		path.posPath[2].push_back(path2.posPath[2].at(t));
		path.dir.push_back(path2.dir.at(t));
	}

	if (CC_closedAlpha==true)
		for (size_t t=0; t<path.dir.size(); ++t)
		{
			if ( ((path.dir.at(t)==0) || (path.dir.at(t)==2)) && (path.posPath[1].at(t)==0))
				path.posPath[1].at(t) = numLines[1]-2;
		}

	return path;
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
