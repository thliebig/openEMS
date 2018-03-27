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

#include <fstream>
#include <algorithm>
#include "operator.h"
#include "engine.h"
#include "extensions/operator_extension.h"
#include "extensions/operator_ext_excitation.h"
#include "Common/processfields.h"
#include "tools/array_ops.h"
#include "tools/vtk_file_writer.h"
#include "fparser.hh"
#include "extensions/operator_ext_excitation.h"

#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkPoints.h"
#include "vtkXMLPolyDataWriter.h"
#include "CSPrimBox.h"
#include "CSPrimCurve.h"

#include "CSPropMaterial.h"
#include "CSPropLumpedElement.h"

Operator* Operator::New()
{
	cout << "Create FDTD operator" << endl;
	Operator* op = new Operator();
	op->Init();
	return op;
}

Operator::Operator() : Operator_Base()
{
	m_Exc = 0;
	m_InvaildTimestep = false;
	m_TimeStepVar = 3;
}

Operator::~Operator()
{
	for (size_t n=0; n<m_Op_exts.size(); ++n)
		delete m_Op_exts.at(n);
	m_Op_exts.clear();

	Delete();
}

Engine* Operator::CreateEngine()
{
	m_Engine = Engine::New(this);
	return m_Engine;
}

void Operator::Init()
{
	CSX = NULL;
	m_Engine = NULL;

	Operator_Base::Init();

	vv=NULL;
	vi=NULL;
	iv=NULL;
	ii=NULL;

	m_epsR=NULL;
	m_kappa=NULL;
	m_mueR=NULL;
	m_sigma=NULL;

	MainOp=NULL;

	for (int n=0; n<3; ++n)
	{
		EC_C[n]=NULL;
		EC_G[n]=NULL;
		EC_L[n]=NULL;
		EC_R[n]=NULL;
	}

	m_Exc = 0;
	m_TimeStepFactor = 1;
	SetMaterialAvgMethod(QuarterCell);
}

void Operator::Delete()
{
	CSX = NULL;

	Delete_N_3DArray(vv,numLines);
	Delete_N_3DArray(vi,numLines);
	Delete_N_3DArray(iv,numLines);
	Delete_N_3DArray(ii,numLines);
	vv=vi=iv=ii=0;
	delete MainOp; MainOp=0;
	for (int n=0; n<3; ++n)
	{
		delete[] EC_C[n];EC_C[n]=0;
		delete[] EC_G[n];EC_G[n]=0;
		delete[] EC_L[n];EC_L[n]=0;
		delete[] EC_R[n];EC_R[n]=0;
	}

	Delete_N_3DArray(m_epsR,numLines);
	m_epsR=0;
	Delete_N_3DArray(m_kappa,numLines);
	m_kappa=0;
	Delete_N_3DArray(m_mueR,numLines);
	m_mueR=0;
	Delete_N_3DArray(m_sigma,numLines);
	m_sigma=0;
}

void Operator::Reset()
{
	Delete();
	Operator_Base::Reset();
}

double Operator::GetDiscLine(int n, unsigned int pos, bool dualMesh) const
{
	if ((n<0) || (n>2)) return 0.0;
	if (pos>=numLines[n]) return 0.0;
	if (dualMesh==false)
		return discLines[n][pos];

	// return dual mesh node
	if (pos<numLines[n]-1)
		return 0.5*(discLines[n][pos] + discLines[n][pos+1]);

	// dual node for the last line (outside the field domain)
	return discLines[n][pos] + 0.5*(discLines[n][pos] - discLines[n][pos-1]);
}

double Operator::GetDiscDelta(int n, unsigned int pos, bool dualMesh) const
{
	if ((n<0) || (n>2)) return 0.0;
	if (pos>=numLines[n]) return 0.0;
	double delta=0;
	if (dualMesh==false)
	{
		if (pos<numLines[n]-1)
			delta = GetDiscLine(n,pos+1,false) - GetDiscLine(n,pos,false);
		else
			delta = GetDiscLine(n,pos,false) - GetDiscLine(n,pos-1,false);
		return delta;
	}
	else
	{
		if (pos>0)
			delta = GetDiscLine(n,pos,true) - GetDiscLine(n,pos-1,true);
		else
			delta = GetDiscLine(n,1,false) - GetDiscLine(n,0,false);
		return delta;
	}
}

bool Operator::GetYeeCoords(int ny, unsigned int pos[3], double* coords, bool dualMesh) const
{
	for (int n=0;n<3;++n)
		coords[n]=GetDiscLine(n,pos[n],dualMesh);
	coords[ny]=GetDiscLine(ny,pos[ny],!dualMesh);

	//check if position is inside the FDTD domain
	if (dualMesh==false) //main grid
	{
		if (pos[ny]>=numLines[ny]-1)
			return false;
	}
	else	//dual grid
	{
		int nP = (ny+1)%3;
		int nPP = (ny+2)%3;
		if ((pos[nP]>=numLines[nP]-1) || (pos[nPP]>=numLines[nPP]-1))
			return false;
	}
	return true;
}

bool Operator::GetNodeCoords(const unsigned int pos[3], double* coords, bool dualMesh, CoordinateSystem c_system) const
{
	for (int n=0;n<3;++n)
		coords[n]=GetDiscLine(n,pos[n],dualMesh);
	TransformCoordSystem(coords,coords,m_MeshType,c_system);
	return true;
}

double Operator::GetEdgeLength(int n, const unsigned int* pos, bool dualMesh) const
{
	return GetDiscDelta(n,pos[n],dualMesh)*gridDelta;
}

double Operator::GetCellVolume(const unsigned int pos[3], bool dualMesh) const
{
	double vol=1;
	for (int n=0;n<3;++n)
		vol*=GetEdgeLength(n,pos,dualMesh);
	return vol;
}

double Operator::GetNodeWidth(int ny, const int pos[3], bool dualMesh) const
{
	if ( (pos[0]<0) || (pos[1]<0) || (pos[2]<0) )
		return 0.0;

	//call the unsigned int version of GetNodeWidth
	unsigned int uiPos[]={(unsigned int)pos[0],(unsigned int)pos[1],(unsigned int)pos[2]};
	return GetNodeWidth(ny, uiPos, dualMesh);
}

double Operator::GetNodeArea(int ny, const unsigned int pos[3], bool dualMesh) const
{
	int nyP = (ny+1)%3;
	int nyPP = (ny+2)%3;
	return GetNodeWidth(nyP,pos,dualMesh) * GetNodeWidth(nyPP,pos,dualMesh);
}

double Operator::GetNodeArea(int ny, const int pos[3], bool dualMesh) const
{
	if ( (pos[0]<0) || (pos[1]<0) || (pos[2]<0) )
		return 0.0;

	//call the unsigned int version of GetNodeArea
	unsigned int uiPos[]={(unsigned int)pos[0],(unsigned int)pos[1],(unsigned int)pos[2]};
	return GetNodeArea(ny, uiPos, dualMesh);
}

unsigned int Operator::SnapToMeshLine(int ny, double coord, bool &inside, bool dualMesh, bool fullMesh) const
{
	inside = false;
	if ((ny<0) || (ny>2))
		return 0;
	if (coord<GetDiscLine(ny,0))
		return 0;
	unsigned int numLines = GetNumberOfLines(ny, fullMesh);
	if (coord>GetDiscLine(ny,numLines-1))
		return numLines-1;
	inside=true;
	if (dualMesh==false)
	{
		for (unsigned int n=0;n<numLines;++n)
		{
			if (coord<=GetDiscLine(ny,n,true))
				return n;
		}
	}
	else
	{
		for (unsigned int n=1;n<numLines;++n)
		{
			if (coord<=GetDiscLine(ny,n,false))
				return n-1;
		}
	}
	//should not happen
	return 0;
}

bool Operator::SnapToMesh(const double* dcoord, unsigned int* uicoord, bool dualMesh, bool fullMesh, bool* inside) const
{
	bool meshInside=false;
	bool ok=true;
	for (int n=0; n<3; ++n)
	{
		uicoord[n] = SnapToMeshLine(n,dcoord[n],meshInside,dualMesh,fullMesh);
		ok &= meshInside;
		if (inside)
			inside[n]=meshInside;
	}
	return ok;
}

int Operator::SnapBox2Mesh(const double* start, const double* stop, unsigned int* uiStart, unsigned int* uiStop, bool dualMesh, bool fullMesh, int SnapMethod, bool* bStartIn, bool* bStopIn) const
{
	double l_start[3], l_stop[3];
	for (int n=0;n<3;++n)
	{
		l_start[n] = fmin(start[n],stop[n]);
		l_stop[n] = fmax(start[n], stop[n]);
		double min = GetDiscLine(n,0);
		double max = GetDiscLine(n,GetNumberOfLines(n, fullMesh)-1);
		if ( ((l_start[n]<min) && (l_stop[n]<min)) || ((l_start[n]>max) && (l_stop[n]>max)) )
		{
			return -2;
		}
	}

	SnapToMesh(l_start, uiStart, dualMesh, fullMesh, bStartIn);
	SnapToMesh(l_stop, uiStop, dualMesh, fullMesh, bStopIn);
	int iDim = 0;

	if (SnapMethod==0)
	{
		for (int n=0;n<3;++n)
			if (uiStop[n]>uiStart[n])
				++iDim;
		return iDim;
	}
	else if (SnapMethod==1)
	{
		for (int n=0;n<3;++n)
		{
			if (uiStop[n]>uiStart[n])
			{
				if ((GetDiscLine( n, uiStart[n], dualMesh ) > l_start[n]) && (uiStart[n]>0))
					--uiStart[n];
				if ((GetDiscLine( n, uiStop[n], dualMesh ) < l_stop[n]) && (uiStop[n]<GetNumberOfLines(n, fullMesh)-1))
					++uiStop[n];
			}
			if (uiStop[n]>uiStart[n])
				++iDim;
		}
		return iDim;
	}
	else if (SnapMethod==2)
	{
		for (int n=0;n<3;++n)
		{
			if (uiStop[n]>uiStart[n])
			{
				if ((GetDiscLine( n, uiStart[n], dualMesh ) < l_start[n]) && (uiStart[n]<GetNumberOfLines(n, fullMesh)-1))
					++uiStart[n];
				if ((GetDiscLine( n, uiStop[n], dualMesh ) > l_stop[n]) && (uiStop[n]>0))
					--uiStop[n];
			}
			if (uiStop[n]>uiStart[n])
				++iDim;
		}
		return iDim;
	}
	else
		cerr << "Operator::SnapBox2Mesh: Unknown snapping method!" << endl;
	return -1;
}

int Operator::SnapLine2Mesh(const double* start, const double* stop, unsigned int* uiStart, unsigned int* uiStop, bool dualMesh, bool fullMesh) const
{
	bool bStartIn[3];
	bool bStopIn[3];
	SnapToMesh(start, uiStart, dualMesh, fullMesh, bStartIn);
	SnapToMesh(stop, uiStop, dualMesh, fullMesh, bStopIn);

	for (int n=0;n<3;++n)
	{
		if ((start[n]<GetDiscLine(n,0)) && (stop[n]<GetDiscLine(n,0)))
			return -1; //lower bound violation
		if ((start[n]>GetDiscLine(n,GetNumberOfLines(n,true)-1)) && (stop[n]>GetDiscLine(n,GetNumberOfLines(n,true)-1)))
			return -1; //upper bound violation
	}

	int ret = 0;
	if (!(bStartIn[0] && bStartIn[1] && bStartIn[2]))
		ret = ret + 1;
	if (!(bStopIn[0] && bStopIn[1] && bStopIn[2]))
		ret = ret + 2;
	if (ret==0)
		return ret;

	//fixme, do we need to do something about start or stop being outside the field domain?
	//maybe caclulate the intersection point and snap to that?
	//it seems to work like this as well...

	return ret;
}


Grid_Path Operator::FindPath(double start[], double stop[])
{
	Grid_Path path;
	unsigned int uiStart[3],uiStop[3],currPos[3];

	int ret = SnapLine2Mesh(start, stop, uiStart, uiStop, false, true);
	if (ret<0)
		return path;

	currPos[0]=uiStart[0];
	currPos[1]=uiStart[1];
	currPos[2]=uiStart[2];
	double meshStart[3] = {discLines[0][uiStart[0]], discLines[1][uiStart[1]], discLines[2][uiStart[2]]};
	double meshStop[3] = {discLines[0][uiStop[0]], discLines[1][uiStop[1]], discLines[2][uiStop[2]]};
	bool UpDir = false;
	double foot=0,dist=0,minFoot=0,minDist=0;
	int minDir=0;
	unsigned int minPos[3];
	double startFoot,stopFoot,currFoot;
	Point_Line_Distance(meshStart,start,stop,startFoot,dist, m_MeshType);
	Point_Line_Distance(meshStop,start,stop,stopFoot,dist, m_MeshType);
	currFoot=startFoot;
	minFoot=startFoot;
	double P[3];

	while (minFoot<stopFoot)
	{
		minDist=1e300;
		for (int n=0; n<3; ++n) //check all 6 surrounding points
		{
			P[0] = discLines[0][currPos[0]];
			P[1] = discLines[1][currPos[1]];
			P[2] = discLines[2][currPos[2]];
			if (((int)currPos[n]-1)>=0)
			{
				P[n] = discLines[n][currPos[n]-1];
				Point_Line_Distance(P,start,stop,foot,dist, m_MeshType);
				if ((foot>currFoot) && (dist<minDist))
				{
					minFoot=foot;
					minDist=dist;
					minDir = n;
					UpDir = false;
				}
			}
			if ((currPos[n]+1)<numLines[n])
			{
				P[n] = discLines[n][currPos[n]+1];
				Point_Line_Distance(P,start,stop,foot,dist, m_MeshType);
				if ((foot>currFoot) && (dist<minDist))
				{
					minFoot=foot;
					minDist=dist;
					minDir = n;
					UpDir = true;
				}
			}
		}
		minPos[0]=currPos[0];
		minPos[1]=currPos[1];
		minPos[2]=currPos[2];
		if (UpDir)
		{
			currPos[minDir]+=1;
		}
		else
		{
			currPos[minDir]+=-1;
			minPos[minDir]-=1;
		}
		//check validity of current postion
		for (int n=0;n<3;++n)
			if (currPos[n]>=numLines[n])
			{
				cerr << __func__ << ": Error, path went out of simulation domain, skipping path!" << endl;
				Grid_Path empty;
				return empty;
			}
		path.posPath[0].push_back(minPos[0]);
		path.posPath[1].push_back(minPos[1]);
		path.posPath[2].push_back(minPos[2]);
		currFoot=minFoot;
		path.dir.push_back(minDir);
	}

	//close missing edges, if currPos is not equal to uiStopPos
	for (int n=0; n<3; ++n)
	{
		if (currPos[n]>uiStop[n])
		{
			--currPos[n];
			path.posPath[0].push_back(currPos[0]);
			path.posPath[1].push_back(currPos[1]);
			path.posPath[2].push_back(currPos[2]);
			path.dir.push_back(n);
		}
		else if (currPos[n]<uiStop[n])
		{
			path.posPath[0].push_back(currPos[0]);
			path.posPath[1].push_back(currPos[1]);
			path.posPath[2].push_back(currPos[2]);
			path.dir.push_back(n);
		}
	}

	return path;
}

void Operator::SetMaterialAvgMethod(MatAverageMethods method)
{
	switch (method)
	{
	default:
	case QuarterCell:
		return SetQuarterCellMaterialAvg();
	case CentralCell:
		return SetCellConstantMaterial();
	}
}

double Operator::GetNumberCells() const
{
	if (numLines)
		return (numLines[0])*(numLines[1])*(numLines[2]); //it's more like number of nodes???
	return 0;
}

void Operator::ShowStat() const
{
	unsigned int OpSize = 12*numLines[0]*numLines[1]*numLines[2]*sizeof(FDTD_FLOAT);
	unsigned int FieldSize = 6*numLines[0]*numLines[1]*numLines[2]*sizeof(FDTD_FLOAT);
	double MBdiff = 1024*1024;

	cout << "------- Stat: FDTD Operator -------" << endl;
	cout << "Dimensions\t\t: " << numLines[0] << "x" << numLines[1] << "x" << numLines[2] << " = " <<  numLines[0]*numLines[1]*numLines[2] << " Cells (" << numLines[0]*numLines[1]*numLines[2]/1e6 << " MCells)" << endl;
	cout << "Size of Operator\t: " << OpSize << " Byte (" << (double)OpSize/MBdiff << " MiB) " << endl;
	cout << "Size of Field-Data\t: " << FieldSize << " Byte (" << (double)FieldSize/MBdiff << " MiB) " << endl;
	cout << "-----------------------------------" << endl;
	cout << "Background materials (epsR/mueR/kappa/sigma): " << GetBackgroundEpsR() << "/" << GetBackgroundMueR() << "/" << GetBackgroundKappa() << "/" << GetBackgroundSigma() << endl;
	cout << "-----------------------------------" << endl;
	cout << "Number of PEC edges\t: " << m_Nr_PEC[0]+m_Nr_PEC[1]+m_Nr_PEC[2] << endl;
	cout << "in " << GetDirName(0) << " direction\t\t: " << m_Nr_PEC[0] << endl;
	cout << "in " << GetDirName(1) << " direction\t\t: " << m_Nr_PEC[1] << endl;
	cout << "in " << GetDirName(2) << " direction\t\t: " << m_Nr_PEC[2] << endl;
	cout << "-----------------------------------" << endl;
	cout << "Timestep (s)\t\t: " << dT ;
	if (opt_dT)
		cout <<"\t(" << opt_dT << ")";
	cout << endl;
	cout << "Timestep method name\t: " << m_Used_TS_Name << endl;
	cout << "Nyquist criteria (TS)\t: " << m_Exc->GetNyquistNum() << endl;
	cout << "Nyquist criteria (s)\t: " << m_Exc->GetNyquistNum()*dT << endl;
	cout << "-----------------------------------" << endl;
}

void Operator::ShowExtStat() const
{
	if (m_Op_exts.size()==0) return;
	cout << "-----------------------------------" << endl;
	for (size_t n=0; n<m_Op_exts.size(); ++n)
		m_Op_exts.at(n)->ShowStat(cout);
	cout << "-----------------------------------" << endl;
}

void Operator::DumpOperator2File(string filename)
{
#ifdef OUTPUT_IN_DRAWINGUNITS
	double discLines_scaling = 1;
#else
	double discLines_scaling = GetGridDelta();
#endif

	cout << "Operator: Dumping FDTD operator information to vtk file: " << filename << " ..." << flush;

	VTK_File_Writer* vtk_Writer = new VTK_File_Writer(filename.c_str(), m_MeshType);
	vtk_Writer->SetMeshLines(discLines,numLines,discLines_scaling);
	vtk_Writer->SetHeader("openEMS - Operator dump");

	vtk_Writer->SetNativeDump(true);

	//find excitation extension
	Operator_Ext_Excitation* Op_Ext_Exc=GetExcitationExtension();

	if (Op_Ext_Exc)
	{
		FDTD_FLOAT**** exc = NULL;
		if (Op_Ext_Exc->Volt_Count>0)
		{
			exc = Create_N_3DArray<FDTD_FLOAT>(numLines);
			for (unsigned int n=0; n<  Op_Ext_Exc->Volt_Count; ++n)
				exc[  Op_Ext_Exc->Volt_dir[n]][  Op_Ext_Exc->Volt_index[0][n]][  Op_Ext_Exc->Volt_index[1][n]][  Op_Ext_Exc->Volt_index[2][n]] =   Op_Ext_Exc->Volt_amp[n];
			vtk_Writer->AddVectorField("exc_volt",exc);
			Delete_N_3DArray(exc,numLines);
		}

		if (  Op_Ext_Exc->Curr_Count>0)
		{
			exc = Create_N_3DArray<FDTD_FLOAT>(numLines);
			for (unsigned int n=0; n<  Op_Ext_Exc->Curr_Count; ++n)
				exc[  Op_Ext_Exc->Curr_dir[n]][  Op_Ext_Exc->Curr_index[0][n]][  Op_Ext_Exc->Curr_index[1][n]][  Op_Ext_Exc->Curr_index[2][n]] =   Op_Ext_Exc->Curr_amp[n];
			vtk_Writer->AddVectorField("exc_curr",exc);
			Delete_N_3DArray(exc,numLines);
		}
	}

	FDTD_FLOAT**** vv_temp = Create_N_3DArray<FDTD_FLOAT>(numLines);
	FDTD_FLOAT**** vi_temp = Create_N_3DArray<FDTD_FLOAT>(numLines);
	FDTD_FLOAT**** iv_temp = Create_N_3DArray<FDTD_FLOAT>(numLines);
	FDTD_FLOAT**** ii_temp = Create_N_3DArray<FDTD_FLOAT>(numLines);

	unsigned int pos[3], n;
	for (n=0; n<3; n++)
		for (pos[0]=0; pos[0]<numLines[0]; pos[0]++)
			for (pos[1]=0; pos[1]<numLines[1]; pos[1]++)
				for (pos[2]=0; pos[2]<numLines[2]; pos[2]++)
				{
					vv_temp[n][pos[0]][pos[1]][pos[2]] = GetVV(n,pos);
					vi_temp[n][pos[0]][pos[1]][pos[2]] = GetVI(n,pos);
					iv_temp[n][pos[0]][pos[1]][pos[2]] = GetIV(n,pos);
					ii_temp[n][pos[0]][pos[1]][pos[2]] = GetII(n,pos);
				}


	vtk_Writer->AddVectorField("vv",vv_temp);
	Delete_N_3DArray(vv_temp,numLines);
	vtk_Writer->AddVectorField("vi",vi_temp);
	Delete_N_3DArray(vi_temp,numLines);
	vtk_Writer->AddVectorField("iv",iv_temp);
	Delete_N_3DArray(iv_temp,numLines);
	vtk_Writer->AddVectorField("ii",ii_temp);
	Delete_N_3DArray(ii_temp,numLines);

	if (vtk_Writer->Write()==false)
		cerr << "Operator::DumpOperator2File: Error: Can't write file... skipping!" << endl;

	delete vtk_Writer;
}

//! \brief dump PEC (perfect electric conductor) information (into VTK-file)
//! visualization via paraview
//! visualize only one component (x, y or z)
void Operator::DumpPEC2File(string filename , unsigned int *range)
{
	cout << "Operator: Dumping PEC information to vtk file: " << filename << " ..." << flush;

#ifdef OUTPUT_IN_DRAWINGUNITS
	double scaling = 1.0;
#else
	double scaling = GetGridDelta();;
#endif

	unsigned int start[3] = {0, 0, 0};
	unsigned int stop[3] = {numLines[0]-1,numLines[1]-1,numLines[2]-1};

	if (range!=NULL)
		for (int n=0;n<3;++n)
		{
			start[n] = range[2*n];
			stop[n] = range[2*n+1];
		}

	vtkPolyData* polydata = vtkPolyData::New();
	vtkCellArray *poly = vtkCellArray::New();
	vtkPoints *points = vtkPoints::New();

	int* pointIdx[2];
	pointIdx[0] = new int[numLines[0]*numLines[1]];
	pointIdx[1] = new int[numLines[0]*numLines[1]];
	// init point idx
	for (unsigned int n=0;n<numLines[0]*numLines[1];++n)
	{
		pointIdx[0][n]=-1;
		pointIdx[1][n]=-1;
	}

	int nP,nPP;
	double coord[3];
	unsigned int pos[3],rpos[3];
	unsigned int mesh_idx=0;
	for (pos[2]=start[2];pos[2]<stop[2];++pos[2])
	{ // each xy-plane
		for (unsigned int n=0;n<numLines[0]*numLines[1];++n)
		{
			pointIdx[0][n]=pointIdx[1][n];
			pointIdx[1][n]=-1;
		}
		for (pos[0]=start[0];pos[0]<stop[0];++pos[0])
			for (pos[1]=start[1];pos[1]<stop[1];++pos[1])
			{
				for (int n=0;n<3;++n)
				{
					nP = (n+1)%3;
					nPP = (n+2)%3;
					if ((GetVV(n,pos) == 0) && (GetVI(n,pos) == 0) && (pos[nP]>0) && (pos[nPP]>0))
					{
						rpos[0]=pos[0];
						rpos[1]=pos[1];
						rpos[2]=pos[2];

						poly->InsertNextCell(2);

						mesh_idx = rpos[0] + rpos[1]*numLines[0];
						if (pointIdx[0][mesh_idx]<0)
						{
							for (int m=0;m<3;++m)
								coord[m] = discLines[m][rpos[m]];
							TransformCoordSystem(coord, coord, m_MeshType, CARTESIAN);
							for (int m=0;m<3;++m)
								coord[m] *= scaling;
							pointIdx[0][mesh_idx] = (int)points->InsertNextPoint(coord);
						}
						poly->InsertCellPoint(pointIdx[0][mesh_idx]);

						++rpos[n];
						mesh_idx = rpos[0] + rpos[1]*numLines[0];
						if (pointIdx[n==2][mesh_idx]<0)
						{
							for (int m=0;m<3;++m)
								coord[m] = discLines[m][rpos[m]];
							TransformCoordSystem(coord, coord, m_MeshType, CARTESIAN);
							for (int m=0;m<3;++m)
								coord[m] *= scaling;
							pointIdx[n==2][mesh_idx] = (int)points->InsertNextPoint(coord);
						}
						poly->InsertCellPoint(pointIdx[n==2][mesh_idx]);
					}
				}
			}
	}
	delete[] pointIdx[0];
	delete[] pointIdx[1];

	polydata->SetPoints(points);
	points->Delete();
	polydata->SetLines(poly);
	poly->Delete();

	vtkXMLPolyDataWriter* writer  = vtkXMLPolyDataWriter::New();
	filename += ".vtp";
	writer->SetFileName(filename.c_str());

#if VTK_MAJOR_VERSION>=6
	writer->SetInputData(polydata);
#else
	writer->SetInput(polydata);
#endif
	writer->Write();

	writer->Delete();
	polydata->Delete();
	cout << " done." << endl;
}

void Operator::DumpMaterial2File(string filename)
{
#ifdef OUTPUT_IN_DRAWINGUNITS
	double discLines_scaling = 1;
#else
	double discLines_scaling = GetGridDelta();
#endif

	cout << "Operator: Dumping material information to vtk file: " << filename << " ..."  << flush;

	FDTD_FLOAT**** epsilon = Create_N_3DArray<FDTD_FLOAT>(numLines);
	FDTD_FLOAT**** mue     = Create_N_3DArray<FDTD_FLOAT>(numLines);
	FDTD_FLOAT**** kappa   = Create_N_3DArray<FDTD_FLOAT>(numLines);
	FDTD_FLOAT**** sigma   = Create_N_3DArray<FDTD_FLOAT>(numLines);

	unsigned int pos[3];
	for (pos[0]=0; pos[0]<numLines[0]; ++pos[0])
	{
		for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
		{
			vector<CSPrimitives*> vPrims = this->GetPrimitivesBoundBox(pos[0], pos[1], -1, CSProperties::MATERIAL);
			for (pos[2]=0; pos[2]<numLines[2]; ++pos[2])
			{
				for (int n=0; n<3; ++n)
				{
					double inMat[4];
					Calc_EffMatPos(n, pos, inMat, vPrims);
					epsilon[n][pos[0]][pos[1]][pos[2]] = inMat[0]/__EPS0__;
					mue[n][pos[0]][pos[1]][pos[2]]     = inMat[2]/__MUE0__;
					kappa[n][pos[0]][pos[1]][pos[2]]   = inMat[1];
					sigma[n][pos[0]][pos[1]][pos[2]]   = inMat[3];
				}
			}
		}
	}

	VTK_File_Writer* vtk_Writer = new VTK_File_Writer(filename.c_str(), m_MeshType);
	vtk_Writer->SetMeshLines(discLines,numLines,discLines_scaling);
	vtk_Writer->SetHeader("openEMS - material dump");

	vtk_Writer->SetNativeDump(true);

	vtk_Writer->AddVectorField("epsilon",epsilon);
	Delete_N_3DArray(epsilon,numLines);
	vtk_Writer->AddVectorField("mue",mue);
	Delete_N_3DArray(mue,numLines);
	vtk_Writer->AddVectorField("kappa",kappa);
	Delete_N_3DArray(kappa,numLines);
	vtk_Writer->AddVectorField("sigma",sigma);
	Delete_N_3DArray(sigma,numLines);

	if (vtk_Writer->Write()==false)
		cerr << "Operator::DumpMaterial2File: Error: Can't write file... skipping!" << endl;

	delete vtk_Writer;
}

 bool Operator::SetupCSXGrid(CSRectGrid* grid)
 {
	 for (int n=0; n<3; ++n)
	 {
		 discLines[n] = grid->GetLines(n,discLines[n],numLines[n],true);
		 if (numLines[n]<3)
		 {
			 cerr << "CartOperator::SetupCSXGrid: you need at least 3 disc-lines in every direction (3D!)!!!" << endl;
			 Reset();
			 return false;
		 }
	 }
	 MainOp = new AdrOp(numLines[0],numLines[1],numLines[2]);
	 MainOp->SetGrid(discLines[0],discLines[1],discLines[2]);
	 if (grid->GetDeltaUnit()<=0)
	 {
		 cerr << "CartOperator::SetupCSXGrid: grid delta unit must not be <=0 !!!" << endl;
		 Reset();
		 return false;
	 }
	 else gridDelta=grid->GetDeltaUnit();
	 MainOp->SetGridDelta(1);
	 MainOp->AddCellAdrOp();

	 //delete the grid clone...
	 delete grid;
	 return true;
 }

bool Operator::SetGeometryCSX(ContinuousStructure* geo)
{
	if (geo==NULL) return false;

	CSX = geo;

	CSBackgroundMaterial* bg_mat=CSX->GetBackgroundMaterial();
	SetBackgroundEpsR(bg_mat->GetEpsilon());
	SetBackgroundMueR(bg_mat->GetMue());
	SetBackgroundKappa(bg_mat->GetKappa());
	SetBackgroundSigma(bg_mat->GetSigma());
	SetBackgroundDensity(0);

	CSRectGrid* grid=CSX->GetGrid();
	return SetupCSXGrid(CSRectGrid::Clone(grid));
}

void Operator::InitOperator()
{
	Delete_N_3DArray(vv,numLines);
	Delete_N_3DArray(vi,numLines);
	Delete_N_3DArray(iv,numLines);
	Delete_N_3DArray(ii,numLines);
	vv = Create_N_3DArray<FDTD_FLOAT>(numLines);
	vi = Create_N_3DArray<FDTD_FLOAT>(numLines);
	iv = Create_N_3DArray<FDTD_FLOAT>(numLines);
	ii = Create_N_3DArray<FDTD_FLOAT>(numLines);
}

void Operator::InitDataStorage()
{
	if (m_StoreMaterial[0])
	{
		if (g_settings.GetVerboseLevel()>0)
			cerr << "Operator::InitDataStorage(): Storing epsR material data..." << endl;
		Delete_N_3DArray(m_epsR,numLines);
		m_epsR = Create_N_3DArray<float>(numLines);
	}
	if (m_StoreMaterial[1])
	{
		if (g_settings.GetVerboseLevel()>0)
			cerr << "Operator::InitDataStorage(): Storing kappa material data..." << endl;
		Delete_N_3DArray(m_kappa,numLines);
		m_kappa = Create_N_3DArray<float>(numLines);
	}
	if (m_StoreMaterial[2])
	{
		if (g_settings.GetVerboseLevel()>0)
			cerr << "Operator::InitDataStorage(): Storing muR material data..." << endl;
		Delete_N_3DArray(m_mueR,numLines);
		m_mueR = Create_N_3DArray<float>(numLines);
	}
	if (m_StoreMaterial[3])
	{
		if (g_settings.GetVerboseLevel()>0)
			cerr << "Operator::InitDataStorage(): Storing sigma material data..." << endl;
		Delete_N_3DArray(m_sigma,numLines);
		m_sigma = Create_N_3DArray<float>(numLines);
	}
}

void Operator::CleanupMaterialStorage()
{
	if (!m_StoreMaterial[0] && m_epsR)
	{
		if (g_settings.GetVerboseLevel()>0)
			cerr << "Operator::CleanupMaterialStorage(): Delete epsR material data..." << endl;
		Delete_N_3DArray(m_epsR,numLines);
		m_epsR = NULL;
	}
	if (!m_StoreMaterial[1] && m_kappa)
	{
		if (g_settings.GetVerboseLevel()>0)
			cerr << "Operator::CleanupMaterialStorage(): Delete kappa material data..." << endl;
		Delete_N_3DArray(m_kappa,numLines);
		m_kappa = NULL;
	}
	if (!m_StoreMaterial[2] && m_mueR)
	{
		if (g_settings.GetVerboseLevel()>0)
			cerr << "Operator::CleanupMaterialStorage(): Delete mueR material data..." << endl;
		Delete_N_3DArray(m_mueR,numLines);
		m_mueR = NULL;
	}
	if (!m_StoreMaterial[3] && m_sigma)
	{
		if (g_settings.GetVerboseLevel()>0)
			cerr << "Operator::CleanupMaterialStorage(): Delete sigma material data..." << endl;
		Delete_N_3DArray(m_sigma,numLines);
		m_sigma = NULL;
	}
}

double Operator::GetDiscMaterial(int type, int n, const unsigned int pos[3]) const
{
	switch (type)
	{
	case 0:
		if (m_epsR==0)
			return 0;
		return m_epsR[n][pos[0]][pos[1]][pos[2]];
	case 1:
		if (m_kappa==0)
			return 0;
		return m_kappa[n][pos[0]][pos[1]][pos[2]];
	case 2:
		if (m_mueR==0)
			return 0;
		return m_mueR[n][pos[0]][pos[1]][pos[2]];
	case 3:
		if (m_sigma==0)
			return 0;
		return m_sigma[n][pos[0]][pos[1]][pos[2]];
	}
	return 0;
}

void Operator::SetExcitationSignal(Excitation* exc)
{
	m_Exc=exc;
}

void Operator::Calc_ECOperatorPos(int n, unsigned int* pos)
{
	unsigned int i = MainOp->SetPos(pos[0],pos[1],pos[2]);
	double C = EC_C[n][i];
	double G = EC_G[n][i];
	if (C>0)
	{
		SetVV(n,pos[0],pos[1],pos[2], (1.0-dT*G/2.0/C)/(1.0+dT*G/2.0/C) );
		SetVI(n,pos[0],pos[1],pos[2], (dT/C)/(1.0+dT*G/2.0/C) );
	}
	else
	{
		SetVV(n,pos[0],pos[1],pos[2], 0 );
		SetVI(n,pos[0],pos[1],pos[2], 0 );
	}

	double L = EC_L[n][i];
	double R = EC_R[n][i];
	if (L>0)
	{
		SetII(n,pos[0],pos[1],pos[2], (1.0-dT*R/2.0/L)/(1.0+dT*R/2.0/L) );
		SetIV(n,pos[0],pos[1],pos[2], (dT/L)/(1.0+dT*R/2.0/L) );
	}
	else
	{
		SetII(n,pos[0],pos[1],pos[2], 0 );
		SetIV(n,pos[0],pos[1],pos[2], 0 );
	}
}

int Operator::CalcECOperator( DebugFlags debugFlags )
{
	Init_EC();
	InitDataStorage();

	if (Calc_EC()==0)
		return -1;

	m_InvaildTimestep = false;
	opt_dT = 0;
	if (dT>0)
	{
		double save_dT = dT;
		CalcTimestep();
		opt_dT = dT;
		if (dT<save_dT)
		{
			cerr << "Operator::CalcECOperator: Warning, forced timestep: " << save_dT << "s is larger than calculated timestep: " << dT << "s! It is not recommended using this timestep!! " << endl;
			m_InvaildTimestep = true;
		}

		dT = save_dT;
	}
	else
		CalcTimestep();

	dT*=m_TimeStepFactor;

	if (m_Exc->GetSignalPeriod()>0)
	{
		unsigned int TS = ceil(m_Exc->GetSignalPeriod()/dT);
		double new_dT = m_Exc->GetSignalPeriod()/TS;
		cout << "Operartor::CalcECOperator: Decreasing timestep by " << round((dT-new_dT)/dT*1000)/10.0 << "% to " << new_dT << " (" << dT << ") to match periodic signal" << endl;
		dT = new_dT;
	}

	m_Exc->Reset(dT);

	InitOperator();

	unsigned int pos[3];

	for (int n=0; n<3; ++n)
	{
		for (pos[0]=0; pos[0]<numLines[0]; ++pos[0])
		{
			for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
			{
				for (pos[2]=0; pos[2]<numLines[2]; ++pos[2])
				{
					Calc_ECOperatorPos(n,pos);
				}
			}
		}
	}

	//Apply PEC to all boundary's
	bool PEC[6]={1,1,1,1,1,1};
	//make an exception for BC == -1
	for (int n=0; n<6; ++n)
		if ((m_BC[n]==-1))
			PEC[n] = false;
	ApplyElectricBC(PEC);

	CalcPEC();

	Calc_LumpedElements();

	bool PMC[6];
	for (int n=0; n<6; ++n)
		PMC[n] = m_BC[n]==1;
	ApplyMagneticBC(PMC);

	//all information available for extension... create now...
	for (size_t n=0; n<m_Op_exts.size(); ++n)
		m_Op_exts.at(n)->BuildExtension();

	//remove inactive extensions
	vector<Operator_Extension*>::iterator it = m_Op_exts.begin();
	while (it!=m_Op_exts.end())
	{
		if ( (*it)->IsActive() == false)
		{
			DeleteExtension((*it));
			it = m_Op_exts.begin(); //restart search for inactive extension
		}
		else
			++it;
	}

	if (debugFlags & debugMaterial)
		DumpMaterial2File( "material_dump" );
	if (debugFlags & debugOperator)
		DumpOperator2File( "operator_dump" );
	if (debugFlags & debugPEC)
		DumpPEC2File( "PEC_dump" );

	//cleanup
	for (int n=0; n<3; ++n)
	{
		delete[] EC_C[n];
		EC_C[n]=NULL;
		delete[] EC_G[n];
		EC_G[n]=NULL;
		delete[] EC_L[n];
		EC_L[n]=NULL;
		delete[] EC_R[n];
		EC_R[n]=NULL;
	}

	return 0;
}

void Operator::ApplyElectricBC(bool* dirs)
{
	if (!dirs)
		return;

	unsigned int pos[3];
	for (int n=0; n<3; ++n)
	{
		int nP = (n+1)%3;
		int nPP = (n+2)%3;
		for (pos[nP]=0; pos[nP]<numLines[nP]; ++pos[nP])
		{
			for (pos[nPP]=0; pos[nPP]<numLines[nPP]; ++pos[nPP])
			{
				if (dirs[2*n])
				{
					// set to PEC
					pos[n] = 0;
					SetVV(nP, pos[0],pos[1],pos[2], 0 );
					SetVI(nP, pos[0],pos[1],pos[2], 0 );
					SetVV(nPP,pos[0],pos[1],pos[2], 0 );
					SetVI(nPP,pos[0],pos[1],pos[2], 0 );
				}

				if (dirs[2*n+1])
				{
					// set to PEC
					pos[n] = numLines[n]-1;
					SetVV(n,  pos[0],pos[1],pos[2], 0 ); // these are outside the FDTD-domain as defined by the main disc
					SetVI(n,  pos[0],pos[1],pos[2], 0 ); // these are outside the FDTD-domain as defined by the main disc

					SetVV(nP, pos[0],pos[1],pos[2], 0 );
					SetVI(nP, pos[0],pos[1],pos[2], 0 );
					SetVV(nPP,pos[0],pos[1],pos[2], 0 );
					SetVI(nPP,pos[0],pos[1],pos[2], 0 );
				}
			}
		}
	}
}

void Operator::ApplyMagneticBC(bool* dirs)
{
	if (!dirs)
		return;

	unsigned int pos[3];
	for (int n=0; n<3; ++n)
	{
		int nP = (n+1)%3;
		int nPP = (n+2)%3;
		for (pos[nP]=0; pos[nP]<numLines[nP]; ++pos[nP])
		{
			for (pos[nPP]=0; pos[nPP]<numLines[nPP]; ++pos[nPP])
			{
				if (dirs[2*n])
				{
					// set to PMC
					pos[n] = 0;
					SetII(n,  pos[0],pos[1],pos[2], 0 );
					SetIV(n,  pos[0],pos[1],pos[2], 0 );
					SetII(nP, pos[0],pos[1],pos[2], 0 );
					SetIV(nP, pos[0],pos[1],pos[2], 0 );
					SetII(nPP,pos[0],pos[1],pos[2], 0 );
					SetIV(nPP,pos[0],pos[1],pos[2], 0 );
				}

				if (dirs[2*n+1])
				{
					// set to PMC
					pos[n] = numLines[n]-2;
					SetII(nP, pos[0],pos[1],pos[2], 0 );
					SetIV(nP, pos[0],pos[1],pos[2], 0 );
					SetII(nPP,pos[0],pos[1],pos[2], 0 );
					SetIV(nPP,pos[0],pos[1],pos[2], 0 );
				}

				// the last current lines are outside the FDTD domain and cannot be iterated by the FDTD engine
				pos[n] = numLines[n]-1;
				SetII(n,  pos[0],pos[1],pos[2], 0 );
				SetIV(n,  pos[0],pos[1],pos[2], 0 );
				SetII(nP, pos[0],pos[1],pos[2], 0 );
				SetIV(nP, pos[0],pos[1],pos[2], 0 );
				SetII(nPP,pos[0],pos[1],pos[2], 0 );
				SetIV(nPP,pos[0],pos[1],pos[2], 0 );
			}
		}
	}
}

bool Operator::Calc_ECPos(int ny, const unsigned int* pos, double* EC, vector<CSPrimitives*> vPrims) const
{
	double EffMat[4];
	Calc_EffMatPos(ny,pos,EffMat, vPrims);

	if (m_epsR)
		m_epsR[ny][pos[0]][pos[1]][pos[2]] =  EffMat[0];
	if (m_kappa)
		m_kappa[ny][pos[0]][pos[1]][pos[2]] =  EffMat[1];
	if (m_mueR)
		m_mueR[ny][pos[0]][pos[1]][pos[2]] =  EffMat[2];
	if (m_sigma)
		m_sigma[ny][pos[0]][pos[1]][pos[2]] =  EffMat[3];

	double delta = GetEdgeLength(ny,pos);
	double area  = GetEdgeArea(ny,pos);

//	if (isnan(EffMat[0]))
//	{
//		cerr << ny << " " << pos[0] << " " << pos[1] << " " << pos[2] << " : " << EffMat[0] << endl;
//	}

	if (delta)
	{
		EC[0] = EffMat[0] * area/delta;
		EC[1] = EffMat[1] * area/delta;
	}
	else
	{
		EC[0] = 0;
		EC[1] = 0;
	}

	delta = GetEdgeLength(ny,pos,true);
	area  = GetEdgeArea(ny,pos,true);

	if (delta)
	{
		EC[2] = EffMat[2] * area/delta;
		EC[3] = EffMat[3] * area/delta;
	}
	else
	{
		EC[2] = 0;
		EC[3] = 0;
	}

	return true;
}

double Operator::GetRawDiscDelta(int ny, const int pos) const
{
	//numLines[ny] is expected to be larger then 1 !

	if (pos<0)
		return (discLines[ny][0] - discLines[ny][1]);
	if (pos>=(int)numLines[ny]-1)
		return (discLines[ny][numLines[ny]-2] - discLines[ny][numLines[ny]-1]);

	return (discLines[ny][pos+1] - discLines[ny][pos]);
}

bool Operator::GetCellCenterMaterialAvgCoord(const int pos[], double coord[3]) const
{
	unsigned int ui_pos[3];
	for (int n=0;n<3;++n)
	{
		if ((pos[n]<0) || (pos[n]>=(int)numLines[n]))
			return false;
		ui_pos[n] = pos[n];
	}
	GetNodeCoords(ui_pos, coord, true);
	return true;
}

double Operator::GetMaterial(int ny, const double* coords, int MatType, vector<CSPrimitives*> vPrims, bool markAsUsed) const
{
	CSProperties* prop = CSX->GetPropertyByCoordPriority(coords,vPrims,markAsUsed);
//	CSProperties* old_prop = CSX->GetPropertyByCoordPriority(coords,CSProperties::MATERIAL,markAsUsed);
//	if (old_prop!=prop)
//	{
//		cerr << "ERROR: Unequal properties!" << endl;
//		exit(-1);
//	}

	CSPropMaterial* mat = dynamic_cast<CSPropMaterial*>(prop);
	if (mat)
	{
		switch (MatType)
		{
		case 0:
			return mat->GetEpsilonWeighted(ny,coords);
		case 1:
			return mat->GetKappaWeighted(ny,coords);
		case 2:
			return mat->GetMueWeighted(ny,coords);
		case 3:
			return mat->GetSigmaWeighted(ny,coords);
		case 4:
			return mat->GetDensityWeighted(coords);
		default:
			cerr << "Operator::GetMaterial: Error: unknown material type" << endl;
			return 0;
		}
	}

	switch (MatType)
	{
	case 0:
		return GetBackgroundEpsR();
	case 1:
		return GetBackgroundKappa();
	case 2:
		return GetBackgroundMueR();
	case 3:
		return GetBackgroundSigma();
	case 4:
		return GetBackgroundDensity();
	default:
		cerr << "Operator::GetMaterial: Error: unknown material type" << endl;
		return 0;
	}
}

bool Operator::AverageMatCellCenter(int ny, const unsigned int* pos, double* EffMat, vector<CSPrimitives *> vPrims) const
{
	int n=ny;
	double coord[3];
	int nP = (n+1)%3;
	int nPP = (n+2)%3;

	int loc_pos[3] = {(int)pos[0],(int)pos[1],(int)pos[2]};
	double A_n;
	double area = 0;
	EffMat[0] = 0;
	EffMat[1] = 0;
	EffMat[2] = 0;
	EffMat[3] = 0;

	//******************************* epsilon,kappa averaging *****************************//
	//shift up-right
	if (GetCellCenterMaterialAvgCoord(loc_pos,coord))
	{
		A_n = GetNodeArea(ny,loc_pos,true);
		EffMat[0] += GetMaterial(n, coord, 0, vPrims)*A_n;
		EffMat[1] += GetMaterial(n, coord, 1, vPrims)*A_n;
		area+=A_n;
	}

	//shift up-left
	--loc_pos[nP];
	if (GetCellCenterMaterialAvgCoord(loc_pos,coord))
	{
		A_n = GetNodeArea(ny,loc_pos,true);
		EffMat[0] += GetMaterial(n, coord, 0, vPrims)*A_n;
		EffMat[1] += GetMaterial(n, coord, 1, vPrims)*A_n;
		area+=A_n;
	}

	//shift down-right
	++loc_pos[nP];
	--loc_pos[nPP];
	if (GetCellCenterMaterialAvgCoord(loc_pos,coord))
	{
		A_n = GetNodeArea(ny,loc_pos,true);
		EffMat[0] += GetMaterial(n, coord, 0, vPrims)*A_n;
		EffMat[1] += GetMaterial(n, coord, 1, vPrims)*A_n;
		area+=A_n;
	}

	//shift down-left
	--loc_pos[nP];
	if (GetCellCenterMaterialAvgCoord(loc_pos,coord))
	{
		A_n = GetNodeArea(ny,loc_pos,true);
		EffMat[0] += GetMaterial(n, coord, 0, vPrims)*A_n;
		EffMat[1] += GetMaterial(n, coord, 1, vPrims)*A_n;
		area+=A_n;
	}

	EffMat[0]*=__EPS0__/area;
	EffMat[1]/=area;

	//******************************* mu,sigma averaging *****************************//
	loc_pos[0]=pos[0];
	loc_pos[1]=pos[1];
	loc_pos[2]=pos[2];
	double length=0;
	double delta_ny,sigma;
	//shift down
	--loc_pos[n];
	if (GetCellCenterMaterialAvgCoord(loc_pos,coord))
	{
		delta_ny = GetNodeWidth(n,loc_pos,true);
		EffMat[2] += delta_ny / GetMaterial(n, coord, 2, vPrims);
		sigma = GetMaterial(n, coord, 3, vPrims);
		if (sigma)
			EffMat[3] += delta_ny / sigma;
		else
			EffMat[3] = 0;
		length+=delta_ny;
	}

	//shift up
	++loc_pos[n];
	if (GetCellCenterMaterialAvgCoord(loc_pos,coord))
	{
		delta_ny = GetNodeWidth(n,loc_pos,true);
		EffMat[2] += delta_ny / GetMaterial(n, coord, 2, vPrims);
		sigma = GetMaterial(n, coord, 3, vPrims);
		if (sigma)
			EffMat[3] += delta_ny / sigma;
		else
			EffMat[3] = 0;
		length+=delta_ny;
	}

	EffMat[2] = length * __MUE0__ / EffMat[2];
	if (EffMat[3]) EffMat[3]=length / EffMat[3];

	for (int n=0; n<4; ++n)
		if (std::isnan(EffMat[n]) || std::isinf(EffMat[n]))
		{
			cerr << "Operator::" << __func__ << ": Error, an effective material parameter is not a valid result, this should NOT have happend... exit..." << endl;
			cerr << ny << "@" << n << " : " << pos[0] << "," << pos[1] << ","  << pos[2] << endl;
			exit(0);
		}
	return true;
}

bool Operator::AverageMatQuarterCell(int ny, const unsigned int* pos, double* EffMat, vector<CSPrimitives*> vPrims) const
{
	int n=ny;
	double coord[3];
	double shiftCoord[3];
	int nP = (n+1)%3;
	int nPP = (n+2)%3;
	coord[0] = discLines[0][pos[0]];
	coord[1] = discLines[1][pos[1]];
	coord[2] = discLines[2][pos[2]];
	double delta=GetRawDiscDelta(n,pos[n]);
	double deltaP=GetRawDiscDelta(nP,pos[nP]);
	double deltaPP=GetRawDiscDelta(nPP,pos[nPP]);
	double delta_M=GetRawDiscDelta(n,pos[n]-1);
	double deltaP_M=GetRawDiscDelta(nP,pos[nP]-1);
	double deltaPP_M=GetRawDiscDelta(nPP,pos[nPP]-1);

	int loc_pos[3] = {(int)pos[0],(int)pos[1],(int)pos[2]};
	double A_n;
	double area = 0;

	//******************************* epsilon,kappa averaging *****************************//
	//shift up-right
	shiftCoord[n] = coord[n]+delta*0.5;
	shiftCoord[nP] = coord[nP]+deltaP*0.25;
	shiftCoord[nPP] = coord[nPP]+deltaPP*0.25;
	A_n = GetNodeArea(ny,loc_pos,true);
	EffMat[0] = GetMaterial(n, shiftCoord, 0, vPrims)*A_n;
	EffMat[1] = GetMaterial(n, shiftCoord, 1, vPrims)*A_n;
	area+=A_n;

	//shift up-left
	shiftCoord[n] = coord[n]+delta*0.5;
	shiftCoord[nP] = coord[nP]-deltaP_M*0.25;
	shiftCoord[nPP] = coord[nPP]+deltaPP*0.25;

	--loc_pos[nP];
	A_n = GetNodeArea(ny,loc_pos,true);
	EffMat[0] += GetMaterial(n, shiftCoord, 0, vPrims)*A_n;
	EffMat[1] += GetMaterial(n, shiftCoord, 1, vPrims)*A_n;
	area+=A_n;

	//shift down-right
	shiftCoord[n] = coord[n]+delta*0.5;
	shiftCoord[nP] = coord[nP]+deltaP*0.25;
	shiftCoord[nPP] = coord[nPP]-deltaPP_M*0.25;
	++loc_pos[nP];
	--loc_pos[nPP];
	A_n = GetNodeArea(ny,loc_pos,true);
	EffMat[0] += GetMaterial(n, shiftCoord, 0, vPrims)*A_n;
	EffMat[1] += GetMaterial(n, shiftCoord, 1, vPrims)*A_n;
	area+=A_n;

	//shift down-left
	shiftCoord[n] = coord[n]+delta*0.5;
	shiftCoord[nP] = coord[nP]-deltaP_M*0.25;
	shiftCoord[nPP] = coord[nPP]-deltaPP_M*0.25;
	--loc_pos[nP];
	A_n = GetNodeArea(ny,loc_pos,true);
	EffMat[0] += GetMaterial(n, shiftCoord, 0, vPrims)*A_n;
	EffMat[1] += GetMaterial(n, shiftCoord, 1, vPrims)*A_n;
	area+=A_n;

	EffMat[0]*=__EPS0__/area;
	EffMat[1]/=area;

	//******************************* mu,sigma averaging *****************************//
	loc_pos[0]=pos[0];
	loc_pos[1]=pos[1];
	loc_pos[2]=pos[2];
	double length=0;

	//shift down
	shiftCoord[n] = coord[n]-delta_M*0.25;
	shiftCoord[nP] = coord[nP]+deltaP*0.5;
	shiftCoord[nPP] = coord[nPP]+deltaPP*0.5;
	--loc_pos[n];
	double delta_ny = GetNodeWidth(n,loc_pos,true);
	EffMat[2] = delta_ny / GetMaterial(n, shiftCoord, 2, vPrims);
	double sigma = GetMaterial(n, shiftCoord, 3, vPrims);
	if (sigma)
		EffMat[3] = delta_ny / sigma;
	else
		EffMat[3] = 0;
	length=delta_ny;

	//shift up
	shiftCoord[n] = coord[n]+delta*0.25;
	shiftCoord[nP] = coord[nP]+deltaP*0.5;
	shiftCoord[nPP] = coord[nPP]+deltaPP*0.5;
	++loc_pos[n];
	delta_ny = GetNodeWidth(n,loc_pos,true);
	EffMat[2] += delta_ny / GetMaterial(n, shiftCoord, 2, vPrims);
	sigma = GetMaterial(n, shiftCoord, 3, vPrims);
	if (sigma)
		EffMat[3] += delta_ny / sigma;
	else
		EffMat[3] = 0;
	length+=delta_ny;

	EffMat[2] = length * __MUE0__ / EffMat[2];
	if (EffMat[3]) EffMat[3]=length / EffMat[3];

	for (int n=0; n<4; ++n)
		if (std::isnan(EffMat[n]) || std::isinf(EffMat[n]))
		{
			cerr << "Operator::" << __func__ << ": Error, An effective material parameter is not a valid result, this should NOT have happend... exit..." << endl;
			cerr << ny << "@" << n << " : " << pos[0] << "," << pos[1] << ","  << pos[2] << endl;
			exit(0);
		}

	return true;
}

bool Operator::Calc_EffMatPos(int ny, const unsigned int* pos, double* EffMat, vector<CSPrimitives *> vPrims) const
{
	switch (m_MatAverageMethod)
	{
	case QuarterCell:
		return AverageMatQuarterCell(ny, pos, EffMat, vPrims);
	case CentralCell:
		return AverageMatCellCenter(ny, pos, EffMat, vPrims);
	default:
		cerr << "Operator:: " << __func__ << ":  Error, unknown material averaging method... exit" << endl;
		exit(1);
	}
	return false;
}

bool Operator::Calc_LumpedElements()
{
	vector<CSProperties*> props = CSX->GetPropertyByType(CSProperties::LUMPED_ELEMENT);
	for (size_t i=0;i<props.size();++i)
	{
		CSPropLumpedElement* PLE = dynamic_cast<CSPropLumpedElement*>(props.at(i));
		if (PLE==NULL)
			return false; //sanity check: this should never happen!
		vector<CSPrimitives*> prims = PLE->GetAllPrimitives();
		for (size_t bn=0;bn<prims.size();++bn)
		{
			CSPrimBox* box = dynamic_cast<CSPrimBox*>(prims.at(bn));
			if (box)
			{	//calculate lumped element parameter

				double C = PLE->GetCapacity();
				if (C<=0)
					C = NAN;
				double R = PLE->GetResistance();
				if (R<0)
					R = NAN;

				if ((std::isnan(R)) && (std::isnan(C)))
				{
					cerr << "Operator::Calc_LumpedElements(): Warning: Lumped Element R or C not specified! skipping. "
							<< " ID: " << prims.at(bn)->GetID() << " @ Property: " << PLE->GetName() << endl;
					continue;
				}

				int ny = PLE->GetDirection();
				if ((ny<0) || (ny>2))
				{
					cerr << "Operator::Calc_LumpedElements(): Warning: Lumped Element direction is invalid! skipping. "
							<< " ID: " << prims.at(bn)->GetID() << " @ Property: " << PLE->GetName() << endl;
					continue;
				}
				int nyP = (ny+1)%3;
				int nyPP = (ny+2)%3;

				unsigned int uiStart[3];
				unsigned int uiStop[3];
				// snap to the native coordinate system
				int Snap_Dimension = Operator::SnapBox2Mesh(box->GetStartCoord()->GetCoords(m_MeshType), box->GetStopCoord()->GetCoords(m_MeshType), uiStart, uiStop, false, true);
				if (Snap_Dimension<=0)
				{
					if (Snap_Dimension>=-1)
						cerr << "Operator::Calc_LumpedElements(): Warning: Lumped Element snapping failed! Dimension is: " << Snap_Dimension << " skipping. "
								<< " ID: " << prims.at(bn)->GetID() << " @ Property: " << PLE->GetName() << endl;
					// Snap_Dimension == -2 means outside the simulation domain --> no special warning, but box probably marked as unused!
					continue;
				}

				if (uiStart[ny]==uiStop[ny])
				{
					cerr << "Operator::Calc_LumpedElements(): Warning: Lumped Element with zero (snapped) length is invalid! skipping. "
							<< " ID: " << prims.at(bn)->GetID() << " @ Property: " << PLE->GetName() << endl;
					continue;
				}

				//calculate geometric property for this lumped element
				unsigned int pos[3];
				double unitGC=0;
				int ipos=0;
				for (pos[ny]=uiStart[ny];pos[ny]<uiStop[ny];++pos[ny])
				{
					double unitGC_Plane=0;
					for (pos[nyP]=uiStart[nyP];pos[nyP]<=uiStop[nyP];++pos[nyP])
					{
						for (pos[nyPP]=uiStart[nyPP];pos[nyPP]<=uiStop[nyPP];++pos[nyPP])
						{
							// capacity/conductivity in parallel: add values
							unitGC_Plane += GetEdgeArea(ny,pos)/GetEdgeLength(ny,pos);
						}
					}

					//capacity/conductivity in series: add reciprocal values
					unitGC += 1/unitGC_Plane;
				}
				unitGC = 1/unitGC;

				bool caps = PLE->GetCaps();
				double kappa = 0;
				double epsilon = 0;
				if (R>0)
					kappa = 1 / R / unitGC;
				if (C>0)
				{
					epsilon =  C / unitGC;

					if (epsilon< __EPS0__)
					{
						cerr << "Operator::Calc_LumpedElements(): Warning: Lumped Element capacity is too small for its size! skipping. "
								<< " ID: " << prims.at(bn)->GetID() << " @ Property: " << PLE->GetName() << endl;
						C = 0;
					}
				}

				for (pos[ny]=uiStart[ny];pos[ny]<uiStop[ny];++pos[ny])
				{
					for (pos[nyP]=uiStart[nyP];pos[nyP]<=uiStop[nyP];++pos[nyP])
					{
						for (pos[nyPP]=uiStart[nyPP];pos[nyPP]<=uiStop[nyPP];++pos[nyPP])
						{
							ipos = MainOp->SetPos(pos[0],pos[1],pos[2]);
							if (C>0)
								EC_C[ny][ipos] = epsilon * GetEdgeArea(ny,pos)/GetEdgeLength(ny,pos);
							if (R>0)
								EC_G[ny][ipos] = kappa * GetEdgeArea(ny,pos)/GetEdgeLength(ny,pos);

							if (R==0) //make lumped element a PEC if resistance is zero
							{
								SetVV(ny,pos[0],pos[1],pos[2], 0 );
								SetVI(ny,pos[0],pos[1],pos[2], 0 );
							}
							else //recalculate operator inside the lumped element
								Calc_ECOperatorPos(ny,pos);
						}
					}
				}

				// setup metal caps
				if (caps)
				{
					for (pos[nyP]=uiStart[nyP];pos[nyP]<=uiStop[nyP];++pos[nyP])
					{
						for (pos[nyPP]=uiStart[nyPP];pos[nyPP]<=uiStop[nyPP];++pos[nyPP])
						{
							pos[ny]=uiStart[ny];
							if (pos[nyP]<uiStop[nyP])
							{
								SetVV(nyP,pos[0],pos[1],pos[2], 0 );
								SetVI(nyP,pos[0],pos[1],pos[2], 0 );
								++m_Nr_PEC[nyP];
							}

							if (pos[nyPP]<uiStop[nyPP])
							{
								SetVV(nyPP,pos[0],pos[1],pos[2], 0 );
								SetVI(nyPP,pos[0],pos[1],pos[2], 0 );
								++m_Nr_PEC[nyPP];
							}

							pos[ny]=uiStop[ny];
							if (pos[nyP]<uiStop[nyP])
							{
								SetVV(nyP,pos[0],pos[1],pos[2], 0 );
								SetVI(nyP,pos[0],pos[1],pos[2], 0 );
								++m_Nr_PEC[nyP];
							}

							if (pos[nyPP]<uiStop[nyPP])
							{
								SetVV(nyPP,pos[0],pos[1],pos[2], 0 );
								SetVI(nyPP,pos[0],pos[1],pos[2], 0 );
								++m_Nr_PEC[nyPP];
							}
						}
					}
				}
				box->SetPrimitiveUsed(true);

			}
			else
				cerr << "Operator::Calc_LumpedElements(): Warning: Primitves other than boxes are not supported for lumped elements! skipping "
						<< prims.at(bn)->GetTypeName() << " ID: " << prims.at(bn)->GetID() << " @ Property: " << PLE->GetName() << endl;
		}
	}
	return true;
}

void Operator::Init_EC()
{
	for (int n=0; n<3; ++n)
	{
		//init x-cell-array
		delete[] EC_C[n];
		delete[] EC_G[n];
		delete[] EC_L[n];
		delete[] EC_R[n];
		EC_C[n] = new FDTD_FLOAT[MainOp->GetSize()];
		EC_G[n] = new FDTD_FLOAT[MainOp->GetSize()];
		EC_L[n] = new FDTD_FLOAT[MainOp->GetSize()];
		EC_R[n] = new FDTD_FLOAT[MainOp->GetSize()];
		for (unsigned int i=0; i<MainOp->GetSize(); i++) //init all
		{
			EC_C[n][i]=0;
			EC_G[n][i]=0;
			EC_L[n][i]=0;
			EC_R[n][i]=0;
		}
	}
}

bool Operator::Calc_EC()
{
	if (CSX==NULL)
	{
		cerr << "CartOperator::Calc_EC: CSX not given or invalid!!!" << endl;
		return false;
	}
	
	MainOp->SetPos(0,0,0);
	Calc_EC_Range(0,numLines[0]-1);
	return true;
}

vector<CSPrimitives*> Operator::GetPrimitivesBoundBox(int posX, int posY, int posZ, CSProperties::PropertyType type) const
{
	double boundBox[6];
	int BBpos[3] = {posX, posY, posZ};
	for (int n=0;n<3;++n)
	{
		if (BBpos[n]<0)
		{
			boundBox[2*n]   = this->GetDiscLine(n,0);
			boundBox[2*n+1] = this->GetDiscLine(n,numLines[n]-1);
		}
		else
		{
			boundBox[2*n]   = this->GetDiscLine(n, max(0, BBpos[n]-1));
			boundBox[2*n+1] = this->GetDiscLine(n, min(int(numLines[n])-1, BBpos[n]+1));
		}
	}

	vector<CSPrimitives*> vPrim = this->CSX->GetPrimitivesByBoundBox(boundBox, true, type);
	return vPrim;
}

void Operator::Calc_EC_Range(unsigned int xStart, unsigned int xStop)
{
//	vector<CSPrimitives*> vPrims = this->CSX->GetAllPrimitives(true, CSProperties::MATERIAL);
	unsigned int ipos;
	unsigned int pos[3];
	double inEC[4];
	for (pos[0]=xStart; pos[0]<=xStop; ++pos[0])
	{
		for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
		{
			vector<CSPrimitives*> vPrims = this->GetPrimitivesBoundBox(pos[0], pos[1], -1, CSProperties::MATERIAL);
			for (pos[2]=0; pos[2]<numLines[2]; ++pos[2])
			{
				ipos = MainOp->GetPos(pos[0],pos[1],pos[2]);
				for (int n=0; n<3; ++n)
				{
					Calc_ECPos(n,pos,inEC,vPrims);
					EC_C[n][ipos]=inEC[0];
					EC_G[n][ipos]=inEC[1];
					EC_L[n][ipos]=inEC[2];
					EC_R[n][ipos]=inEC[3];
				}
			}
		}
	}
}

void Operator::SetTimestepFactor(double factor)
{
	if ((factor<=0) || (factor>1))
	{
		cerr << "Operator::SetTimestepFactor: Warning, invalid timestep factor, skipping!" << endl;
		return;
	}

	cout << "Operator::SetTimestepFactor: Setting timestep factor to " << factor << endl;
	m_TimeStepFactor=factor;
}

double Operator::CalcTimestep()
{
	if (m_TimeStepVar==3)
		return CalcTimestep_Var3(); //the biggest one for cartesian meshes

	//variant 1 is default
	return CalcTimestep_Var1();
}

////Berechnung nach Andreas Rennings Dissertation 2008, Seite 66, Formel 4.52
double Operator::CalcTimestep_Var1()
{
	m_Used_TS_Name = string("Rennings_1");
//	cout << "Operator::CalcTimestep(): Using timestep algorithm by Andreas Rennings, Dissertation @ University Duisburg-Essen, 2008, pp. 66, eq. 4.52" << endl;
	dT=1e200;
	double newT;
	unsigned int pos[3];
	unsigned int smallest_pos[3] = {0, 0, 0};
	unsigned int smallest_n = 0;
	unsigned int ipos;
	unsigned int ipos_PM;
	unsigned int ipos_PPM;
	MainOp->SetReflection2Cell();
	for (int n=0; n<3; ++n)
	{
		int nP = (n+1)%3;
		int nPP = (n+2)%3;

		for (pos[2]=0; pos[2]<numLines[2]; ++pos[2])
		{
			for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
			{
				for (pos[0]=0; pos[0]<numLines[0]; ++pos[0])
				{
					ipos = MainOp->SetPos(pos[0],pos[1],pos[2]);
					ipos_PM = MainOp->Shift(nP,-1);
					MainOp->ResetShift();
					ipos_PPM= MainOp->Shift(nPP,-1);
					MainOp->ResetShift();
					newT = 2/sqrt( ( 4/EC_L[nP][ipos] + 4/EC_L[nP][ipos_PPM] + 4/EC_L[nPP][ipos] + 4/EC_L[nPP][ipos_PM]) / EC_C[n][ipos] );
					if ((newT<dT) && (newT>0.0))
					{
						dT=newT;
						smallest_pos[0]=pos[0];smallest_pos[1]=pos[1];smallest_pos[2]=pos[2];
						smallest_n = n;
					}
				}
			}
		}
	}
	if (dT==0)
	{
		cerr << "Operator::CalcTimestep: Timestep is zero... this is not supposed to happen!!! exit!" << endl;
		exit(3);
	}
	if (g_settings.GetVerboseLevel()>1)
	{
		cout << "Operator::CalcTimestep_Var1: Smallest timestep (" << dT << "s) found at position: " <<  smallest_n << " : " << smallest_pos[0] << ";" <<  smallest_pos[1] << ";" <<  smallest_pos[2] << endl;
	}
	return 0;
}

double min(double* val, unsigned int count)
{
	if (count==0)
		return 0.0;
	double min = val[0];
	for (unsigned int n=1; n<count; ++n)
		if (val[n]<min)
			min = val[n];
	return min;
}

//Berechnung nach Andreas Rennings Dissertation 2008, Seite 76 ff, Formel 4.77 ff
double Operator::CalcTimestep_Var3()
{
	dT=1e200;
	m_Used_TS_Name = string("Rennings_2");
//	cout << "Operator::CalcTimestep(): Using timestep algorithm by Andreas Rennings, Dissertation @ University Duisburg-Essen, 2008, pp. 76, eq. 4.77 ff." << endl;
	double newT;
	unsigned int pos[3];
	unsigned int smallest_pos[3] = {0, 0, 0};
	unsigned int smallest_n = 0;
	unsigned int ipos;
	double w_total=0;
	double wqp=0,wt1=0,wt2=0;
	double wt_4[4]={0,0,0,0};
	MainOp->SetReflection2Cell();
	for (int n=0; n<3; ++n)
	{
		int nP = (n+1)%3;
		int nPP = (n+2)%3;

		for (pos[2]=0; pos[2]<numLines[2]; ++pos[2])
		{
			for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
			{
				for (pos[0]=0; pos[0]<numLines[0]; ++pos[0])
				{
					MainOp->ResetShift();
					ipos = MainOp->SetPos(pos[0],pos[1],pos[2]);
					wqp  = 1/(EC_L[nPP][ipos]*EC_C[n][MainOp->GetShiftedPos(nP ,1)]) + 1/(EC_L[nPP][ipos]*EC_C[n][ipos]);
					wqp += 1/(EC_L[nP ][ipos]*EC_C[n][MainOp->GetShiftedPos(nPP,1)]) + 1/(EC_L[nP ][ipos]*EC_C[n][ipos]);
					ipos = MainOp->Shift(nP,-1);
					wqp += 1/(EC_L[nPP][ipos]*EC_C[n][MainOp->GetShiftedPos(nP ,1)]) + 1/(EC_L[nPP][ipos]*EC_C[n][ipos]);
					ipos = MainOp->Shift(nPP,-1);
					wqp += 1/(EC_L[nP ][ipos]*EC_C[n][MainOp->GetShiftedPos(nPP,1)]) + 1/(EC_L[nP ][ipos]*EC_C[n][ipos]);

					MainOp->ResetShift();
					ipos = MainOp->SetPos(pos[0],pos[1],pos[2]);
					wt_4[0] = 1/(EC_L[nPP][ipos]						  *EC_C[nP ][ipos]);
					wt_4[1] = 1/(EC_L[nPP][MainOp->GetShiftedPos(nP ,-1)] *EC_C[nP ][ipos]);
					wt_4[2] = 1/(EC_L[nP ][ipos]						  *EC_C[nPP][ipos]);
					wt_4[3] = 1/(EC_L[nP ][MainOp->GetShiftedPos(nPP,-1)] *EC_C[nPP][ipos]);

					wt1 = wt_4[0]+wt_4[1]+wt_4[2]+wt_4[3] - 2*min(wt_4,4);

					MainOp->ResetShift();
					ipos = MainOp->SetPos(pos[0],pos[1],pos[2]);
					wt_4[0] = 1/(EC_L[nPP][ipos]						  *EC_C[nP ][MainOp->GetShiftedPos(n,1)]);
					wt_4[1] = 1/(EC_L[nPP][MainOp->GetShiftedPos(nP ,-1)] *EC_C[nP ][MainOp->GetShiftedPos(n,1)]);
					wt_4[2] = 1/(EC_L[nP ][ipos]						  *EC_C[nPP][MainOp->GetShiftedPos(n,1)]);
					wt_4[3] = 1/(EC_L[nP ][MainOp->GetShiftedPos(nPP,-1)] *EC_C[nPP][MainOp->GetShiftedPos(n,1)]);

					wt2 = wt_4[0]+wt_4[1]+wt_4[2]+wt_4[3] - 2*min(wt_4,4);

					w_total = wqp + wt1 + wt2;
					newT = 2/sqrt( w_total );
					if ((newT<dT) && (newT>0.0))
					{
						dT=newT;
						smallest_pos[0]=pos[0];smallest_pos[1]=pos[1];smallest_pos[2]=pos[2];
						smallest_n = n;
					}
				}
			}
		}
	}
	if (dT==0)
	{
		cerr << "Operator::CalcTimestep: Timestep is zero... this is not supposed to happen!!! exit!" << endl;
		exit(3);
	}
	if (g_settings.GetVerboseLevel()>1)
	{
		cout << "Operator::CalcTimestep_Var3: Smallest timestep (" << dT << "s) found at position: " <<  smallest_n << " : " << smallest_pos[0] << ";" <<  smallest_pos[1] << ";" <<  smallest_pos[2] << endl;
	}
	return 0;
}

bool Operator::CalcPEC()
{
	m_Nr_PEC[0]=0;
	m_Nr_PEC[1]=0;
	m_Nr_PEC[2]=0;

	CalcPEC_Range(0,numLines[0]-1,m_Nr_PEC);

	CalcPEC_Curves();

	return true;
}

void Operator::CalcPEC_Range(unsigned int startX, unsigned int stopX, unsigned int* counter)
{
	double coord[3];
	unsigned int pos[3];
	for (pos[0]=startX; pos[0]<=stopX; ++pos[0])
	{
		for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
		{
			vector<CSPrimitives*> vPrims = this->GetPrimitivesBoundBox(pos[0], pos[1], -1, (CSProperties::PropertyType)(CSProperties::MATERIAL | CSProperties::METAL));
			for (pos[2]=0; pos[2]<numLines[2]; ++pos[2])
			{
				for (int n=0; n<3; ++n)
				{
					GetYeeCoords(n,pos,coord,false);
					CSProperties* prop = CSX->GetPropertyByCoordPriority(coord, vPrims, true);
//					CSProperties* old_prop = CSX->GetPropertyByCoordPriority(coord, (CSProperties::PropertyType)(CSProperties::MATERIAL | CSProperties::METAL), true);
//					if (old_prop!=prop)
//					{
//						cerr << "CalcPEC_Range: " << old_prop << " vs " << prop << endl;
//						exit(-1);
//					}
					if (prop)
					{
						if (prop->GetType()==CSProperties::METAL) //set to PEC
						{
							SetVV(n,pos[0],pos[1],pos[2], 0 );
							SetVI(n,pos[0],pos[1],pos[2], 0 );
							++counter[n];
						}
					}
				}
			}
		}
	}
}

void Operator::CalcPEC_Curves()
{
	//special treatment for primitives of type curve (treated as wires)
	double p1[3];
	double p2[3];
	Grid_Path path;
	vector<CSProperties*> vec_prop = CSX->GetPropertyByType(CSProperties::METAL);
	for (size_t p=0; p<vec_prop.size(); ++p)
	{
		CSProperties* prop = vec_prop.at(p);
		for (size_t n=0; n<prop->GetQtyPrimitives(); ++n)
		{
			CSPrimitives* prim = prop->GetPrimitive(n);
			CSPrimCurve* curv = prim->ToCurve();
			if (curv)
			{
				for (size_t i=1; i<curv->GetNumberOfPoints(); ++i)
				{
					curv->GetPoint(i-1,p1,m_MeshType);
					curv->GetPoint(i,p2,m_MeshType);
					path = FindPath(p1,p2);
					if (path.dir.size()>0)
						prim->SetPrimitiveUsed(true);
					for (size_t t=0; t<path.dir.size(); ++t)
					{
						SetVV(path.dir.at(t),path.posPath[0].at(t),path.posPath[1].at(t),path.posPath[2].at(t), 0 );
						SetVI(path.dir.at(t),path.posPath[0].at(t),path.posPath[1].at(t),path.posPath[2].at(t), 0 );
						++m_Nr_PEC[path.dir.at(t)];
					}
				}
			}
		}
	}
}

Operator_Ext_Excitation* Operator::GetExcitationExtension() const
{
	//search for excitation extension
	Operator_Ext_Excitation* Op_Ext_Exc=0;
	for (size_t n=0; n<m_Op_exts.size(); ++n)
	{
		Op_Ext_Exc = dynamic_cast<Operator_Ext_Excitation*>(m_Op_exts.at(n));
		if (Op_Ext_Exc)
			break;
	}
	return Op_Ext_Exc;
}

void Operator::AddExtension(Operator_Extension* op_ext)
{
	m_Op_exts.push_back(op_ext);
}

void Operator::DeleteExtension(Operator_Extension* op_ext)
{
	for (size_t n=0;n<m_Op_exts.size();++n)
	{
		if (m_Op_exts.at(n)==op_ext)
		{
			m_Op_exts.erase(m_Op_exts.begin()+n);
			return;
		}
	}
}

double Operator::CalcNumericPhaseVelocity(unsigned int start[3], unsigned int stop[3], double propDir[3], float freq) const
{
	double average_mesh_disc[3];
	double c0 = __C0__/sqrt(GetBackgroundEpsR()*GetBackgroundMueR());

	//calculate average mesh deltas
	for (int n=0;n<3;++n)
	{
		average_mesh_disc[n] = fabs(GetDiscLine(n,start[n])-GetDiscLine(n,stop[n]))*GetGridDelta() / (fabs(stop[n]-start[n]));
	}

	// if propagation is in a Cartesian direction, return analytic solution
	for (int n=0;n<3;++n)
	{
		int nP = (n+1)%3;
		int nPP = (n+2)%3;
		if ((fabs(propDir[n])==1) && (propDir[nP]==0) && (propDir[nPP]==0))
		{
			double kx = asin(average_mesh_disc[0]/c0/dT*sin(2*PI*freq*dT/2))*2/average_mesh_disc[0];
			return 2*PI*freq/kx;
		}
	}

	// else, do an newton iterative estimation
	double k0=2*PI*freq/c0;
	double k=k0;
	double RHS = pow(sin(2*PI*freq*dT/2)/c0/dT,2);
	double fk=1,fdk=0;
	double old_phv=0;
	double phv=c0;
	double err_est = 1e-6;
	int it_count=0;
	while (fabs(old_phv-phv)>err_est)
	{
		++it_count;
		old_phv=phv;
		fk=0;
		fdk=0;
		for (int n=0;n<3;++n)
		{
			fk+= pow(sin(propDir[n]*k*average_mesh_disc[n]/2)/average_mesh_disc[n],2);
			fdk+= propDir[n]*sin(propDir[n]*k*average_mesh_disc[n]/2)*cos(propDir[n]*k*average_mesh_disc[n]/2)/average_mesh_disc[n];
		}
		fk -= RHS;
		k-=fk/fdk;

		// do not allow a speed greater than c0 due to a numerical inaccuracy
		if (k<k0)
			k=k0;

		phv=2*PI*freq/k;

		//abort if iteration count is getting to high!
		if (it_count>99)
		{
			cerr << "Operator::CalcNumericPhaseVelocity: Error, newton iteration estimation can't find a solution!!" << endl;
			break;
		}
	}

	if (g_settings.GetVerboseLevel()>1)
		cerr << "Operator::CalcNumericPhaseVelocity: Newton iteration estimated solution: " << phv/__C0__ << "*c0 in " <<  it_count << " iterations." << endl;

	return phv;
}
