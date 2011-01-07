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
#include "operator.h"
#include "engine.h"
#include "extensions/operator_extension.h"
#include "Common/processfields.h"
#include "tools/array_ops.h"
#include "fparser.hh"

Operator* Operator::New()
{
	cout << "Create FDTD operator" << endl;
	Operator* op = new Operator();
	op->Init();
	return op;
}

Operator::Operator() : Operator_Base()
{
	Exc = 0;
	m_InvaildTimestep = false;
}

Operator::~Operator()
{
	for (size_t n=0; n<m_Op_exts.size(); ++n)
		delete m_Op_exts.at(n);
	m_Op_exts.clear();
	Reset();
}

Engine* Operator::CreateEngine() const
{
	Engine* eng = Engine::New(this);
	return eng;
}

void Operator::Init()
{
	Operator_Base::Init();

	CSX = NULL;

	vv=NULL;
	vi=NULL;
	iv=NULL;
	ii=NULL;

	m_epsR=NULL;
	m_kappa=NULL;
	m_mueR=NULL;
	m_sigma=NULL;

	MainOp=NULL;
	DualOp=NULL;

	for (int n=0; n<3; ++n)
	{
		EC_C[n]=NULL;
		EC_G[n]=NULL;
		EC_L[n]=NULL;
		EC_R[n]=NULL;
	}

	Exc = 0;
}

void Operator::Reset()
{
	Delete_N_3DArray(vv,numLines);
	Delete_N_3DArray(vi,numLines);
	Delete_N_3DArray(iv,numLines);
	Delete_N_3DArray(ii,numLines);
	delete MainOp;
	delete DualOp;
	for (int n=0; n<3; ++n)
	{
		delete[] EC_C[n];
		delete[] EC_G[n];
		delete[] EC_L[n];
		delete[] EC_R[n];
	}

	delete Exc;

	Delete_N_3DArray(m_epsR,numLines);
	Delete_N_3DArray(m_kappa,numLines);
	Delete_N_3DArray(m_mueR,numLines);
	Delete_N_3DArray(m_sigma,numLines);

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

double Operator::GetEdgeLength(int n, const unsigned int* pos, bool dualMesh) const
{
	if ((n<0) || (n>2)) return 0.0;
	if (pos[n]>=numLines[n]) return 0.0;
	double delta=0;
	if (dualMesh==false)
	{
		if (pos[n]<numLines[n]-1)
			delta = GetDiscLine(n,pos[n]+1,false) - GetDiscLine(n,pos[n],false);
		else
			delta = GetDiscLine(n,pos[n],false) - GetDiscLine(n,pos[n]-1,false);
		return delta*gridDelta;
	}
	else
	{
		if (pos[n]>0)
			delta = GetDiscLine(n,pos[n],true) - GetDiscLine(n,pos[n]-1,true);
		else
			delta = GetDiscLine(n,1,false) - GetDiscLine(n,0,false);
		return delta*gridDelta;
	}
}

double Operator::GetNodeArea(int ny, const unsigned int pos[3], bool dualMesh) const
{
	int nyP = (ny+1)%3;
	int nyPP = (ny+2)%3;
	return GetNodeWidth(nyP,pos,dualMesh) * GetNodeWidth(nyPP,pos,dualMesh);
}

bool Operator::SnapToMesh(const double* dcoord, unsigned int* uicoord, bool lower, bool* inside) const
{
	bool ok=true;
	unsigned int numLines[3];
	for (int n=0; n<3; ++n)
	{
		numLines[n] = GetNumberOfLines(n);
		if (inside) //set defaults
			inside[n] = true;
		uicoord[n]=0;
		if (dcoord[n]<discLines[n][0])
		{
			ok=false;
			uicoord[n]=0;
			if (inside) inside[n] = false;
		}
		else if (dcoord[n]==discLines[n][0])
			uicoord[n]=0;
		else if (dcoord[n]>discLines[n][numLines[n]-1])
		{
			ok=false;
			uicoord[n]=numLines[n]-1;
			if (lower) uicoord[n]=numLines[n]-2;
			if (inside) inside[n] = false;
		}
		else if (dcoord[n]==discLines[n][numLines[n]-1])
		{
			uicoord[n]=numLines[n]-1;
			if (lower) uicoord[n]=numLines[n]-2;
		}
		else
			for (unsigned int i=1; i<numLines[n]; ++i)
			{
				if (dcoord[n]<discLines[n][i])
				{
					if (fabs(dcoord[n]-discLines[n][i])<(fabs(dcoord[n]-discLines[n][i-1])))
						uicoord[n]=i;
					else
						uicoord[n]=i-1;
					if (lower) uicoord[n]=i-1;
					i = numLines[n];
				}
			}
	}
//	cerr << "Operator::SnapToMesh Wish: " << dcoord[0] << " " << dcoord[1] << " " << dcoord[2] << endl;
//	cerr << "Operator::SnapToMesh Found: " << discLines[0][uicoord[0]] << " " << discLines[1][uicoord[1]] << " " << discLines[2][uicoord[2]] << endl;
//	cerr << "Operator::SnapToMesh Index: " << uicoord[0] << " " << uicoord[1] << " " << uicoord[2] << endl;
	return ok;
}

struct Operator::Grid_Path Operator::FindPath(double start[], double stop[])
{
	struct Grid_Path path;
//	double dV[] = {stop[0]-start[0],stop[1]-start[1],stop[2]-start[2]};

	unsigned int uiStart[3],uiStop[3],currPos[3];
	SnapToMesh(start,uiStart);
	SnapToMesh(stop,uiStop);
	currPos[0]=uiStart[0];
	currPos[1]=uiStart[1];
	currPos[2]=uiStart[2];
	double meshStart[] = {discLines[0][uiStart[0]], discLines[1][uiStart[1]], discLines[2][uiStart[2]]};
	double meshStop[] = {discLines[0][uiStop[0]], discLines[1][uiStop[1]], discLines[2][uiStop[2]]};

	bool UpDir = false;
	double foot=0,dist=0,minFoot=0,minDist=0;
	int minDir=0;
	unsigned int minPos[3];
	double startFoot,stopFoot,currFoot;
	Point_Line_Distance(meshStart,start,stop,startFoot,dist);
	Point_Line_Distance(meshStop,start,stop,stopFoot,dist);
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
				Point_Line_Distance(P,start,stop,foot,dist);
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
				Point_Line_Distance(P,start,stop,foot,dist);
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
	cout << "Voltage excitations\t: " << Exc->Volt_Count    << "\t (" << Exc->Volt_Count_Dir[0] << ", " << Exc->Volt_Count_Dir[1] << ", " << Exc->Volt_Count_Dir[2] << ")" << endl;
	cout << "Current excitations\t: " << Exc->Curr_Count << "\t (" << Exc->Curr_Count_Dir[0] << ", " << Exc->Curr_Count_Dir[1] << ", " << Exc->Curr_Count_Dir[2] << ")" << endl;
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
	cout << "Nyquist criteria (TS)\t: " << Exc->GetNyquistNum() << endl;
	cout << "Nyquist criteria (s)\t: " << Exc->GetNyquistNum()*dT << endl;
	cout << "Excitation Length (TS)\t: " << Exc->Length << endl;
	cout << "Excitation Length (s)\t: " << Exc->Length*dT << endl;
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

	ofstream file(filename.c_str(),ios_base::out);
	if (!file.is_open())
	{
		cerr << "Operator::DumpOperator2File(): Can't open file: " << filename << endl;
		return;
	}

	cout << "Operator: Dumping FDTD operator information to vtk file: " << filename << " ..." << flush;

	FDTD_FLOAT**** exc = Create_N_3DArray<FDTD_FLOAT>(numLines);
	if (Exc)
	{
		for (unsigned int n=0; n<Exc->Volt_Count; ++n)
			exc[Exc->Volt_dir[n]][Exc->Volt_index[0][n]][Exc->Volt_index[1][n]][Exc->Volt_index[2][n]] = Exc->Volt_amp[n];
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

	string names[] = {"vv", "vi", "iv" , "ii", "exc"};
	FDTD_FLOAT**** array[] = {vv_temp,vi_temp,iv_temp,ii_temp,exc};

	ProcessFields::DumpMultiVectorArray2VTK(file, names , array , 5, discLines, numLines, 6, "Operator dump" , (ProcessFields::MeshType)m_MeshType, discLines_scaling);

	Delete_N_3DArray(ii_temp,numLines);
	Delete_N_3DArray(iv_temp,numLines);
	Delete_N_3DArray(vi_temp,numLines);
	Delete_N_3DArray(vv_temp,numLines);
	Delete_N_3DArray(exc,numLines);

	file.close();

	cout << " done!" << endl;
}

//! \brief dump PEC (perfect electric conductor) information (into VTK-file)
//! visualization via paraview
//! visualize only one component (x, y or z)
void Operator::DumpPEC2File( string filename )
{
	ofstream file( filename.c_str() );
	if (!file.is_open())
	{
		cerr << "Operator::DumpPEC2File(): Can't open file: " << filename << endl;
		return;
	}

	cout << "Operator: Dumping PEC information to vtk file: " << filename << " ..." << flush;

	FDTD_FLOAT**** pec = Create_N_3DArray<FDTD_FLOAT>( numLines );
	unsigned int pos[3];

#ifdef OUTPUT_IN_DRAWINGUNITS
	double scaling = 1.0/GetGridDelta();
#else
	double scaling = 1;
#endif

	for (pos[0]=0; pos[0]<numLines[0]-1; pos[0]++)
	{
		for (pos[1]=0; pos[1]<numLines[1]-1; pos[1]++)
		{
			for (pos[2]=0; pos[2]<numLines[2]-1; pos[2]++)
			{
				if ((pos[1] != 0) && (pos[2] != 0))
				{
					// PEC surrounds the computational area; do not output this
					if ((GetVV(0,pos) == 0) && (GetVI(0,pos) == 0))
						pec[0][pos[0]][pos[1]][pos[2]] = GetEdgeLength( 0, pos ) * scaling; // PEC-x found
				}
				if ((pos[0] != 0) && (pos[2] != 0))
				{
					// PEC surrounds the computational area; do not output this
					if ((GetVV(1,pos) == 0) && (GetVI(1,pos) == 0))
						pec[1][pos[0]][pos[1]][pos[2]] = GetEdgeLength( 1, pos ) * scaling; // PEC-y found
				}
				if ((pos[0] != 0) && (pos[1] != 0))
				{
					// PEC surrounds the computational area; do not output this
					if ((GetVV(2,pos) == 0) && (GetVI(2,pos) == 0))
						pec[2][pos[0]][pos[1]][pos[2]] = GetEdgeLength( 2, pos ) * scaling; // PEC-z found
				}
			}
		}
	}

	// evaluate boundary conditions
	for (int n=0; n<3; n++)
	{
		int nP = (n+1)%3;
		int nPP = (n+2)%3;
		for (pos[nP]=0; pos[nP]<numLines[nP]; pos[nP]++)
		{
			for (pos[nPP]=0; pos[nPP]<numLines[nPP]; pos[nPP]++)
			{
				pos[n] = 0;
				if ((pos[nP] != numLines[nP]-1) && (m_BC[2*n] == 0))
					pec[nP ][pos[0]][pos[1]][pos[2]] = GetEdgeLength( nP,  pos ) * scaling;
				if ((pos[nPP] != numLines[nPP]-1) && (m_BC[2*n] == 0))
					pec[nPP][pos[0]][pos[1]][pos[2]] = GetEdgeLength( nPP, pos ) * scaling;

				pos[n] = numLines[n]-1;
				if ((pos[nP] != numLines[nP]-1) && (m_BC[2*n+1] == 0))
					pec[nP ][pos[0]][pos[1]][pos[2]] = GetEdgeLength( nP,  pos ) * scaling;
				if ((pos[nPP] != numLines[nPP]-1) && (m_BC[2*n+1] == 0))
					pec[nPP][pos[0]][pos[1]][pos[2]] = GetEdgeLength( nPP, pos ) * scaling;
			}
		}
	}

#ifdef OUTPUT_IN_DRAWINGUNITS
	scaling = 1;
#else
	scaling = GetGridDelta();
#endif
	ProcessFields::DumpVectorArray2VTK( file, "PEC", pec, discLines, numLines, 6, "PEC dump" , (ProcessFields::MeshType)m_MeshType, scaling );

	file.close();

	cout << " done!" << endl;
}

void Operator::DumpMaterial2File(string filename)
{
#ifdef OUTPUT_IN_DRAWINGUNITS
	double discLines_scaling = 1;
#else
	double discLines_scaling = GetGridDelta();
#endif

	ofstream file(filename.c_str(),ios_base::out);
	if (!file.is_open())
	{
		cerr << "Operator::DumpMaterial2File(): Can't open file: " << filename << endl;
		return;
	}

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
			for (pos[2]=0; pos[2]<numLines[2]; ++pos[2])
			{
				for (int n=0; n<3; ++n)
				{
					double inMat[4];
					Calc_EffMatPos(n, pos, inMat);
					epsilon[n][pos[0]][pos[1]][pos[2]] = inMat[0]/__EPS0__;
					mue[n][pos[0]][pos[1]][pos[2]]     = inMat[2]/__MUE0__;
					kappa[n][pos[0]][pos[1]][pos[2]]   = inMat[1];
					sigma[n][pos[0]][pos[1]][pos[2]]   = inMat[3];
				}
			}
		}
	}

	string names[] = {"epsilon","mue","kappa","sigma"};
	FDTD_FLOAT**** array[] = {epsilon,mue,kappa,sigma};

	ProcessFields::DumpMultiVectorArray2VTK(file, names, array, 4, discLines, numLines,  6, "Material dump" , (ProcessFields::MeshType)m_MeshType, discLines_scaling);

	Delete_N_3DArray(epsilon,numLines);
	Delete_N_3DArray(mue,numLines);
	Delete_N_3DArray(kappa,numLines);
	Delete_N_3DArray(sigma,numLines);

	file.close();

	cout << " done!" << endl;
}

bool Operator::SetGeometryCSX(ContinuousStructure* geo)
{
	if (geo==NULL) return false;

	CSX = geo;

	CSRectGrid* grid=CSX->GetGrid();
	for (int n=0; n<3; ++n)
	{
		discLines[n] = grid->GetLines(n,discLines[n],numLines[n],true);
		if (n==1)
			if (numLines[n]<3)
			{
				cerr << "CartOperator::SetGeometryCSX: you need at least 3 disc-lines in every direction (3D!)!!!" << endl;
				Reset();
				return false;
			}
	}
	MainOp = new AdrOp(numLines[0],numLines[1],numLines[2]);
	MainOp->SetGrid(discLines[0],discLines[1],discLines[2]);
	if (grid->GetDeltaUnit()<=0)
	{
		cerr << "CartOperator::SetGeometryCSX: grid delta unit must not be <=0 !!!" << endl;
		Reset();
		return false;
	}
	else gridDelta=grid->GetDeltaUnit();
	MainOp->SetGridDelta(1);
	MainOp->AddCellAdrOp();
	return true;
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

void Operator::InitExcitation()
{
	delete Exc;
	Exc = new Excitation( dT );
}

void Operator::Calc_ECOperatorPos(int n, unsigned int* pos)
{
	unsigned int i = MainOp->SetPos(pos[0],pos[1],pos[2]);
	if (EC_C[n][i]>0)
	{
		SetVV(n,pos[0],pos[1],pos[2], (1-dT*EC_G[n][i]/2/EC_C[n][i])/(1+dT*EC_G[n][i]/2/EC_C[n][i]) );
		SetVI(n,pos[0],pos[1],pos[2], (dT/EC_C[n][i])/(1+dT*EC_G[n][i]/2/EC_C[n][i]) );
	}
	else
	{
		SetVV(n,pos[0],pos[1],pos[2], 0 );
		SetVI(n,pos[0],pos[1],pos[2], 0 );
	}
	if (EC_L[n][i]>0)
	{
		SetII(n,pos[0],pos[1],pos[2], (1-dT*EC_R[n][i]/2/EC_L[n][i])/(1+dT*EC_R[n][i]/2/EC_L[n][i]) );
		SetIV(n,pos[0],pos[1],pos[2], (dT/EC_L[n][i])/(1+dT*EC_R[n][i]/2/EC_L[n][i]) );
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

	bool PMC[6];
	for (int n=0; n<6; ++n)
		PMC[n] = m_BC[n]==1;
	ApplyMagneticBC(PMC);

	InitExcitation();

	if (CalcFieldExcitation()==false) return -1;

	//all information available for extension... create now...
	for (size_t n=0; n<m_Op_exts.size(); ++n)
		m_Op_exts.at(n)->BuildExtension();

	if (debugFlags & debugMaterial)
		DumpMaterial2File( "material_dump.vtk" );
	if (debugFlags & debugOperator)
		DumpOperator2File( "operator_dump.vtk" );
	if (debugFlags & debugPEC)
		DumpPEC2File( "PEC_dump.vtk" );

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

bool Operator::Calc_ECPos(int ny, const unsigned int* pos, double* EC) const
{
	double EffMat[4];
	Calc_EffMatPos(ny,pos,EffMat);

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

bool Operator::Calc_EffMatPos(int ny, const unsigned int* pos, double* EffMat) const
{
	int n=ny;
	double coord[3];
	double shiftCoord[3];
	int nP = (n+1)%3;
	int nPP = (n+2)%3;
	coord[0] = discLines[0][pos[0]];
	coord[1] = discLines[1][pos[1]];
	coord[2] = discLines[2][pos[2]];
	double delta=MainOp->GetIndexDelta(n,pos[n]);
	double deltaP=MainOp->GetIndexDelta(nP,pos[nP]);
	double deltaPP=MainOp->GetIndexDelta(nPP,pos[nPP]);
	double delta_M=MainOp->GetIndexDelta(n,pos[n]-1);
	double deltaP_M=MainOp->GetIndexDelta(nP,pos[nP]-1);
	double deltaPP_M=MainOp->GetIndexDelta(nPP,pos[nPP]-1);

	int loc_pos[3]={pos[0],pos[1],pos[2]};
	double A_n;
	double area = 0;

	//******************************* epsilon,kappa averaging *****************************//
	//shift up-right
	shiftCoord[n] = coord[n]+delta*0.5;
	shiftCoord[nP] = coord[nP]+deltaP*0.25;
	shiftCoord[nPP] = coord[nPP]+deltaPP*0.25;
	A_n = GetNodeArea(ny,(unsigned int*)loc_pos,true);
//	{
//		cerr << ny << " " << pos[0] << " " <<  pos[1] << " " <<  pos[2] << ": " << A_n << endl;
//		exit(0);
//	}
	CSProperties* prop = CSX->GetPropertyByCoordPriority(shiftCoord,CSProperties::MATERIAL,true);
	if (prop)
	{
		CSPropMaterial* mat = prop->ToMaterial();
		EffMat[0] = mat->GetEpsilonWeighted(n,shiftCoord)*A_n;
		EffMat[1] = mat->GetKappaWeighted(n,shiftCoord)*A_n;
	}
	else
	{
		EffMat[0] = 1*A_n;
		EffMat[1] = 0;
	}
	area+=A_n;

	//shift up-left
	shiftCoord[n] = coord[n]+delta*0.5;
	shiftCoord[nP] = coord[nP]-deltaP_M*0.25;
	shiftCoord[nPP] = coord[nPP]+deltaPP*0.25;

	--loc_pos[nP];
	A_n = GetNodeArea(ny,(unsigned int*)loc_pos,true);
//	cerr << A_n << endl;
	prop = CSX->GetPropertyByCoordPriority(shiftCoord,CSProperties::MATERIAL,true);
	if (prop)
	{
		CSPropMaterial* mat = prop->ToMaterial();
		EffMat[0] += mat->GetEpsilonWeighted(n,shiftCoord)*A_n;
		EffMat[1] += mat->GetKappaWeighted(n,shiftCoord)*A_n;
	}
	else
	{
		EffMat[0] += 1*A_n;
		EffMat[1] += 0;
	}
	area+=A_n;

	//shift down-right
	shiftCoord[n] = coord[n]+delta*0.5;
	shiftCoord[nP] = coord[nP]+deltaP*0.25;
	shiftCoord[nPP] = coord[nPP]-deltaPP_M*0.25;
	++loc_pos[nP];
	--loc_pos[nPP];
	A_n = GetNodeArea(ny,(unsigned int*)loc_pos,true);
	prop = CSX->GetPropertyByCoordPriority(shiftCoord,CSProperties::MATERIAL,true);
	if (prop)
	{
		CSPropMaterial* mat = prop->ToMaterial();
		EffMat[0] += mat->GetEpsilonWeighted(n,shiftCoord)*A_n;
		EffMat[1] += mat->GetKappaWeighted(n,shiftCoord)*A_n;
	}
	else
	{
		EffMat[0] += 1*A_n;
		EffMat[1] += 0;
	}
	area+=A_n;

	//shift down-left
	shiftCoord[n] = coord[n]+delta*0.5;
	shiftCoord[nP] = coord[nP]-deltaP_M*0.25;
	shiftCoord[nPP] = coord[nPP]-deltaPP_M*0.25;
	--loc_pos[nP];
	A_n = GetNodeArea(ny,(unsigned int*)loc_pos,true);
	prop = CSX->GetPropertyByCoordPriority(shiftCoord,CSProperties::MATERIAL,true);
	if (prop)
	{
		CSPropMaterial* mat = prop->ToMaterial();
		EffMat[0] += mat->GetEpsilonWeighted(n,shiftCoord)*A_n;
		EffMat[1] += mat->GetKappaWeighted(n,shiftCoord)*A_n;
	}
	else
	{
		EffMat[0] += 1*A_n;
		EffMat[1] += 0;
	}
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
	double delta_ny = GetNodeWidth(n,(unsigned int*)loc_pos,true);
	prop = CSX->GetPropertyByCoordPriority(shiftCoord,CSProperties::MATERIAL,true);
	if (prop)
	{
		CSPropMaterial* mat = prop->ToMaterial();
		EffMat[2] = delta_ny / mat->GetMueWeighted(n,shiftCoord);
		if (mat->GetSigmaWeighted(n,shiftCoord))
			EffMat[3] = delta_ny / mat->GetSigmaWeighted(n,shiftCoord);
		else
			EffMat[3] = 0;
	}
	else
	{
		EffMat[2] = delta_ny;
		EffMat[3] = 0;
	}
	length=delta_ny;

	//shift up
	shiftCoord[n] = coord[n]+delta*0.25;
	shiftCoord[nP] = coord[nP]+deltaP*0.5;
	shiftCoord[nPP] = coord[nPP]+deltaPP*0.5;
	++loc_pos[n];
	delta_ny = GetNodeWidth(n,(unsigned int*)loc_pos,true);
	prop = CSX->GetPropertyByCoordPriority(shiftCoord,CSProperties::MATERIAL,true);
	if (prop)
	{
		CSPropMaterial* mat = prop->ToMaterial();
		EffMat[2] += delta_ny / mat->GetMueWeighted(n,shiftCoord);
		if (mat->GetSigmaWeighted(n,shiftCoord))
			EffMat[3] += delta_ny/mat->GetSigmaWeighted(n,shiftCoord);
		else
			EffMat[3] = 0;
	}
	else
	{
		EffMat[2] += 1*delta_ny;
		EffMat[3] = 0;
	}
	length+=delta_ny;

	EffMat[2] = length * __MUE0__ / EffMat[2];
	if (EffMat[3]) EffMat[3]=length / EffMat[3];

	for (int n=0; n<4; ++n)
		if (isnan(EffMat[n]) || isinf(EffMat[n]))
		{
			cerr << "Operator::Calc_EffMatPos: An effective material parameter is not a valid result, this should NOT have happend... exit..." << endl;
			exit(0);
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
		EC_C[n] = new double[MainOp->GetSize()];
		EC_G[n] = new double[MainOp->GetSize()];
		EC_L[n] = new double[MainOp->GetSize()];
		EC_R[n] = new double[MainOp->GetSize()];
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
	unsigned int ipos;
	unsigned int pos[3];
	double inEC[4];
	for (int n=0; n<3; ++n)
	{
		for (pos[2]=0; pos[2]<numLines[2]; ++pos[2])
		{
			for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
			{
				for (pos[0]=0; pos[0]<numLines[0]; ++pos[0])
				{
					Calc_ECPos(n,pos,inEC);
					ipos = MainOp->SetPos(pos[0],pos[1],pos[2]);
					EC_C[n][ipos]=inEC[0];
					EC_G[n][ipos]=inEC[1];
					EC_L[n][ipos]=inEC[2];
					EC_R[n][ipos]=inEC[3];
				}
			}
		}
	}

	return true;
}


double Operator::CalcTimestep()
{
#if 1  //use the new timestep-calc (1) or the old one (0)
	return CalcTimestep_Var3(); //the biggest one for cartesian meshes
#else
	return CalcTimestep_Var1();
#endif
}

////Berechnung nach Andreas Rennings Dissertation 2008, Seite 66, Formel 4.52
double Operator::CalcTimestep_Var1()
{
	m_Used_TS_Name = string("Rennings_1");
//	cout << "Operator::CalcTimestep(): Using timestep algorithm by Andreas Rennings, Dissertation @ University Duisburg-Essen, 2008, pp. 66, eq. 4.52" << endl;
	dT=1e200;
	double newT;
	unsigned int pos[3];
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
					if ((newT<dT) && (newT>0.0)) dT=newT;
				}
			}
		}
	}
	if (dT==0)
	{
		cerr << "Operator::CalcTimestep: Timestep is zero... this is not supposed to happen!!! exit!" << endl;
		exit(3);
	}
//	cerr << "Operator Timestep: " << dT << endl;
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
						dT=newT;
				}
			}
		}
	}
	if (dT==0)
	{
		cerr << "Operator::CalcTimestep: Timestep is zero... this is not supposed to happen!!! exit!" << endl;
		exit(3);
	}
//	cerr << "Operator Timestep: " << dT << endl;
	return 0;
}

bool Operator::CalcFieldExcitation()
{
	if (dT==0)
		return false;
	if (Exc==0)
		return false;

	unsigned int pos[3];
	double delta[3];
	double amp=0;

	vector<unsigned int> volt_vIndex[3];
	vector<FDTD_FLOAT> volt_vExcit;
	vector<unsigned int> volt_vDelay;
	vector<unsigned int> volt_vDir;
	double volt_coord[3];

	vector<unsigned int> curr_vIndex[3];
	vector<FDTD_FLOAT> curr_vExcit;
	vector<unsigned int> curr_vDelay;
	vector<unsigned int> curr_vDir;
	double curr_coord[3];

	vector<CSProperties*> vec_prop = CSX->GetPropertyByType(CSProperties::ELECTRODE);

	if (vec_prop.size()==0)
	{
		cerr << "Operator::CalcFieldExcitation: Warning, no excitation properties found" << endl;
		return false;
	}

	CSPropElectrode* elec=NULL;
	CSProperties* prop=NULL;
	int priority=0;

	for (pos[2]=0; pos[2]<numLines[2]; ++pos[2])
	{
		delta[2]=fabs(MainOp->GetIndexDelta(2,pos[2]));
		for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
		{
			delta[1]=fabs(MainOp->GetIndexDelta(1,pos[1]));
			for (pos[0]=0; pos[0]<numLines[0]; ++pos[0])
			{
				delta[0]=fabs(MainOp->GetIndexDelta(0,pos[0]));

				//electric field excite
				for (int n=0; n<3; ++n)
				{
					volt_coord[0] = discLines[0][pos[0]];
					volt_coord[1] = discLines[1][pos[1]];
					volt_coord[2] = discLines[2][pos[2]];
					volt_coord[n]+=delta[n]*0.5;
					for (size_t p=0; p<vec_prop.size(); ++p)
					{
						prop = vec_prop.at(p);
						elec = prop->ToElectrode();
						if (elec==NULL)
							continue;
						if (prop->CheckCoordInPrimitive(volt_coord,priority,true))
						{
							if ((elec->GetActiveDir(n)) && ( (elec->GetExcitType()==0) || (elec->GetExcitType()==1) ))//&& (pos[n]<numLines[n]-1))
							{
								amp = elec->GetWeightedExcitation(n,volt_coord)*GetEdgeLength(n,pos);// delta[n]*gridDelta;
								if (amp!=0)
								{
									volt_vExcit.push_back(amp);
									volt_vDelay.push_back((unsigned int)(elec->GetDelay()/dT));
									volt_vDir.push_back(n);
									volt_vIndex[0].push_back(pos[0]);
									volt_vIndex[1].push_back(pos[1]);
									volt_vIndex[2].push_back(pos[2]);
								}
								if (elec->GetExcitType()==1) //hard excite
								{
									SetVV(n,pos[0],pos[1],pos[2], 0 );
									SetVI(n,pos[0],pos[1],pos[2], 0 );
								}
							}
						}
					}
				}

				//magnetic field excite
				for (int n=0; n<3; ++n)
				{
					if ((pos[0]>=numLines[0]-1) || (pos[1]>=numLines[1]-1) || (pos[2]>=numLines[2]-1))
						continue;  //skip the last H-Line which is outside the FDTD-domain
					int nP = (n+1)%3;
					int nPP = (n+2)%3;
					curr_coord[0] = discLines[0][pos[0]];
					curr_coord[1] = discLines[1][pos[1]];
					curr_coord[2] = discLines[2][pos[2]];
					curr_coord[nP] +=delta[nP]*0.5;
					curr_coord[nPP] +=delta[nPP]*0.5;
					for (size_t p=0; p<vec_prop.size(); ++p)
					{
						prop = vec_prop.at(p);
						elec = prop->ToElectrode();
						if (elec==NULL)
							continue;
						if (prop->CheckCoordInPrimitive(curr_coord,priority,true))
						{
							if ((elec->GetActiveDir(n)) && ( (elec->GetExcitType()==2) || (elec->GetExcitType()==3) ))
							{
								amp = elec->GetWeightedExcitation(n,curr_coord)*GetEdgeLength(n,pos,true);// delta[n]*gridDelta;
								if (amp!=0)
								{
									curr_vExcit.push_back(amp);
									curr_vDelay.push_back((unsigned int)(elec->GetDelay()/dT));
									curr_vDir.push_back(n);
									curr_vIndex[0].push_back(pos[0]);
									curr_vIndex[1].push_back(pos[1]);
									curr_vIndex[2].push_back(pos[2]);
								}
								if (elec->GetExcitType()==3) //hard excite
								{
									SetII(n,pos[0],pos[1],pos[2], 0 );
									SetIV(n,pos[0],pos[1],pos[2], 0 );
								}
							}
						}
					}
				}

			}
		}
	}

	//special treatment for primitives of type curve (treated as wires) see also Calc_PEC
	double p1[3];
	double p2[3];
	double deltaN=0.0;
	struct Grid_Path path;
	for (size_t p=0; p<vec_prop.size(); ++p)
	{
		prop = vec_prop.at(p);
		elec = prop->ToElectrode();
		for (size_t n=0; n<prop->GetQtyPrimitives(); ++n)
		{
			CSPrimitives* prim = prop->GetPrimitive(n);
			CSPrimCurve* curv = prim->ToCurve();
			if (curv)
			{
				for (size_t i=1; i<curv->GetNumberOfPoints(); ++i)
				{
					curv->GetPoint(i-1,p1);
					curv->GetPoint(i,p2);
					path = FindPath(p1,p2);
					if (path.dir.size()>0)
						prim->SetPrimitiveUsed(true);
					for (size_t t=0; t<path.dir.size(); ++t)
					{
						n = path.dir.at(t);
						pos[0] = path.posPath[0].at(t);
						pos[1] = path.posPath[1].at(t);
						pos[2] = path.posPath[2].at(t);
						MainOp->SetPos(pos[0],pos[1],pos[2]);
						deltaN=fabs(MainOp->GetIndexDelta(n,pos[n]));
						volt_coord[0] = discLines[0][pos[0]];
						volt_coord[1] = discLines[1][pos[1]];
						volt_coord[2] = discLines[2][pos[2]];
						volt_coord[n] += 0.5*deltaN;
//						cerr << n << " " << coord[0] << " " << coord[1] << " " << coord[2] << endl;
						if (elec!=NULL)
						{
							if ((elec->GetActiveDir(n)) && (pos[n]<numLines[n]-1) && ( (elec->GetExcitType()==0) || (elec->GetExcitType()==1) ))
							{
								amp = elec->GetWeightedExcitation(n,volt_coord)*deltaN*gridDelta;
								if (amp!=0)
								{
									volt_vExcit.push_back(amp);
									volt_vDelay.push_back((unsigned int)(elec->GetDelay()/dT));
									volt_vDir.push_back(n);
									volt_vIndex[0].push_back(pos[0]);
									volt_vIndex[1].push_back(pos[1]);
									volt_vIndex[2].push_back(pos[2]);
								}
								if (elec->GetExcitType()==1) //hard excite
								{
									SetVV(n,pos[0],pos[1],pos[2], 0 );
									SetVI(n,pos[0],pos[1],pos[2], 0 );
								}
							}
						}
					}
				}
			}
		}
	}


	// set voltage excitations
	Exc->setupVoltageExcitation( volt_vIndex, volt_vExcit, volt_vDelay, volt_vDir );

	// set current excitations
	Exc->setupCurrentExcitation( curr_vIndex, curr_vExcit, curr_vDelay, curr_vDir );

	return true;
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
	double delta;
	unsigned int pos[3];
	for (pos[0]=startX; pos[0]<=stopX; ++pos[0])
	{
		for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
		{
			for (pos[2]=0; pos[2]<numLines[2]; ++pos[2])
			{
				for (int n=0; n<3; ++n)
				{
					coord[0] = discLines[0][pos[0]];
					coord[1] = discLines[1][pos[1]];
					coord[2] = discLines[2][pos[2]];
					delta=MainOp->GetIndexDelta(n,pos[n]);
					coord[n]= discLines[n][pos[n]] + delta*0.5;
					CSProperties* prop = CSX->GetPropertyByCoordPriority(coord, (CSProperties::PropertyType)(CSProperties::MATERIAL | CSProperties::METAL), true);
					if (prop)
					{
						if (prop->GetType()==CSProperties::METAL) //set to PEC
						{
							SetVV(n,pos[0],pos[1],pos[2], 0 );
							SetVI(n,pos[0],pos[1],pos[2], 0 );
							++counter[n];
							//							cerr << "CartOperator::CalcPEC: PEC found at " << pos[0] << " ; "  << pos[1] << " ; " << pos[2] << endl;
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
	struct Grid_Path path;
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
					curv->GetPoint(i-1,p1);
					curv->GetPoint(i,p2);
					path = FindPath(p1,p2);
					if (path.dir.size()>0)
						prim->SetPrimitiveUsed(true);
					for (size_t t=0; t<path.dir.size(); ++t)
					{
//						cerr << path.dir.at(t) << " " << path.posPath[0].at(t) << " " << path.posPath[1].at(t) << " " << path.posPath[2].at(t) << endl;
						SetVV(path.dir.at(t),path.posPath[0].at(t),path.posPath[1].at(t),path.posPath[2].at(t), 0 );
						SetVI(path.dir.at(t),path.posPath[0].at(t),path.posPath[1].at(t),path.posPath[2].at(t), 0 );
						++m_Nr_PEC[path.dir.at(t)];
					}
//					cerr << "found path size: " << path.dir.size() << endl;
				}
			}
		}
	}
}

void Operator::AddExtension(Operator_Extension* op_ext)
{
	m_Op_exts.push_back(op_ext);
}
