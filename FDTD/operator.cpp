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
#include "operator_extension.h"
#include "processfields.h"
#include "tools/array_ops.h"
#include "fparser.hh"

Operator* Operator::New()
{
	cout << "Create FDTD operator" << endl;
	Operator* op = new Operator();
	op->Init();
	return op;
}

Operator::Operator()
{
	m_MeshType = ProcessFields::CARTESIAN_MESH;
	Exc = 0;
}

Operator::~Operator()
{
	for (size_t n=0;n<m_Op_exts.size();++n)
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
	CSX = NULL;

	vv=NULL;
	vi=NULL;
	iv=NULL;
	ii=NULL;
	for (int n=0;n<3;++n)
		discLines[n]=NULL;

	MainOp=NULL;
	DualOp=NULL;

	for (int n=0;n<3;++n)
	{
		EC_C[n]=NULL;
		EC_G[n]=NULL;
		EC_L[n]=NULL;
		EC_R[n]=NULL;
	}

	for (int n=0;n<6;++n)
		m_BC[n]=0;

	Exc = 0;
}

void Operator::Reset()
{
	Delete_N_3DArray(vv,numLines);
	Delete_N_3DArray(vi,numLines);
	Delete_N_3DArray(iv,numLines);
	Delete_N_3DArray(ii,numLines);
	for (int n=0;n<3;++n)
		delete[] discLines[n];
	delete MainOp;
	delete DualOp;
	for (int n=0;n<3;++n)
	{
		delete[] EC_C[n];
		delete[] EC_G[n];
		delete[] EC_L[n];
		delete[] EC_R[n];
	}

	delete Exc;

	Init();
}

string Operator::GetDirName(int ny) const
{
	if (ny==0) return "x";
	if (ny==1) return "y";
	if (ny==2) return "z";
	return "";
}

double Operator::GetMeshDelta(int n, const unsigned int* pos, bool dualMesh) const
{
	if ((n<0) || (n>2)) return 0.0;
	int i_pos[] = {pos[0],pos[1],pos[2]};
	return GetMeshDelta(n,i_pos,dualMesh);
}

double Operator::GetMeshDelta(int n, const int* pos, bool dualMesh) const
{
	if ((n<0) || (n>2)) return 0.0;
	if (dualMesh==false)
		return fabs(MainOp->GetIndexDelta(n,pos[n]))*gridDelta;
	else
		return fabs(MainOp->GetIndexWidth(n,pos[n]))*gridDelta;
}

double Operator::GetDiscLine(int n, unsigned int pos, bool dualMesh) const
{
	return GetDiscLine(n,(int)pos,dualMesh);
}

double Operator::GetDiscLine(int n, int pos, bool dualMesh) const
{
	if ((n<0) || (n>2)) return 0.0;
	if ((pos<0) || (pos>=(int)numLines[n])) return 0.0;
	if (dualMesh==false)
		return discLines[n][pos];
	else
	{
		if (pos<(int)numLines[n]-1)
			return 0.5*(discLines[n][pos+1]+discLines[n][pos]);
		else
			return 0.5*(discLines[n][pos]+discLines[n][pos-1]);
	}
}

bool Operator::SnapToMesh(double* dcoord, unsigned int* uicoord, bool lower, bool* inside)
{
	bool ok=true;
	unsigned int numLines[3];
	for (int n=0;n<3;++n)
	{
		numLines[n] = GetNumberOfLines(n);
		if (inside) //set defaults
			inside[n] = true;
		uicoord[n]=0;
		if (dcoord[n]<discLines[n][0])
		{
			ok=false;uicoord[n]=0;
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
			for (unsigned int i=1;i<numLines[n];++i)
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
		for (int n=0;n<3;++n) //check all 6 surrounding points
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
	cout << "Voltage excitations\t: " << Exc->E_Count << endl;
	cout << "Current excitations\t: " << Exc->Curr_Count << endl;
	cout << "-----------------------------------" << endl;
	cout << "Number of PEC edges\t: " << m_Nr_PEC[0]+m_Nr_PEC[1]+m_Nr_PEC[2] << endl;
	cout << "in " << GetDirName(0) << " direction\t\t: " << m_Nr_PEC[0] << endl;
	cout << "in " << GetDirName(1) << " direction\t\t: " << m_Nr_PEC[1] << endl;
	cout << "in " << GetDirName(2) << " direction\t\t: " << m_Nr_PEC[2] << endl;
	cout << "-----------------------------------" << endl;
	cout << "Timestep (s)\t\t: " << dT << endl;
	cout << "Nyquist criteria (TS)\t: " << Exc->GetNyquistNum() << endl;
	cout << "Nyquist criteria (s)\t: " << Exc->GetNyquistNum()*dT << endl;
	cout << "Excitation Length (TS)\t: " << Exc->Length << endl;
	cout << "Excitation Length (s)\t: " << Exc->Length*dT << endl;
	cout << "-----------------------------------" << endl;
}


void Operator::DumpOperator2File(string filename)
{
	ofstream file(filename.c_str(),ios_base::out);
	if (file.is_open()==false)
	{
		cerr << "Operator::DumpOperator2File: Can't open file: " << filename << endl;
		return;
	}

	cout << "Dumping FDTD operator information to vtk file: " << filename << " ..." << flush ;

	FDTD_FLOAT**** exc = Create_N_3DArray<FDTD_FLOAT>(numLines);
	if (Exc) {
		for (unsigned int n=0;n<Exc->E_Count;++n)
			exc[Exc->E_dir[n]][Exc->E_index[0][n]][Exc->E_index[1][n]][Exc->E_index[2][n]] = Exc->E_amp[n];
	}

	string names[] = {"vv", "vi", "iv" , "ii", "exc"};
	FDTD_FLOAT**** array[] = {vv,vi,iv,ii,exc};

	ProcessFields::DumpMultiVectorArray2VTK(file, names , array , 5, discLines, numLines, 6, "Operator dump" , (ProcessFields::MeshType)m_MeshType);

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
	if (!file.is_open()) {
		cerr << "Operator::DumpPEC2File: Can't open file: " << filename << endl;
		return;
	}

	cout << "Dumping PEC information to vtk file: " << filename << " ..." << flush;

	FDTD_FLOAT**** pec = Create_N_3DArray<FDTD_FLOAT>( numLines );
	unsigned int pos[3];

	for (pos[0]=0; pos[0]<numLines[0]; pos[0]++) {
		for (pos[1]=0; pos[1]<numLines[1]; pos[1]++) {
			for (pos[2]=0; pos[2]<numLines[2]; pos[2]++) {
				if ((GetVV(0,pos[0],pos[1],pos[2]) == 0) && (GetVI(0,pos[0],pos[1],pos[2]) == 0))
					pec[0][pos[0]][pos[1]][pos[2]] = MainOp->GetIndexDelta( 0, pos[0] ); // PEC-x found
				if ((GetVV(1,pos[0],pos[1],pos[2]) == 0) && (GetVI(1,pos[0],pos[1],pos[2]) == 0))
					pec[1][pos[0]][pos[1]][pos[2]] = MainOp->GetIndexDelta( 1, pos[1] ); // PEC-y found
				if ((GetVV(2,pos[0],pos[1],pos[2]) == 0) && (GetVI(2,pos[0],pos[1],pos[2]) == 0))
					pec[2][pos[0]][pos[1]][pos[2]] = MainOp->GetIndexDelta( 2, pos[2] ); // PEC-z found
			}
		}
	}

	ProcessFields::DumpVectorArray2VTK( file, "PEC", pec, discLines, numLines, 6, "PEC dump" , (ProcessFields::MeshType)m_MeshType );

	file.close();

	cout << " done!" << endl;
}

void Operator::DumpMaterial2File(string filename)
{
	ofstream file(filename.c_str(),ios_base::out);
	if (file.is_open()==false)
	{
		cerr << "Operator::DumpMaterial2File: Can't open file: " << filename << endl;
		return;
	}

	cout << "Dumping material information to vtk file: " << filename << " ..."  << flush;

	FDTD_FLOAT*** epsilon;
	FDTD_FLOAT*** mue;
	FDTD_FLOAT*** kappa;
	FDTD_FLOAT*** sigma;
	unsigned int pos[3];
	double inMat[4];

	epsilon = Create3DArray<FDTD_FLOAT>( numLines);
	mue = Create3DArray<FDTD_FLOAT>( numLines);
	kappa = Create3DArray<FDTD_FLOAT>( numLines);
	sigma = Create3DArray<FDTD_FLOAT>( numLines);
	for (pos[0]=0;pos[0]<numLines[0];++pos[0])
	{
		for (pos[1]=0;pos[1]<numLines[1];++pos[1])
		{
			for (pos[2]=0;pos[2]<numLines[2];++pos[2])
			{
				for (int n=0;n<3;++n)
				{
					Calc_EffMatPos(n, pos, inMat);
					epsilon[pos[0]][pos[1]][pos[2]]+=inMat[0]/__EPS0__;
					mue[pos[0]][pos[1]][pos[2]]+=inMat[2]/__MUE0__;
					kappa[pos[0]][pos[1]][pos[2]]+=inMat[1];
					sigma[pos[0]][pos[1]][pos[2]]+=inMat[3];
				}
				epsilon[pos[0]][pos[1]][pos[2]]/=3;
				mue[pos[0]][pos[1]][pos[2]]/=3;
				kappa[pos[0]][pos[1]][pos[2]]/=3;
				sigma[pos[0]][pos[1]][pos[2]]/=3;
			}
		}
	}

	string names[] = {"epsilon","mue","kappa","sigma"};
	FDTD_FLOAT*** array[] = {epsilon,mue,kappa,sigma};
	ProcessFields::DumpMultiScalarArray2VTK(file, names, array, 4, discLines, numLines,  6, "Material dump" , (ProcessFields::MeshType)m_MeshType);
	Delete3DArray(epsilon,numLines);
	Delete3DArray(mue,numLines);
	Delete3DArray(kappa,numLines);
	Delete3DArray(sigma,numLines);
	file.close();

	cout << " done!" << endl;
}

bool Operator::SetGeometryCSX(ContinuousStructure* geo)
{
	if (geo==NULL) return false;

	Reset();
	CSX = geo;

	CSRectGrid* grid=CSX->GetGrid();
	for (int n=0;n<3;++n)
	{
		discLines[n] = grid->GetLines(n,discLines[n],numLines[n],true);
		if (n==1)
		if (numLines[n]<3) {cerr << "CartOperator::SetGeometryCSX: you need at least 3 disc-lines in every direction (3D!)!!!" << endl; Reset(); return false;}
	}
	MainOp = new AdrOp(numLines[0],numLines[1],numLines[2]);
	MainOp->SetGrid(discLines[0],discLines[1],discLines[2]);
	if (grid->GetDeltaUnit()<=0)  {cerr << "CartOperator::SetGeometryCSX: grid delta unit must not be <=0 !!!" << endl; Reset(); return false;}
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
		GetVV(n,pos[0],pos[1],pos[2]) = (1-dT*EC_G[n][i]/2/EC_C[n][i])/(1+dT*EC_G[n][i]/2/EC_C[n][i]);
		GetVI(n,pos[0],pos[1],pos[2]) = (dT/EC_C[n][i])/(1+dT*EC_G[n][i]/2/EC_C[n][i]);
	}
	else
	{
		GetVV(n,pos[0],pos[1],pos[2]) = 0;
		GetVI(n,pos[0],pos[1],pos[2]) = 0;
	}
	if (EC_L[n][i]>0)
	{
		GetII(n,pos[0],pos[1],pos[2]) = (1-dT*EC_R[n][i]/2/EC_L[n][i])/(1+dT*EC_R[n][i]/2/EC_L[n][i]);
		GetIV(n,pos[0],pos[1],pos[2]) = (dT/EC_L[n][i])/(1+dT*EC_R[n][i]/2/EC_L[n][i]);
	}
	else
	{
		GetII(n,pos[0],pos[1],pos[2]) = 0;
		GetIV(n,pos[0],pos[1],pos[2]) = 0;
	}
}

int Operator::CalcECOperator()
{
	Init_EC();

	if (Calc_EC()==0)
		return -1;

	CalcTimestep();

	InitOperator();

	unsigned int pos[3];

	for (int n=0;n<3;++n)
	{
		for (pos[0]=0;pos[0]<numLines[0];++pos[0])
		{
			for (pos[1]=0;pos[1]<numLines[1];++pos[1])
			{
				for (pos[2]=0;pos[2]<numLines[2];++pos[2])
				{
					Calc_ECOperatorPos(n,pos);
				}
			}
		}
	}

	//all information available for extension... create now...
	for (size_t n=0;n<m_Op_exts.size();++n)
		m_Op_exts.at(n)->BuildExtension();

	//cleanup
	for (int n=0;n<3;++n)
	{
		delete[] EC_C[n];EC_C[n]=NULL;
		delete[] EC_G[n];EC_G[n]=NULL;
		delete[] EC_L[n];EC_L[n]=NULL;
		delete[] EC_R[n];EC_R[n]=NULL;
	}

	//Always apply PEC to all boundary's
	bool PEC[6]={1,1,1,1,1,1};
	ApplyElectricBC(PEC);

	InitExcitation();

	if (CalcFieldExcitation()==false) return -1;

	CalcPEC();

	bool PMC[6];
	for (int n=0;n<6;++n)
		PMC[n] = m_BC[n]==1;
	ApplyMagneticBC(PMC);

	return 0;
}

void Operator::ApplyElectricBC(bool* dirs)
{
	if (dirs==NULL) return;
	unsigned int pos[3];
	for (int n=0;n<3;++n)
	{
		int nP = (n+1)%3;
		int nPP = (n+2)%3;
		for (pos[nP]=0;pos[nP]<numLines[nP];++pos[nP])
		{
			for (pos[nPP]=0;pos[nPP]<numLines[nPP];++pos[nPP])
			{
				pos[n]=0;
				GetVV(nP,pos[0],pos[1],pos[2]) *= (FDTD_FLOAT)!dirs[2*n];
				GetVI(nP,pos[0],pos[1],pos[2]) *= (FDTD_FLOAT)!dirs[2*n];
				GetVV(nPP,pos[0],pos[1],pos[2]) *= (FDTD_FLOAT)!dirs[2*n];
				GetVI(nPP,pos[0],pos[1],pos[2]) *= (FDTD_FLOAT)!dirs[2*n];
				pos[n]=numLines[n]-1;
				GetVV(n,pos[0],pos[1],pos[2]) *= (FDTD_FLOAT)!dirs[2*n+1];
				GetVI(n,pos[0],pos[1],pos[2]) *= (FDTD_FLOAT)!dirs[2*n+1];
				GetVV(nP,pos[0],pos[1],pos[2]) *= (FDTD_FLOAT)!dirs[2*n+1];
				GetVI(nP,pos[0],pos[1],pos[2]) *= (FDTD_FLOAT)!dirs[2*n+1];
				GetVV(nPP,pos[0],pos[1],pos[2]) *= (FDTD_FLOAT)!dirs[2*n+1];
				GetVI(nPP,pos[0],pos[1],pos[2]) *= (FDTD_FLOAT)!dirs[2*n+1];
			}
		}
	}
}

void Operator::ApplyMagneticBC(bool* dirs)
{
	if (dirs==NULL) return;
	unsigned int pos[3];
	for (int n=0;n<3;++n)
	{
		int nP = (n+1)%3;
		int nPP = (n+2)%3;
		for (pos[nP]=0;pos[nP]<numLines[nP];++pos[nP])
		{
			for (pos[nPP]=0;pos[nPP]<numLines[nPP];++pos[nPP])
			{
				pos[n]=0;
				GetII(n,pos[0],pos[1],pos[2]) *= (FDTD_FLOAT)!dirs[2*n];
				GetIV(n,pos[0],pos[1],pos[2]) *= (FDTD_FLOAT)!dirs[2*n];
				GetII(nP,pos[0],pos[1],pos[2]) *= (FDTD_FLOAT)!dirs[2*n];
				GetIV(nP,pos[0],pos[1],pos[2]) *= (FDTD_FLOAT)!dirs[2*n];
				GetII(nPP,pos[0],pos[1],pos[2]) *= (FDTD_FLOAT)!dirs[2*n];
				GetIV(nPP,pos[0],pos[1],pos[2]) *= (FDTD_FLOAT)!dirs[2*n];

				pos[n]=numLines[n]-2;
				GetII(nP,pos[0],pos[1],pos[2]) *= (FDTD_FLOAT)!dirs[2*n+1];
				GetIV(nP,pos[0],pos[1],pos[2]) *= (FDTD_FLOAT)!dirs[2*n+1];
				GetII(nPP,pos[0],pos[1],pos[2]) *= (FDTD_FLOAT)!dirs[2*n+1];
				GetIV(nPP,pos[0],pos[1],pos[2]) *= (FDTD_FLOAT)!dirs[2*n+1];
			}
		}
	}
}


bool Operator::Calc_ECPos(int n, unsigned int* pos, double* inEC)
{
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

	//******************************* epsilon,kappa averaging *****************************//
	//shift up-right
	shiftCoord[n] = coord[n]+delta*0.5;
	shiftCoord[nP] = coord[nP]+deltaP*0.25;
	shiftCoord[nPP] = coord[nPP]+deltaPP*0.25;
	CSProperties* prop = CSX->GetPropertyByCoordPriority(shiftCoord,CSProperties::MATERIAL);
	if (prop)
	{
		CSPropMaterial* mat = prop->ToMaterial();
		inEC[0] = mat->GetEpsilonWeighted(n,shiftCoord)*fabs(deltaP*deltaPP);
		inEC[1] = mat->GetKappaWeighted(n,shiftCoord)*fabs(deltaP*deltaPP);
	}
	else
	{
		inEC[0] = 1*fabs(deltaP*deltaPP);
		inEC[1] = 0;
	}
	//shift up-left
	shiftCoord[n] = coord[n]+delta*0.5;
	shiftCoord[nP] = coord[nP]-deltaP_M*0.25;
	shiftCoord[nPP] = coord[nPP]+deltaPP*0.25;
	prop = CSX->GetPropertyByCoordPriority(shiftCoord,CSProperties::MATERIAL);
	if (prop)
	{
		CSPropMaterial* mat = prop->ToMaterial();
		inEC[0] += mat->GetEpsilonWeighted(n,shiftCoord)*fabs(deltaP_M*deltaPP);
		inEC[1] += mat->GetKappaWeighted(n,shiftCoord)*fabs(deltaP_M*deltaPP);
	}
	else
	{
		inEC[0] += 1*fabs(deltaP_M*deltaPP);
		inEC[1] += 0;
	}

	//shift down-right
	shiftCoord[n] = coord[n]+delta*0.5;
	shiftCoord[nP] = coord[nP]+deltaP*0.25;
	shiftCoord[nPP] = coord[nPP]-deltaPP_M*0.25;
	prop = CSX->GetPropertyByCoordPriority(shiftCoord,CSProperties::MATERIAL);
	if (prop)
	{
		CSPropMaterial* mat = prop->ToMaterial();
		inEC[0] += mat->GetEpsilonWeighted(n,shiftCoord)*fabs(deltaP*deltaPP_M);
		inEC[1] += mat->GetKappaWeighted(n,shiftCoord)*fabs(deltaP*deltaPP_M);
	}
	else
	{
		inEC[0] += 1*fabs(deltaP*deltaPP_M);
		inEC[1] += 0;
	}

	//shift down-left
	shiftCoord[n] = coord[n]+delta*0.5;
	shiftCoord[nP] = coord[nP]-deltaP_M*0.25;
	shiftCoord[nPP] = coord[nPP]-deltaPP_M*0.25;
	prop = CSX->GetPropertyByCoordPriority(shiftCoord,CSProperties::MATERIAL);
	if (prop)
	{
		CSPropMaterial* mat = prop->ToMaterial();
		inEC[0] += mat->GetEpsilonWeighted(n,shiftCoord)*fabs(deltaP_M*deltaPP_M);
		inEC[1] += mat->GetKappaWeighted(n,shiftCoord)*fabs(deltaP_M*deltaPP_M);
	}
	else
	{
		inEC[0] += 1*fabs(deltaP_M*deltaPP_M);
		inEC[1] += 0;
	}

	inEC[0]*=gridDelta/fabs(delta)/4.0*__EPS0__;
	inEC[1]*=gridDelta/fabs(delta)/4.0;

	//******************************* mu,sigma averaging *****************************//
	//shift down
	shiftCoord[n] = coord[n]-delta_M*0.25;
	shiftCoord[nP] = coord[nP]+deltaP*0.5;
	shiftCoord[nPP] = coord[nPP]+deltaPP*0.5;
	prop = CSX->GetPropertyByCoordPriority(shiftCoord,CSProperties::MATERIAL);
	if (prop)
	{
		CSPropMaterial* mat = prop->ToMaterial();
		inEC[2] = fabs(delta_M) / mat->GetMueWeighted(n,shiftCoord);
		if (mat->GetSigma(n))
			inEC[3] = fabs(delta_M) / mat->GetSigmaWeighted(n,shiftCoord);
		else
			inEC[3] = 0;
	}
	else
	{
		inEC[2] = fabs(delta_M);
		inEC[3] = 0;
	}
	//shift up
	shiftCoord[n] = coord[n]+delta*0.25;
	shiftCoord[nP] = coord[nP]+deltaP*0.5;
	shiftCoord[nPP] = coord[nPP]+deltaPP*0.5;
	prop = CSX->GetPropertyByCoordPriority(shiftCoord,CSProperties::MATERIAL);
	if (prop)
	{
		CSPropMaterial* mat = prop->ToMaterial();
		inEC[2] += mat->GetMue(n)*fabs(delta);
		if (mat->GetSigmaWeighted(n,shiftCoord))
			inEC[3] += fabs(delta)/mat->GetSigmaWeighted(n,shiftCoord);
		else
			inEC[3] = 0;
	}
	else
	{
		inEC[2] += 1*fabs(delta);
		inEC[3] = 0;
	}

	inEC[2] = gridDelta * fabs(deltaP*deltaPP) * 2.0 * __MUE0__ / inEC[2];
	if (inEC[3]) inEC[3]=gridDelta*fabs(deltaP*deltaPP) * 2.0 / inEC[3];

	return true;
}

bool Operator::Calc_EffMatPos(int n, unsigned int* pos, double* inMat)
{
	int nP = (n+1)%3;
	int nPP = (n+2)%3;

	double delta=MainOp->GetIndexDelta(n,pos[n]);
	double deltaP=MainOp->GetIndexDelta(nP,pos[nP]);
	double deltaPP=MainOp->GetIndexDelta(nPP,pos[nPP]);

	double delta_M=MainOp->GetIndexDelta(n,pos[n]-1);
	double deltaP_M=MainOp->GetIndexDelta(nP,pos[nP]-1);
	double deltaPP_M=MainOp->GetIndexDelta(nPP,pos[nPP]-1);

	this->Calc_ECPos(n,pos,inMat);

	inMat[0] *= fabs(delta)/(0.25*(fabs(deltaP_M) + fabs(deltaP))*(fabs(deltaPP_M) + fabs(deltaPP)))/gridDelta;
	inMat[1] *= fabs(delta)/(0.25*(fabs(deltaP_M) + fabs(deltaP))*(fabs(deltaPP_M) + fabs(deltaPP)))/gridDelta;

	inMat[2] *= 0.5*(fabs(delta_M) + fabs(delta)) / fabs(deltaP*deltaPP) / gridDelta;
	inMat[3] *= 0.5*(fabs(delta_M) + fabs(delta)) / fabs(deltaP*deltaPP) / gridDelta;

	return true;
}

void Operator::Init_EC()
{
	for (int n=0;n<3;++n)
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
		for (unsigned int i=0;i<MainOp->GetSize();i++) //init all
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
	if (CSX==NULL) {cerr << "CartOperator::Calc_EC: CSX not given or invalid!!!" << endl; return false;}
	unsigned int ipos;
	unsigned int pos[3];
	double inEC[4];
	for (int n=0;n<3;++n)
	{
		for (pos[2]=0;pos[2]<numLines[2];++pos[2])
		{
			for (pos[1]=0;pos[1]<numLines[1];++pos[1])
			{
				for (pos[0]=0;pos[0]<numLines[0];++pos[0])
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
#if 1  //use the old timestep-calc (1) or the new one (0)
	return CalcTimestep_Var3();
#else
	return CalcTimestep_Var1();
#endif
}

////Berechnung nach Andreas Rennings Dissertation 2008, Seite 66, Formel 4.52
double Operator::CalcTimestep_Var1()
{
	cout << "Operator::CalcTimestep(): Using timestep algorithm by Andreas Rennings, Dissertation @ University Duisburg-Essen, 2008, pp. 66, eq. 4.52" << endl;
	dT=1e200;
	double newT;
	unsigned int pos[3];
	unsigned int ipos;
	unsigned int ipos_PM;
	unsigned int ipos_PPM;
	MainOp->SetReflection2Cell();
	for (int n=0;n<3;++n)
	{
		int nP = (n+1)%3;
		int nPP = (n+2)%3;

		for (pos[2]=0;pos[2]<numLines[2];++pos[2])
		{
			for (pos[1]=0;pos[1]<numLines[1];++pos[1])
			{
				for (pos[0]=0;pos[0]<numLines[0];++pos[0])
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
	for (unsigned int n=1;n<count;++n)
		if (val[n]<min)
			min = val[n];
	return min;
}

//Berechnung nach Andreas Rennings Dissertation 2008, Seite 76 ff, Formel 4.77 ff
double Operator::CalcTimestep_Var3()
{
	dT=1e200;
	cout << "Operator::CalcTimestep(): Using timestep algorithm by Andreas Rennings, Dissertation @ University Duisburg-Essen, 2008, pp. 76, eq. 4.77 ff." << endl;
	double newT;
	unsigned int pos[3];
	unsigned int ipos;
	double w_total=0;
	double wqp=0,wt1=0,wt2=0;
	double wt_4[4]={0,0,0,0};
	MainOp->SetReflection2Cell();
	for (int n=0;n<3;++n)
	{
		int nP = (n+1)%3;
		int nPP = (n+2)%3;

		for (pos[2]=0;pos[2]<numLines[2];++pos[2])
		{
			for (pos[1]=0;pos[1]<numLines[1];++pos[1])
			{
				for (pos[0]=0;pos[0]<numLines[0];++pos[0])
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

	for (pos[2]=0;pos[2]<numLines[2];++pos[2])
	{
		delta[2]=fabs(MainOp->GetIndexDelta(2,pos[2]));
		for (pos[1]=0;pos[1]<numLines[1];++pos[1])
		{
			delta[1]=fabs(MainOp->GetIndexDelta(1,pos[1]));
			for (pos[0]=0;pos[0]<numLines[0];++pos[0])
			{
				delta[0]=fabs(MainOp->GetIndexDelta(0,pos[0]));

				//electric field excite
				for (int n=0;n<3;++n)
				{
					volt_coord[0] = discLines[0][pos[0]];
					volt_coord[1] = discLines[1][pos[1]];
					volt_coord[2] = discLines[2][pos[2]];
					volt_coord[n]+=delta[n]*0.5;
					for (size_t p=0;p<vec_prop.size();++p)
					{
						prop = vec_prop.at(p);
						elec = prop->ToElectrode();
						if (prop->CheckCoordInPrimitive(volt_coord,priority)==false)
							elec=NULL;
						if (elec!=NULL)
						{
							if ((elec->GetActiveDir(n)) && ( (elec->GetExcitType()==0) || (elec->GetExcitType()==1) ))//&& (pos[n]<numLines[n]-1))
							{
								amp = elec->GetWeightedExcitation(n,volt_coord)*GetMeshDelta(n,pos);// delta[n]*gridDelta;
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
									GetVV(n,pos[0],pos[1],pos[2]) = 0;
									GetVI(n,pos[0],pos[1],pos[2]) = 0;
								}
							}
						}
					}
				}

				//magnetic field excite
				for (int n=0;n<3;++n)
				{
					int nP = (n+1)%3;
					int nPP = (n+2)%3;
					curr_coord[0] = discLines[0][pos[0]];
					curr_coord[1] = discLines[1][pos[1]];
					curr_coord[2] = discLines[2][pos[2]];
					curr_coord[nP] +=delta[nP]*0.5;
					curr_coord[nPP] +=delta[nPP]*0.5;
					for (size_t p=0;p<vec_prop.size();++p)
					{
						prop = vec_prop.at(p);
						elec = prop->ToElectrode();
						if (prop->CheckCoordInPrimitive(curr_coord,priority)==false)
							elec=NULL;
						if (elec!=NULL)
						{
							if ((elec->GetActiveDir(n)) && ( (elec->GetExcitType()==2) || (elec->GetExcitType()==3) ))//&& (pos[n]<numLines[n]-1))
							{
								amp = elec->GetWeightedExcitation(n,curr_coord)*GetMeshDelta(n,pos,true);// delta[n]*gridDelta;
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
									GetII(n,pos[0],pos[1],pos[2]) = 0;
									GetIV(n,pos[0],pos[1],pos[2]) = 0;
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
	for (size_t p=0;p<vec_prop.size();++p)
	{
		prop = vec_prop.at(p);
		elec = prop->ToElectrode();
		for (size_t n=0;n<prop->GetQtyPrimitives();++n)
		{
			CSPrimitives* prim = prop->GetPrimitive(n);
			CSPrimCurve* curv = prim->ToCurve();
			if (curv)
			{
				for (size_t i=1;i<curv->GetNumberOfPoints();++i)
				{
					curv->GetPoint(i-1,p1);
					curv->GetPoint(i,p2);
					path = FindPath(p1,p2);
					if (path.dir.size()>0)
						prim->SetPrimitiveUsed(true);
					for (size_t t=0;t<path.dir.size();++t)
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
									GetVV(n,pos[0],pos[1],pos[2]) = 0;
									GetVI(n,pos[0],pos[1],pos[2]) = 0;
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
	m_Nr_PEC[0]=0;	m_Nr_PEC[1]=0;	m_Nr_PEC[2]=0;

	CalcPEC_Range(0,numLines[0]-1,m_Nr_PEC);

	CalcPEC_Curves();

	return true;
}

void Operator::CalcPEC_Range(unsigned int startX, unsigned int stopX, unsigned int* counter)
{
	double coord[3];
	double delta;
	unsigned int pos[3];
	for (pos[0]=startX;pos[0]<=stopX;++pos[0])
	{
		for (pos[1]=0;pos[1]<numLines[1];++pos[1])
		{
			for (pos[2]=0;pos[2]<numLines[2];++pos[2])
			{
				for (int n=0;n<3;++n)
				{
					coord[0] = discLines[0][pos[0]];
					coord[1] = discLines[1][pos[1]];
					coord[2] = discLines[2][pos[2]];
					delta=MainOp->GetIndexDelta(n,pos[n]);
					coord[n]= discLines[n][pos[n]] + delta*0.5;
					CSProperties* prop = CSX->GetPropertyByCoordPriority(coord, (CSProperties::PropertyType)(CSProperties::MATERIAL | CSProperties::METAL));
					if (prop)
					{
						if (prop->GetType()==CSProperties::METAL) //set to PEC
						{
							GetVV(n,pos[0],pos[1],pos[2]) = 0;
							GetVI(n,pos[0],pos[1],pos[2]) = 0;
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
	for (size_t p=0;p<vec_prop.size();++p)
	{
		CSProperties* prop = vec_prop.at(p);
		for (size_t n=0;n<prop->GetQtyPrimitives();++n)
		{
			CSPrimitives* prim = prop->GetPrimitive(n);
			CSPrimCurve* curv = prim->ToCurve();
			if (curv)
			{
				for (size_t i=1;i<curv->GetNumberOfPoints();++i)
				{
					curv->GetPoint(i-1,p1);
					curv->GetPoint(i,p2);
					path = FindPath(p1,p2);
					if (path.dir.size()>0)
						prim->SetPrimitiveUsed(true);
					for (size_t t=0;t<path.dir.size();++t)
					{
//						cerr << path.dir.at(t) << " " << path.posPath[0].at(t) << " " << path.posPath[1].at(t) << " " << path.posPath[2].at(t) << endl;
						GetVV(path.dir.at(t),path.posPath[0].at(t),path.posPath[1].at(t),path.posPath[2].at(t)) = 0;
						GetVI(path.dir.at(t),path.posPath[0].at(t),path.posPath[1].at(t),path.posPath[2].at(t)) = 0;
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
