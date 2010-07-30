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

#include "operator_ext_pml_sf.h"
#include "engine_ext_pml_sf.h"
#include "operator_cylinder.h"

#include "tools/array_ops.h"

bool Build_Split_Field_PML(Operator* op, int BC[6], int size[6])
{
	for (int n=0;n<6;++n)
	{
		if (BC[n]==3) //split field PML
		{
			cerr << "Build_Split_Field_PML:: Warning, currently only pml planes are implemented... edges and corner coming soon..." << endl;
			Operator_Ext_PML_SF_Plane* op_pml_sf = new Operator_Ext_PML_SF_Plane(op);
			op_pml_sf->SetDirection(n/2,n%2);
			op_pml_sf->SetPMLLength(size[n]);
			op_pml_sf->SetBoundaryCondition(BC);
			op->AddExtension(op_pml_sf);
		}
	}
	return true;
}

/************************************************  Operator_Ext_PML_SF  **************************************************************************/
Operator_Ext_PML_SF::Operator_Ext_PML_SF(Operator* op) : Operator_Extension(op)
{
	m_SetupDone = false;

	m_numLines[0]=0;
	m_numLines[1]=0;
	m_numLines[2]=0;

	vv[0] = NULL;	vv[1] = NULL;
	vi[0] = NULL;	vi[1] = NULL;
	ii[0] = NULL;	ii[1] = NULL;
	iv[0] = NULL;	iv[1] = NULL;

	for (int n=0;n<6;++n)
		m_BC[n]=0;
}

Operator_Ext_PML_SF::~Operator_Ext_PML_SF()
{
	DeleteOP();
}

void Operator_Ext_PML_SF::InitOP()
{
	if (!m_SetupDone)
		return;
	vv[0] = Create_N_3DArray<FDTD_FLOAT>(m_numLines);
	vv[1] = Create_N_3DArray<FDTD_FLOAT>(m_numLines);

	vi[0] = Create_N_3DArray<FDTD_FLOAT>(m_numLines);
	vi[1] = Create_N_3DArray<FDTD_FLOAT>(m_numLines);

	ii[0] = Create_N_3DArray<FDTD_FLOAT>(m_numLines);
	ii[1] = Create_N_3DArray<FDTD_FLOAT>(m_numLines);

	iv[0] = Create_N_3DArray<FDTD_FLOAT>(m_numLines);
	iv[1] = Create_N_3DArray<FDTD_FLOAT>(m_numLines);
}


void Operator_Ext_PML_SF::DeleteOP()
{
	if (!m_SetupDone)
		return;
	Delete_N_3DArray<FDTD_FLOAT>(vv[0],m_numLines);
	vv[0] = NULL;
	Delete_N_3DArray<FDTD_FLOAT>(vv[1],m_numLines);
	vv[1] = NULL;

	Delete_N_3DArray<FDTD_FLOAT>(vi[0],m_numLines);
	vi[0] = NULL;
	Delete_N_3DArray<FDTD_FLOAT>(vi[1],m_numLines);
	vi[1] = NULL;

	Delete_N_3DArray<FDTD_FLOAT>(ii[0],m_numLines);
	ii[0] = NULL;
	Delete_N_3DArray<FDTD_FLOAT>(ii[1],m_numLines);
	ii[1] = NULL;

	Delete_N_3DArray<FDTD_FLOAT>(iv[0],m_numLines);
	iv[0] = NULL;
	Delete_N_3DArray<FDTD_FLOAT>(iv[1],m_numLines);
	iv[1] = NULL;
}

bool Operator_Ext_PML_SF::BuildExtension()
{
	if (!m_SetupDone)
	{
		cerr << "Operator_Ext_PML_SF::BuildExtension: Warning, Extension not initialized! Abort build!!" << endl;
		return false;
	}
	double dT = m_Op->GetTimestep();
	unsigned int pos[] = {0,0,0};
	DeleteOP();
	InitOP();

	double inEC[4];

	for (int n=0;n<3;++n)
	{
		for (pos[0]=0;pos[0]<m_numLines[0];++pos[0])
		{
			for (pos[1]=0;pos[1]<m_numLines[1];++pos[1])
			{
				for (pos[2]=0;pos[2]<m_numLines[2];++pos[2])
				{
					Calc_ECPos(0,n,pos,inEC);
					if (inEC[0]>0)
						GetVV(0,n,pos[0],pos[1],pos[2]) = (1-dT*inEC[1]/2/inEC[0])/(1+dT*inEC[1]/2/inEC[0]);
					if (inEC[2]>0)
						GetII(0,n,pos[0],pos[1],pos[2]) = (1-dT*inEC[3]/2/inEC[2])/(1+dT*inEC[3]/2/inEC[2]);

					if (inEC[0]>0)
						GetVI(0,n,pos[0],pos[1],pos[2]) = (dT/inEC[0])/(1+dT*inEC[1]/2/inEC[0]);
					if (inEC[2]>0)
						GetIV(0,n,pos[0],pos[1],pos[2]) = (dT/inEC[2])/(1+dT*inEC[3]/2/inEC[2]);

//					if (n==0)
//						cerr << pos[0] << " " << pos[1] << " " << pos[2] << " " << inEC[1] << endl;

					Calc_ECPos(1,n,pos,inEC);
					if (inEC[0]>0)
						GetVV(1,n,pos[0],pos[1],pos[2]) = (1-dT*inEC[1]/2/inEC[0])/(1+dT*inEC[1]/2/inEC[0]);
					if (inEC[2]>0)
						GetII(1,n,pos[0],pos[1],pos[2]) = (1-dT*inEC[3]/2/inEC[2])/(1+dT*inEC[3]/2/inEC[2]);

					if (inEC[0]>0)
						GetVI(1,n,pos[0],pos[1],pos[2]) = (dT/inEC[0])/(1+dT*inEC[1]/2/inEC[0]);
					if (inEC[2]>0)
						GetIV(1,n,pos[0],pos[1],pos[2]) = (dT/inEC[2])/(1+dT*inEC[3]/2/inEC[2]);

//					if (n==0)
//						cerr << pos[0] << " " << pos[1] << " " << pos[2] << " " << inEC[1] << endl;
				}
			}
		}
	}
	ApplyBC();
	return true;
}

/************************************************  Operator_Ext_PML_SF_Plane  **************************************************************************/
Operator_Ext_PML_SF_Plane::Operator_Ext_PML_SF_Plane(Operator* op) : Operator_Ext_PML_SF(op)
{
}

Operator_Ext_PML_SF_Plane::~Operator_Ext_PML_SF_Plane()
{
}

void Operator_Ext_PML_SF_Plane::SetDirection(int ny, bool top_ny)
{
	if ((ny<0) || (ny>2))
		return;

	m_ny = ny;
	m_nyP = (ny+1)%3;
	m_nyPP = (ny+2)%3;

	m_top = top_ny;

	m_numLines[m_ny] = 8; //default width of the pml plane
	m_numLines[m_nyP] = m_Op->GetNumberOfLines(m_nyP);
	m_numLines[m_nyPP] = m_Op->GetNumberOfLines(m_nyPP);

	unsigned int pos[] = {0,0,0};
	m_LineNr = (unsigned int)((int)m_top * (int)(m_Op->GetNumberOfLines(m_ny)-1));
	pos[m_ny] = m_LineNr;

	m_pml_delta = m_Op->GetMeshDelta(m_ny,pos);
}

void Operator_Ext_PML_SF_Plane::SetPMLLength(int width)
{
	if (m_ny<0)
	{
		cerr << "Operator_Ext_PML_SF_Plane::SetPMLLength: Warning, Direction not set! Use SetDirection first!!" << endl;
		return;
	}

	if (width<4)
	{
		cerr << "Operator_Ext_PML_SF_Plane::SetPMLLength: Warning: A pml width smaller than 4 lines is not allowed, skipping..." << endl;
		return;
	}
	if (width>50)
	{
		cerr << "Operator_Ext_PML_SF_Plane::SetPMLLength: Warning: A pml width greater than 20 lines is not allowed, skipping..." << endl;
		return;
	}
	m_SetupDone = true;
	m_numLines[m_ny] = width;

	m_pml_width = (width - 1.5) * m_pml_delta;

}

double Operator_Ext_PML_SF_Plane::GetNodeArea(int ny, unsigned int pos[3], bool dualMesh) const
{
	unsigned int l_pos[] = {pos[0],pos[1],pos[2]};
	l_pos[m_ny] = m_LineNr;

	return m_Op->GetNodeArea(ny,l_pos,dualMesh);
}

double Operator_Ext_PML_SF_Plane::GetNodeLength(int ny, unsigned int pos[3], bool dualMesh) const
{
	if (ny==m_ny)
		return m_pml_delta;

	unsigned int l_pos[] = {pos[0],pos[1],pos[2]};
	l_pos[m_ny] = m_LineNr;
	return m_Op->GetMeshDelta(ny,l_pos,dualMesh);
}

double Operator_Ext_PML_SF_Plane::GetKappaGraded(double depth, double Zm) const
{
	if (depth<0)
		return 0.0;

	//todo: use fparser to allow arbitrary, user-defined profiles and parameter
	double g = 2.5;
	double R0 = 1e-6;
	double kappa0 = -log(R0)*log(g)/(2*m_pml_delta * pow(g,m_pml_width/m_pml_delta) -1);
	return pow(g,depth/m_pml_delta)*kappa0 / Zm;
}

bool Operator_Ext_PML_SF_Plane::Calc_ECPos(int nP, int n, unsigned int* pos, double* inEC) const
{
	unsigned int l_pos[] = {pos[0],pos[1],pos[2]};
	l_pos[m_ny] = m_LineNr;

	double inMat[4];
	m_Op->Calc_EffMatPos(n,l_pos,inMat);

	double Zm2 = inMat[2] / inMat[0]; // Zm^2 = mue/eps
	double Zm = sqrt(Zm2);			  // Zm   = sqrt(Zm^2) = sqrt(mue/eps)
	double kappa = 0;
	double sigma = 0;
	double depth = 0;
	if ( (n + nP + 1)%3 == m_ny )
	{
		if (m_top)
		{
			depth = pos[m_ny]*m_pml_delta - 0.5*m_pml_delta;
			kappa = GetKappaGraded(depth, Zm);
			sigma = GetKappaGraded(depth + 0.5*m_pml_delta, Zm) * Zm2;
		}
		else
		{
			depth = m_pml_width - (pos[m_ny])*m_pml_delta;
			kappa = GetKappaGraded(depth, Zm) ;
			sigma = GetKappaGraded(depth-0.5*m_pml_delta, Zm) * Zm2;
		}
		if ((inMat[0]<=0) || (inMat[2]<=0)) //check if material properties are valid (necessary for cylindrical coords)
		{
			kappa = sigma = 0;
		}
	}

	double geomFactor = GetNodeArea(n,pos) / GetNodeLength(n,pos);
	if (geomFactor<=0 || isnan(geomFactor) || isinf(geomFactor)) //check if geomFactor is positive, not zero and a valid number (necessary for cylindrical coords)
		geomFactor = 0;
	inEC[0] = inMat[0] * geomFactor;
	inEC[1] = (inMat[1]+kappa) * geomFactor;

	geomFactor = GetNodeArea(n,pos) / GetNodeLength(n,pos);
	if (geomFactor<=0 || isnan(geomFactor) || isinf(geomFactor)) //check if geomFactor is positive, not zero and a valid number (necessary for cylindrical coords)
		geomFactor = 0;
	inEC[2] = inMat[2] * geomFactor;
	inEC[3] = (inMat[3]+sigma) * geomFactor;

	return true;
}

void Operator_Ext_PML_SF_Plane::ApplyBC()
{
	bool PEC[6] = {1,1,1,1,1,1};
	bool PMC[6] = {0,0,0,0,0,0};

	if (m_top==false)
		PEC[2*m_ny+1] = 0;

	for (int n=0;n<6;++n)
	{
		PMC[n] = (m_BC[n] == 1);
		if (n/2 == m_ny)
			PMC[n] = false;
	}

	//apply BC
	unsigned int pos[3] = {0,0,0};
	for (int n=0;n<3;++n)
	{
		int nP = (n+1)%3;
		int nPP = (n+2)%3;
		for (pos[nP]=0;pos[nP]<m_numLines[nP];++pos[nP])
		{
			for (pos[nPP]=0;pos[nPP]<m_numLines[nPP];++pos[nPP])
			{
				for (int m=0;m<2;++m)
				{
					pos[n]=0;
					GetVV(m,nP,pos[0],pos[1],pos[2]) *= (FDTD_FLOAT)!PEC[2*n];
					GetVI(m,nP,pos[0],pos[1],pos[2]) *= (FDTD_FLOAT)!PEC[2*n];
					GetVV(m,nPP,pos[0],pos[1],pos[2]) *= (FDTD_FLOAT)!PEC[2*n];
					GetVI(m,nPP,pos[0],pos[1],pos[2]) *= (FDTD_FLOAT)!PEC[2*n];

					GetII(m,n,pos[0],pos[1],pos[2]) *= (FDTD_FLOAT)!PMC[2*n];
					GetIV(m,n,pos[0],pos[1],pos[2]) *= (FDTD_FLOAT)!PMC[2*n];
					GetII(m,nP,pos[0],pos[1],pos[2]) *= (FDTD_FLOAT)!PMC[2*n];
					GetIV(m,nP,pos[0],pos[1],pos[2]) *= (FDTD_FLOAT)!PMC[2*n];
					GetII(m,nPP,pos[0],pos[1],pos[2]) *= (FDTD_FLOAT)!PMC[2*n];
					GetIV(m,nPP,pos[0],pos[1],pos[2]) *= (FDTD_FLOAT)!PMC[2*n];

					pos[n]=m_numLines[n]-1;
					GetVV(m,n,pos[0],pos[1],pos[2]) = 0; // these are outside the FDTD-domain as defined by the main disc
					GetVI(m,n,pos[0],pos[1],pos[2]) = 0; // these are outside the FDTD-domain as defined by the main disc

					GetVV(m,nP,pos[0],pos[1],pos[2]) *= (FDTD_FLOAT)!PEC[2*n+1];
					GetVI(m,nP,pos[0],pos[1],pos[2]) *= (FDTD_FLOAT)!PEC[2*n+1];
					GetVV(m,nPP,pos[0],pos[1],pos[2]) *= (FDTD_FLOAT)!PEC[2*n+1];
					GetVI(m,nPP,pos[0],pos[1],pos[2]) *= (FDTD_FLOAT)!PEC[2*n+1];

					pos[n]=m_numLines[n]-2;
					GetII(m,nP,pos[0],pos[1],pos[2]) *= (FDTD_FLOAT)!PMC[2*n+1];
					GetIV(m,nP,pos[0],pos[1],pos[2]) *= (FDTD_FLOAT)!PMC[2*n+1];
					GetII(m,nPP,pos[0],pos[1],pos[2]) *= (FDTD_FLOAT)!PMC[2*n+1];
					GetIV(m,nPP,pos[0],pos[1],pos[2]) *= (FDTD_FLOAT)!PMC[2*n+1];
				}
			}
		}
	}
}

Engine_Extension* Operator_Ext_PML_SF_Plane::CreateEngineExtention()
{
	Engine_Ext_PML_SF_Plane* eng_ext = new Engine_Ext_PML_SF_Plane(this);
	return eng_ext;
}

bool Operator_Ext_PML_SF_Plane::IsCylinderCoordsSave()  const
{
	if (m_ny==2)
	{
		Operator_Cylinder* op_cyl = dynamic_cast<Operator_Cylinder*>(m_Op);
		if (op_cyl==NULL)
		{
			cerr << "Operator_Ext_PML_SF_Plane::IsCylinderCoordsSave(): Error!!! Sanity check failed!!! ==> Developer is not sane.... this should never have happend.. exit..." << endl;
			exit(0);
		}
		if (op_cyl->GetClosedAlpha())
		{
			cerr << "Operator_Ext_PML_SF_Plane::IsCylinderCoordsSave(): Warning... this extension can not handle a closed alpha cylinder operator... " << endl;
			return false;
		}
		return true;
	}
	return false;
}

void Operator_Ext_PML_SF_Plane::ShowStat(ostream &ostr)  const
{
	Operator_Extension::ShowStat(ostr);
	string XYZ[3] = {"x","y","z"};
	string top_bot[2] = {"bottom", "top"};
	ostr << " Active direction\t: " << XYZ[m_ny] << " (" << top_bot[m_top] << ")" << endl;
	ostr << " PML width (cells)\t: " << m_numLines[m_ny] << endl;
}
