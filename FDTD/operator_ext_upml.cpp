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

#include "operator_ext_upml.h"
#include "operator_cylinder.h"
#include "engine_ext_upml.h"
#include "tools/array_ops.h"
#include "fparser.hh"

Operator_Ext_UPML::Operator_Ext_UPML(Operator* op) : Operator_Extension(op)
{
	m_GradingFunction = new FunctionParser();
	//default grading function
	SetGradingFunction(" -log(1e-6)*log(2.5)/(2*dl*pow(2.5,W/dl)-1) * pow(2.5, D/dl) / Z ");

	for (int n=0;n<6;++n)
	{
		m_BC[n]=0;
		m_Size[n]=0;
	}
	for (int n=0;n<3;++n)
	{
		m_StartPos[n]=0;
		m_numLines[n]=0;
	}

	vv = NULL;
	vvfo = NULL;
	vvfn = NULL;
	ii = NULL;
	iifo = NULL;
	iifn = NULL;
}

Operator_Ext_UPML::~Operator_Ext_UPML()
{
	delete m_GradingFunction;
	m_GradingFunction = NULL;
	DeleteOp();
}

bool Operator_Ext_UPML::Create_UPML(Operator* op, int BC[6], unsigned int size[6], string gradFunc)
{
	Operator_Ext_UPML* op_ext_upml=NULL;
	unsigned int start[3]={0 ,0 ,0};
	unsigned int stop[3] ={op->GetNumberOfLines(0)-1,op->GetNumberOfLines(1)-1,op->GetNumberOfLines(2)-1};

	//create a pml in x-direction over the full width of yz-space
	if (BC[0]==3)
	{
		op_ext_upml = new Operator_Ext_UPML(op);
		op_ext_upml->SetGradingFunction(gradFunc);
		start[0]=0;
		stop[0] =size[0];
		op_ext_upml->SetBoundaryCondition(BC, size);
		op_ext_upml->SetRange(start,stop);
		op->AddExtension(op_ext_upml);
	}
	if (BC[1]==3)
	{
		op_ext_upml = new Operator_Ext_UPML(op);
		op_ext_upml->SetGradingFunction(gradFunc);
		start[0]=op->GetNumberOfLines(0)-1-size[1];
		stop[0] =op->GetNumberOfLines(0)-1;
		op_ext_upml->SetBoundaryCondition(BC, size);
		op_ext_upml->SetRange(start,stop);
		op->AddExtension(op_ext_upml);
	}

	//create a pml in y-direction over the xz-space (if a pml in x-direction already exists, skip that corner regions)
	start[0]=(size[0]+1)*(BC[0]==3);
	stop[0] =op->GetNumberOfLines(0)-1-(size[0]+1)*(BC[1]==3);

	if (BC[2]==3)
	{
		op_ext_upml = new Operator_Ext_UPML(op);
		op_ext_upml->SetGradingFunction(gradFunc);
		start[1]=0;
		stop[1] =size[2];
		op_ext_upml->SetBoundaryCondition(BC, size);
		op_ext_upml->SetRange(start,stop);
		op->AddExtension(op_ext_upml);
	}
	if (BC[3]==3)
	{
		op_ext_upml = new Operator_Ext_UPML(op);
		op_ext_upml->SetGradingFunction(gradFunc);
		start[1]=op->GetNumberOfLines(1)-1-size[3];
		stop[1] =op->GetNumberOfLines(1)-1;
		op_ext_upml->SetBoundaryCondition(BC, size);
		op_ext_upml->SetRange(start,stop);
		op->AddExtension(op_ext_upml);
	}

	//create a pml in z-direction over the xy-space (if a pml in x- and/or y-direction already exists, skip that corner/edge regions)
	start[1]=(size[2]+1)*(BC[2]==3);
	stop[1] =op->GetNumberOfLines(1)-1-(size[3]+1)*(BC[3]==3);

	if (BC[4]==3)
	{
		op_ext_upml = new Operator_Ext_UPML(op);
		op_ext_upml->SetGradingFunction(gradFunc);
		start[2]=0;
		stop[2] =size[4];
		op_ext_upml->SetBoundaryCondition(BC, size);
		op_ext_upml->SetRange(start,stop);
		op->AddExtension(op_ext_upml);
	}
	if (BC[5]==3)
	{
		op_ext_upml = new Operator_Ext_UPML(op);
		op_ext_upml->SetGradingFunction(gradFunc);
		start[2]=op->GetNumberOfLines(2)-1-size[5];
		stop[2] =op->GetNumberOfLines(2)-1;
		op_ext_upml->SetBoundaryCondition(BC, size);
		op_ext_upml->SetRange(start,stop);
		op->AddExtension(op_ext_upml);
	}

	return true;
}


void Operator_Ext_UPML::DeleteOp()
{
	Delete_N_3DArray<FDTD_FLOAT>(vv,m_numLines);
	vv = NULL;
	Delete_N_3DArray<FDTD_FLOAT>(vvfo,m_numLines);
	vvfo = NULL;
	Delete_N_3DArray<FDTD_FLOAT>(vvfn,m_numLines);
	vvfn = NULL;
	Delete_N_3DArray<FDTD_FLOAT>(ii,m_numLines);
	ii = NULL;
	Delete_N_3DArray<FDTD_FLOAT>(iifo,m_numLines);
	iifo = NULL;
	Delete_N_3DArray<FDTD_FLOAT>(iifn,m_numLines);
	iifn = NULL;
}


bool Operator_Ext_UPML::SetGradingFunction(string func)
{
	if (func.empty())
		return true;

	m_GradFunc = func;
	int res = m_GradingFunction->Parse(m_GradFunc.c_str(), "D,dl,W,Z,N");
	if(res < 0) return true;

	cerr << "Operator_Ext_UPML::SetGradingFunction: Warning, an error occured parsing the pml grading function (see below) ..." << endl;
	cerr << func << "\n" << string(res, ' ') << "^\n" << m_GradingFunction->ErrorMsg() << "\n";
	return false;
}

void Operator_Ext_UPML::CalcGradingKappa(int ny, unsigned int pos[3], double Zm[3], double kappa_v[3], double kappa_i[3])
{
	double depth=0;
	double width=0;
	for (int n=0;n<3;++n)
	{
		if ((pos[n] <= m_Size[2*n]) && (m_BC[2*n]==3))  //lower pml in n-dir
		{
			width = (m_Op->GetDiscLine(n,m_Size[2*n]) - m_Op->GetDiscLine(n,0))*m_Op->GetGridDelta();
			depth = width - (m_Op->GetDiscLine(n,pos[n]) - m_Op->GetDiscLine(n,0))*m_Op->GetGridDelta();

			if (n==ny)
				depth-=m_Op->GetMeshDelta(n,pos)/2;
			double vars[5] = {depth, width/m_Size[2*n], width, Zm[n], m_Size[2*n]};
			if (depth>0)
				kappa_v[n] = m_GradingFunction->Eval(vars);
			else
				kappa_v[n]=0;
			if (n==ny)
				depth+=m_Op->GetMeshDelta(n,pos)/2;

			if (n!=ny)
				depth-=m_Op->GetMeshDelta(n,pos)/2;
			if (depth<0)
				depth=0;
			vars[0]=depth;
			if (depth>0)
				kappa_i[n] = m_GradingFunction->Eval(vars);
			else
				kappa_i[n] = 0;
		}
		else if ((pos[n] >= m_Op->GetNumberOfLines(n) -1 -m_Size[2*n+1]) && (m_BC[2*n+1]==3))  //upper pml in n-dir
		{
			width = (m_Op->GetDiscLine(n,m_Op->GetNumberOfLines(n)-1) - m_Op->GetDiscLine(n,m_Op->GetNumberOfLines(n)-m_Size[2*n+1]-1))*m_Op->GetGridDelta();
			depth = width - (m_Op->GetDiscLine(n,m_Op->GetNumberOfLines(n)-1) - m_Op->GetDiscLine(n,pos[n]))*m_Op->GetGridDelta();

			if (n==ny)
				depth+=m_Op->GetMeshDelta(n,pos)/2;
			double vars[5] = {depth, width/(m_Size[2*n]), width, Zm[n], m_Size[2*n]};
			if (depth>0)
				kappa_v[n] = m_GradingFunction->Eval(vars);
			else
				kappa_v[n]=0;
			if (n==ny)
				depth-=m_Op->GetMeshDelta(n,pos)/2;

			if (n!=ny)
				depth+=m_Op->GetMeshDelta(n,pos)/2;
			if (depth>width)
				depth=0;
			vars[0]=depth;
			if (depth>0)
				kappa_i[n] = m_GradingFunction->Eval(vars);
			else
				kappa_i[n]=0;
		}
		else
		{
			kappa_v[n] = 0;
			kappa_i[n] = 0;
		}
	}
}

bool Operator_Ext_UPML::BuildExtension()
{
	/*Calculate the upml coefficients as defined in:
	  Allen Taflove, computational electrodynamics - the FDTD method, third edition, chapter 7.8, pages 297-300
	  - modified by Thorsten Liebig to match the equivalent circuit (EC) FDTD method
	  - kappa is used for conductivities (instead of sigma)
	*/
	if (m_Op==NULL)
		return false;

	DeleteOp();
	vv = Create_N_3DArray<FDTD_FLOAT>(m_numLines);
	vvfo = Create_N_3DArray<FDTD_FLOAT>(m_numLines);
	vvfn = Create_N_3DArray<FDTD_FLOAT>(m_numLines);
	ii = Create_N_3DArray<FDTD_FLOAT>(m_numLines);
	iifo = Create_N_3DArray<FDTD_FLOAT>(m_numLines);
	iifn = Create_N_3DArray<FDTD_FLOAT>(m_numLines);

	unsigned int pos[3];
	unsigned int loc_pos[3];
	int nP,nPP;
	double kappa_v[3]={0,0,0};
	double kappa_i[3]={0,0,0};
	double eff_Mat[3][4];
	double Zm[3];
	double dT = m_Op->GetTimestep();

	for (loc_pos[0]=0;loc_pos[0]<m_numLines[0];++loc_pos[0])
	{
		pos[0] = loc_pos[0] + m_StartPos[0];
		for (loc_pos[1]=0;loc_pos[1]<m_numLines[1];++loc_pos[1])
		{
			pos[1] = loc_pos[1] + m_StartPos[1];
			for (loc_pos[2]=0;loc_pos[2]<m_numLines[2];++loc_pos[2])
			{
				pos[2] = loc_pos[2] + m_StartPos[2];
				//calc kappa here for all dir ...
				m_Op->Calc_EffMatPos(0,pos,eff_Mat[0]);
				m_Op->Calc_EffMatPos(1,pos,eff_Mat[1]);
				m_Op->Calc_EffMatPos(2,pos,eff_Mat[2]);
				Zm[0] = sqrt(eff_Mat[0][2]/eff_Mat[0][0]);
				Zm[1] = sqrt(eff_Mat[1][2]/eff_Mat[1][0]);
				Zm[2] = sqrt(eff_Mat[2][2]/eff_Mat[2][0]);
				for (int n=0;n<3;++n)
				{
					CalcGradingKappa(n, pos,Zm,kappa_v,kappa_i);
					nP = (n+1)%3;
					nPP = (n+2)%3;
					//check if pos is on PEC
					if ( (m_Op->GetVV(n,pos[0],pos[1],pos[2]) + m_Op->GetVI(n,pos[0],pos[1],pos[2])) != 0 )
					{
						//modify the original operator to perform eq. (7.85) by the main engine (EC-FDTD: equation is multiplied by delta_n)
						//the engine extension will replace the original voltages with the "voltage flux" (volt*eps0) prior to the voltage updates
						//after the updates are done the extension will calculate the new voltages (see below) and place them back into the main field domain
						m_Op->GetVV(n,pos[0],pos[1],pos[2]) = (2*__EPS0__ - kappa_v[nP]*dT) / (2*__EPS0__ + kappa_v[nP]*dT);
						m_Op->GetVI(n,pos[0],pos[1],pos[2]) = (2*__EPS0__*dT) / (2*__EPS0__ + kappa_v[nP]*dT) * m_Op->GetEdgeLength(n,pos) / m_Op->GetEdgeArea(n,pos);

						//operators needed by eq. (7.88) to calculate new voltages from old voltages and old and new "voltage fluxes"
						vv  [n][loc_pos[0]][loc_pos[1]][loc_pos[2]] = (2*__EPS0__ - kappa_v[nPP]*dT) / (2*__EPS0__ + kappa_v[nPP]*dT);
						vvfn[n][loc_pos[0]][loc_pos[1]][loc_pos[2]] = (2*__EPS0__ + kappa_v[n]*dT)   / (2*__EPS0__ + kappa_v[nPP]*dT)/eff_Mat[n][0];
						vvfo[n][loc_pos[0]][loc_pos[1]][loc_pos[2]] = (2*__EPS0__ - kappa_v[n]*dT)   / (2*__EPS0__ + kappa_v[nPP]*dT)/eff_Mat[n][0];
					}

					//check if pos is on PMC
					if ( (m_Op->GetII(n,pos[0],pos[1],pos[2]) + m_Op->GetIV(n,pos[0],pos[1],pos[2])) != 0 )
					{
						//modify the original operator to perform eq. (7.89) by the main engine (EC-FDTD: equation is multiplied by delta_n)
						//the engine extension will replace the original currents with the "current flux" (curr*mu0) prior to the current updates
						//after the updates are done the extension will calculate the new currents (see below) and place them back into the main field domain
						m_Op->GetII(n,pos[0],pos[1],pos[2]) = (2*__EPS0__ - kappa_i[nP]*dT) / (2*__EPS0__ + kappa_i[nP]*dT);
						m_Op->GetIV(n,pos[0],pos[1],pos[2]) = (2*__EPS0__*dT) / (2*__EPS0__ + kappa_i[nP]*dT) * m_Op->GetEdgeLength(n,pos,true) / m_Op->GetEdgeArea(n,pos,true);

						//operators needed by eq. (7.90) to calculate new currents from old currents and old and new "current fluxes"
						ii  [n][loc_pos[0]][loc_pos[1]][loc_pos[2]] = (2*__EPS0__ - kappa_i[nPP]*dT) / (2*__EPS0__ + kappa_i[nPP]*dT);
						iifn[n][loc_pos[0]][loc_pos[1]][loc_pos[2]] = (2*__EPS0__ + kappa_i[n]*dT)   / (2*__EPS0__ + kappa_i[nPP]*dT)/eff_Mat[n][2];
						iifo[n][loc_pos[0]][loc_pos[1]][loc_pos[2]] = (2*__EPS0__ - kappa_i[n]*dT)   / (2*__EPS0__ + kappa_i[nPP]*dT)/eff_Mat[n][2];
					}
				}
			}
		}
	}
	return true;
}

Engine_Extension* Operator_Ext_UPML::CreateEngineExtention()
{
	Engine_Ext_UPML* eng_ext = new Engine_Ext_UPML(this);
	return eng_ext;
}

void Operator_Ext_UPML::ShowStat(ostream &ostr)  const
{
	Operator_Extension::ShowStat(ostr);

	ostr << " PML range\t\t: " << "[" << m_StartPos[0]<< "," << m_StartPos[1]<< "," << m_StartPos[2]<< "] to ["
			<<  m_StartPos[0]+m_numLines[0]-1 << "," << m_StartPos[1]+m_numLines[1]-1 << "," << m_StartPos[2]+m_numLines[2]-1 << "]" << endl;
	ostr << " Grading function\t: \"" << m_GradFunc << "\"" << endl;
}
