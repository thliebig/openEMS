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
#include "FDTD/operator_cylindermultigrid.h"
#include "engine_ext_upml.h"
#include "tools/array_ops.h"
#include "fparser.hh"

using namespace std;

Operator_Ext_UPML::Operator_Ext_UPML(Operator* op) : Operator_Extension(op)
{
	m_GradingFunction = new FunctionParser();
	//default grading function
	SetGradingFunction(" -log(1e-6)*log(2.5)/(2*dl*Z*(pow(2.5,W/dl)-1)) * pow(2.5, D/dl) ");

	for (int n=0; n<6; ++n)
	{
		m_BC[n]=0;
		m_Size[n]=0;
	}
	for (int n=0; n<3; ++n)
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

void Operator_Ext_UPML::SetBoundaryCondition(const int* BCs, const unsigned int size[6])
{
	for (int n=0; n<6; ++n)
	{
		m_BC[n]=BCs[n];
		m_Size[n]=size[n];
	}
}

void Operator_Ext_UPML::SetRange(const unsigned int start[3], const unsigned int stop[3])
{
	for (int n=0; n<3; ++n)
	{
		m_StartPos[n]=start[n];
		m_numLines[n]=stop[n]-start[n]+1;
	}
}

bool Operator_Ext_UPML::Create_UPML(Operator* op, const int ui_BC[6], const unsigned int ui_size[6], string gradFunc)
{
	int BC[6]={ui_BC[0],ui_BC[1],ui_BC[2],ui_BC[3],ui_BC[4],ui_BC[5]};
	unsigned int size[6]={ui_size[0],ui_size[1],ui_size[2],ui_size[3],ui_size[4],ui_size[5]};

	//check if mesh is large enough to support the pml
	for (int n=0; n<3; ++n)
		if ( (size[2*n]*(BC[2*n]==3)+size[2*n+1]*(BC[2*n+1]==3)) >= op->GetNumberOfLines(n,true) )
		{
			cerr << "Operator_Ext_UPML::Create_UPML: Warning: Not enough lines in direction: " << n << ", resetting to PEC" << endl;
			BC[2*n]=0;
			size[2*n]=0;
			BC[2*n+1]=0;
			size[2*n+1]=0;
		}

	//check cylindrical coord compatiblility
	Operator_Cylinder* op_cyl = dynamic_cast<Operator_Cylinder*>(op);
	if (op_cyl)
	{
		if ((BC[0]==3) && (op_cyl->GetClosedAlpha() || op_cyl->GetR0Included()))
		{
			BC[0]=0;
			size[0]=0;
			cerr << "Operator_Ext_UPML::Create_UPML: Warning: An upml in r-min direction is not possible, resetting to PEC..." << endl;
		}
		if ( (BC[2]==3) && (op_cyl->GetClosedAlpha()) )
		{
			BC[2]=0;
			size[2]=0;
			cerr << "Operator_Ext_UPML::Create_UPML: Warning: An upml in alpha-min direction is not possible, resetting to PEC..." << endl;
		}
		if ( (BC[3]==3) && (op_cyl->GetClosedAlpha()) )
		{
			BC[3]=0;
			size[3]=0;
			cerr << "Operator_Ext_UPML::Create_UPML: Warning: An upml in alpha-max direction is not possible, resetting to PEC..." << endl;
		}
	}

	//check cylindrical coord compatiblility
	if (dynamic_cast<Operator_CylinderMultiGrid*>(op))
	{
		if (BC[2]==3)
		{
			BC[2]=0;
			size[2]=0;
			cerr << "Operator_Ext_UPML::Create_UPML: Warning: An upml in alpha direction is not possible for a cylindrical multi-grid, resetting to PEC..." << endl;
		}
		if (BC[3]==3)
		{
			BC[3]=0;
			size[3]=0;
			cerr << "Operator_Ext_UPML::Create_UPML: Warning: An upml in alpha direction is not possible for a cylindrical multi-grid, resetting to PEC..." << endl;
		}
	}


	Operator_Ext_UPML* op_ext_upml=NULL;
	unsigned int start[3]={0 ,0 ,0};
	unsigned int stop[3] ={op->GetNumberOfLines(0,true)-1,op->GetNumberOfLines(1,true)-1,op->GetNumberOfLines(2,true)-1};

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
		start[0]=op->GetNumberOfLines(0,true)-1-size[1];
		stop[0] =op->GetNumberOfLines(0,true)-1;
		op_ext_upml->SetBoundaryCondition(BC, size);
		op_ext_upml->SetRange(start,stop);
		op->AddExtension(op_ext_upml);
	}

	//create a pml in y-direction over the xz-space (if a pml in x-direction already exists, skip that corner regions)
	start[0]=(size[0]+1)*(BC[0]==3);
	stop[0] =op->GetNumberOfLines(0,true)-1-(size[0]+1)*(BC[1]==3);

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
		start[1]=op->GetNumberOfLines(1,true)-1-size[3];
		stop[1] =op->GetNumberOfLines(1,true)-1;
		op_ext_upml->SetBoundaryCondition(BC, size);
		op_ext_upml->SetRange(start,stop);
		op->AddExtension(op_ext_upml);
	}

	//create a pml in z-direction over the xy-space (if a pml in x- and/or y-direction already exists, skip that corner/edge regions)
	start[1]=(size[2]+1)*(BC[2]==3);
	stop[1] =op->GetNumberOfLines(1,true)-1-(size[3]+1)*(BC[3]==3);

	//exclude x-lines that does not belong to the base multi-grid operator;
	Operator_CylinderMultiGrid* op_cyl_MG = dynamic_cast<Operator_CylinderMultiGrid*>(op);
	if (op_cyl_MG)
		start[0] = op_cyl_MG->GetSplitPos()-1;

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
		start[2]=op->GetNumberOfLines(2,true)-1-size[5];
		stop[2] =op->GetNumberOfLines(2,true)-1;
		op_ext_upml->SetBoundaryCondition(BC, size);
		op_ext_upml->SetRange(start,stop);
		op->AddExtension(op_ext_upml);
	}

	BC[1]=0;
	size[1]=0;
	//create pml extensions (in z-direction only) for child operators in cylindrical multigrid operators
	while (op_cyl_MG)
	{
		Operator_Cylinder* op_child = op_cyl_MG->GetInnerOperator();
		op_cyl_MG = dynamic_cast<Operator_CylinderMultiGrid*>(op_child);
		for (int n=0; n<2; ++n)
		{
			start[n]=0;
			stop[n]=op_child->GetNumberOfLines(n,true)-1;
		}

		if (op_cyl_MG)
			start[0] = op_cyl_MG->GetSplitPos()-1;

		if (BC[4]==3)
		{
			op_ext_upml = new Operator_Ext_UPML(op_child);
			op_ext_upml->SetGradingFunction(gradFunc);
			start[2]=0;
			stop[2] =size[4];
			op_ext_upml->SetBoundaryCondition(BC, size);
			op_ext_upml->SetRange(start,stop);
			op_child->AddExtension(op_ext_upml);
		}
		if (BC[5]==3)
		{
			op_ext_upml = new Operator_Ext_UPML(op_child);
			op_ext_upml->SetGradingFunction(gradFunc);
			start[2]=op->GetNumberOfLines(2,true)-1-size[5];
			stop[2] =op->GetNumberOfLines(2,true)-1;
			op_ext_upml->SetBoundaryCondition(BC, size);
			op_ext_upml->SetRange(start,stop);
			op_child->AddExtension(op_ext_upml);
		}
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
	if (res < 0) return true;

	cerr << "Operator_Ext_UPML::SetGradingFunction: Warning, an error occured parsing the pml grading function (see below) ..." << endl;
	cerr << func << "\n" << string(res, ' ') << "^\n" << m_GradingFunction->ErrorMsg() << "\n";
	return false;
}

void Operator_Ext_UPML::CalcGradingKappa(int ny, unsigned int pos[3], double Zm, double kappa_v[3], double kappa_i[3])
{
	double depth=0;
	double width=0;
	for (int n=0; n<3; ++n)
	{
		if ((pos[n] <= m_Size[2*n]) && (m_BC[2*n]==3))  //lower pml in n-dir
		{
			width = (m_Op->GetDiscLine(n,m_Size[2*n]) - m_Op->GetDiscLine(n,0))*m_Op->GetGridDelta();
			depth = width - (m_Op->GetDiscLine(n,pos[n]) - m_Op->GetDiscLine(n,0))*m_Op->GetGridDelta();

			if ((m_Op_Cyl) && (n==1))
			{
				width *= m_Op_Cyl->GetDiscLine(0,pos[0]);
				depth *= m_Op_Cyl->GetDiscLine(0,pos[0]);
			}

			if (n==ny)
				depth-=m_Op->GetEdgeLength(n,pos)/2;
			double vars[5] = {depth, width/m_Size[2*n], width, Zm, (double)m_Size[2*n]};
			if (depth>0)
				kappa_v[n] = m_GradingFunction->Eval(vars);
			else
				kappa_v[n]=0;
			if (n==ny)
				depth+=m_Op->GetEdgeLength(n,pos)/2;

			if (n!=ny)
				depth-=m_Op->GetEdgeLength(n,pos)/2;
			if (depth<0)
				depth=0;
			vars[0]=depth;
			if (depth>0)
				kappa_i[n] = m_GradingFunction->Eval(vars);
			else
				kappa_i[n] = 0;
		}
		else if ((pos[n] >= m_Op->GetNumberOfLines(n,true) -1 -m_Size[2*n+1]) && (m_BC[2*n+1]==3))  //upper pml in n-dir
		{
			width = (m_Op->GetDiscLine(n,m_Op->GetNumberOfLines(n,true)-1) - m_Op->GetDiscLine(n,m_Op->GetNumberOfLines(n,true)-m_Size[2*n+1]-1))*m_Op->GetGridDelta();
			depth = width - (m_Op->GetDiscLine(n,m_Op->GetNumberOfLines(n,true)-1) - m_Op->GetDiscLine(n,pos[n]))*m_Op->GetGridDelta();

			if ((m_Op_Cyl) && (n==1))
			{
				width *= m_Op_Cyl->GetDiscLine(0,pos[0]);
				depth *= m_Op_Cyl->GetDiscLine(0,pos[0]);
			}

			if (n==ny)
				depth+=m_Op->GetEdgeLength(n,pos)/2;
			double vars[5] = {depth, width/(m_Size[2*n]), width, Zm, (double)m_Size[2*n]};
			if (depth>0)
				kappa_v[n] = m_GradingFunction->Eval(vars);
			else
				kappa_v[n]=0;
			if (n==ny)
				depth-=m_Op->GetEdgeLength(n,pos)/2;

			if (n!=ny)
				depth+=m_Op->GetEdgeLength(n,pos)/2;
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
	double eff_Mat[4];
	double dT = m_Op->GetTimestep();

	for (loc_pos[0]=0; loc_pos[0]<m_numLines[0]; ++loc_pos[0])
	{
		pos[0] = loc_pos[0] + m_StartPos[0];
		for (loc_pos[1]=0; loc_pos[1]<m_numLines[1]; ++loc_pos[1])
		{
			pos[1] = loc_pos[1] + m_StartPos[1];
			vector<CSPrimitives*> vPrims = m_Op->GetPrimitivesBoundBox(pos[0], pos[1], -1, CSProperties::MATERIAL);
			for (loc_pos[2]=0; loc_pos[2]<m_numLines[2]; ++loc_pos[2])
			{
				pos[2] = loc_pos[2] + m_StartPos[2];
				for (int n=0; n<3; ++n)
				{
					m_Op->Calc_EffMatPos(n,pos,eff_Mat,vPrims);
					CalcGradingKappa(n, pos,__Z0__ ,kappa_v ,kappa_i);
					nP = (n+1)%3;
					nPP = (n+2)%3;
					if ((kappa_v[0]+kappa_v[1]+kappa_v[2])!=0)
					{
						//check if pos is on PEC
						if ( (m_Op->GetVV(n,pos[0],pos[1],pos[2]) + m_Op->GetVI(n,pos[0],pos[1],pos[2])) != 0 )
						{
							//modify the original operator to perform eq. (7.85) by the main engine (EC-FDTD: equation is multiplied by delta_n)
							//the engine extension will replace the original voltages with the "voltage flux" (volt*eps0) prior to the voltage updates
							//after the updates are done the extension will calculate the new voltages (see below) and place them back into the main field domain
							m_Op->SetVV(n,pos[0],pos[1],pos[2], (2*__EPS0__ - kappa_v[nP]*dT) / (2*__EPS0__ + kappa_v[nP]*dT) );
							m_Op->SetVI(n,pos[0],pos[1],pos[2], (2*__EPS0__*dT) / (2*__EPS0__ + kappa_v[nP]*dT) * m_Op->GetEdgeLength(n,pos) / m_Op->GetEdgeArea(n,pos) );


							//operators needed by eq. (7.88) to calculate new voltages from old voltages and old and new "voltage fluxes"
							GetVV(n,loc_pos)   = (2*__EPS0__ - kappa_v[nPP]*dT) / (2*__EPS0__ + kappa_v[nPP]*dT);
							GetVVFN(n,loc_pos) = (2*__EPS0__ + kappa_v[n]*dT)   / (2*__EPS0__ + kappa_v[nPP]*dT)/eff_Mat[0];
							GetVVFO(n,loc_pos) = (2*__EPS0__ - kappa_v[n]*dT)   / (2*__EPS0__ + kappa_v[nPP]*dT)/eff_Mat[0];
						}
					}
					else
					{
						//disable upml
						GetVV(n,loc_pos) = m_Op->GetVV(n,pos[0],pos[1],pos[2]);
						m_Op->SetVV(n,pos[0],pos[1],pos[2], 0 );
						GetVVFO(n,loc_pos) = 0;
						GetVVFN(n,loc_pos) = 1;
					}

					if ((kappa_i[0]+kappa_i[1]+kappa_i[2])!=0)
					{
						//check if pos is on PMC
						if ( (m_Op->GetII(n,pos[0],pos[1],pos[2]) + m_Op->GetIV(n,pos[0],pos[1],pos[2])) != 0 )
						{
							//modify the original operator to perform eq. (7.89) by the main engine (EC-FDTD: equation is multiplied by delta_n)
							//the engine extension will replace the original currents with the "current flux" (curr*mu0) prior to the current updates
							//after the updates are done the extension will calculate the new currents (see below) and place them back into the main field domain
							m_Op->SetII(n,pos[0],pos[1],pos[2], (2*__EPS0__ - kappa_i[nP]*dT) / (2*__EPS0__ + kappa_i[nP]*dT) );
							m_Op->SetIV(n,pos[0],pos[1],pos[2], (2*__EPS0__*dT) / (2*__EPS0__ + kappa_i[nP]*dT) * m_Op->GetEdgeLength(n,pos,true) / m_Op->GetEdgeArea(n,pos,true) );

							//operators needed by eq. (7.90) to calculate new currents from old currents and old and new "current fluxes"
							GetII(n,loc_pos)   = (2*__EPS0__ - kappa_i[nPP]*dT) / (2*__EPS0__ + kappa_i[nPP]*dT);
							GetIIFN(n,loc_pos) = (2*__EPS0__ + kappa_i[n]*dT)   / (2*__EPS0__ + kappa_i[nPP]*dT)/eff_Mat[2];
							GetIIFO(n,loc_pos) = (2*__EPS0__ - kappa_i[n]*dT)   / (2*__EPS0__ + kappa_i[nPP]*dT)/eff_Mat[2];
						}
					}
					else
					{
						//disable upml
						GetII(n,loc_pos) = m_Op->GetII(n,pos[0],pos[1],pos[2]);
						m_Op->SetII(n,pos[0],pos[1],pos[2], 0 );
						GetIIFO(n,loc_pos) = 0;
						GetIIFN(n,loc_pos) = 1;
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
