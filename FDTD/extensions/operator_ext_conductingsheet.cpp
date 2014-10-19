/*
*	Copyright (C) 2012 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#include "operator_ext_conductingsheet.h"
#include "tools/array_ops.h"
#include "tools/constants.h"
#include "cond_sheet_parameter.h"

#include "CSPropConductingSheet.h"

Operator_Ext_ConductingSheet::Operator_Ext_ConductingSheet(Operator* op, double f_max) : Operator_Ext_LorentzMaterial(op)
{
	m_f_max = f_max;
}

Operator_Ext_ConductingSheet::Operator_Ext_ConductingSheet(Operator* op, Operator_Ext_ConductingSheet* op_ext) : Operator_Ext_LorentzMaterial(op, op_ext)
{
	m_f_max = op_ext->m_f_max;
}

Operator_Extension* Operator_Ext_ConductingSheet::Clone(Operator* op)
{
	if (dynamic_cast<Operator_Ext_ConductingSheet*>(this)==NULL)
		return NULL;
	return new Operator_Ext_ConductingSheet(op, this);
}

bool Operator_Ext_ConductingSheet::BuildExtension()
{
	double dT = m_Op->GetTimestep();
	unsigned int pos[] = {0,0,0};
	double coord[3];
	unsigned int numLines[3] = {m_Op->GetNumberOfLines(0,true),m_Op->GetNumberOfLines(1,true),m_Op->GetNumberOfLines(2,true)};

	m_Order = 0;
	vector<unsigned int> v_pos[3];
	int ****tanDir = Create_N_3DArray<int>(numLines);
	float ****Conductivity = Create_N_3DArray<float>(numLines);
	float ****Thickness = Create_N_3DArray<float>(numLines);

	CSPrimitives* cs_sheet = NULL;
	double box[6];
	int nP, nPP;
	bool b_pos_on;
	bool disable_pos;
	for (pos[0]=0; pos[0]<numLines[0]; ++pos[0])
	{
		for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
		{
			vector<CSPrimitives*> vPrims = m_Op->GetPrimitivesBoundBox(pos[0], pos[1], -1, (CSProperties::PropertyType)(CSProperties::MATERIAL | CSProperties::METAL));
			for (pos[2]=0; pos[2]<numLines[2]; ++pos[2])
			{
				b_pos_on = false;
				disable_pos = false;
				// disable conducting sheet model inside the boundary conditions, especially inside a pml
				for (int m=0;m<3;++m)
					if ((pos[m]<=(unsigned int)m_Op->GetBCSize(2*m)) || (pos[m]>=(numLines[m]-m_Op->GetBCSize(2*m+1)-1)))
						disable_pos = true;

				for (int n=0; n<3; ++n)
				{
					nP = (n+1)%3;
					nPP = (n+2)%3;

					tanDir[n][pos[0]][pos[1]][pos[2]] = -1; //deactivate by default
					Conductivity[n][pos[0]][pos[1]][pos[2]] = 0; //deactivate by default
					Thickness[n][pos[0]][pos[1]][pos[2]] = 0; //deactivate by default

					if (m_Op->GetYeeCoords(n,pos,coord,false)==false)
						continue;

					// Ez at r==0 not supported --> set to PEC
					if (m_CC_R0_included && (n==2) && (pos[0]==0))
						disable_pos = true;

//					CSProperties* prop = m_Op->GetGeometryCSX()->GetPropertyByCoordPriority(coord,(CSProperties::PropertyType)(CSProperties::METAL | CSProperties::MATERIAL), false, &cs_sheet);
					CSProperties* prop = m_Op->GetGeometryCSX()->GetPropertyByCoordPriority(coord, vPrims, false, &cs_sheet);
					CSPropConductingSheet* cs_prop = dynamic_cast<CSPropConductingSheet*>(prop);
					if (cs_prop)
					{
						if (cs_sheet==NULL)
							return false; //sanity check, this should never happen
						if (cs_sheet->GetDimension()!=2)
						{
							cerr << "Operator_Ext_ConductingSheet::BuildExtension: A conducting sheet primitive (ID: " << cs_sheet->GetID() << ") with dimension: " << cs_sheet->GetDimension() << " found, fallback to PEC!" << endl;
							m_Op->SetVV(n,pos[0],pos[1],pos[2], 0 );
							m_Op->SetVI(n,pos[0],pos[1],pos[2], 0 );
							++m_Op->m_Nr_PEC[n];
							continue;
						}
						cs_sheet->SetPrimitiveUsed(true);

						if (disable_pos)
						{
							m_Op->SetVV(n,pos[0],pos[1],pos[2], 0 );
							m_Op->SetVI(n,pos[0],pos[1],pos[2], 0 );
							++m_Op->m_Nr_PEC[n];
							continue;
						}

						Conductivity[n][pos[0]][pos[1]][pos[2]] = cs_prop->GetConductivity();
						Thickness[n][pos[0]][pos[1]][pos[2]] = cs_prop->GetThickness();

						if ((Conductivity[n][pos[0]][pos[1]][pos[2]]<=0) || (Thickness[n][pos[0]][pos[1]][pos[2]]<=0))
						{
							cerr << "Operator_Ext_ConductingSheet::BuildExtension: Warning: Zero conductivity or thickness detected... fallback to PEC!" << endl;
							m_Op->SetVV(n,pos[0],pos[1],pos[2], 0 );
							m_Op->SetVI(n,pos[0],pos[1],pos[2], 0 );
							++m_Op->m_Nr_PEC[n];
							continue;
						}

						cs_sheet->GetBoundBox(box);
						if (box[2*nP]!=box[2*nP+1])
							tanDir[n][pos[0]][pos[1]][pos[2]] = nP;
						if (box[2*nPP]!=box[2*nPP+1])
							tanDir[n][pos[0]][pos[1]][pos[2]] = nPP;
						b_pos_on = true;
					}
				}
				if (b_pos_on)
				{
					for (int n=0; n<3; ++n)
						v_pos[n].push_back(pos[n]);
				}
			}
		}
	}

	size_t numCS = v_pos[0].size();
	if (numCS==0)
		return false;

	m_LM_Count.push_back(numCS);
	m_LM_Count.push_back(numCS);

	m_Order	= 2;
	m_volt_ADE_On = new bool[m_Order];
	m_volt_ADE_On[0] = m_volt_ADE_On[1]=true;
	m_curr_ADE_On = new bool[m_Order];
	m_curr_ADE_On[0] = m_curr_ADE_On[1]=false;

	m_volt_Lor_ADE_On = new bool[m_Order];
	m_volt_Lor_ADE_On[0] = m_volt_Lor_ADE_On[1]=false;
	m_curr_Lor_ADE_On = new bool[m_Order];
	m_curr_Lor_ADE_On[0] = m_curr_Lor_ADE_On[1]=false;

	m_LM_pos = new unsigned int**[m_Order];
	m_LM_pos[0] = new unsigned int*[3];
	m_LM_pos[1] = new unsigned int*[3];

	v_int_ADE = new FDTD_FLOAT**[m_Order];
	v_ext_ADE = new FDTD_FLOAT**[m_Order];

	v_int_ADE[0] = new FDTD_FLOAT*[3];
	v_ext_ADE[0] = new FDTD_FLOAT*[3];
	v_int_ADE[1] = new FDTD_FLOAT*[3];
	v_ext_ADE[1] = new FDTD_FLOAT*[3];

	for (int n=0; n<3; ++n)
	{
		m_LM_pos[0][n] = new unsigned int[numCS];
		m_LM_pos[1][n] = new unsigned int[numCS];
		for (unsigned int i=0; i<numCS; ++i)
		{
			m_LM_pos[0][n][i] = v_pos[n].at(i);
			m_LM_pos[1][n][i] = v_pos[n].at(i);
		}
		v_int_ADE[0][n]  = new FDTD_FLOAT[numCS];
		v_int_ADE[1][n]  = new FDTD_FLOAT[numCS];
		v_ext_ADE[0][n]  = new FDTD_FLOAT[numCS];
		v_ext_ADE[1][n]  = new FDTD_FLOAT[numCS];
	}

	unsigned int index;
	float w_stop = m_f_max*2*PI;
	float Omega_max=0;
	float G,L1,L2,R1,R2,Lmin;
	float G0, w0;
	float wtl; //width to length factor
	float factor=1;
	int t_dir=0; //tangential sheet direction
	unsigned int tpos[] = {0,0,0};
	unsigned int optParaPos;
	for (unsigned int i=0;i<numCS;++i)
	{
		pos[0]=m_LM_pos[0][0][i];pos[1]=m_LM_pos[0][1][i];pos[2]=m_LM_pos[0][2][i];
		tpos[0]=pos[0];tpos[1]=pos[1];tpos[2]=pos[2];
		index = m_Op->MainOp->SetPos(pos[0],pos[1],pos[2]);
		for (int n=0;n<3;++n)
		{
			tpos[0]=pos[0];tpos[1]=pos[1];tpos[2]=pos[2];
			t_dir = tanDir[n][pos[0]][pos[1]][pos[2]];
			G0 = Conductivity[n][pos[0]][pos[1]][pos[2]]*Thickness[n][pos[0]][pos[1]][pos[2]];
			w0 = 8.0/ G0 / Thickness[n][pos[0]][pos[1]][pos[2]]/__MUE0__;
			Omega_max = w_stop/w0;
			for (optParaPos=0;optParaPos<numOptPara;++optParaPos)
				if (omega_stop[optParaPos]>Omega_max)
					break;
			if (optParaPos>=numOptPara)
			{
				cerr << "Operator_Ext_ConductingSheet::BuildExtension(): Error, conductor thickness, conductivity or max. simulation frequency of interest is too high! Check parameter!" << endl;
				cerr << " --> max f: " << m_f_max << "Hz,  Conductivity: " << Conductivity[n][pos[0]][pos[1]][pos[2]] << "S/m, Thickness " << Thickness[n][pos[0]][pos[1]][pos[2]]*1e6 << "um" << endl;
				optParaPos = numOptPara-1;
			}
			v_int_ADE[0][n][i]=0;
			v_ext_ADE[0][n][i]=0;
			v_int_ADE[1][n][i]=0;
			v_ext_ADE[1][n][i]=0;
			if (t_dir>=0)
			{
				wtl = m_Op->GetEdgeLength(n,pos)/m_Op->GetNodeWidth(t_dir,pos);
				factor = 1;
				if (tanDir[t_dir][tpos[0]][tpos[1]][tpos[2]]<0)
					factor = 2;
				--tpos[t_dir];
				if (tanDir[t_dir][tpos[0]][tpos[1]][tpos[2]]<0)
					factor = 2;

				L1 = l1[optParaPos]/G0/w0*factor;
				L2 = l2[optParaPos]/G0/w0*factor;
				R1 = r1[optParaPos]/G0*factor;
				R2 = r2[optParaPos]/G0*factor;
				G = G0*g[optParaPos]/factor;

				L1*=wtl;
				L2*=wtl;
				R1*=wtl;
				R2*=wtl;
				G/=wtl;

				Lmin = L1;
				if (L2<L1)
					Lmin = L2;
				m_Op->EC_G[n][index]= G;
				m_Op->EC_C[n][index]= dT*dT/4.0*(16.0/Lmin + 1/L1 + 1/L2);
				m_Op->Calc_ECOperatorPos(n,pos);

				v_int_ADE[0][n][i]=(2.0*L1-dT*R1)/(2.0*L1+dT*R1);
				v_ext_ADE[0][n][i]=dT/(L1+dT*R1/2.0)*m_Op->GetVI(n,pos[0],pos[1],pos[2]);
				v_int_ADE[1][n][i]=(2.0*L2-dT*R2)/(2.0*L2+dT*R2);
				v_ext_ADE[1][n][i]=dT/(L2+dT*R2/2.0)*m_Op->GetVI(n,pos[0],pos[1],pos[2]);
			}
		}
	}
	return true;
}
