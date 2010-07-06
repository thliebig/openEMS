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

#include "operator_ext_lorentzmaterial.h"
#include "engine_ext_lorentzmaterial.h"

Operator_Ext_LorentzMaterial::Operator_Ext_LorentzMaterial(Operator* op) : Operator_Ext_Dispersive(op)
{
	for (int n=0;n<3;++n)
	{
		v_int_ADE[n] = NULL;
		v_ext_ADE[n] = NULL;
		i_int_ADE[n] = NULL;
		i_ext_ADE[n] = NULL;
	}
}

Operator_Ext_LorentzMaterial::~Operator_Ext_LorentzMaterial()
{
	for (int n=0;n<3;++n)
	{
		delete[] v_int_ADE[n];
		v_int_ADE[n] = NULL;
		delete[] v_ext_ADE[n];
		v_ext_ADE[n] = NULL;
		delete[] i_int_ADE[n];
		i_int_ADE[n] = NULL;
		delete[] i_ext_ADE[n];
		i_ext_ADE[n] = NULL;
	}
}

bool Operator_Ext_LorentzMaterial::BuildExtension()
{
	double dT = m_Op->GetTimestep();
	unsigned int pos[] = {0,0,0};
	double coord[3];
	unsigned int numLines[3] = {m_Op->GetNumberOfLines(0),m_Op->GetNumberOfLines(1),m_Op->GetNumberOfLines(2)};
	CSPropLorentzMaterial* mat = NULL;

	double w_plasma;
	bool b_pos_on;
	double L_D[3], C_D[3];
	vector<double> v_int[3];
	vector<double> v_ext[3];
	vector<double> i_int[3];
	vector<double> i_ext[3];
	vector<unsigned int> v_pos[3];

	for (pos[0]=0;pos[0]<numLines[0];++pos[0])
	{
		for (pos[1]=0;pos[1]<numLines[1];++pos[1])
		{
			for (pos[2]=0;pos[2]<numLines[2];++pos[2])
			{
				unsigned int index = m_Op->MainOp->SetPos(pos[0],pos[1],pos[2]);
				//calc epsilon lorentz material
				b_pos_on = false;
				for (int n=0;n<3;++n)
				{
					L_D[n]=0;
					coord[0] = m_Op->GetDiscLine(0,pos[0]);
					coord[1] = m_Op->GetDiscLine(1,pos[1]);
					coord[2] = m_Op->GetDiscLine(2,pos[2]);
					coord[n] = m_Op->GetDiscLine(n,pos[n],true); //pos of E_n

					CSProperties* prop = m_Op->GetGeometryCSX()->GetPropertyByCoordPriority(coord,CSProperties::LORENTZMATERIAL);
					if ((mat = prop->ToLorentzMaterial()))
					{
						w_plasma = mat->GetEpsPlasmaFreqWeighted(n,coord) * 2 * PI;
						if (w_plasma)
						{
							b_pos_on = true;
							m_volt_ADE_On = true;
							L_D[n] = 1/(w_plasma*w_plasma*m_Op->EC_C[n][index]);
						}
					}
				}

				for (int n=0;n<3;++n)
				{
					C_D[n]=0;
					coord[0] = m_Op->GetDiscLine(0,pos[0],true);
					coord[1] = m_Op->GetDiscLine(1,pos[1],true);
					coord[2] = m_Op->GetDiscLine(2,pos[2],true);
					coord[n] = m_Op->GetDiscLine(n,pos[n]); //pos of H_n

					CSProperties* prop = m_Op->GetGeometryCSX()->GetPropertyByCoordPriority(coord,CSProperties::LORENTZMATERIAL);
					if ((mat = prop->ToLorentzMaterial()))
					{
						w_plasma = mat->GetMuePlasmaFreqWeighted(n,coord) * 2 * PI;
						if (w_plasma)
						{
							b_pos_on = true;
							m_curr_ADE_On = true;
							C_D[n] = 1/(w_plasma*w_plasma*m_Op->EC_L[n][index]);
						}
					}
				}

				if (b_pos_on) //this position has active lorentz material
				{
					for (unsigned int n=0;n<3;++n)
					{
						v_pos[n].push_back(pos[n]);
						v_int[n].push_back(1);
						i_int[n].push_back(1);
						if (L_D[n]>0)
							v_ext[n].push_back(dT/L_D[n]	*m_Op->GetVI(n,pos[0],pos[1],pos[2]));
						else
							v_ext[n].push_back(0);
						if (C_D[n]>0)
							i_ext[n].push_back(dT/C_D[n]	*m_Op->GetIV(n,pos[0],pos[1],pos[2]));
						else
							i_ext[n].push_back(0);
//						cerr << v_int[n].back() << " " << v_ext[n].back() << " " << i_int[n].back() << " " << i_ext[n].back() << endl;
					}
				}
			}
		}
	}

	//copy all vectors into the array's
	m_LM_Count = v_pos[0].size();

	for (int n=0;n<3;++n)
	{
		m_LM_pos[n] = new unsigned int[m_LM_Count];
		for (unsigned int i=0;i<m_LM_Count;++i)
			m_LM_pos[n][i] = v_pos[n].at(i);
		if (m_volt_ADE_On)
		{
			v_int_ADE[n]  = new FDTD_FLOAT[m_LM_Count];
			v_ext_ADE[n]  = new FDTD_FLOAT[m_LM_Count];

			for (unsigned int i=0;i<m_LM_Count;++i)
			{
				v_int_ADE[n][i] = v_int[n].at(i);
				v_ext_ADE[n][i] = v_ext[n].at(i);
			}
		}
		if (m_curr_ADE_On)
		{
			i_int_ADE[n]  = new FDTD_FLOAT[m_LM_Count];
			i_ext_ADE[n]  = new FDTD_FLOAT[m_LM_Count];

			for (unsigned int i=0;i<m_LM_Count;++i)
			{
				i_int_ADE[n][i] = i_int[n].at(i);
				i_ext_ADE[n][i] = i_ext[n].at(i);
			}
		}

	}

	return true;
}

Engine_Extension* Operator_Ext_LorentzMaterial::CreateEngineExtention()
{
	Engine_Ext_LorentzMaterial* eng_ext_lor = new Engine_Ext_LorentzMaterial(this);
	return eng_ext_lor;
}
