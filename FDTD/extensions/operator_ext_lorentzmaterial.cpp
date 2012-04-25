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
	v_int_ADE = NULL;
	v_ext_ADE = NULL;
	i_int_ADE = NULL;
	i_ext_ADE = NULL;
}

Operator_Ext_LorentzMaterial::~Operator_Ext_LorentzMaterial()
{
	for (int i=0;i<m_Order;++i)
	{
		for (int n=0; n<3; ++n)
		{
			if (m_volt_ADE_On[i])
			{
				delete[] v_int_ADE[i][n];
				delete[] v_ext_ADE[i][n];
			}
			if (m_curr_ADE_On[i])
			{
				delete[] i_int_ADE[i][n];
				delete[] i_ext_ADE[i][n];
			}
		}
		if (m_volt_ADE_On[i])
		{
			delete[] v_int_ADE[i];
			delete[] v_ext_ADE[i];
		}
		if (m_curr_ADE_On[i])
		{
			delete[] i_int_ADE[i];
			delete[] i_ext_ADE[i];
		}
	}
	delete[] v_int_ADE;
	delete[] v_ext_ADE;
	delete[] i_int_ADE;
	delete[] i_ext_ADE;
	v_int_ADE = NULL;
	v_ext_ADE = NULL;
	i_int_ADE = NULL;
	i_ext_ADE = NULL;
}

bool Operator_Ext_LorentzMaterial::BuildExtension()
{
	double dT = m_Op->GetTimestep();
	unsigned int pos[] = {0,0,0};
	double coord[3];
	unsigned int numLines[3] = {m_Op->GetNumberOfLines(0),m_Op->GetNumberOfLines(1),m_Op->GetNumberOfLines(2)};
	CSPropLorentzMaterial* mat = NULL;

	double w_plasma,t_relax;
	bool b_pos_on;
	double L_D[3], C_D[3];
	double R_D[3], G_D[3];
	vector<double> v_int[3];
	vector<double> v_ext[3];
	vector<double> i_int[3];
	vector<double> i_ext[3];
	vector<unsigned int> v_pos[3];

	m_Order = 0;
	vector<CSProperties*> LD_props = m_Op->CSX->GetPropertyByType(CSProperties::LORENTZMATERIAL);
	for (size_t n=0;n<LD_props.size();++n)
	{
		CSPropLorentzMaterial* LorMat = dynamic_cast<CSPropLorentzMaterial*>(LD_props.at(n));
		if (LorMat==NULL)
			return false; //sanity check, this should not happen
		if (LorMat->GetDispersionOrder()>m_Order)
			m_Order=LorMat->GetDispersionOrder();
	}

	m_volt_ADE_On = new bool[m_Order];
	m_curr_ADE_On = new bool[m_Order];
	m_LM_pos = new unsigned int**[m_Order];

	v_int_ADE = new FDTD_FLOAT**[m_Order];
	v_ext_ADE = new FDTD_FLOAT**[m_Order];
	i_int_ADE = new FDTD_FLOAT**[m_Order];
	i_ext_ADE = new FDTD_FLOAT**[m_Order];

	for (int order=0;order<m_Order;++order)
	{
		m_volt_ADE_On[order]=false;
		m_curr_ADE_On[order]=false;

		for (int n=0;n<3;++n)
		{
			v_int[n].clear();
			v_ext[n].clear();
			i_int[n].clear();
			i_ext[n].clear();
			v_pos[n].clear();
		}

		for (pos[0]=0; pos[0]<numLines[0]; ++pos[0])
		{
			for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
			{
				for (pos[2]=0; pos[2]<numLines[2]; ++pos[2])
				{
					unsigned int index = m_Op->MainOp->SetPos(pos[0],pos[1],pos[2]);
					//calc epsilon lorentz material
					b_pos_on = false;
					for (int n=0; n<3; ++n)
					{
						L_D[n]=0;
						R_D[n]=0;
						coord[0] = m_Op->GetDiscLine(0,pos[0]);
						coord[1] = m_Op->GetDiscLine(1,pos[1]);
						coord[2] = m_Op->GetDiscLine(2,pos[2]);
						coord[n] = m_Op->GetDiscLine(n,pos[n],true); //pos of E_n

						CSProperties* prop = m_Op->GetGeometryCSX()->GetPropertyByCoordPriority(coord,CSProperties::LORENTZMATERIAL, true);
						if ((mat = prop->ToLorentzMaterial()))
						{
							w_plasma = mat->GetEpsPlasmaFreqWeighted(order,n,coord) * 2 * PI;
							if (w_plasma>0)
							{
								b_pos_on = true;
								m_volt_ADE_On[order] = true;
								L_D[n] = 1/(w_plasma*w_plasma*m_Op->EC_C[n][index]);
							}
							t_relax = mat->GetEpsRelaxTimeWeighted(order,n,coord);
							if ((t_relax>0) && m_volt_ADE_On[order])
							{
								R_D[n] = L_D[n]/t_relax;
							}
						}
					}

					for (int n=0; n<3; ++n)
					{
						C_D[n]=0;
						G_D[n]=0;
						coord[0] = m_Op->GetDiscLine(0,pos[0],true);
						coord[1] = m_Op->GetDiscLine(1,pos[1],true);
						coord[2] = m_Op->GetDiscLine(2,pos[2],true);
						coord[n] = m_Op->GetDiscLine(n,pos[n]); //pos of H_n

						CSProperties* prop = m_Op->GetGeometryCSX()->GetPropertyByCoordPriority(coord,CSProperties::LORENTZMATERIAL, true);
						if ((mat = prop->ToLorentzMaterial()))
						{
							w_plasma = mat->GetMuePlasmaFreqWeighted(order,n,coord) * 2 * PI;
							if (w_plasma>0)
							{
								b_pos_on = true;
								m_curr_ADE_On[order] = true;
								C_D[n] = 1/(w_plasma*w_plasma*m_Op->EC_L[n][index]);
							}
							t_relax = mat->GetMueRelaxTimeWeighted(order,n,coord);
							if ((t_relax>0) && m_curr_ADE_On[order])
							{
								G_D[n] = C_D[n]/t_relax;
							}
						}
					}

					if (b_pos_on) //this position has active lorentz material
					{
						for (unsigned int n=0; n<3; ++n)
						{
							v_pos[n].push_back(pos[n]);
							if (L_D[n]>0)
							{
								v_int[n].push_back((2*L_D[n]-dT*R_D[n])/(2*L_D[n]+dT*R_D[n]));
								v_ext[n].push_back(dT/(L_D[n]+dT*R_D[n]/2)*m_Op->GetVI(n,pos[0],pos[1],pos[2]));
							}
							else
							{
								v_int[n].push_back(1);
								v_ext[n].push_back(0);
							}
							if (C_D[n]>0)
							{
								i_int[n].push_back((2*C_D[n]-dT*G_D[n])/(2*C_D[n]+dT*G_D[n]));
								i_ext[n].push_back(dT/(C_D[n]+dT*G_D[n]/2)*m_Op->GetIV(n,pos[0],pos[1],pos[2]));
							}
							else
							{
								i_int[n].push_back(1);
								i_ext[n].push_back(0);
							}
							//						cerr << v_int[n].back() << " " << v_ext[n].back() << " " << i_int[n].back() << " " << i_ext[n].back() << endl;
						}
					}
				}
			}
		}

		//copy all vectors into the array's
		m_LM_Count.push_back(v_pos[0].size());

		m_LM_pos[order] = new unsigned int*[3];

		v_int_ADE[order] = new FDTD_FLOAT*[3];
		v_ext_ADE[order] = new FDTD_FLOAT*[3];
		i_int_ADE[order] = new FDTD_FLOAT*[3];
		i_ext_ADE[order] = new FDTD_FLOAT*[3];

		for (int n=0; n<3; ++n)
		{
			m_LM_pos[order][n] = new unsigned int[m_LM_Count.at(order)];
			for (unsigned int i=0; i<m_LM_Count.at(order); ++i)
				m_LM_pos[order][n][i] = v_pos[n].at(i);
			if (m_volt_ADE_On)
			{
				v_int_ADE[order][n]  = new FDTD_FLOAT[m_LM_Count.at(order)];
				v_ext_ADE[order][n]  = new FDTD_FLOAT[m_LM_Count.at(order)];

				for (unsigned int i=0; i<m_LM_Count.at(order); ++i)
				{
					v_int_ADE[order][n][i] = v_int[n].at(i);
					v_ext_ADE[order][n][i] = v_ext[n].at(i);
				}
			}
			if (m_curr_ADE_On)
			{
				i_int_ADE[order][n]  = new FDTD_FLOAT[m_LM_Count.at(order)];
				i_ext_ADE[order][n]  = new FDTD_FLOAT[m_LM_Count.at(order)];

				for (unsigned int i=0; i<m_LM_Count.at(order); ++i)
				{
					i_int_ADE[order][n][i] = i_int[n].at(i);
					i_ext_ADE[order][n][i] = i_ext[n].at(i);
				}
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

void Operator_Ext_LorentzMaterial::ShowStat(ostream &ostr)  const
{
	Operator_Ext_Dispersive::ShowStat(ostr);
}
