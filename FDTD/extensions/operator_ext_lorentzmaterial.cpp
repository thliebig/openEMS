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

#include <algorithm>
#include <utility>
#include "operator_ext_lorentzmaterial.h"
#include "engine_ext_lorentzmaterial.h"
#include "operator_ext_cylinder.h"
#include "../operator_cylinder.h"

#include "CSPropLorentzMaterial.h"

using std::endl;

Operator_Ext_LorentzMaterial::Operator_Ext_LorentzMaterial(Operator* op) : Operator_Ext_Dispersive(op)
{
	v_int_ADE = NULL;
	v_ext_ADE = NULL;
	i_int_ADE = NULL;
	i_ext_ADE = NULL;

	v_Lor_ADE = NULL;
	i_Lor_ADE = NULL;

	m_curr_Lor_ADE_On = NULL;
	m_volt_Lor_ADE_On = NULL;
}

Operator_Ext_LorentzMaterial::Operator_Ext_LorentzMaterial(Operator* op, Operator_Ext_LorentzMaterial* op_ext) : Operator_Ext_Dispersive(op,op_ext)
{
	v_int_ADE = NULL;
	v_ext_ADE = NULL;
	i_int_ADE = NULL;
	i_ext_ADE = NULL;

	v_Lor_ADE = NULL;
	i_Lor_ADE = NULL;

	m_curr_Lor_ADE_On = NULL;
	m_volt_Lor_ADE_On = NULL;
}

Operator_Ext_LorentzMaterial::~Operator_Ext_LorentzMaterial()
{
	if ((m_curr_Lor_ADE_On==NULL) || (m_volt_Lor_ADE_On==NULL)) // all data already cleaned up or never used
		return;
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
			if (m_volt_Lor_ADE_On[i])
				delete[] v_Lor_ADE[i][n];
			if (m_curr_Lor_ADE_On[i])
				delete[] i_Lor_ADE[i][n];
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
		if (m_volt_Lor_ADE_On[i])
			delete[] v_Lor_ADE[i];
		if (m_curr_Lor_ADE_On[i])
			delete[] i_Lor_ADE[i];
	}
	delete[] v_int_ADE;
	delete[] v_ext_ADE;
	delete[] i_int_ADE;
	delete[] i_ext_ADE;
	v_int_ADE = NULL;
	v_ext_ADE = NULL;
	i_int_ADE = NULL;
	i_ext_ADE = NULL;

	delete[] v_Lor_ADE;
	delete[] i_Lor_ADE;
	v_Lor_ADE = NULL;
	i_Lor_ADE = NULL;

	delete[] m_curr_Lor_ADE_On;
	delete[] m_volt_Lor_ADE_On;
	m_curr_Lor_ADE_On = NULL;
	m_volt_Lor_ADE_On = NULL;
}

Operator_Extension* Operator_Ext_LorentzMaterial::Clone(Operator* op)
{
	return new Operator_Ext_LorentzMaterial(op, this);
}

bool Operator_Ext_LorentzMaterial::BuildExtension()
{
	double dT = m_Op->GetTimestep();
	unsigned int numLines[3] = {m_Op->GetNumberOfLines(0,true),m_Op->GetNumberOfLines(1,true),m_Op->GetNumberOfLines(2,true)};

	//! the cells a thread found, and whether they need the ADEs
	struct Cells
	{
		std::vector<unsigned int> pos[3];
		std::vector<FDTD_FLOAT> v_int[3], v_ext[3], i_int[3], i_ext[3], v_Lor[3], i_Lor[3];   //!< FDTD_FLOAT as the arrays: the same conversion, earlier
		bool volt = false, curr = false, volt_Lor = false, curr_Lor = false;
	};

	m_Order = 0;
	std::vector<CSProperties*> LD_props = m_Op->CSX->GetPropertyByType(CSProperties::LORENTZMATERIAL);
	for (size_t n=0;n<LD_props.size();++n)
	{
		CSPropLorentzMaterial* LorMat = dynamic_cast<CSPropLorentzMaterial*>(LD_props.at(n));
		if (LorMat==NULL)
			return false; //sanity check, this should not happen
		if (LorMat->GetDispersionOrder()>m_Order)
			m_Order=LorMat->GetDispersionOrder();
	}
	m_LM_pos = new unsigned int**[m_Order];

	m_volt_ADE_On = new bool[m_Order];
	m_curr_ADE_On = new bool[m_Order];
	m_volt_Lor_ADE_On = new bool[m_Order];
	m_curr_Lor_ADE_On = new bool[m_Order];

	v_int_ADE = new FDTD_FLOAT**[m_Order];
	v_ext_ADE = new FDTD_FLOAT**[m_Order];
	i_int_ADE = new FDTD_FLOAT**[m_Order];
	i_ext_ADE = new FDTD_FLOAT**[m_Order];

	v_Lor_ADE = new FDTD_FLOAT**[m_Order];
	i_Lor_ADE = new FDTD_FLOAT**[m_Order];

	for (int order=0;order<m_Order;++order)
	{
		m_volt_ADE_On[order]=false;
		m_curr_ADE_On[order]=false;

		m_volt_Lor_ADE_On[order]=false;
		m_curr_Lor_ADE_On[order]=false;

		// the columns (x, y) on the build threads, their cells joined in the order of the columns;
		// a column without a primitive of a Lorentz material in its bounding box has none
		std::vector<Cells> found(m_Op->GetBuildThreads());
		m_Op->ParallelFor((size_t)numLines[0]*numLines[1], [&](unsigned int thread, size_t begin, size_t end) {
			Cells& c = found[thread];
			unsigned int pos[3];
			double coord[3];
			bool b_pos_on;
			CSPropLorentzMaterial* mat = NULL;
			double w_plasma, t_relax, w_Lor_Pol;
			double L_D[3], C_D[3], R_D[3], G_D[3], C_L[3], L_L[3];
			for (size_t col=begin; col<end; ++col)
			{
				pos[0] = col/numLines[1];
				pos[1] = col%numLines[1];
				std::vector<CSPrimitives*> vPrims = m_Op->GetPrimitivesBoundBox(
					pos[0], pos[1], -1,
					(CSProperties::PropertyType)(CSProperties::MATERIAL | CSProperties::METAL)
				);
				bool any = false;
				for (size_t p=0; p<vPrims.size() && !any; ++p)
					any = vPrims[p]->GetProperty() && vPrims[p]->GetProperty()->ToLorentzMaterial();
				if (!any)
					continue;

				for (pos[2]=0; pos[2]<numLines[2]; ++pos[2])
				{
					const size_t index = m_Op->GetECIndex(pos);
					//calc epsilon lorentz material
					b_pos_on = false;
					for (int n=0; n<3; ++n)
					{
						L_D[n]=0;
						R_D[n]=0;
						C_L[n]=0;
						if (m_Op->GetYeeCoords(n,pos,coord,false)==false)
							continue;
						if (m_CC_R0_included && (n==2) && (pos[0]==0))
							coord[1] = m_Op->GetDiscLine(1,0);

						if (m_Op->GetVI(n,pos[0],pos[1],pos[2])==0)
							continue;

//						CSProperties* prop = m_Op->GetGeometryCSX()->GetPropertyByCoordPriority(coord,(CSProperties::PropertyType)(CSProperties::METAL | CSProperties::MATERIAL), true);
						CSProperties* prop = m_Op->GetGeometryCSX()->GetPropertyByCoordPriority(coord, vPrims, true);

						if (prop==NULL) continue;

						if ((mat = prop->ToLorentzMaterial()))
						{
							w_plasma = mat->GetEpsPlasmaFreqWeighted(order,n,coord) * 2 * PI;
							if ((w_plasma>0) && (m_Op->EC_C[n][index]>0))
							{
								b_pos_on = true;
								c.volt = true;
								L_D[n] = 1/(w_plasma*w_plasma*m_Op->EC_C[n][index]);
							}
							t_relax = mat->GetEpsRelaxTimeWeighted(order,n,coord);
							if ((t_relax>0) && c.volt)
							{
								R_D[n] = L_D[n]/t_relax;
							}
							w_Lor_Pol = mat->GetEpsLorPoleFreqWeighted(order,n,coord) * 2 * PI;
							if ((w_Lor_Pol>0) && (L_D[n]>0))
							{
								c.volt_Lor = true;
								C_L[n] = 1/(w_Lor_Pol*w_Lor_Pol*L_D[n]);
							}
						}
					}

					for (int n=0; n<3; ++n)
					{
						C_D[n]=0;
						G_D[n]=0;
						L_L[n]=0;
						if (m_Op->GetYeeCoords(n,pos,coord,true)==false)
							continue;
						if (m_Op->GetIV(n,pos[0],pos[1],pos[2])==0)
							continue;

//						CSProperties* prop = m_Op->GetGeometryCSX()->GetPropertyByCoordPriority(coord,(CSProperties::PropertyType)(CSProperties::METAL | CSProperties::MATERIAL), true);
						CSProperties* prop = m_Op->GetGeometryCSX()->GetPropertyByCoordPriority(coord, vPrims, true);

						if (prop==NULL) continue;

						if ((mat = prop->ToLorentzMaterial()))
						{
							w_plasma = mat->GetMuePlasmaFreqWeighted(order,n,coord) * 2 * PI;
							if ((w_plasma>0) && (m_Op->EC_L[n][index]>0))
							{
								b_pos_on = true;
								c.curr = true;
								C_D[n] = 1/(w_plasma*w_plasma*m_Op->EC_L[n][index]);
							}
							t_relax = mat->GetMueRelaxTimeWeighted(order,n,coord);
							if ((t_relax>0) && c.curr)
							{
								G_D[n] = C_D[n]/t_relax;
							}
							w_Lor_Pol = mat->GetMueLorPoleFreqWeighted(order,n,coord) * 2 * PI;
							if ((w_Lor_Pol>0) && (C_D[n]>0))
							{
								c.curr_Lor = true;
								L_L[n] = 1/(w_Lor_Pol*w_Lor_Pol*C_D[n]);
							}
						}
					}

					if (b_pos_on) //this position has active drude material
					{
						for (unsigned int n=0; n<3; ++n)
						{
							c.pos[n].push_back(pos[n]);
							if (L_D[n]>0)
							{
								c.v_int[n].push_back((2.0*L_D[n]-dT*R_D[n])/(2.0*L_D[n]+dT*R_D[n]));
								// check for r==0 in clyindrical coords and get special VI coefficient
								if (m_CC_R0_included && n==2 && pos[0]==0)
									c.v_ext[n].push_back(dT/(L_D[n]+dT*R_D[n]/2.0)*m_Op_Cyl->m_Cyl_Ext->vi_R0[pos[2]]);
								else
									c.v_ext[n].push_back(dT/(L_D[n]+dT*R_D[n]/2.0)*m_Op->GetVI(n,pos[0],pos[1],pos[2]));
							}
							else
							{
								c.v_int[n].push_back(1);
								c.v_ext[n].push_back(0);
							}
							if (C_D[n]>0)
							{
								c.i_int[n].push_back((2.0*C_D[n]-dT*G_D[n])/(2.0*C_D[n]+dT*G_D[n]));
								c.i_ext[n].push_back(dT/(C_D[n]+dT*G_D[n]/2.0)*m_Op->GetIV(n,pos[0],pos[1],pos[2]));
							}
							else
							{
								c.i_int[n].push_back(1);
								c.i_ext[n].push_back(0);
							}
							if (C_L[n]>0)
								c.v_Lor[n].push_back(dT/C_L[n]/m_Op->GetVI(n,pos[0],pos[1],pos[2]));
							else
								c.v_Lor[n].push_back(0);
							if (L_L[n]>0)
								c.i_Lor[n].push_back(dT/L_L[n]/m_Op->GetIV(n,pos[0],pos[1],pos[2]));
							else
								c.i_Lor[n].push_back(0);
						}
					}
				}
			}
		});
		// into the arrays, thread after thread, each thread's lists freed once copied: no second copy of all cells
		size_t count = 0;
		for (size_t t=0; t<found.size(); ++t)
		{
			const Cells& c = found[t];
			count += c.pos[0].size();
			m_volt_ADE_On[order] = m_volt_ADE_On[order] || c.volt;
			m_curr_ADE_On[order] = m_curr_ADE_On[order] || c.curr;
			m_volt_Lor_ADE_On[order] = m_volt_Lor_ADE_On[order] || c.volt_Lor;
			m_curr_Lor_ADE_On[order] = m_curr_Lor_ADE_On[order] || c.curr_Lor;
		}
		m_LM_Count.push_back(count);

		m_LM_pos[order] = new unsigned int*[3];
		v_int_ADE[order] = m_volt_ADE_On[order] ? new FDTD_FLOAT*[3] : NULL;
		v_ext_ADE[order] = m_volt_ADE_On[order] ? new FDTD_FLOAT*[3] : NULL;
		i_int_ADE[order] = m_curr_ADE_On[order] ? new FDTD_FLOAT*[3] : NULL;
		i_ext_ADE[order] = m_curr_ADE_On[order] ? new FDTD_FLOAT*[3] : NULL;
		v_Lor_ADE[order] = m_volt_Lor_ADE_On[order] ? new FDTD_FLOAT*[3] : NULL;
		i_Lor_ADE[order] = m_curr_Lor_ADE_On[order] ? new FDTD_FLOAT*[3] : NULL;
		for (int n=0; n<3; ++n)
		{
			m_LM_pos[order][n] = new unsigned int[count];
			if (m_volt_ADE_On[order])
			{
				v_int_ADE[order][n] = new FDTD_FLOAT[count];
				v_ext_ADE[order][n] = new FDTD_FLOAT[count];
			}
			if (m_curr_ADE_On[order])
			{
				i_int_ADE[order][n] = new FDTD_FLOAT[count];
				i_ext_ADE[order][n] = new FDTD_FLOAT[count];
			}
			if (m_volt_Lor_ADE_On[order])
				v_Lor_ADE[order][n] = new FDTD_FLOAT[count];
			if (m_curr_Lor_ADE_On[order])
				i_Lor_ADE[order][n] = new FDTD_FLOAT[count];
		}
		size_t offset = 0;
		for (size_t t=0; t<found.size(); ++t)
		{
			Cells c;
			std::swap(c, found[t]);   // freed at the end of this iteration
			const size_t num = c.pos[0].size();
			for (int n=0; n<3; ++n)
			{
				std::copy(c.pos[n].begin(), c.pos[n].end(), m_LM_pos[order][n]+offset);
				if (m_volt_ADE_On[order])
				{
					std::copy(c.v_int[n].begin(), c.v_int[n].end(), v_int_ADE[order][n]+offset);
					std::copy(c.v_ext[n].begin(), c.v_ext[n].end(), v_ext_ADE[order][n]+offset);
				}
				if (m_curr_ADE_On[order])
				{
					std::copy(c.i_int[n].begin(), c.i_int[n].end(), i_int_ADE[order][n]+offset);
					std::copy(c.i_ext[n].begin(), c.i_ext[n].end(), i_ext_ADE[order][n]+offset);
				}
				if (m_volt_Lor_ADE_On[order])
					std::copy(c.v_Lor[n].begin(), c.v_Lor[n].end(), v_Lor_ADE[order][n]+offset);
				if (m_curr_Lor_ADE_On[order])
					std::copy(c.i_Lor[n].begin(), c.i_Lor[n].end(), i_Lor_ADE[order][n]+offset);
			}
			offset += num;
		}
	}

	return true;
}

Engine_Extension* Operator_Ext_LorentzMaterial::CreateEngineExtention()
{
	Engine_Ext_LorentzMaterial* eng_ext_lor = new Engine_Ext_LorentzMaterial(this);
	return eng_ext_lor;
}

void Operator_Ext_LorentzMaterial::ShowStat(std::ostream &ostr)  const
{
	Operator_Extension::ShowStat(ostr);
	std::string On_Off[2] = {"Off", "On"};
	ostr << " Max. Dispersion Order N = " << m_Order << endl;
	for (int i=0;i<m_Order;++i)
	{
		ostr << " N=" << i << ":\t Active cells\t\t: " << 	m_LM_Count.at(i) << endl;
		ostr << " N=" << i << ":\t Voltage ADE is \t: " << On_Off[m_volt_ADE_On[i]] << endl;
		ostr << " N=" << i << ":\t Voltage Lor-ADE is \t: " << On_Off[m_volt_Lor_ADE_On[i]] << endl;
		ostr << " N=" << i << ":\t Current ADE is \t: " << On_Off[m_curr_ADE_On[i]] << endl;
		ostr << " N=" << i << ":\t Current Lor-ADE is \t: " << On_Off[m_curr_Lor_ADE_On[i]] << endl;
	}
}
