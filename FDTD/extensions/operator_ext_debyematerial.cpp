/*
*	Copyright (C) 2026 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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
#include <atomic>
#include <utility>
#include "operator_ext_debyematerial.h"
#include "engine_ext_debyematerial.h"
#include "operator_ext_cylinder.h"
#include "../operator_cylinder.h"

#include "tools/constants.h"
#include "CSPropDebyeMaterial.h"

#include <vector>

using std::cerr;
using std::endl;

//! A Debye pole needs its relaxation time resolved by the timestep. The
//! trapezoidal update stays stable up to dT/tau of about 0.35, independent of
//! d_eps, so refuse a pole below three timesteps per relaxation time.
#define DEBYE_MIN_TAU_PER_DT 3.0

Operator_Ext_DebyeMaterial::Operator_Ext_DebyeMaterial(Operator* op) : Operator_Ext_Dispersive(op)
{
	m_PoleCount = 0;
	v_relax_ADE = NULL;
	v_drive_ADE = NULL;
	v_solve_ADE = NULL;
}

Operator_Ext_DebyeMaterial::Operator_Ext_DebyeMaterial(Operator* op, Operator_Ext_DebyeMaterial* op_ext) : Operator_Ext_Dispersive(op,op_ext)
{
	m_PoleCount = 0;
	v_relax_ADE = NULL;
	v_drive_ADE = NULL;
	v_solve_ADE = NULL;
}

void Operator_Ext_DebyeMaterial::DeleteArrays()
{
	for (int o=0; o<m_PoleCount; ++o)
	{
		if (v_relax_ADE!=NULL)
		{
			for (int n=0; n<3; ++n)
			{
				delete[] v_relax_ADE[o][n];
				delete[] v_drive_ADE[o][n];
			}
			delete[] v_relax_ADE[o];
			delete[] v_drive_ADE[o];
		}
	}
	delete[] v_relax_ADE;
	delete[] v_drive_ADE;
	v_relax_ADE = NULL;
	v_drive_ADE = NULL;

	if (v_solve_ADE!=NULL)
		for (int n=0; n<3; ++n)
			delete[] v_solve_ADE[n];
	delete[] v_solve_ADE;
	v_solve_ADE = NULL;
}

Operator_Ext_DebyeMaterial::~Operator_Ext_DebyeMaterial()
{
	DeleteArrays();
	m_PoleCount = 0;
}

Operator_Extension* Operator_Ext_DebyeMaterial::Clone(Operator* op)
{
	return new Operator_Ext_DebyeMaterial(op, this);
}

bool Operator_Ext_DebyeMaterial::BuildExtension()
{
	double dT = m_Op->GetTimestep();
	unsigned int numLines[3] = {m_Op->GetNumberOfLines(0,true),m_Op->GetNumberOfLines(1,true),m_Op->GetNumberOfLines(2,true)};

	m_PoleCount = 0;
	std::vector<CSProperties*> props = m_Op->GetGeometryCSX()->GetPropertyByType(CSProperties::DEBYEMATERIAL);
	for (size_t n=0; n<props.size(); ++n)
	{
		CSPropDebyeMaterial* mat = dynamic_cast<CSPropDebyeMaterial*>(props.at(n));
		if (mat==NULL)
			return false; //sanity check, this should not happen
		if (mat->GetDispersionOrder()>m_PoleCount)
			m_PoleCount = mat->GetDispersionOrder();
	}
	if (m_PoleCount<1)
		return false; //no Debye material, drop this extension

	//! the cells a thread found: one entry per cell, m_PoleCount coefficient sets per direction, FDTD_FLOAT as the arrays (the same conversion, earlier)
	struct Cells
	{
		std::vector<unsigned int> pos[3];
		std::vector<FDTD_FLOAT> solve[3];
		std::vector<std::vector<FDTD_FLOAT> > relax[3], drive[3];
	};
	std::atomic<bool> warn_once(true);

	// the columns (x, y) on the build threads, their cells joined in the order of the columns;
	// a column without a primitive of a Debye material in its bounding box has none
	// every thread's lists sized, also of those without columns (fewer columns than threads)
	std::vector<Cells> found(m_Op->GetBuildThreads());
	for (Cells& c : found)
		for (int n=0; n<3; ++n)
		{
			c.relax[n].resize(m_PoleCount);
			c.drive[n].resize(m_PoleCount);
		}
	m_Op->ParallelFor((size_t)numLines[0]*numLines[1], [&](unsigned int thread, size_t begin, size_t end) {
		Cells& c = found[thread];
		unsigned int pos[3];
		double coord[3];
		// per cell scratch
		std::vector<std::vector<double> > relax(3), drive(3);
		for (int n=0; n<3; ++n)
		{
			relax[n].resize(m_PoleCount);
			drive[n].resize(m_PoleCount);
		}
		double S[3];

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
				any = vPrims[p]->GetProperty() && vPrims[p]->GetProperty()->ToDebyeMaterial();
			if (!any)
				continue;

			for (pos[2]=0; pos[2]<numLines[2]; ++pos[2])
			{
				bool b_pos_on = false;

				for (int n=0; n<3; ++n)
				{
					S[n] = 0.0;
					for (int o=0; o<m_PoleCount; ++o)
					{
						relax[n][o] = 0.0;
						drive[n][o] = 0.0;
					}

					if (m_Op->GetYeeCoords(n,pos,coord,false)==false)
						continue;
					if (m_CC_R0_included && (n==2) && (pos[0]==0))
						coord[1] = m_Op->GetDiscLine(1,0);

					// the cell's own voltage update coefficient: it carries the cell
					// losses and, at r==0 in cylindrical coords, the special case
					double vi;
					if (m_CC_R0_included && (n==2) && (pos[0]==0))
						vi = m_Op_Cyl->m_Cyl_Ext->vi_R0[pos[2]];
					else
						vi = m_Op->GetVI(n,pos[0],pos[1],pos[2]);
					if (vi==0)
						continue;

					CSProperties* prop = m_Op->GetGeometryCSX()->GetPropertyByCoordPriority(coord, vPrims, true);
					if (prop==NULL)
						continue;
					CSPropDebyeMaterial* mat = prop->ToDebyeMaterial();
					if (mat==NULL)
						continue;

					double A_l = m_Op->GetEdgeArea(n,pos)/m_Op->GetEdgeLength(n,pos);
					int order = mat->GetDispersionOrder();
					if (order>m_PoleCount)
						order = m_PoleCount;

					for (int o=0; o<order; ++o)
					{
						double d_eps   = mat->GetEpsDeltaWeighted(o,n,coord);
						double t_relax = mat->GetEpsRelaxTimeWeighted(o,n,coord);
						if ((d_eps<=0) || (t_relax<=0))
							continue;
						if (t_relax < DEBYE_MIN_TAU_PER_DT*dT)
						{
							if (warn_once.exchange(false))
							{
								cerr << "Operator_Ext_DebyeMaterial::BuildExtension(): Warning, "
								     << "relaxation time (" << t_relax << "s) of material \""
								     << mat->GetName() << "\" needs at least "
								     << DEBYE_MIN_TAU_PER_DT << " timesteps (dT=" << dT
								     << "s), skipping this pole..." << endl;
							}
							continue;
						}
						double C_L = EPS0 * d_eps * A_l;
						relax[n][o] = (2.0*t_relax - dT)/(2.0*t_relax + dT);
						// drive = r*c2, with r = C_L*vi/dT the branch to cell
						// capacitance ratio and c2 = dT/(2 tau + dT)
						drive[n][o] = C_L*vi/(2.0*t_relax + dT);
						S[n] += drive[n][o];
						b_pos_on = true;
					}
				}

				if (b_pos_on==false)
					continue;

				for (int n=0; n<3; ++n)
				{
					c.pos[n].push_back(pos[n]);
					c.solve[n].push_back(1.0/(1.0 + S[n]));
					for (int o=0; o<m_PoleCount; ++o)
					{
						c.relax[n][o].push_back(relax[n][o]);
						c.drive[n][o].push_back(drive[n][o]);
					}
				}
			}
		}
	});
	size_t count = 0;
	for (size_t t=0; t<found.size(); ++t)
		count += found[t].pos[0].size();
	if (count==0)
		return false; //no active cell, drop this extension

	// The base applies one correction per slot; the implicit solve couples all
	// poles into a single joint correction, so there is exactly one slot and the
	// cell list is shared, not replicated per pole.
	m_Order = 1;
	m_LM_Count.push_back(count);
	m_volt_ADE_On = new bool[1];
	m_curr_ADE_On = new bool[1];
	m_volt_ADE_On[0] = true;
	m_curr_ADE_On[0] = false;

	m_LM_pos = new unsigned int**[1];
	m_LM_pos[0] = new unsigned int*[3];

	v_solve_ADE = new FDTD_FLOAT*[3];
	v_relax_ADE  = new FDTD_FLOAT**[m_PoleCount];
	v_drive_ADE   = new FDTD_FLOAT**[m_PoleCount];
	for (int o=0; o<m_PoleCount; ++o)
	{
		v_relax_ADE[o] = new FDTD_FLOAT*[3];
		v_drive_ADE[o]  = new FDTD_FLOAT*[3];
	}
	for (int n=0; n<3; ++n)
	{
		m_LM_pos[0][n] = new unsigned int[count];
		v_solve_ADE[n] = new FDTD_FLOAT[count];
		for (int o=0; o<m_PoleCount; ++o)
		{
			v_relax_ADE[o][n] = new FDTD_FLOAT[count];
			v_drive_ADE[o][n] = new FDTD_FLOAT[count];
		}
	}

	// into the arrays, thread after thread, each thread's lists freed once copied: no second copy of all cells
	size_t offset = 0;
	for (size_t t=0; t<found.size(); ++t)
	{
		Cells c;
		std::swap(c, found[t]);   // freed at the end of this iteration
		for (int n=0; n<3; ++n)
		{
			std::copy(c.pos[n].begin(), c.pos[n].end(), m_LM_pos[0][n]+offset);
			std::copy(c.solve[n].begin(), c.solve[n].end(), v_solve_ADE[n]+offset);
			for (int o=0; o<m_PoleCount; ++o)
			{
				std::copy(c.relax[n][o].begin(), c.relax[n][o].end(), v_relax_ADE[o][n]+offset);
				std::copy(c.drive[n][o].begin(), c.drive[n][o].end(), v_drive_ADE[o][n]+offset);
			}
		}
		offset += c.pos[0].size();
	}

	return true;
}

Engine_Extension* Operator_Ext_DebyeMaterial::CreateEngineExtention()
{
	return new Engine_Ext_DebyeMaterial(this);
}

void Operator_Ext_DebyeMaterial::ShowStat(std::ostream &ostr)  const
{
	Operator_Extension::ShowStat(ostr);
	ostr << " Debye poles N       \t: " << m_PoleCount << endl;
	ostr << " Active cells        \t: " << m_LM_Count.at(0) << endl;
}
