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

#include "operator_ext_mur_abc.h"
#include "engine_ext_mur_abc.h"

#include "tools/array_ops.h"

#include "CSPropMaterial.h"

Operator_Ext_Mur_ABC::Operator_Ext_Mur_ABC(Operator* op) : Operator_Extension(op)
{
	Initialize();
}

Operator_Ext_Mur_ABC::~Operator_Ext_Mur_ABC()
{
	Delete2DArray(m_Mur_Coeff_nyP,m_numLines);
	m_Mur_Coeff_nyP = NULL;
	Delete2DArray(m_Mur_Coeff_nyPP,m_numLines);
	m_Mur_Coeff_nyPP = NULL;
}

Operator_Ext_Mur_ABC::Operator_Ext_Mur_ABC(Operator* op, Operator_Ext_Mur_ABC* op_ext) : Operator_Extension(op, op_ext)
{
	Initialize();
	m_v_phase = op_ext->m_v_phase;
	SetDirection(op_ext->m_ny,op_ext->m_top);
}

Operator_Extension* Operator_Ext_Mur_ABC::Clone(Operator* op)
{
	if (dynamic_cast<Operator_Ext_Mur_ABC*>(this)==NULL)
		return NULL;
	return new Operator_Ext_Mur_ABC(op, this);
}

bool Operator_Ext_Mur_ABC::IsCylinderCoordsSave(bool closedAlpha, bool R0_included) const
{
	if ((m_ny==0) && (!m_top) && (R0_included || closedAlpha))
		return false;
	if ((m_ny==1) && (closedAlpha))
		return false;
	return true;
}

bool Operator_Ext_Mur_ABC::IsCylindricalMultiGridSave(bool child) const
{
	if (m_ny==2) //always allow in z-direction
		return true;
	if ((m_ny==0) && (m_top) && (!child)) //if top r-direction and is not a child grid allow Mur...
		return true;
	//in all other cases this ABC is not save to use in CylindricalMultiGrid
	return false;
}

void Operator_Ext_Mur_ABC::Initialize()
{
	m_ny = -1;
	m_nyP = -1;
	m_nyPP = -1;
	m_LineNr = 0;
	m_LineNr_Shift = 0;

	m_v_phase = 0.0;

	m_Mur_Coeff_nyP = NULL;
	m_Mur_Coeff_nyPP = NULL;

	m_numLines[0]=0;
	m_numLines[1]=0;
}

void Operator_Ext_Mur_ABC::SetDirection(int ny, bool top_ny)
{
	if ((ny<0) || (ny>2))
		return;

	Delete2DArray(m_Mur_Coeff_nyP,m_numLines);
	Delete2DArray(m_Mur_Coeff_nyPP,m_numLines);

	m_ny = ny;
	m_top = top_ny;
	m_nyP = (ny+1)%3;
	m_nyPP = (ny+2)%3;
	if (!top_ny)
	{
		m_LineNr = 0;
		m_LineNr_Shift = 1;
	}
	else
	{
		m_LineNr = m_Op->GetNumberOfLines(m_ny,true)-1;
		m_LineNr_Shift = m_Op->GetNumberOfLines(m_ny,true) - 2;
	}

	m_numLines[0] = m_Op->GetNumberOfLines(m_nyP,true);
	m_numLines[1] = m_Op->GetNumberOfLines(m_nyPP,true);

	m_Mur_Coeff_nyP = Create2DArray<FDTD_FLOAT>(m_numLines);
	m_Mur_Coeff_nyPP = Create2DArray<FDTD_FLOAT>(m_numLines);

}

bool Operator_Ext_Mur_ABC::BuildExtension()
{
	if (m_ny<0)
	{
		cerr << "Operator_Ext_Mur_ABC::BuildExtension: Warning, Extension not initialized! Use SetDirection!! Abort build!!" << endl;
		return false;
	}
	double dT = m_Op->GetTimestep();
	unsigned int pos[] = {0,0,0};
	pos[m_ny] = m_LineNr;
	double delta = fabs(m_Op->GetEdgeLength(m_ny,pos));
	double coord[] = {0,0,0};
	coord[0] = m_Op->GetDiscLine(0,pos[0]);
	coord[1] = m_Op->GetDiscLine(1,pos[1]);
	coord[2] = m_Op->GetDiscLine(2,pos[2]);

	double eps,mue;
	double c0t;

	if (m_LineNr==0)
		coord[m_ny] = m_Op->GetDiscLine(m_ny,pos[m_ny]) + delta/2 / m_Op->GetGridDelta();
	else
		coord[m_ny] = m_Op->GetDiscLine(m_ny,pos[m_ny]) - delta/2 / m_Op->GetGridDelta();

	int posBB[3];
	posBB[m_ny]  =pos[m_ny];
	posBB[m_nyPP]=-1;

	for (pos[m_nyP]=0; pos[m_nyP]<m_numLines[0]; ++pos[m_nyP])
	{
		posBB[m_nyP]=pos[m_nyP];
		vector<CSPrimitives*> vPrims = m_Op->GetPrimitivesBoundBox(posBB[0], posBB[1], posBB[2], CSProperties::MATERIAL);
		coord[m_nyP] = m_Op->GetDiscLine(m_nyP,pos[m_nyP]);
		for (pos[m_nyPP]=0; pos[m_nyPP]<m_numLines[1]; ++pos[m_nyPP])
		{
			coord[m_nyPP] = m_Op->GetDiscLine(m_nyPP,pos[m_nyPP]);
//			CSProperties* prop = m_Op->GetGeometryCSX()->GetPropertyByCoordPriority(coord, CSProperties::MATERIAL, false);
			CSProperties* prop = m_Op->GetGeometryCSX()->GetPropertyByCoordPriority(coord, vPrims, false);
			if (prop)
			{
				CSPropMaterial* mat = prop->ToMaterial();

				//nP
				eps = mat->GetEpsilonWeighted(m_nyP,coord);
				mue = mat->GetMueWeighted(m_nyP,coord);
				if (m_v_phase>0.0)
					c0t = m_v_phase * dT;
				else
					c0t = __C0__ * dT / sqrt(eps*mue);
				m_Mur_Coeff_nyP[pos[m_nyP]][pos[m_nyPP]] = (c0t - delta) / (c0t + delta);

				//nPP
				eps = mat->GetEpsilonWeighted(m_nyPP,coord);
				mue = mat->GetMueWeighted(m_nyPP,coord);
				if (m_v_phase>0.0)
					c0t = m_v_phase * dT;
				else
					c0t = __C0__ * dT / sqrt(eps*mue);
				m_Mur_Coeff_nyPP[pos[m_nyP]][pos[m_nyPP]] = (c0t - delta) / (c0t + delta);

			}
			else
			{
				if (m_v_phase>0.0)
					c0t = m_v_phase * dT;
				else
					c0t = __C0__ / sqrt(m_Op->GetBackgroundEpsR()*m_Op->GetBackgroundMueR()) * dT;
				m_Mur_Coeff_nyP[pos[m_nyP]][pos[m_nyPP]] = (c0t - delta) / (c0t + delta);
				m_Mur_Coeff_nyPP[pos[m_nyP]][pos[m_nyPP]] = m_Mur_Coeff_nyP[pos[m_nyP]][pos[m_nyPP]];
			}
//			cerr << m_Mur_Coeff_nyP[pos[m_nyP]][pos[m_nyPP]] << " : " << m_Mur_Coeff_nyP[pos[m_nyP]][pos[m_nyPP]] << endl;
		}
	}
//	cerr << "Operator_Ext_Mur_ABC::BuildExtension(): " << m_ny << " @ " << m_LineNr << endl;
	return true;
}

Engine_Extension* Operator_Ext_Mur_ABC::CreateEngineExtention()
{
	Engine_Ext_Mur_ABC* eng_ext = new Engine_Ext_Mur_ABC(this);
	return eng_ext;
}


void Operator_Ext_Mur_ABC::ShowStat(ostream &ostr)  const
{
	Operator_Extension::ShowStat(ostr);
	string XYZ[3] = {"x","y","z"};
	ostr << " Active direction\t: " << XYZ[m_ny] << " at line: " << m_LineNr << endl;
	if (m_v_phase>0.0)
		ostr << " Used phase velocity\t: " << m_v_phase << " (" << m_v_phase/__C0__ << " * c_0)" <<endl;
}
