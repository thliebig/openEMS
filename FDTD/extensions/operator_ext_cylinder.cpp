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

#include "operator_ext_cylinder.h"
#include "FDTD/operator_cylinder.h"
#include "engine_ext_cylinder.h"

Operator_Ext_Cylinder::Operator_Ext_Cylinder(Operator_Cylinder* op) : Operator_Extension(op)
{
	m_Op_Cyl = op;

	CC_R0_included=m_Op_Cyl->GetR0Included();
	CC_closedAlpha=m_Op_Cyl->GetClosedAlpha();

	vv_R0 = NULL;
	vi_R0 = NULL;
}

Operator_Ext_Cylinder::~Operator_Ext_Cylinder()
{
	delete[] vv_R0;
	vv_R0=NULL;
	delete[] vi_R0;
	vi_R0=NULL;
}

bool Operator_Ext_Cylinder::BuildExtension()
{
	delete[] vv_R0;
	vv_R0=NULL;
	delete[] vi_R0;
	vi_R0=NULL;

	//if r=0 is not included -> obviously no special treatment for r=0
	//if alpha direction is not closed, PEC-BC at r=0 necessary and already set...
	if (CC_R0_included==false)
		return true;

	vv_R0 = new FDTD_FLOAT[m_Op->GetNumberOfLines(2,true)];
	vi_R0 = new FDTD_FLOAT[m_Op->GetNumberOfLines(2,true)];

	unsigned int pos[3];
	double coord[3];
	double inEC[4];
	double dT = m_Op->GetTimestep();
	pos[0]=0;
	vector<CSPrimitives*> vPrims_metal = m_Op->GetPrimitivesBoundBox(pos[0], -1, -1, (CSProperties::PropertyType)(CSProperties::MATERIAL | CSProperties::METAL));
	for (pos[2]=0; pos[2]<m_Op->GetNumberOfLines(2,true); ++pos[2])
	{
		double C=0;
		double G=0;
		vector<CSPrimitives*> vPrims_mat   = m_Op->GetPrimitivesBoundBox(pos[0], -1, pos[2], CSProperties::MATERIAL);
		for (pos[1]=0; pos[1]<m_Op->GetNumberOfLines(1,true)-2; ++pos[1])
		{
			m_Op_Cyl->Calc_ECPos(2,pos,inEC,vPrims_mat);
			C+=inEC[0];
			G+=inEC[1];
		}
		m_Op->SetVV(2,0,0,pos[2], 1);
		vv_R0[pos[2]] = (1-dT*G/2/C)/(1+dT*G/2/C);
		vi_R0[pos[2]] = (dT/C)/(1+dT*G/2/C);

		for (unsigned int i=0; i<m_Op->GetNumberOfLines(1,true); ++i)
		{
			m_Op->EC_C[2][m_Op->MainOp->SetPos(0,i,pos[2])] = C;
			m_Op->EC_G[2][m_Op->MainOp->SetPos(0,i,pos[2])] = G;
		}

		//search for metal on z-axis
		m_Op_Cyl->GetYeeCoords(2,pos,coord,false);
		CSProperties* prop = m_Op->CSX->GetPropertyByCoordPriority(coord, vPrims_metal, true);
		if (prop)
		{
			if (prop->GetType()==CSProperties::METAL) //set to PEC
			{
				m_Op->SetVV(2,0,0,pos[2], 0);
				vv_R0[pos[2]] = 0;
				vi_R0[pos[2]] = 0;
				m_Op->EC_C[2][m_Op->MainOp->SetPos(0,0,pos[2])] = 0;
				m_Op->EC_G[2][m_Op->MainOp->SetPos(0,0,pos[2])] = 0;
			}
		}
	}
	return true;
}

Engine_Extension* Operator_Ext_Cylinder::CreateEngineExtention()
{
	Engine_Ext_Cylinder* eng_ext = new Engine_Ext_Cylinder(this);
	return eng_ext;
}


void Operator_Ext_Cylinder::ShowStat(ostream &ostr)  const
{
	Operator_Extension::ShowStat(ostr);
	string On_Off[2] = {"Off", "On"};
	ostr << " Zeroth Radius\t\t: "   << On_Off[CC_R0_included] << endl;
	ostr << " Closed Rotation\t: " << On_Off[CC_closedAlpha] << endl;
}
