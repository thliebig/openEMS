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

#include "operator_base.h"

using namespace std;

Operator_Base::Operator_Base()
{
	Init();
	m_MeshType = CARTESIAN;
	m_StoreMaterial[0]=false;
	m_StoreMaterial[1]=false;
	m_StoreMaterial[2]=false;
	m_StoreMaterial[3]=false;
}

Operator_Base::~Operator_Base()
{
	Delete();
}

bool Operator_Base::SetGeometryCSX(ContinuousStructure* geo)
{
	if (geo==NULL) return false;
	CSX = geo;

	return true;
}

std::string Operator_Base::GetDirName(int ny) const
{
	if (ny==0) return "x";
	if (ny==1) return "y";
	if (ny==2) return "z";
	return "";
}

void Operator_Base::Init()
{
	CSX = NULL;

	dT = 0;
	for (int n=0; n<3; ++n)
		discLines[n]=NULL;
	for (int n=0; n<6; ++n)
		m_BC[n]=0;

	SetBackgroundMaterial(1,1,0,0);
}

void Operator_Base::Delete()
{
	for (int n=0; n<3; ++n)
	{
		delete[] discLines[n];
		discLines[n]=0;
	}
	for (int n=0; n<6; ++n)
		m_BC[n]=0;
	dT = 0;
}

void Operator_Base::Reset()
{
	Delete();
}

void Operator_Base::SetMaterialStoreFlags(int type, bool val)
{
	if ((type<0) || (type>4))
		return;
	m_StoreMaterial[type]=val;
}

bool Operator_Base::GetCellCenterMaterialAvgCoord(const unsigned int pos[3], double coord[3]) const
{
	int ipos[3] = {(int)pos[0], (int)pos[1], (int)pos[2]};
	return GetCellCenterMaterialAvgCoord(ipos, coord);
}

void Operator_Base::SetBackgroundMaterial(double epsR, double mueR, double kappa, double sigma, double density)
{
	SetBackgroundEpsR(epsR);
	SetBackgroundMueR(mueR);
	SetBackgroundKappa(kappa);
	SetBackgroundSigma(sigma);
	SetBackgroundDensity(density);
}

void Operator_Base::SetBackgroundEpsR(double val)
{
	if (val<1)
	{
		cerr << __func__ << ": Warning, a relative electric permittivity <1 it not supported, skipping" << endl;
		return;
	}
	m_BG_epsR=val;
}

void Operator_Base::SetBackgroundMueR(double val)
{
	if (val<1)
	{
		cerr << __func__ << ": Warning, a relative magnetic permeability <1 it not supported, skipping" << endl;
		return;
	}
	m_BG_mueR=val;
}

void Operator_Base::SetBackgroundKappa(double val)
{
	if (val<0)
	{
		cerr << __func__ << ": Warning, an electric conductivity <0 it not supported, skipping" << endl;
		return;
	}
	m_BG_kappa=val;
}

void Operator_Base::SetBackgroundSigma(double val)
{
	if (val<0)
	{
		cerr << __func__ << ": Warning, an artifival magnetic conductivity <0 it not supported, skipping" << endl;
		return;
	}
	m_BG_sigma=val;
}


void Operator_Base::SetBackgroundDensity(double val)
{
	if (val<0)
	{
		cerr << __func__ << ": Warning, a mass density <0 it not supported, skipping" << endl;
		return;
	}
	m_BG_density=val;
}
