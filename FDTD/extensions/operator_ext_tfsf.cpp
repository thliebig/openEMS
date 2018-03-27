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

#include "operator_ext_tfsf.h"
#include "engine_ext_tfsf.h"
#include <cmath>

#include "CSPrimBox.h"
#include "CSPropExcitation.h"

Operator_Ext_TFSF::Operator_Ext_TFSF(Operator* op) : Operator_Extension(op)
{
	Init();
}

Operator_Ext_TFSF::~Operator_Ext_TFSF()
{
	Reset();
}

void Operator_Ext_TFSF::Init()
{
	for (int n=0;n<3;++n)
		for (int l=0;l<2;++l)
			for (int c=0;c<2;++c)
			{
				m_VoltDelay[n][l][c]=NULL;
				m_VoltDelayDelta[n][l][c]=NULL;
				m_VoltAmp[n][l][c]=NULL;
				m_CurrDelay[n][l][c]=NULL;
				m_CurrDelayDelta[n][l][c]=NULL;
				m_CurrAmp[n][l][c]=NULL;
			}

	m_Frequency = 0.0;
	m_PhVel = __C0__;
	Operator_Extension::Init();
}

void Operator_Ext_TFSF::Reset()
{
	for (int n=0;n<3;++n)
		for (int l=0;l<2;++l)
			for (int c=0;c<2;++c)
			{
				delete[] m_VoltDelay[n][l][c];
				m_VoltDelay[n][l][c]=NULL;
				delete[] m_VoltDelayDelta[n][l][c];
				m_VoltDelayDelta[n][l][c]=NULL;
				delete[] m_VoltAmp[n][l][c];
				m_VoltAmp[n][l][c]=NULL;
				delete[] m_CurrDelay[n][l][c];
				m_CurrDelay[n][l][c]=NULL;
				delete[] m_CurrDelayDelta[n][l][c];
				m_CurrDelayDelta[n][l][c]=NULL;
				delete[] m_CurrAmp[n][l][c];
				m_CurrAmp[n][l][c]=NULL;
			}
	Operator_Extension::Reset();
}

Operator_Extension* Operator_Ext_TFSF::Clone(Operator* op)
{
	UNUSED(op);
	return NULL;
}

bool Operator_Ext_TFSF::BuildExtension()
{
	m_Exc = m_Op->GetExcitationSignal();
	double dT = m_Op->GetTimestep();
	if (dT==0)
		return false;
	if (m_Exc==0)
		return false;

	Reset();
	ContinuousStructure* CSX = m_Op->GetGeometryCSX();

	vector<CSProperties*> vec_prop = CSX->GetPropertyByType(CSProperties::EXCITATION);

	if (vec_prop.size()==0)
	{
		cerr << "Operator_Ext_TFSF::BuildExtension: Warning, no excitation properties found" << endl;
		SetActive(false);
		return false;
	}

	double ref_index = sqrt(m_Op->GetBackgroundEpsR()*m_Op->GetBackgroundMueR());
	m_PhVel = __C0__/ref_index;
	CSPropExcitation* elec=NULL;
	CSProperties* prop=NULL;
	CSPrimitives* prim=NULL;
	CSPrimBox* box=NULL;
	for (size_t p=0; p<vec_prop.size(); ++p)
	{
		prop = vec_prop.at(p);
		elec = prop->ToExcitation();
		if (elec->GetExcitType()!=10)
			continue;
		if (prop->GetQtyPrimitives()!=1)
		{
			cerr << "Operator_Ext_TFSF::BuildExtension: Warning, plane wave excitation found with more (or less) than one primitive, skipping..." << endl;
			continue;
		}
		prim = prop->GetPrimitive(0);
		if (prim->GetType()!=CSPrimitives::BOX)
		{
			cerr << "Operator_Ext_TFSF::BuildExtension: Warning, plane wave excitation found with false non-Box primitive, skipping..." << endl;
			continue;
		}
		box = prim->ToBox();
		if (box==NULL) //sanity check, should not happen!
		{
			SetActive(false);
			return false;
		}

		// found a plane-wave excite with exactly one box
		for (int n=0;n<3;++n)
			m_PropDir[n] = elec->GetPropagationDir(n);
		double dir_norm = sqrt(m_PropDir[0]*m_PropDir[0]+m_PropDir[1]*m_PropDir[1]+m_PropDir[2]*m_PropDir[2]);
		if (dir_norm==0)
		{
			cerr << "Operator_Ext_TFSF::BuildExtension: Warning, plane wave direction is zero, skipping..." << endl;
			SetActive(false);
			return false;
		}

		//make it a unit vector
		m_PropDir[0]/=dir_norm;m_PropDir[1]/=dir_norm;m_PropDir[2]/=dir_norm;

		if (m_Op->SnapBox2Mesh(box->GetStartCoord()->GetNativeCoords(), box->GetStopCoord()->GetNativeCoords(), m_Start, m_Stop)!=3)
		{
			cerr << "Operator_Ext_TFSF::BuildExtension: Warning, plane wave box dimension is invalid, skipping..." << endl;
			SetActive(false);
			return false;
		}

		m_Frequency = elec->GetFrequency();
//		if (m_Frequency<=0)
//			m_Frequency = m_Op->GetExcitationSignal()->GetFrequencyOfInterest();
		if (m_Frequency<=0)
			m_PhVel=__C0__/ref_index;
		else
			m_PhVel=m_Op->CalcNumericPhaseVelocity(m_Start,m_Stop,m_PropDir,m_Frequency);

		if ((m_PhVel<0) || (m_PhVel>__C0__/ref_index) || std::isnan(m_PhVel))
		{
			cerr << "Operator_Ext_TFSF::BuildExtension: Warning, invalid phase velocity found, resetting to c0! " << endl;
			m_PhVel = __C0__/ref_index;
		}

		double origin[3];
		unsigned int ui_origin[3];
		for (int n=0;n<3;++n)
		{
			m_E_Amp[n] = elec->GetExcitation(n);
			m_numLines[n] = m_Stop[n]-m_Start[n]+1;
			m_IncLow[n] = m_PropDir[n]>=0;

			if (m_Start[n]==0)
				m_ActiveDir[n][0]=false;
			else
				m_ActiveDir[n][0]=true;
			if (m_Stop[n]==m_Op->GetNumberOfLines(n,true)-1)
				m_ActiveDir[n][1]=false;
			else
				m_ActiveDir[n][1]=true;

			if (m_IncLow[n])
			{
				ui_origin[n] = m_Start[n]-1;
			}
			else
			{
				ui_origin[n] = m_Stop[n]+1;
			}
			origin[n] = m_Op->GetDiscLine(n,ui_origin[n]);
		}

		double dotEk = (m_E_Amp[0]*m_PropDir[0] + m_E_Amp[1]*m_PropDir[1] + m_E_Amp[2]*m_PropDir[2]);
		double angle = acos( dotEk / (m_E_Amp[0]*m_E_Amp[0] + m_E_Amp[1]*m_E_Amp[1] + m_E_Amp[2]*m_E_Amp[2]) ) / PI * 180;

		if (angle==0)
		{
			cerr << "Operator_Ext_TFSF::BuildExtension: Warning, plane wave direction and polarization is identical, skipping..." << endl;
			SetActive(false);
			return false;
		}
		if (angle!=90)
		{
			cerr << "Operator_Ext_TFSF::BuildExtension: Warning, angle between propagation direction and polarization is not 90deg, resetting E-polarization to : (";
			for (int n=0;n<3;++n)
				m_E_Amp[n]-=m_PropDir[n]*dotEk;
			cerr << m_E_Amp[0] << "," << m_E_Amp[1] << "," << m_E_Amp[2] << ")" << endl;
		}

		int nP,nPP;
		for (int n=0;n<3;++n)
		{
			nP = (n+1)%3;
			nPP = (n+2)%3;
			m_H_Amp[n] = m_PropDir[nP]*m_E_Amp[nPP] - m_PropDir[nPP]*m_E_Amp[nP];
			m_H_Amp[n] /= __Z0__*sqrt(m_Op->GetBackgroundMueR()/m_Op->GetBackgroundEpsR());
		}

		double coord[3];
		double unit = m_Op->GetGridDelta();
		double delay;
		double dist;
		unsigned int pos[3];
		unsigned int numP, ui_pos;
		m_maxDelay = 0;
		for (int n=0;n<3;++n)
		{
			nP = (n+1)%3;
			nPP = (n+2)%3;
			pos[n] = 0; //unused
			pos[nP] = m_Start[nP];

			numP = m_numLines[nP]*m_numLines[nPP];

			if (!m_ActiveDir[n][0] && !m_ActiveDir[n][1])
				continue;

			for (int l=0;l<2;++l)
				for (int c=0;c<2;++c)
				{
					if (m_ActiveDir[n][l])
					{
						m_VoltDelay[n][l][c]=new unsigned int[numP];
						m_VoltDelayDelta[n][l][c]=new FDTD_FLOAT[numP];
						m_VoltAmp[n][l][c]=new FDTD_FLOAT[numP];

						m_CurrDelay[n][l][c]=new unsigned int[numP];
						m_CurrDelayDelta[n][l][c]=new FDTD_FLOAT[numP];
						m_CurrAmp[n][l][c]=new FDTD_FLOAT[numP];
					}
				}

			ui_pos = 0;
			for (unsigned int i=0;i<m_numLines[nP];++i)
			{
				pos[nPP] = m_Start[nPP];
				for (unsigned int j=0;j<m_numLines[nPP];++j)
				{
					// current updates
					pos[n] = m_Start[n];

					if (m_ActiveDir[n][0])
					{
						m_Op->GetYeeCoords(nP,pos,coord,false);
						dist = fabs((coord[0]-origin[0])*m_PropDir[0])+fabs((coord[1]-origin[1])*m_PropDir[1])+fabs((coord[2]-origin[2])*m_PropDir[2]);
						delay = dist*unit/m_PhVel/dT;
						m_maxDelay = max((unsigned int)delay,m_maxDelay);
						m_CurrDelay[n][0][1][ui_pos] = floor(delay);
						m_CurrDelayDelta[n][0][1][ui_pos] = delay - floor(delay);
						m_CurrAmp[n][0][1][ui_pos] = m_E_Amp[nP]*m_Op->GetEdgeLength(nP,pos);

						m_Op->GetYeeCoords(nPP,pos,coord,false);
						dist = fabs((coord[0]-origin[0])*m_PropDir[0])+fabs((coord[1]-origin[1])*m_PropDir[1])+fabs((coord[2]-origin[2])*m_PropDir[2]);
						delay = dist*unit/m_PhVel/dT;
						m_maxDelay = max((unsigned int)delay,m_maxDelay);
						m_CurrDelay[n][0][0][ui_pos] = floor(delay);
						m_CurrDelayDelta[n][0][0][ui_pos] = delay - floor(delay);
						m_CurrAmp[n][0][0][ui_pos] = m_E_Amp[nPP]*m_Op->GetEdgeLength(nPP,pos);

						--pos[n];
						m_CurrAmp[n][0][0][ui_pos]*=m_Op->GetIV(nP,pos);
						m_CurrAmp[n][0][1][ui_pos]*=m_Op->GetIV(nPP,pos);
					}

					if (m_ActiveDir[n][1])
					{
						pos[n] = m_Stop[n];
						m_Op->GetYeeCoords(nP,pos,coord,false);
						dist = fabs((coord[0]-origin[0])*m_PropDir[0])+fabs((coord[1]-origin[1])*m_PropDir[1])+fabs((coord[2]-origin[2])*m_PropDir[2]);
						delay = dist*unit/m_PhVel/dT;
						m_maxDelay = max((unsigned int)delay,m_maxDelay);
						m_CurrDelay[n][1][1][ui_pos] = floor(delay);
						m_CurrDelayDelta[n][1][1][ui_pos] = delay - floor(delay);
						m_CurrAmp[n][1][1][ui_pos] = m_E_Amp[nP]*m_Op->GetEdgeLength(nP,pos);

						m_Op->GetYeeCoords(nPP,pos,coord,false);
						dist = fabs((coord[0]-origin[0])*m_PropDir[0])+fabs((coord[1]-origin[1])*m_PropDir[1])+fabs((coord[2]-origin[2])*m_PropDir[2]);
						delay = dist*unit/m_PhVel/dT;
						m_maxDelay = max((unsigned int)delay,m_maxDelay);
						m_CurrDelay[n][1][0][ui_pos] = floor(delay);
						m_CurrDelayDelta[n][1][0][ui_pos] = delay - floor(delay);
						m_CurrAmp[n][1][0][ui_pos] = m_E_Amp[nPP]*m_Op->GetEdgeLength(nPP,pos);

						m_CurrAmp[n][1][0][ui_pos]*=m_Op->GetIV(nP,pos);
						m_CurrAmp[n][1][1][ui_pos]*=m_Op->GetIV(nPP,pos);
					}

					if (m_ActiveDir[n][0])
						m_CurrAmp[n][0][0][ui_pos]*=-1;
					if (m_ActiveDir[n][1])
						m_CurrAmp[n][1][1][ui_pos]*=-1;

					if (pos[nP]==m_Stop[nP])
					{
						if (m_ActiveDir[n][0])
							m_CurrAmp[n][0][1][ui_pos]=0;
						if (m_ActiveDir[n][1])
							m_CurrAmp[n][1][1][ui_pos]=0;
					}
					if (pos[nPP]==m_Stop[nPP])
					{
						if (m_ActiveDir[n][0])
							m_CurrAmp[n][0][0][ui_pos]=0;
						if (m_ActiveDir[n][1])
							m_CurrAmp[n][1][0][ui_pos]=0;
					}

					// voltage updates
					pos[n] = m_Start[n]-1;
					if (m_ActiveDir[n][0])
					{
						m_Op->GetYeeCoords(nP,pos,coord,true);
						dist = fabs((coord[0]-origin[0])*m_PropDir[0])+fabs((coord[1]-origin[1])*m_PropDir[1])+fabs((coord[2]-origin[2])*m_PropDir[2]);
						delay = dist*unit/m_PhVel/dT + 1.0;
						m_maxDelay = max((unsigned int)delay,m_maxDelay);
						m_VoltDelay[n][0][1][ui_pos] = floor(delay);
						m_VoltDelayDelta[n][0][1][ui_pos] = delay - floor(delay);
						m_VoltAmp[n][0][1][ui_pos] = m_H_Amp[nP]*m_Op->GetEdgeLength(nP,pos,true);

						m_Op->GetYeeCoords(nPP,pos,coord,true);
						dist = fabs((coord[0]-origin[0])*m_PropDir[0])+fabs((coord[1]-origin[1])*m_PropDir[1])+fabs((coord[2]-origin[2])*m_PropDir[2]);
						delay = dist*unit/m_PhVel/dT + 1.0;
						m_maxDelay = max((unsigned int)delay,m_maxDelay);
						m_VoltDelay[n][0][0][ui_pos] = floor(delay);
						m_VoltDelayDelta[n][0][0][ui_pos] = delay - floor(delay);
						m_VoltAmp[n][0][0][ui_pos] = m_H_Amp[nPP]*m_Op->GetEdgeLength(nPP,pos,true);

						++pos[n];
						m_VoltAmp[n][0][0][ui_pos]*=m_Op->GetVI(nP,pos);
						m_VoltAmp[n][0][1][ui_pos]*=m_Op->GetVI(nPP,pos);
					}

					pos[n] = m_Stop[n];
					if (m_ActiveDir[n][1])
					{
						m_Op->GetYeeCoords(nP,pos,coord,true);
						dist = fabs((coord[0]-origin[0])*m_PropDir[0])+fabs((coord[1]-origin[1])*m_PropDir[1])+fabs((coord[2]-origin[2])*m_PropDir[2]);
						delay = dist*unit/m_PhVel/dT + 1.0;
						m_maxDelay = max((unsigned int)delay,m_maxDelay);
						m_VoltDelay[n][1][1][ui_pos] = floor(delay);
						m_VoltDelayDelta[n][1][1][ui_pos] = delay - floor(delay);
						m_VoltAmp[n][1][1][ui_pos] = m_H_Amp[nP]*m_Op->GetEdgeLength(nP,pos,true);

						m_Op->GetYeeCoords(nPP,pos,coord,true);
						dist = fabs((coord[0]-origin[0])*m_PropDir[0])+fabs((coord[1]-origin[1])*m_PropDir[1])+fabs((coord[2]-origin[2])*m_PropDir[2]);
						delay = dist*unit/m_PhVel/dT + 1.0;
						m_maxDelay = max((unsigned int)delay,m_maxDelay);
						m_VoltDelay[n][1][0][ui_pos] = floor(delay);
						m_VoltDelayDelta[n][1][0][ui_pos] = delay - floor(delay);
						m_VoltAmp[n][1][0][ui_pos] = m_H_Amp[nPP]*m_Op->GetEdgeLength(nPP,pos,true);

						m_VoltAmp[n][1][0][ui_pos]*=m_Op->GetVI(nP,pos);
						m_VoltAmp[n][1][1][ui_pos]*=m_Op->GetVI(nPP,pos);
					}

					if (m_ActiveDir[n][1])
						m_VoltAmp[n][1][0][ui_pos]*=-1;
					if (m_ActiveDir[n][0])
						m_VoltAmp[n][0][1][ui_pos]*=-1;

					if (pos[nP]==m_Stop[nP])
					{
						if (m_ActiveDir[n][0])
							m_VoltAmp[n][0][0][ui_pos]=0;
						if (m_ActiveDir[n][1])
							m_VoltAmp[n][1][0][ui_pos]=0;
					}
					if (pos[nPP]==m_Stop[nPP])
					{
						if (m_ActiveDir[n][0])
							m_VoltAmp[n][0][1][ui_pos]=0;
						if (m_ActiveDir[n][1])
							m_VoltAmp[n][1][1][ui_pos]=0;
					}

					++pos[nPP];
					++ui_pos;
				}
				++pos[nP];
			}
		}
		++m_maxDelay;
		return true;
	}
	SetActive(false);
	return false;
}

Engine_Extension* Operator_Ext_TFSF::CreateEngineExtention()
{
	return new Engine_Ext_TFSF(this);
}

void Operator_Ext_TFSF::ShowStat(ostream &ostr) const
{
	Operator_Extension::ShowStat(ostr);
	cout << "Active directions\t: " << "(" << m_ActiveDir[0][0] << "/" << m_ActiveDir[0][1] << ", " << m_ActiveDir[1][0] << "/" << m_ActiveDir[1][1]  << ", " << m_ActiveDir[2][0] << "/" << m_ActiveDir[2][1]  << ")" << endl;
	cout << "Propagation direction\t: " << "(" << m_PropDir[0] << ", " << m_PropDir[1] << ", " << m_PropDir[2] << ")" << endl;
	cout << "Rel. propagation speed\t: " << m_PhVel/__C0__ << "*c0  @ " << m_Frequency << " Hz" << endl;
	cout << "E-field amplitude (V/m)\t: " << "(" << m_E_Amp[0] << ", " << m_E_Amp[1] << ", " << m_E_Amp[2] << ")" << endl;
	cout << "H-field amplitude (A/m)\t: " << "(" << m_H_Amp[0] << ", " << m_H_Amp[1] << ", " << m_H_Amp[2] << ")" << endl;
	cout << "Box Dimensions\t\t: " << m_numLines[0] << " x " << m_numLines[1] << " x " << m_numLines[2] << endl;
	cout << "Max. Delay (TS)\t\t: " << m_maxDelay << endl;
	int dirs = m_ActiveDir[0][0] + m_ActiveDir[0][1] + m_ActiveDir[1][0] + m_ActiveDir[1][1] + m_ActiveDir[2][0] + m_ActiveDir[2][1] ;
	cout << "Memory usage (est.)\t: ~" << m_numLines[0] * m_numLines[1] * m_numLines[2] * dirs * 4 * 4 / 1024 << " kiB" << endl;
}
