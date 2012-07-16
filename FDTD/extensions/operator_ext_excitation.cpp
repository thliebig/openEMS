/*
*	Copyright (C) 2011 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#include "operator_ext_excitation.h"
#include "engine_ext_excitation.h"
#include "FDTD/excitation.h"
#include "ContinuousStructure.h"

Operator_Ext_Excitation::Operator_Ext_Excitation(Operator* op, Excitation* exc) : Operator_Extension(op)
{
	m_Exc = exc;
}

Operator_Ext_Excitation::~Operator_Ext_Excitation()
{
}

Operator_Ext_Excitation::Operator_Ext_Excitation(Operator* op, Operator_Ext_Excitation* op_ext) : Operator_Extension(op, op_ext)
{
	m_Exc = NULL;
}

bool Operator_Ext_Excitation::BuildExtension()
{
	double dT = m_Exc->GetTimestep();
	if (dT==0)
		return false;
	if (m_Exc==0)
		return false;

	ContinuousStructure* CSX = m_Op->GetGeometryCSX();

	unsigned int pos[3];
	double amp=0;

	vector<unsigned int> volt_vIndex[3];
	vector<FDTD_FLOAT> volt_vExcit;
	vector<unsigned int> volt_vDelay;
	vector<unsigned int> volt_vDir;
	double volt_coord[3];

	vector<unsigned int> curr_vIndex[3];
	vector<FDTD_FLOAT> curr_vExcit;
	vector<unsigned int> curr_vDelay;
	vector<unsigned int> curr_vDir;
	double curr_coord[3];

	vector<CSProperties*> vec_prop = CSX->GetPropertyByType(CSProperties::EXCITATION);

	if (vec_prop.size()==0)
	{
		cerr << "Operator::CalcFieldExcitation: Warning, no excitation properties found" << endl;
		return false;
	}

	CSPropExcitation* elec=NULL;
	CSProperties* prop=NULL;
	int priority=0;

	unsigned int numLines[] = {m_Op->GetOriginalNumLines(0),m_Op->GetOriginalNumLines(1),m_Op->GetOriginalNumLines(2)};

	for (pos[2]=0; pos[2]<numLines[2]; ++pos[2])
	{
		for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
		{
			for (pos[0]=0; pos[0]<numLines[0]; ++pos[0])
			{
				//electric field excite
				for (int n=0; n<3; ++n)
				{
					if (m_Op->GetYeeCoords(n,pos,volt_coord,false)==false)
						continue;
					if (m_CC_R0_included && (n==2) && (pos[0]==0))
						volt_coord[1] = m_Op->GetDiscLine(1,0);

					if (m_CC_R0_included && (n==1) && (pos[0]==0))
						continue;

					for (size_t p=0; p<vec_prop.size(); ++p)
					{
						prop = vec_prop.at(p);
						elec = prop->ToExcitation();
						if (elec==NULL)
							continue;
						if (prop->CheckCoordInPrimitive(volt_coord,priority,true))
						{
							if ((elec->GetActiveDir(n)) && ( (elec->GetExcitType()==0) || (elec->GetExcitType()==1) ))//&& (pos[n]<numLines[n]-1))
							{
								amp = elec->GetWeightedExcitation(n,volt_coord)*m_Op->GetEdgeLength(n,pos);// delta[n]*gridDelta;
								if (amp!=0)
								{
									volt_vExcit.push_back(amp);
									volt_vDelay.push_back((unsigned int)(elec->GetDelay()/dT));
									volt_vDir.push_back(n);
									volt_vIndex[0].push_back(pos[0]);
									volt_vIndex[1].push_back(pos[1]);
									volt_vIndex[2].push_back(pos[2]);
								}
								if (elec->GetExcitType()==1) //hard excite
								{
									m_Op->SetVV(n,pos[0],pos[1],pos[2], 0 );
									m_Op->SetVI(n,pos[0],pos[1],pos[2], 0 );
								}
							}
						}
					}
				}

				//magnetic field excite
				for (int n=0; n<3; ++n)
				{
					if ((pos[0]>=numLines[0]-1) || (pos[1]>=numLines[1]-1) || (pos[2]>=numLines[2]-1))
						continue;  //skip the last H-Line which is outside the FDTD-domain
					if (m_Op->GetYeeCoords(n,pos,curr_coord,true)==false)
						continue;
					for (size_t p=0; p<vec_prop.size(); ++p)
					{
						prop = vec_prop.at(p);
						elec = prop->ToExcitation();
						if (elec==NULL)
							continue;
						if (prop->CheckCoordInPrimitive(curr_coord,priority,true))
						{
							if ((elec->GetActiveDir(n)) && ( (elec->GetExcitType()==2) || (elec->GetExcitType()==3) ))
							{
								amp = elec->GetWeightedExcitation(n,curr_coord)*m_Op->GetEdgeLength(n,pos,true);// delta[n]*gridDelta;
								if (amp!=0)
								{
									curr_vExcit.push_back(amp);
									curr_vDelay.push_back((unsigned int)(elec->GetDelay()/dT));
									curr_vDir.push_back(n);
									curr_vIndex[0].push_back(pos[0]);
									curr_vIndex[1].push_back(pos[1]);
									curr_vIndex[2].push_back(pos[2]);
								}
								if (elec->GetExcitType()==3) //hard excite
								{
									m_Op->SetII(n,pos[0],pos[1],pos[2], 0 );
									m_Op->SetIV(n,pos[0],pos[1],pos[2], 0 );
								}
							}
						}
					}
				}

			}
		}
	}

	//special treatment for primitives of type curve (treated as wires) see also Calc_PEC
	double p1[3];
	double p2[3];
	struct Operator::Grid_Path path;
	for (size_t p=0; p<vec_prop.size(); ++p)
	{
		prop = vec_prop.at(p);
		elec = prop->ToExcitation();
		for (size_t n=0; n<prop->GetQtyPrimitives(); ++n)
		{
			CSPrimitives* prim = prop->GetPrimitive(n);
			CSPrimCurve* curv = prim->ToCurve();
			if (curv)
			{
				for (size_t i=1; i<curv->GetNumberOfPoints(); ++i)
				{
					curv->GetPoint(i-1,p1);
					curv->GetPoint(i,p2);
					path = m_Op->FindPath(p1,p2);
					if (path.dir.size()>0)
						prim->SetPrimitiveUsed(true);
					for (size_t t=0; t<path.dir.size(); ++t)
					{
						n = path.dir.at(t);
						pos[0] = path.posPath[0].at(t);
						pos[1] = path.posPath[1].at(t);
						pos[2] = path.posPath[2].at(t);
						m_Op->GetYeeCoords(n,pos,volt_coord,false);
						if (elec!=NULL)
						{
							if ((elec->GetActiveDir(n)) && (pos[n]<numLines[n]-1) && ( (elec->GetExcitType()==0) || (elec->GetExcitType()==1) ))
							{
								amp = elec->GetWeightedExcitation(n,volt_coord)*m_Op->GetEdgeLength(n,pos);
								if (amp!=0)
								{
									volt_vExcit.push_back(amp);
									volt_vDelay.push_back((unsigned int)(elec->GetDelay()/dT));
									volt_vDir.push_back(n);
									volt_vIndex[0].push_back(pos[0]);
									volt_vIndex[1].push_back(pos[1]);
									volt_vIndex[2].push_back(pos[2]);
								}
								if (elec->GetExcitType()==1) //hard excite
								{
									m_Op->SetVV(n,pos[0],pos[1],pos[2], 0 );
									m_Op->SetVI(n,pos[0],pos[1],pos[2], 0 );
								}
							}
						}
					}
				}
			}
		}
	}

	// set voltage excitations
	m_Exc->setupVoltageExcitation( volt_vIndex, volt_vExcit, volt_vDelay, volt_vDir );

	// set current excitations
	m_Exc->setupCurrentExcitation( curr_vIndex, curr_vExcit, curr_vDelay, curr_vDir );

	return true;
}

Engine_Extension* Operator_Ext_Excitation::CreateEngineExtention()
{
	return new Engine_Ext_Excitation(this);
}

void Operator_Ext_Excitation::ShowStat(ostream &ostr)  const
{
	Operator_Extension::ShowStat(ostr);
	cout << "Voltage excitations\t: " << m_Exc->Volt_Count    << "\t (" << m_Exc->Volt_Count_Dir[0] << ", " << m_Exc->Volt_Count_Dir[1] << ", " << m_Exc->Volt_Count_Dir[2] << ")" << endl;
	cout << "Current excitations\t: " << m_Exc->Curr_Count << "\t (" << m_Exc->Curr_Count_Dir[0] << ", " << m_Exc->Curr_Count_Dir[1] << ", " << m_Exc->Curr_Count_Dir[2] << ")" << endl;
	cout << "Excitation Length (TS)\t: " << m_Exc->Length << endl;
	cout << "Excitation Length (s)\t: " << m_Exc->Length*m_Op->GetTimestep() << endl;
}

