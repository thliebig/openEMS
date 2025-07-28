/*
*	Copyright (C) 2023 Gadi Lahav (gadi@rfwithcare.com)
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

#include "../operator.h"
#include "tools/array_ops.h"
#include "tools/constants.h"
//#include "cond_sheet_parameter.h"
#include "tools/AdrOp.h"

#include "operator_ext_lumpedRLC.h"
#include "engine_ext_lumpedRLC.h"

#include "CSPrimBox.h"
#include "CSProperties.h"
#include "CSPropLumpedElement.h"

#define COPY_V2A(V,A) std::copy(V.begin(),V.end(),A)

Operator_Ext_LumpedRLC::Operator_Ext_LumpedRLC(Operator* op) : Operator_Extension(op)
{
	// Parallel circuit coefficients
	v_RLC_ilv = NULL;
	v_RLC_i2v = NULL;

	// Series circuit coefficients
	v_RLC_vv2	= NULL;	// Coefficient for [n-2] time of Vd update in Vd equation
	v_RLC_vj1	= NULL;	// Coefficient for [n-1] time of J update in Vd equation
	v_RLC_vj2	= NULL;	// Coefficient for [n-2] time of J update in Vd equation
	v_RLC_vvd	= NULL;	// Coefficient to multiply all Vd in the Vd update equation
	v_RLC_ib0	= NULL;	// Inverse of beta_0
	v_RLC_b1	= NULL;	// beta_1
	v_RLC_b2	= NULL;	// beta_2

	// Additional containers
	v_RLC_dir = NULL;
	v_RLC_pos = NULL;

	RLC_count = 0;
}

Operator_Ext_LumpedRLC::Operator_Ext_LumpedRLC(Operator* op, Operator_Ext_LumpedRLC* op_ext) : Operator_Extension(op,op_ext)
{
	// Parallel circuit coefficients
	v_RLC_ilv = NULL;
	v_RLC_i2v = NULL;

	// Series circuit coefficients
	v_RLC_vv2	= NULL;	// Coefficient for [n-2] time of Vd update in Vd equation
	v_RLC_vj1	= NULL;	// Coefficient for [n-1] time of J update in Vd equation
	v_RLC_vj2	= NULL;	// Coefficient for [n-2] time of J update in Vd equation
	v_RLC_vvd	= NULL;	// Coefficient to multiply all Vd in the Vd update equation
	v_RLC_ib0	= NULL;	// Inverse of beta_0
	v_RLC_b1	= NULL;	// beta_1
	v_RLC_b2	= NULL;	// beta_2

	// Additional containers
	v_RLC_dir = NULL;
	v_RLC_pos = NULL;

	RLC_count = 0;
}

Operator_Ext_LumpedRLC::~Operator_Ext_LumpedRLC()
{
	if (this->RLC_count)
	{
		// Parallel circuit coefficients
		delete[] v_RLC_ilv;
		delete[] v_RLC_i2v;

		// Series circuit coefficients
		delete[] v_RLC_vv2;
		delete[] v_RLC_vj1;
		delete[] v_RLC_vj2;
		delete[] v_RLC_vvd;
		delete[] v_RLC_ib0;
		delete[] v_RLC_b1;
		delete[] v_RLC_b2;

		// Additional containers
		delete[] v_RLC_dir;

		for (unsigned int dIdx = 0 ; dIdx < 3 ; dIdx++)
			delete[] v_RLC_pos[dIdx];

		delete[] v_RLC_pos;

	}
}

Operator_Extension* Operator_Ext_LumpedRLC::Clone(Operator* op)
{
	if (dynamic_cast<Operator_Ext_LumpedRLC*>(this)==NULL)
		return NULL;
	return new Operator_Ext_LumpedRLC(op, this);
}

bool Operator_Ext_LumpedRLC::BuildExtension()
{
	double 	dT = m_Op->GetTimestep();

	double	fMax = 	m_Op->GetExcitationSignal()->GetCenterFreq()
					+
					m_Op->GetExcitationSignal()->GetCutOffFreq();

	unsigned int 	pos[] = {0,0,0};

	vector<CSProperties*> cs_props;

	int 			dir;
	CSPropLumpedElement::LEtype lumpedType;

	vector<unsigned int> 	v_pos[3];

	vector<int>		v_dir;

	vector<double>  v_ilv;
	vector<double>	v_i2v;

	vector<double>	v_vv2;
	vector<double>	v_vj1;
	vector<double>	v_vj2;
	vector<double>	v_vvd;
	vector<double>	v_ib0;
	vector<double>	v_b1;
	vector<double>	v_b2;

	// Lumped RLC parameters
	double R, L, C;

	// clear all vectors to initialize them
	for (unsigned int dIdx = 0 ; dIdx < 3 ; dIdx++)
		v_pos[dIdx].clear();

	v_dir.clear();

	v_ilv.clear();
	v_i2v.clear();

	v_vv2.clear();
	v_vj1.clear();
	v_vj2.clear();
	v_vvd.clear();
	v_ib0.clear();
	v_b1.clear();
	v_b2.clear();

	// Obtain from CSX (continuous structure) all the lumped RLC properties
	// Properties are material properties, not the objects themselves
	cs_props = m_Op->CSX->GetPropertyByType(CSProperties::LUMPED_ELEMENT);

	// Iterate through various properties. In theory, there should be a property set per-
	// primitive, as each lumped element should have it's own unique properties.
	for(size_t n = 0 ; n < cs_props.size() ; ++n)
	{
		// Cast current property to lumped RLC property continuous structure properties
		CSPropLumpedElement* cs_RLC_props = dynamic_cast<CSPropLumpedElement*>(cs_props.at(n));
		if (cs_RLC_props==NULL)
			return false; //sanity check: this should never happen!

		// Store direction and type
		dir = cs_RLC_props->GetDirection();
		lumpedType = cs_RLC_props->GetLEtype();

		//		if (lumpedType == LEtype::INVALID
		if (lumpedType == CSPropLumpedElement::INVALID)
		{
			cerr << "Operator_Ext_LumpedRLC::BuildExtension(): Warning: RLCtype is invalid! considering as parallel. "
					<< " ID: " << cs_RLC_props->GetID() << " @ Property: " << cs_RLC_props->GetName() << endl;
			lumpedType = CSPropLumpedElement::PARALLEL;
		}

		// Extract R, L and C from property class
		C = cs_RLC_props->GetCapacity();
		if (C < 0.0)
		{
			cerr << "Operator_Ext_LumpedRLC::BuildExtension(): Warning:";
			if (C < 0.0)
				cerr << " Value of C is smaller than zero, automatically set to 0.";
			if (lumpedType == CSPropLumpedElement::SERIES)
				cerr << " In a series RLC, C = 0 will be considered as Short-Circuit.";

			cerr << " ID: " << cs_RLC_props->GetID() << " @ Property: " << cs_RLC_props->GetName() << endl;
			C = 0.0;
		}
		R = cs_RLC_props->GetResistance();
		if (R < 0.0)
		{
			cerr << "Operator_Ext_LumpedRLC::BuildExtension(): Warning:";
			if (R < 0.0)
				cerr << "Value of R is smaller than zero, automatically set to 0.";
			if (lumpedType == CSPropLumpedElement::PARALLEL)
				cerr << " In a parallel RLC, R = 0 will be considered as Open-Circuit.";
			cerr << " ID: " << cs_RLC_props->GetID() << " @ Property: " << cs_RLC_props->GetName() << endl;
			R = 0;

		}
		L = cs_RLC_props->GetInductance();
		if (L < 0.0)
		{
			cerr << "Operator_Ext_LumpedRLC::BuildExtension(): Warning:";
			if (L < 0.0)
				cerr << " Value of L is smaller than zero, automatically set to zero.";
			if (lumpedType == CSPropLumpedElement::PARALLEL)
				cerr << " In a parallel RLC, L = 0 will be considered as Open-Circuit.";
			cerr << " ID: " << cs_RLC_props->GetID() << " @ Property: " << cs_RLC_props->GetName() << endl;
			L = 0;
		}

		// Validate they are non NaNs
		if (std::isnan(C)) C = 0.0;
		if (std::isnan(L)) L = 0.0;
		if (std::isnan(R)) R = 0.0;

		// Check that this is a lumped RLC. If this is a regular parallel LC, this is handled by the regular engine.
		if (!(this->IsLElumpedRLC(cs_RLC_props)))
			continue;

		if ((dir < 0) || (dir > 2))
		{
			cerr << "Operator_Ext_LumpedRLC::BuildExtension(): Warning: Lumped Element direction is invalid! skipping. "
					<< " ID: " << cs_RLC_props->GetID() << " @ Property: " << cs_RLC_props->GetName() << endl;
			continue;
		}

		// Initialize other two direction containers
		int dir_p1 = (dir + 1) % 3;
		int dir_p2 = (dir + 2) % 3;

		// Now iterate through primitive(s). I still think there should be only one per-
		// material definition, but maybe I'm wrong...
		vector<CSPrimitives*> cs_RLC_prims = cs_RLC_props->GetAllPrimitives();

		for (size_t boxIdx = 0 ; boxIdx < cs_RLC_prims.size() ; ++boxIdx)
		{
			CSPrimBox* cBox = dynamic_cast<CSPrimBox*>(cs_RLC_prims.at(boxIdx));

			if (cBox)
			{

				// Get box start and stop positions
				unsigned int 	uiStart[3],
								uiStop[3];


				// snap to the native coordinate system
				int Snap_Dimension =
						m_Op->SnapBox2Mesh(
								cBox->GetStartCoord()->GetCoords(m_Op->m_MeshType),	// Start Coord
								cBox->GetStopCoord()->GetCoords(m_Op->m_MeshType),	// Stop Coord
								uiStart,	// Start Index
								uiStop,		// Stop Index
								false,		// Dual (doublet) Grid?
								true);		// Full mesh?

				// Verify that snapped dimension is correct
				if (Snap_Dimension<=0)
				{
					if (Snap_Dimension>=-1)
						cerr << "Operator_Ext_LumpedRLC::BuildExtension(): Warning: Lumped RLC snapping failed! Dimension is: " << Snap_Dimension << " skipping. "
								<< " ID: " << cs_RLC_prims.at(boxIdx)->GetID() << " @ Property: " << cs_RLC_props->GetName() << endl;
					// Snap_Dimension == -2 means outside the simulation domain --> no special warning, but box probably marked as unused!
					continue;
				}

				// Verify that in the direction of the current propagation, the size isn't zero.
				if (uiStart[dir]==uiStop[dir])
				{
					cerr << "Operator_Ext_LumpedRLC::BuildExtension(): Warning: Lumped RLC with zero (snapped) length is invalid! skipping. "
							<< " ID: " << cs_RLC_prims.at(boxIdx)->GetID() << " @ Property: " << cs_RLC_props->GetName() << endl;
					continue;
				}

				// Calculate number of cells per-direction
				unsigned int 	Ncells_0 = uiStop[dir] - uiStart[dir],
								Ncells_1 = uiStop[dir_p1] - uiStart[dir_p1] + 1,
								Ncells_2 = uiStop[dir_p2] - uiStart[dir_p2] + 1;

				// All cells in directions 1 and 2 are considered parallel connection
				unsigned int Npar = Ncells_1*Ncells_2;

				// Separate elements such that individual elements can be calculated.
				double	dL = L*Npar/Ncells_0,
						dR = R*Npar/Ncells_0,
						dG = (1/R)*Ncells_0/Npar,
						dC = C*Ncells_0/Npar;

				// Elements for series RLC
				double	ib0	= 2.0*dT*dC/(4.0*dL*dC + 2.0*dT*dR*dC + dT*dT),
						b1	= (dT*dT - 4.0*dL*dC)/(dT*dC),
						b2	= (4.0*dL*dC - 2.0*dT*dR*dC + dT*dT)/(2.0*dT*dC);

				// In a series RLC setup, re-set the coefficients in case they are zero
				if ((lumpedType == CSPropLumpedElement::SERIES)  && (C == 0))
				{
					ib0 = dT/(2.0*dL + dT*dR);
					b1 = -4.0*dL/dT;
					b2 = (2.0*dL - dT*dR)/dT;
				}

				// Special case: If this is a parallel resonant circuit, and there is no
				// parallel resistor, use zero conductivity. May be risky when low-loss
				// simulations are involved
				if (lumpedType == CSPropLumpedElement::PARALLEL)
					if (R == 0.0)
						dG = 0.0;

				int iPos = 0;

				double Zmin,Zcd_min;

				// In the various positions, update the capacitors and "inverse" resistors
				for (pos[dir] = uiStart[dir] ; pos[dir] < uiStop[dir] ; ++pos[dir])
				{
					for (pos[dir_p1] = uiStart[dir_p1] ; pos[dir_p1] <= uiStop[dir_p1] ; ++pos[dir_p1])
					{
						for (pos[dir_p2] = uiStart[dir_p2] ; pos[dir_p2] <= uiStop[dir_p2] ; ++pos[dir_p2])
						{
							iPos = m_Op->MainOp->SetPos(pos[0],pos[1],pos[2]);


							// Separate to two different cases. Parallel and series
							switch (lumpedType)
							{
								case CSPropLumpedElement::PARALLEL:
									// If a capacitor was set, use it. Otherwise, use the existing node capacitor
									if (dC > 0)
										m_Op->EC_C[dir][iPos] = dC;
									else
										// This case takes the "natural" capacitor into account.
										dC = m_Op->EC_C[dir][iPos];

									v_i2v.push_back((dT/dC)/(1.0 + dT*dG/(2.0*dC)));

									// Update conductivity
									if (R >= 0)
										m_Op->EC_G[dir][iPos] = dG;

									// Update coefficients with respect to the parallel inductance
									if (L > 0)
										v_ilv.push_back(dT/dL);
									else
										v_ilv.push_back(0.0);

									// Take into account the case that the "natural" capacitor is too small
									// with respect to the inductor or the resistor, and add a warning.
									if (dC == 0)
									{
										double Cd = m_Op->EC_C[dir][iPos];
										Zmin = max(dR,2*PI*fMax*dL);
										Zcd_min = 1.0/(2.0*PI*fMax*Cd);

										// Check if the "parasitic" capcitance is not small enough
										if (Zcd_min < LUMPED_RLC_Z_FACT*Zmin)
										{
											Cd = 1.0/(2*PI*fMax*Zmin*LUMPED_RLC_Z_FACT);
											m_Op->EC_C[dir][iPos] = Cd;
										}
									}

									v_vv2.push_back(0.0);
									v_vj1.push_back(0.0);
									v_vj2.push_back(0.0);
									v_vvd.push_back(1.0);
									v_ib0.push_back(0.0);
									v_b1.push_back(0.0);
									v_b2.push_back(0.0);

									// Update with discrete component values of
									m_Op->Calc_ECOperatorPos(dir,pos);

									v_dir.push_back(dir);

									break;

								case CSPropLumpedElement::SERIES:
									m_Op->EC_G[dir][iPos] = 0.0;

									// is a series inductor, modeled separately.
									FDTD_FLOAT Cd = m_Op->EC_C[dir][iPos];

									// Calculate minimum impedance, at maximum frequency
									if (dC)
										Zmin = sqrt(pow(dR,2) + pow(2*PI*fMax*dL - 1.0/(dC*2*PI*fMax),2));
									else
										Zmin = sqrt(pow(dR,2) + pow(2*PI*fMax*dL,2));

									Zcd_min = 1.0/(2.0*PI*fMax*Cd);

									// Check if the "parasitic" capcitance is not small enough
									if (Zcd_min < LUMPED_RLC_Z_FACT*Zmin)
									{
										Cd = 1.0/(2*PI*fMax*Zmin*LUMPED_RLC_Z_FACT);
										m_Op->EC_C[dir][iPos] = Cd;
									}

									// No contribution from parallel inductor
									v_ilv.push_back(0.0);
									v_i2v.push_back(0.0);

									// Contributions from series resistor and inductor
									v_vv2.push_back(0.5*dT*ib0/Cd);
									v_vj1.push_back(0.5*dT*(b1*ib0 - 1.0)/Cd);
									v_vj2.push_back(0.5*dT*b2*ib0/Cd);
									v_vvd.push_back(1.0/(1.0 + 0.5*dT*ib0/Cd));
									v_ib0.push_back(ib0);
									v_b1.push_back(b1);
									v_b2.push_back(b2);

									m_Op->Calc_ECOperatorPos(dir,pos);

									v_dir.push_back(dir);

									break;
							}

							// Store position and direction
							for (unsigned int dIdx = 0 ; dIdx < 3 ; ++dIdx)
								v_pos[dIdx].push_back(pos[dIdx]);

						}
					}
				}

				// Build metallic caps
				if (cs_RLC_props->GetCaps())
					for (pos[dir_p1] = uiStart[dir_p1] ; pos[dir_p1] <= uiStop[dir_p1] ; ++pos[dir_p1])
					{
						for (pos[dir_p2] = uiStart[dir_p2] ; pos[dir_p2] <= uiStop[dir_p2] ; ++pos[dir_p2])
						{
							pos[dir]=uiStart[dir];
							if (pos[dir_p1]<uiStop[dir_p1])
							{
								m_Op->SetVV(dir_p1,pos[0],pos[1],pos[2], 0 );
								m_Op->SetVI(dir_p1,pos[0],pos[1],pos[2], 0 );
								++(m_Op->m_Nr_PEC[dir_p1]);
							}

							if (pos[dir_p2]<uiStop[dir_p2])
							{
								m_Op->SetVV(dir_p2,pos[0],pos[1],pos[2], 0 );
								m_Op->SetVI(dir_p2,pos[0],pos[1],pos[2], 0 );
								++(m_Op->m_Nr_PEC[dir_p2]);
							}

							pos[dir]=uiStop[dir];
							if (pos[dir_p1]<uiStop[dir_p1])
							{
								m_Op->SetVV(dir_p1,pos[0],pos[1],pos[2], 0 );
								m_Op->SetVI(dir_p1,pos[0],pos[1],pos[2], 0 );
								++(m_Op->m_Nr_PEC[dir_p1]);
							}

							if (pos[dir_p2]<uiStop[dir_p2])
							{
								m_Op->SetVV(dir_p2,pos[0],pos[1],pos[2], 0 );
								m_Op->SetVI(dir_p2,pos[0],pos[1],pos[2], 0 );
								++(m_Op->m_Nr_PEC[dir_p2]);
							}
						}
					}


				// Mark as used
				cBox->SetPrimitiveUsed(true);
			}
		}

	}

	// Start data storage
	RLC_count = v_dir.size();

	// values
	if (RLC_count)
	{
		// Allocate space to all variables
		v_RLC_dir 	= new int[RLC_count];

		// Parallel circuit coefficients
		v_RLC_ilv 	= new FDTD_FLOAT[RLC_count];
		v_RLC_i2v 	= new FDTD_FLOAT[RLC_count];

		// Series circuit coefficients
		v_RLC_vv2 = new FDTD_FLOAT[RLC_count];
		v_RLC_vj1 = new FDTD_FLOAT[RLC_count];
		v_RLC_vj2 = new FDTD_FLOAT[RLC_count];
		v_RLC_vvd = new FDTD_FLOAT[RLC_count];
		v_RLC_ib0 = new FDTD_FLOAT[RLC_count];
		v_RLC_b1 = new FDTD_FLOAT[RLC_count];
		v_RLC_b2 = new FDTD_FLOAT[RLC_count];

		v_RLC_pos = new unsigned int*[3];
		for (unsigned int dIdx = 0 ; dIdx < 3 ; ++dIdx)
			v_RLC_pos[dIdx] = new unsigned int[RLC_count];

		// Copy all vectors to arrays
		COPY_V2A(v_dir, v_RLC_dir);

		COPY_V2A(v_ilv, v_RLC_ilv);
		COPY_V2A(v_i2v, v_RLC_i2v);

		COPY_V2A(v_vv2,v_RLC_vv2);
		COPY_V2A(v_vj1,v_RLC_vj1);
		COPY_V2A(v_vj2,v_RLC_vj2);
		COPY_V2A(v_vvd,v_RLC_vvd);
		COPY_V2A(v_ib0,v_RLC_ib0);
		COPY_V2A(v_b1,v_RLC_b1);
		COPY_V2A(v_b2,v_RLC_b2);

		for (unsigned int dIdx = 0 ; dIdx < 3 ; ++dIdx)
			COPY_V2A(v_pos[dIdx],v_RLC_pos[dIdx]);
	}

	return true;
}

Engine_Extension* Operator_Ext_LumpedRLC::CreateEngineExtention()
{
	Engine_Ext_LumpedRLC* eng_ext_RLC = new Engine_Ext_LumpedRLC(this);
	return eng_ext_RLC;
}

void Operator_Ext_LumpedRLC::ShowStat(ostream &ostr)  const
{
	Operator_Extension::ShowStat(ostr);
	string On_Off[2] = {"Off", "On"};

	ostr << "Active cells\t\t: " << RLC_count << endl;
}

bool Operator_Ext_LumpedRLC::IsLElumpedRLC(const CSPropLumpedElement* const p_prop)
{
	CSPropLumpedElement::LEtype lumpedType = p_prop->GetLEtype();

	double L = p_prop->GetInductance();

	bool isParallelRLC = (lumpedType == CSPropLumpedElement::PARALLEL) && (L > 0.0);
	bool isSeriesRLC = lumpedType == CSPropLumpedElement::SERIES;

	// This needs to be something that isn't a parallel RC circuit to add data to this extension.
	return isParallelRLC || isSeriesRLC;
}

