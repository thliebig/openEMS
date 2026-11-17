/*
*	Copyright (C) 2023 Gadi Lahav (gadi@rfwithcare.com)
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

#include "../operator.h"
#include "tools/array_ops.h"
#include "tools/constants.h"
//#include "cond_sheet_parameter.h"

#include "operator_ext_lumpedRLC.h"
#include "engine_ext_lumpedRLC.h"

#include "CSPrimBox.h"
#include "CSProperties.h"
#include "CSPropLumpedElement.h"

#define COPY_V2A(V,A) std::copy(V.begin(),V.end(),A)

using std::cerr;
using std::endl;

Operator_Ext_LumpedRLC::Operator_Ext_LumpedRLC(Operator* op) : Operator_Extension(op)
{
	// Parallel circuit coefficients
	v_RLC_ilv = NULL;
	v_RLC_i2v = NULL;

	// Series circuit coefficients
	v_RLC_vvd	= NULL;
	v_RLC_vv1	= NULL;
	v_RLC_vv2	= NULL;
	v_RLC_vj1	= NULL;
	v_RLC_vj2	= NULL;
	v_RLC_jv0	= NULL;
	v_RLC_jv1	= NULL;
	v_RLC_jv2	= NULL;
	v_RLC_jj1	= NULL;
	v_RLC_jj2	= NULL;

	// Additional containers
	v_RLC_dir = NULL;
	v_RLC_pos = NULL;

	RLC_count = 0;
}

Operator_Ext_LumpedRLC::Operator_Ext_LumpedRLC(Operator* op, Operator_Ext_LumpedRLC* op_ext) : Operator_Extension(op,op_ext)
{
	RLC_count = 0;
	v_RLC_ilv = NULL; v_RLC_i2v = NULL;
	v_RLC_vvd = NULL; v_RLC_vv1 = NULL; v_RLC_vv2 = NULL; v_RLC_vj1 = NULL; v_RLC_vj2 = NULL;
	v_RLC_jv0 = NULL; v_RLC_jv1 = NULL; v_RLC_jv2 = NULL; v_RLC_jj1 = NULL; v_RLC_jj2 = NULL;
	v_RLC_dir = NULL;
	v_RLC_pos = NULL;
}

Operator_Ext_LumpedRLC::~Operator_Ext_LumpedRLC()
{
	if (this->RLC_count)
	{
		// Parallel circuit coefficients
		delete[] v_RLC_ilv;
		delete[] v_RLC_i2v;

		// Series circuit coefficients
		delete[] v_RLC_vvd;
		delete[] v_RLC_vv1;
		delete[] v_RLC_vv2;
		delete[] v_RLC_vj1;
		delete[] v_RLC_vj2;
		delete[] v_RLC_jv0;
		delete[] v_RLC_jv1;
		delete[] v_RLC_jv2;
		delete[] v_RLC_jj1;
		delete[] v_RLC_jj2;

		// Additional containers
		delete[] v_RLC_dir;

		for (unsigned int dIdx = 0 ; dIdx < 3 ; dIdx++)
			delete[] v_RLC_pos[dIdx];

		delete[] v_RLC_pos;

	}
}

Operator_Extension* Operator_Ext_LumpedRLC::Clone(Operator* op)
{
	return new Operator_Ext_LumpedRLC(op, this);
}

bool Operator_Ext_LumpedRLC::BuildExtension()
{
	double 	dT = m_Op->GetTimestep();

	double	fMax = 	m_Op->GetExcitationSignal()->GetMaxFreq();

	unsigned int 	pos[] = {0,0,0};

	std::vector<CSProperties*> cs_props;

	int 			dir;
	CSPropLumpedElement::LEtype lumpedType;

	std::vector<unsigned int> 	v_pos[3];

	std::vector<int>		v_dir;

	std::vector<double>  v_ilv;
	std::vector<double>	v_i2v;

	std::vector<double>	v_vvd;
	std::vector<double>	v_vv1;
	std::vector<double>	v_vv2;
	std::vector<double>	v_vj1;
	std::vector<double>	v_vj2;
	std::vector<double>	v_jv0;
	std::vector<double>	v_jv1;
	std::vector<double>	v_jv2;
	std::vector<double>	v_jj1;
	std::vector<double>	v_jj2;

	// Lumped RLC parameters
	double R, L, C;

	// clear all vectors to initialize them
	for (unsigned int dIdx = 0 ; dIdx < 3 ; dIdx++)
		v_pos[dIdx].clear();

	v_dir.clear();

	v_ilv.clear();
	v_i2v.clear();

	v_vvd.clear();
	v_vv1.clear();
	v_vv2.clear();
	v_vj1.clear();
	v_vj2.clear();
	v_jv0.clear();
	v_jv1.clear();
	v_jv2.clear();
	v_jj1.clear();
	v_jj2.clear();

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

		// Extract R, L and C from property class
		C = cs_RLC_props->GetCapacity();
		R = cs_RLC_props->GetResistance();
		L = cs_RLC_props->GetInductance();

		// NaN means "not present" — silent and valid.  Negative values are clamped to
		// zero with a warning.  Explicit zero on a component that is non-physical in
		// the chosen topology also gets a dedicated warning.
		if (C < 0.0)
		{
			cerr 	<< "Operator_Ext_LumpedRLC::BuildExtension(): Warning: Value of C is smaller than zero, automatically set to 0. ID:"
					<< cs_RLC_props->GetID() << " @ Property: " << cs_RLC_props->GetName() << endl;
			C = 0.0;
		}
		else if (!std::isnan(C) && C == 0.0 && lumpedType == CSPropLumpedElement::SERIES)
			cerr	<< "Operator_Ext_LumpedRLC::BuildExtension(): Warning: C = 0 in a series circuit is non-physical; treating as absent (series RL only). ID: "
					<< cs_RLC_props->GetID() << " @ Property: " << cs_RLC_props->GetName() << endl;

		if (R < 0.0)
		{
			cerr 	<< "Operator_Ext_LumpedRLC::BuildExtension(): Warning: Value of R is smaller than zero, automatically set to 0. ID:"
					<< cs_RLC_props->GetID() << " @ Property: " << cs_RLC_props->GetName() << endl;
			R = 0.0;
		}
		else if (!std::isnan(R) && R == 0.0 && lumpedType == CSPropLumpedElement::PARALLEL)
			cerr	<< "Operator_Ext_LumpedRLC::BuildExtension(): Warning: R = 0 in a parallel circuit would be a short circuit; use AddMetal() for an intentional short. Treating as absent (open). ID: "
					<< cs_RLC_props->GetID() << " @ Property: " << cs_RLC_props->GetName() << endl;

		if (L < 0.0)
		{
			cerr 	<< "Operator_Ext_LumpedRLC::BuildExtension(): Warning: Value of L is smaller than zero, automatically set to zero. ID: "
					<< cs_RLC_props->GetID() << " @ Property: " << cs_RLC_props->GetName() << endl;
			L = 0.0;
		}
		else if (!std::isnan(L) && L == 0.0 && lumpedType == CSPropLumpedElement::PARALLEL)
			cerr	<< "Operator_Ext_LumpedRLC::BuildExtension(): Warning: L = 0 in a parallel circuit would be a short circuit; treating as absent. ID: "
					<< cs_RLC_props->GetID() << " @ Property: " << cs_RLC_props->GetName() << endl;

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

		std::vector<CSPrimitives*> cs_RLC_prims = cs_RLC_props->GetAllPrimitives();

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

				// Scale to per-cell values.  NaN means the component is absent:
				//   series  R/L absent → wire (0)
				//   parallel R absent  → open path (G = 0)
				//   C absent           → use natural grid capacitance (handled below)
				double	dL = std::isnan(L) ? 0.0 : L*Npar/Ncells_0,
						dR = std::isnan(R) ? 0.0 : R*Npar/Ncells_0,
						dG = (std::isnan(R) || R == 0.0) ? 0.0 : (1.0/R)*Ncells_0/Npar,
						dC = std::isnan(C) ? 0.0 : C*Ncells_0/Npar;

				// Series ADE, the current J of the element (trapezoidal rule):
				//   J[n] = jv0*Vd[n] + jv1*Vd[n-1] + jv2*Vd[n-2] + jj1*J[n-1] + jj2*J[n-2]
				// computed only for SERIES to avoid a division by zero when dC == 0.
				double jv0 = 0.0, jv1 = 0.0, jv2 = 0.0, jj1 = 0.0, jj2 = 0.0;
				if (lumpedType == CSPropLumpedElement::SERIES)
				{
					if (dC == 0.0)
					{
						// RL: (2L/dT + R)*J[n] - (2L/dT - R)*J[n-1] = Vd[n] + Vd[n-1]. Not the
						// difference of two timesteps of it (the RLC form below with C -> inf):
						// that has a second root at z = 1, which the float coefficients and
						// their rounding errors turn into a slow growth of the fields, with any
						// R > 0 into an exponential one (#241).
						const double b0 = 2.0*dL/dT + dR;
						jv0 = 1.0/b0;
						jv1 = 1.0/b0;
						jj1 = (2.0*dL/dT - dR)/b0;
					}
					else
					{
						// RLC: b0*J[n] + b1*J[n-1] + b2*J[n-2] = Vd[n] - Vd[n-2]
						const double b0 = 2.0*dL/dT + dR + dT/(2.0*dC);
						const double b1 = dT/dC - 4.0*dL/dT;
						const double b2 = 2.0*dL/dT - dR + dT/(2.0*dC);
						jv0 = 1.0/b0;
						jv2 = -1.0/b0;
						jj1 = -b1/b0;
						jj2 = -b2/b0;
					}
				}


				double Zmin, Zcd_min;

				// In the various positions, update the capacitors and "inverse" resistors
				for (pos[dir] = uiStart[dir] ; pos[dir] < uiStop[dir] ; ++pos[dir])
				{
					for (pos[dir_p1] = uiStart[dir_p1] ; pos[dir_p1] <= uiStop[dir_p1] ; ++pos[dir_p1])
					{
						for (pos[dir_p2] = uiStart[dir_p2] ; pos[dir_p2] <= uiStop[dir_p2] ; ++pos[dir_p2])
						{

							switch (lumpedType)
							{
								case CSPropLumpedElement::PARALLEL:
								{
									// Use the explicit capacitor if set; otherwise fall back to the
									// natural grid capacitance, adjusting upward if it is too small
									// to keep the ADE stable relative to R and L.
									double dC_cell;
									double ec[4] = {0, 0, 0, 0};
									unsigned int mask = 0;
									if (dC > 0)
									{
										ec[Operator::EC_CAP] = dC;
										mask |= 1u<<Operator::EC_CAP;
										dC_cell = dC;
									}
									else
									{
										dC_cell = m_Op->GetEC(Operator::EC_CAP, dir, pos);
										Zmin    = std::max(dR, 2*PI*fMax*dL);
										Zcd_min = 1.0/(2.0*PI*fMax*dC_cell);
										if (Zcd_min < LUMPED_RLC_Z_FACT*Zmin)
										{
											dC_cell = 1.0/(2*PI*fMax*Zmin*LUMPED_RLC_Z_FACT);
											ec[Operator::EC_CAP] = dC_cell;
											mask |= 1u<<Operator::EC_CAP;
										}
									}

									v_i2v.push_back((dT/dC_cell)/(1.0 + dT*dG/(2.0*dC_cell)));

									// Update conductivity
									if (R >= 0)
									{
										ec[Operator::EC_CON] = dG;
										mask |= 1u<<Operator::EC_CON;
									}
									if (mask)
										m_Op->SetEC(dir, pos, mask, ec);

									// Update coefficients with respect to the parallel inductance
									v_ilv.push_back((L > 0) ? dT/dL : 0.0);

									v_vvd.push_back(1.0);
									v_vv1.push_back(0.0);
									v_vv2.push_back(0.0);
									v_vj1.push_back(0.0);
									v_vj2.push_back(0.0);
									v_jv0.push_back(0.0);
									v_jv1.push_back(0.0);
									v_jv2.push_back(0.0);
									v_jj1.push_back(0.0);
									v_jj2.push_back(0.0);

									m_Op->Calc_ECOperatorPos(dir,pos);

									v_dir.push_back(dir);

									break;
								}

								case CSPropLumpedElement::SERIES:
								{
									FDTD_FLOAT Cd = m_Op->GetEC(Operator::EC_CAP, dir, pos);
									unsigned int mask = 1u<<Operator::EC_CON;   // G is 0

									// Check if the natural grid capacitance is too small
									// relative to the series impedance at fMax.
									if (dC)
										Zmin = sqrt(pow(dR,2) + pow(2*PI*fMax*dL - 1.0/(dC*2*PI*fMax),2));
									else
										Zmin = sqrt(pow(dR,2) + pow(2*PI*fMax*dL,2));

									Zcd_min = 1.0/(2.0*PI*fMax*Cd);

									if (Zcd_min < LUMPED_RLC_Z_FACT*Zmin)
									{
										Cd = 1.0/(2*PI*fMax*Zmin*LUMPED_RLC_Z_FACT);
										mask |= 1u<<Operator::EC_CAP;
									}
									const double ec[4] = {Cd, 0, 0, 0};
									m_Op->SetEC(dir, pos, mask, ec);

									// No contribution from parallel inductor or capacitor update
									v_ilv.push_back(0.0);
									v_i2v.push_back(0.0);

									// Contributions from series R, L, C via ADE coefficients: the
									// node takes Cd*(Vd[n]-Vd[n-1])/dT = I - (J[n]+J[n-1])/2, with
									// J[n] from above
									const double k = 0.5*dT/Cd;
									v_vvd.push_back(1.0/(1.0 + k*jv0));
									v_vv1.push_back(-k*jv1);
									v_vv2.push_back(-k*jv2);
									v_vj1.push_back(-k*(jj1 + 1.0));
									v_vj2.push_back(-k*jj2);
									v_jv0.push_back(jv0);
									v_jv1.push_back(jv1);
									v_jv2.push_back(jv2);
									v_jj1.push_back(jj1);
									v_jj2.push_back(jj2);

									m_Op->Calc_ECOperatorPos(dir,pos);

									v_dir.push_back(dir);

									break;
								}
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
		v_RLC_vvd = new FDTD_FLOAT[RLC_count];
		v_RLC_vv1 = new FDTD_FLOAT[RLC_count];
		v_RLC_vv2 = new FDTD_FLOAT[RLC_count];
		v_RLC_vj1 = new FDTD_FLOAT[RLC_count];
		v_RLC_vj2 = new FDTD_FLOAT[RLC_count];
		v_RLC_jv0 = new FDTD_FLOAT[RLC_count];
		v_RLC_jv1 = new FDTD_FLOAT[RLC_count];
		v_RLC_jv2 = new FDTD_FLOAT[RLC_count];
		v_RLC_jj1 = new FDTD_FLOAT[RLC_count];
		v_RLC_jj2 = new FDTD_FLOAT[RLC_count];

		v_RLC_pos = new unsigned int*[3];
		for (unsigned int dIdx = 0 ; dIdx < 3 ; ++dIdx)
			v_RLC_pos[dIdx] = new unsigned int[RLC_count];

		// Copy all vectors to arrays
		COPY_V2A(v_dir, v_RLC_dir);

		COPY_V2A(v_ilv, v_RLC_ilv);
		COPY_V2A(v_i2v, v_RLC_i2v);

		COPY_V2A(v_vvd,v_RLC_vvd);
		COPY_V2A(v_vv1,v_RLC_vv1);
		COPY_V2A(v_vv2,v_RLC_vv2);
		COPY_V2A(v_vj1,v_RLC_vj1);
		COPY_V2A(v_vj2,v_RLC_vj2);
		COPY_V2A(v_jv0,v_RLC_jv0);
		COPY_V2A(v_jv1,v_RLC_jv1);
		COPY_V2A(v_jv2,v_RLC_jv2);
		COPY_V2A(v_jj1,v_RLC_jj1);
		COPY_V2A(v_jj2,v_RLC_jj2);

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

void Operator_Ext_LumpedRLC::ShowStat(std::ostream &ostr)  const
{
	Operator_Extension::ShowStat(ostr);
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

