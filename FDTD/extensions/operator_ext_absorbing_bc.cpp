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

#include "operator_ext_absorbing_bc.h"
#include "engine_ext_absorbing_bc.h"

#include "tools/array_ops.h"

#include "CSPrimBox.h"

#define COPY_V2A(V,A) std::copy(V.begin(),V.end(),A)

#define SAFE_DELETE(v) 	\
	if (v != NULL) 		\
	{ 					\
		delete[] v;		\
		v = NULL;		\
	}					\

Operator_Ext_Absorbing_BC::Operator_Ext_Absorbing_BC(Operator* op) : Operator_Extension(op)
{
	Initialize();
}

Operator_Ext_Absorbing_BC::Operator_Ext_Absorbing_BC(Operator* op, uint numThreads) : Operator_Extension(op)
{
	m_Nthreads = numThreads;
	Initialize();
}


Operator_Ext_Absorbing_BC::~Operator_Ext_Absorbing_BC()
{
	for (uint primIdx = 0 ; primIdx < m_numPrims ; primIdx++)
	{
		SAFE_DELETE(m_sheetX0[primIdx]);
		SAFE_DELETE(m_sheetX1[primIdx]);

		SAFE_DELETE(m_dir[primIdx]);

		for (uint dimIdx = 0 ; dimIdx < 3 ; dimIdx++)
		{
			SAFE_DELETE(m_K1[primIdx][dimIdx]);
			SAFE_DELETE(m_K2[primIdx][dimIdx]);
		}

		SAFE_DELETE(m_K1[primIdx]);
		SAFE_DELETE(m_K2[primIdx]);
	}

	// Delete all 2D arrays
	SAFE_DELETE(m_sheetX0);
	SAFE_DELETE(m_sheetX1);
	SAFE_DELETE(m_dir);
	SAFE_DELETE(m_K1);
	SAFE_DELETE(m_K2);

	// Delete all 1D arrays
	SAFE_DELETE(m_numCells);
	SAFE_DELETE(m_normDir);
	SAFE_DELETE(m_boundaryTypes);

}

Operator_Ext_Absorbing_BC::Operator_Ext_Absorbing_BC(Operator* op, Operator_Ext_Absorbing_BC* op_ext) : Operator_Extension(op, op_ext)
{
	Initialize();
}

Operator_Extension* Operator_Ext_Absorbing_BC::Clone(Operator* op)
{
	if (dynamic_cast<Operator_Ext_Absorbing_BC*>(this)==NULL)
		return NULL;
	return new Operator_Ext_Absorbing_BC(op, this);
}

void Operator_Ext_Absorbing_BC::Initialize()
{
	m_sheetX0		= NULL;
	m_sheetX1		= NULL;
	m_dir			= NULL;
	m_normDir		= NULL;
	m_K1			= NULL;
	m_K2			= NULL;
	m_boundaryTypes = NULL;
	m_numCells		= NULL;
	m_numPrims		= 0;

}


bool Operator_Ext_Absorbing_BC::BuildExtension()
{

	double dT	= m_Op->GetTimestep();
	//double fMax	= m_Op->GetExcitationSignal()->GetCenterFreq();

	uint pos[] = {0,0,0};

	// Containers for primitive initialization
	int							normDir;
	double						phaseVelocity;
	CSPropAbsorbingBC::BCtype	boundaryType;

	vector<CSProperties*>		cs_props;

	// Containers for data storage
	vector<uint>			v_numCells;
	vector<uint*>			v_dir;
	vector<int>				v_normDir;
	vector<uint*>			v_sheetX0;
	vector<uint*>			v_sheetX1;
	vector<FDTD_FLOAT**>	v_K1;
	vector<FDTD_FLOAT**>	v_K2;
	vector<int>				v_boundaryTypes;

	// Initialize all vectors
	v_numCells.clear();
	v_dir.clear();
	v_normDir.clear();
	v_sheetX0.clear();
	v_sheetX1.clear();
	v_K1.clear();
	v_K2.clear();
	v_boundaryTypes.clear();

	// Obtain from CSX (continuous structure) all the lumped RLC properties
	// Properties are material properties, not the objects themselves
	cs_props = m_Op->CSX->GetPropertyByType(CSProperties::ABSORBING_BC);

	// Iterate through various properties. In theory, there should be a property set per-
	// primitive, as each "port" should have it's own unique properties.
	for(size_t n = 0 ; n < cs_props.size() ; ++n)
	{
		// Cast current property to absorbing BC sheet (Absorbing Boundary-condition Sheet : ABS) property continuous structure properties
		CSPropAbsorbingBC * cs_ABS_props = dynamic_cast<CSPropAbsorbingBC*>(cs_props.at(n));
		if (cs_ABS_props == NULL)
			return false; //sanity check: this should never happen!

		// Store direction and type
		normDir			= cs_ABS_props->GetNormDir();
		boundaryType	= cs_ABS_props->GetBoundaryType();
		phaseVelocity	= cs_ABS_props->GetPhaseVelocity();

		// Initialize other two direction containers. Keep in mind that 'normDir'
		// has the values 1,2,3 or -1,-2,-3. 0 is considered illegal.
		uint dir[3] = {	 abs(normDir) - 1		,
						 abs(normDir)      % 3	,
						(abs(normDir) + 1) % 3};

		// Now iterate through primitive(s). I still think there should be only one per-
		// material definition, but maybe I'm wrong...
		vector<CSPrimitives*> cs_ABS_prims = cs_ABS_props->GetAllPrimitives();
		for (size_t sheetIdx = 0 ; sheetIdx < cs_ABS_prims.size() ; ++sheetIdx)
		{
			CSPrimBox* cSheet = dynamic_cast<CSPrimBox*>(cs_ABS_prims.at(sheetIdx));

			// If this just so happen to not be a sheet, ignore this
			if (!cSheet)
			{
				cerr << "Operator_Ext_Absorbing_BC::BuildExtension(): Warning: Absorbing sheet validation failed, skipping. "
											<< " ID: " << cs_ABS_prims.at(sheetIdx)->GetID() << " @ Property: " << cs_ABS_props->GetName() << endl;
				continue;
			}


			// Check that this is actually a sheet
			// Get box start and stop positions
			uint 	uiStart[3],
					uiStop[3];

			// snap to the native coordinate system
			int Snap_Dimension =
				m_Op->SnapBox2Mesh(
						cSheet->GetStartCoord()->GetCoords(m_Op->m_MeshType),	// Start Coord
						cSheet->GetStopCoord()->GetCoords(m_Op->m_MeshType),	// Stop Coord
						uiStart,	// Start Index
						uiStop,		// Stop Index
						false,		// Dual (doublet) Grid?
						true);		// Full mesh?

			// Verify that snapped dimension is correct
			if (Snap_Dimension <= 0)
			{
				if (Snap_Dimension>=-1)
					cerr << "Operator_Ext_Absorbing_BC::BuildExtension(): Warning: Absorbing sheet snapping failed! Dimension is: " << Snap_Dimension << " skipping. "
							<< " ID: " << cs_ABS_prims.at(sheetIdx)->GetID() << " @ Property: " << cs_ABS_props->GetName() << endl;
				// Snap_Dimension == -2 means outside the simulation domain --> no special warning, but box probably marked as unused!
				continue;
			}

			uint	Ncells[3] = {0,0,0},
					totCells = 1;
			for (uint dimIdx = 0 ; dimIdx < 3 ; dimIdx++)
			{
				Ncells[dimIdx] = uiStop[dimIdx] - uiStart[dimIdx] + 1;
				totCells *= Ncells[dimIdx];
			}


			// check that one of the dimensions is 1, and only one
			uint sheetCheck = 0;
			for (uint dimIdx = 0 ; dimIdx < 3 ; dimIdx++)
				sheetCheck += Ncells[dimIdx] == 1;
			if (sheetCheck != 1)
			{
				cerr 	<< "Operator_Ext_Absorbing_BC::BuildExtension(): Error: Absorbing sheet is not a sheet! Skipping. "
						<< " ID: " << cs_ABS_prims.at(sheetIdx)->GetID() << " @ Property: " << cs_ABS_props->GetName() << endl;
				continue;
			}

			// Now check that the direction is compatible with the sheet
			if (Ncells[dir[0]] != 1)
			{
				cerr 	<< "Operator_Ext_Absorbing_BC::BuildExtension(): Error: Absorbing sheet is not a sheet! Skipping. "
						<< " ID: " << cs_ABS_prims.at(sheetIdx)->GetID() << " @ Property: " << cs_ABS_props->GetName() << endl;
				continue;
			}

			// The position starts and stops in the same place.
			pos[dir[0]] = uiStart[dir[0]];

			// The position is considered in the middle of the mesh cell.
			double delta = fabs(m_Op->GetEdgeLength(dir[0],pos));


			// Initialize coefficients for this engine extension
			FDTD_FLOAT vt = phaseVelocity*dT;


			uint 		*temp_dir = new uint[3];
			uint		*temp_sheetX0 = new uint[3];
			uint		*temp_sheetX1 = new uint[3];


			for (uint dimIdx = 0 ; dimIdx < 3 ; dimIdx++)
			{
				temp_dir[dimIdx] = dir[dimIdx];
				temp_sheetX0[dimIdx] = uiStart[dimIdx];
				temp_sheetX1[dimIdx] = uiStop[dimIdx];
			}

			if ((boundaryType == CSPropAbsorbingBC::MUR_1ST_1PV) || (boundaryType == CSPropAbsorbingBC::MUR_1ST_1PV_SA))
			{
				FDTD_FLOAT	**temp_K1 = new FDTD_FLOAT*[3];

				for (int dimIdx = 0 ; dimIdx < 3 ; dimIdx++)
				{
					// In this revision, initialize all coefficients
					vector<FDTD_FLOAT> temp_v_K1(totCells,(vt - delta) / (vt + delta));
					temp_K1[dimIdx] = new FDTD_FLOAT[totCells];
					COPY_V2A(temp_v_K1,temp_K1[dimIdx]);
				}

				v_K1.push_back(temp_K1);
			}

			if (boundaryType == CSPropAbsorbingBC::MUR_1ST_1PV_SA)
			{
				FDTD_FLOAT	**temp_K2 = new FDTD_FLOAT*[3];

				for (int dimIdx = 0 ; dimIdx < 3 ; dimIdx++)
				{
					// In this revision, initialize all coefficients
					vector<FDTD_FLOAT> temp_v_K2(totCells,vt/delta);
					temp_K2[dimIdx] = new FDTD_FLOAT[totCells];
					COPY_V2A(temp_v_K2,temp_K2[dimIdx]);
				}

				v_K2.push_back(temp_K2);
			}

			v_numCells.push_back(totCells);
			v_dir.push_back(temp_dir);
			v_normDir.push_back(normDir);
			v_sheetX0.push_back(temp_sheetX0);
			v_sheetX1.push_back(temp_sheetX1);
			v_boundaryTypes.push_back((int)boundaryType);


			// Mark as used
			cSheet->SetPrimitiveUsed(true);
		}
		// Clear before next iteration
		cs_ABS_prims.clear();
	}

	// Count number of primitives that are absorbers
	m_numPrims = v_numCells.size();

	// Copy all vectors to arrays
	m_numCells	= new uint[m_numPrims];
	m_dir		= new uint*[m_numPrims];
	m_sheetX0	= new uint*[m_numPrims];
	m_sheetX1	= new uint*[m_numPrims];
	m_normDir	= new int[m_numPrims];
	m_boundaryTypes = new int[m_numPrims];

	m_K1 = new FDTD_FLOAT**[m_numPrims];

	COPY_V2A(v_numCells,m_numCells);
	COPY_V2A(v_dir,m_dir);
	COPY_V2A(v_sheetX0,m_sheetX0);
	COPY_V2A(v_sheetX1,m_sheetX1);
	COPY_V2A(v_normDir,m_normDir);
	COPY_V2A(v_boundaryTypes,m_boundaryTypes);

	COPY_V2A(v_K1,m_K1);

	if (boundaryType == CSPropAbsorbingBC::MUR_1ST_1PV_SA)
	{
		m_K2 = new FDTD_FLOAT**[m_numPrims];
		COPY_V2A(v_K2,m_K2);
	}



	return true;
}


Engine_Extension* Operator_Ext_Absorbing_BC::CreateEngineExtention()
{
	Engine_Ext_Absorbing_BC* eng_ext = new Engine_Ext_Absorbing_BC(this);
	return eng_ext;
}

void Operator_Ext_Absorbing_BC::ShowStat(ostream &ostr) const
{
	Operator_Extension::ShowStat(ostr);
//	string XYZ[3] = {"x","y","z"};
//	char pmSign = (m_normDir > 0) ? '+' : '-';
//
//
//	ostr << " Active direction\t: " << pmSign << XYZ[abs(m_ny)] << " at line: " << m_LineNr << endl;
//	if (m_v_phase)
//		ostr << " Used phase velocity\t: " << m_v_phase << " (" << m_v_phase/__C0__ << " * c_0)" << endl;

	uint totCells = 0;
	for (uint primIdx ; primIdx < m_numPrims ; primIdx++)
	{
		totCells += m_numCells[primIdx];
	}

	ostr << " Number of absorbing BCs\t: " << m_numPrims << " total cells: " << totCells << endl;

}









