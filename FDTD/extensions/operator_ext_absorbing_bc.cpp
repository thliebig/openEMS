/*
*	Copyright (C) 2023-2025 Gadi Lahav (gadi@rfwithcare.com), Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

using std::cerr;
using std::endl;

Operator_Ext_Absorbing_BC::Operator_Ext_Absorbing_BC(Operator* op) : Operator_Extension(op)
{
	Initialize();
}

Operator_Ext_Absorbing_BC::~Operator_Ext_Absorbing_BC()
{
}

Operator_Ext_Absorbing_BC::Operator_Ext_Absorbing_BC(Operator* op, Operator_Ext_Absorbing_BC* op_ext) : Operator_Extension(op, op_ext)
{
	Initialize();
}

Operator_Extension* Operator_Ext_Absorbing_BC::Clone(Operator* op)
{
	return new Operator_Ext_Absorbing_BC(op, this);
}

void Operator_Ext_Absorbing_BC::Initialize()
{
	m_numCells = 0;

	m_ny = m_nyP = m_nyPP = -1;

	for (unsigned int dirIdx = 0 ; dirIdx < 3 ; m_sheetX0[dirIdx++] = 0);
	for (unsigned int dirIdx = 0 ; dirIdx < 3 ; m_sheetX1[dirIdx++] = 0);

	m_normalSignPositive = true;

	m_ABCtype = ABCtype::UNDEFINED;

	m_phaseVelocity = 0.0;
}

bool Operator_Ext_Absorbing_BC::SetInitParams(CSPrimitives* prim, CSPropAbsorbingBC* abc_prop)
{
	CSPrimBox* cSheet = dynamic_cast<CSPrimBox*>(prim);
	// If this just so happen to not be a sheet, ignore this
	if (!cSheet)
	{
		cerr << "Operator_Ext_Absorbing_BC::SetInitParams(): Warning: Absorbing sheet validation failed, skipping. "
									<< " ID: " << prim->GetID() << " @ Property: " << abc_prop->GetName() << endl;
		return false;
	}

	// Check that this is actually a sheet
	// Get box start and stop positions


	// snap to the native coordinate system
	int Snap_Dimension =
		m_Op->SnapBox2Mesh(
				cSheet->GetStartCoord()->GetCoords(m_Op->m_MeshType),	// Start Coord
				cSheet->GetStopCoord()->GetCoords(m_Op->m_MeshType),	// Stop Coord
				m_sheetX0,	// Start Index
				m_sheetX1,	// Stop Index
				false,		// Dual H-mesh
				true);		// Full mesh?

	// Verify that snapped dimension is correct
	if (Snap_Dimension <= 0)
	{
		if (Snap_Dimension >= -1)
			cerr << "Operator_Ext_Absorbing_BC::SetInitParams(): Warning: Absorbing sheet snapping failed! Dimension is: " << Snap_Dimension << " skipping. "
					<< " ID: " << prim->GetID() << " @ Property: " << abc_prop->GetName() << endl;
		// Snap_Dimension == -2 means outside the simulation domain --> no special warning, but box probably marked as unused!
		return false;
	}

	// check that one of the dimensions is 1, and only one
	unsigned int sheetCheck = 0;
	m_ny = -1;
	for (int dimIdx = 0 ; dimIdx < 3 ; dimIdx++)
	{
		unsigned int Ncells = m_sheetX1[dimIdx] - m_sheetX0[dimIdx] + 1;

		// Count the number of dimensions that are single-cell length (start == stop)
		sheetCheck += Ncells == 1;

		// Update the normal direction
		m_ny = (Ncells == 1) ? dimIdx : m_ny;
	}

	if (sheetCheck != 1)
	{
		cerr 	<< "Operator_Ext_Absorbing_BC::SetInitParams(): Warning: Absorbing sheet is not a sheet! Skipping. "
				<< " ID: " << prim->GetID() << " @ Property: " << abc_prop->GetName() << endl;
		return false;
	}

	m_phaseVelocity = abc_prop->GetPhaseVelocity();
	// Currently not handling arbitrary phase velocity
	if (m_phaseVelocity == 0.0)
	{
		cerr 	<< "Operator_Ext_Absorbing_BC::SetInitParams(): Warning: Absorbing sheet Currently does not support per-material velocity. Setting to C0 "
				<< " ID: " << prim->GetID() << " @ Property: " << abc_prop->GetName() << endl;

		m_phaseVelocity = __C0__;
	}

	// Copy all of the relevant data, so BuildExtension ca0n work
	m_ABCtype = (ABCtype)(abc_prop->GetAbsorbingBoundaryType());
	m_normalSignPositive = abc_prop->GetNormalSignPositive();

	prim->SetPrimitiveUsed(true);

	return true;
}

bool Operator_Ext_Absorbing_BC::BuildExtension()
{

	double dT	= m_Op->GetTimestep();

	unsigned int	pos[] = {0,0,0};
	double			coord[] = {0.0,0.0,0.0};

	double 			delta;

	unsigned int	Ncells[] = {0,0,0},
					totCells = 1;

	// Count the number of cells in each direction
	for (unsigned int dimIdx = 0 ; dimIdx < 3 ; dimIdx++)
	{
		Ncells[dimIdx] = m_sheetX1[dimIdx] - m_sheetX0[dimIdx] + 1;
		totCells *= Ncells[dimIdx];
	}

	m_numCells = totCells;

	// Update remaining directions
	m_nyP 	= (m_ny + 1) % 3;
	m_nyPP 	= (m_ny + 2) % 3;

	// The position starts and stops in the same place.
	pos[m_ny] = m_sheetX0[m_ny];

	// The position is considered in the middle of the mesh cell.
	delta = fabs(m_Op->GetEdgeLength(m_ny,pos));

	// Initialize coefficients for this engine extension
	FDTD_FLOAT	vt_nyP  = m_phaseVelocity*dT,
				vt_nyPP = m_phaseVelocity*dT;

	// Initialize number of lines to be used in the arrayIJ
	m_numLines[0] = Ncells[m_nyP];
	m_numLines[1] = Ncells[m_nyPP];

	// Initialize containers. If there are more BCs in the future, this needs to be updated with the respective conditions.
	m_K1_nyP.Init("K1_Coeff_nyP", m_numLines);
	m_K1_nyPP.Init("K1_Coeff_nyPP", m_numLines);

	if (m_ABCtype == ABCtype::MUR_1ST_SA)
	{
		m_K2_nyP.Init("K2_Coeff_nyP", m_numLines);
		m_K2_nyPP.Init("K2_Coeff_nyPP", m_numLines);
	}

	// Prepare containers for per-material assignment
	coord[m_ny] = m_Op->GetDiscLine(m_ny,pos[m_ny]);
	// Check if this is the first or last cell
	if (m_sheetX0[m_ny] == 0)
		coord[m_ny] +=  delta/2 / m_Op->GetGridDelta();
	if (m_sheetX0[m_ny] == (m_Op->GetNumberOfLines(m_ny,true)-1))
		coord[m_ny] += -delta/2 / m_Op->GetGridDelta();

	// Initialize array coefficients
	unsigned int	arrI = 0,
					arrJ;

	for (pos[m_nyP] = m_sheetX0[m_nyP] ; pos[m_nyP] <= m_sheetX1[m_nyP] ; ++pos[m_nyP])
	{
		coord[m_nyP] = m_Op->GetDiscLine(m_nyP,pos[m_nyP]);
		arrJ = 0;
		for (pos[m_nyPP] = m_sheetX0[m_nyPP] ; pos[m_nyPP] <= m_sheetX1[m_nyPP] ; ++pos[m_nyPP])
		{
			coord[m_nyPP] = m_Op->GetDiscLine(m_nyPP,pos[m_nyPP]);

			if (m_phaseVelocity == 0.0)
			{
				double eps,mue;
				CSProperties** prop = m_Op->GetGeometryCSX()->GetPropertiesByCoordsPriority(coord, CSProperties::MATERIAL, false);
				if(prop != NULL)
				{
					cerr << "Operator_Ext_Absorbing_BC::BuildExtension(): Warning: This shouldn't happen...";
					/*CSPropMaterial* mat = (*prop)->ToMaterial();

					eps = mat->GetEpsilonWeighted(int(m_dir[1]),coord);
					mue = mat->GetMueWeighted(int(m_dir[1]),coord);
					vt_nyP = __C0__ * dT / sqrt(eps*mue);

					eps = mat->GetEpsilonWeighted(int(m_dir[2]),coord);
					mue = mat->GetMueWeighted(int(m_dir[2]),coord);
					vt_nyPP = __C0__ * dT / sqrt(eps*mue);*/
				}
				else
				{
					eps = m_Op->GetBackgroundEpsR();
					mue = m_Op->GetBackgroundMueR();

					vt_nyP  = __C0__ * dT / sqrt(eps*mue);
					vt_nyPP = __C0__ * dT / sqrt(eps*mue);
				}

			}

			// If more boundary types of boundary conditions are added in the future,
			// this needs to be in a condition, as well.
			m_K1_nyP(arrI,arrJ)  = (vt_nyP  - delta) / (vt_nyP  + delta);
			m_K1_nyPP(arrI,arrJ) = (vt_nyPP - delta) / (vt_nyPP + delta);

			if (m_ABCtype == ABCtype::MUR_1ST_SA)
			{
				m_K2_nyP(arrI,arrJ)  = vt_nyP /delta;
				m_K2_nyPP(arrI,arrJ) = vt_nyPP/delta;
			}

			arrJ++;
		}
		arrI++;
	}

	return true;
}

Engine_Extension* Operator_Ext_Absorbing_BC::CreateEngineExtention()
{
	Engine_Ext_Absorbing_BC* eng_ext = new Engine_Ext_Absorbing_BC(this);
	return eng_ext;
}

void Operator_Ext_Absorbing_BC::ShowStat(std::ostream &ostr) const
{
	Operator_Extension::ShowStat(ostr);

	ostr << " Total cells: " << m_numCells << endl;
}









