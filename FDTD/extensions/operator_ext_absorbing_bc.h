/*
*	Copyright (C) 2023-2025 Gadi Lahav (gadi@rfwithcare.com)
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

/*
 * The localized boundary conditions currently implemented are denoted:
 * - First order Mur BC, single phased velocity
 * - First order Mur BC with super absorption, single phase velocity.
 * These were chosen after experimentation with the various boundary
 * conditions suggested in [1]. The amount of post processing and
 * absorption performance pointed towards the Single PV 1st order Mur.
 * Since this method also requires an extra mesh cell in the doublet grid,
 * the option to use the regular 1st order Mur was left available.
 *
 * References:
 * [1] Betz, Vaughn Timothy, and R. Mittra. "Absorbing boundary conditions for the finite-difference time-domain analysis of guided-wave structures." Coordinated Science Laboratory Report no. UILU-ENG-93-2243 (1993).‏
 */


#ifndef OPERATOR_EXT_ABSORBING_BC_H
#define OPERATOR_EXT_ABSORBING_BC_H

#include "FDTD/operator.h"
#include "operator_extension.h"
#include "tools/arraylib/array_ij.h"

#include "CSPropAbsorbingBC.h"

class Operator_Ext_Absorbing_BC : public Operator_Extension
{
	friend class Engine_Ext_Absorbing_BC;
public:

	// This should be a replica of the CSXCAD property, but can also be something
	// else, later
	enum ABCtype
	{
		UNDEFINED	= 0,
		MUR_1ST 	= 1,	// Mur's BC, 1st order
		MUR_1ST_SA 	= 2		// Mur's BC, 1st order, with Super Absorption
	};

	Operator_Ext_Absorbing_BC(Operator* op);
	~Operator_Ext_Absorbing_BC();

	virtual Operator_Extension* Clone(Operator* op);

	virtual bool BuildExtension();

	virtual Engine_Extension* CreateEngineExtention();

	//virtual bool IsMPISave() const {return true;}

	virtual string GetExtensionName() const {return string("Local absorbing boundary condition sheet");}

	virtual void ShowStat(ostream &ostr) const;

	//! Initialize all parameters, so the extension can be built later
	virtual bool SetInitParams(CSPrimitives* prim, CSPropAbsorbingBC* abc_prop);

protected:

	Operator_Ext_Absorbing_BC(Operator* op, Operator_Ext_Absorbing_BC* op_ext);
	void Initialize();

	unsigned int			m_numCells;		// Number of cells in each primitive

	// Storage for the directions. ny is the normal direction to the sheet
	int				m_ny,
					m_nyP,
					m_nyPP;

	// Storage for sheet bounding box start and stop.
	unsigned int	m_sheetX0[3];
	unsigned int	m_sheetX1[3];

	unsigned int 	m_numLines[2];

	bool 			m_normalSignPositive;


	ABCtype			m_ABCtype;

	double			m_phaseVelocity;

	// Coefficients, to be initialized on-demand.
	ArrayLib::ArrayIJ<FDTD_FLOAT>	m_K1_nyP;
	ArrayLib::ArrayIJ<FDTD_FLOAT>	m_K1_nyPP;
	ArrayLib::ArrayIJ<FDTD_FLOAT> 	m_K2_nyP;
	ArrayLib::ArrayIJ<FDTD_FLOAT>	m_K2_nyPP;



};

#endif // OPERATOR_EXT_ABSORBING_BC_H
