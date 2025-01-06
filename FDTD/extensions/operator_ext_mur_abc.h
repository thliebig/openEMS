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

// Boundary condition "MUR_ABC" implements an absorbing boundary condition
// first proposed by Gerrit Mur [1]. Better explanations are available in
// standard FDTD textbooks, so readers are recommended to use textbooks (see
// Learn More) on the openEMS website instead of Mur's original paper.
//
// G. Mur, "Absorbing Boundary Conditions for the Finite-Difference
// Approximation of the Time-Domain Electromagnetic-Field Equations," in
// IEEE Transactions on Electromagnetic Compatibility, vol. EMC-23, no. 4,
// pp. 377-382, Nov. 1981, doi: 10.1109/TEMC.1981.303970.

#ifndef OPERATOR_EXT_MUR_ABC_H
#define OPERATOR_EXT_MUR_ABC_H

#include "FDTD/operator.h"
#include "operator_extension.h"
#include "tools/arraylib/array_ij.h"

class Operator_Ext_Mur_ABC : public Operator_Extension
{
	friend class Engine_Ext_Mur_ABC;
public:
	Operator_Ext_Mur_ABC(Operator* op);
	~Operator_Ext_Mur_ABC();

	virtual Operator_Extension* Clone(Operator* op);

	//! Define the direction of this ABC: \a ny=0,1,2 -> x,y,z and if at bottom_ny -> e.g. x=0 or x=end
	void SetDirection(int ny, bool top_ny);

	//! Set (override) the expected phase velocity of the incoming wave
	void SetPhaseVelocity(double c_phase) {m_v_phase=c_phase;};

	virtual bool BuildExtension();

	virtual Engine_Extension* CreateEngineExtention();

	virtual bool IsCylinderCoordsSave(bool closedAlpha, bool R0_included) const;
	virtual bool IsCylindricalMultiGridSave(bool child) const;
	virtual bool IsMPISave() const {return true;}

	virtual string GetExtensionName() const {return string("Mur ABC Extension");}

	virtual void ShowStat(ostream &ostr) const;

protected:
	Operator_Ext_Mur_ABC(Operator* op, Operator_Ext_Mur_ABC* op_ext);
	void Initialize();

	// The simulation box contains 6 faces, m_ny gives the normal direction
	// of the surface on which Mur's ABC is active, can be 0, 1, 2 for x, y, z.
	// If ABC is enabled on more than one face, multiple instances are created.
	int m_ny; // ny = axis (n) of the Yee (y) mesh

	// The next and second next orthogonal direction with respect to ABC's
	// normal direction m_ny. If the ABC's normal direction is Y (1),
	// m_nyP and m_nyPP would be 2 (Z) and 0 (X). This value wraps back to X.
	int m_nyP, m_nyPP; // P = Plus 1, PP = Plus 2.

	// Whether the current face is on the "top" (or leftmost, or front)
	// side of the simulation box.
	bool m_top;

	// Size of the Mur ABC 2D surface controlled by us. The physical significance
	// is that Mur's ABC is installed at a single face of the simulation box,
	// thus from our point of view, we only need to keep track of a 2D surface
	// m_numLines[2], not m_numLines[3]. m_numLines[0] and m_numLines[1] are the
	// number of Yee mesh lines on the direction m_nyP and m_nyPP.
	unsigned int m_numLines[2];

	// The line at which the Mur ABC is active, with respect to m_ny.
	unsigned int m_LineNr;

	// The nearest neighboring line from Mur ABC's active face, with respect to
	// m_ny. If the ABC is enabled on the "last" face of the plane (m_top == true),
	// m_LineNr is set to 0, and m_LineNr_Shift is set to 1. If the ABC is
	// enabled on the "first" face of the plane (!m_top), m_LineNr is set to
	// the last line of the mesh, and m_LineNr_Shift is set to the second last
	// line.
	int m_LineNr_Shift;

	double m_v_phase;
	ArrayLib::ArrayIJ<FDTD_FLOAT> m_Mur_Coeff_nyP;
	ArrayLib::ArrayIJ<FDTD_FLOAT> m_Mur_Coeff_nyPP;
};

#endif // OPERATOR_EXT_MUR_ABC_H
