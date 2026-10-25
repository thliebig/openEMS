/*
*	Copyright (C) 2026 Roel Budé (roel.bude@proton.me)
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

#ifndef OPERATOR_EXT_LOSSYMETAL_H
#define OPERATOR_EXT_LOSSYMETAL_H

#include "FDTD/operator.h"
#include "operator_ext_lorentzmaterial.h"
#include "surface_impedance.h"

#include <map>

class CSPropLossyMetal;

/*!
  Surface impedance of lossy (and rough) metal bodies, see CSPropLossyMetal.

  The inside of a LossyMetal body stays PEC. Every E edge on its surface carries the
  surface admittance of the faces it touches, Ys*width/length, as a passive rational fit
  G + sum_k 1/(R_k + j w L_k) (SurfaceImpedance::FitFace) on the Lorentz/ADE machinery,
  so all engines that run Operator_Ext_LorentzMaterial run this too.

  The width of an edge is the length of metal surface inside its dual face (see
  AddContributions): a flat face gives the node width, a convex or concave corner the
  two half-widths on the two faces. A body thinner than the mesh (no cell centre inside)
  is a sheet with both faces open, i.e. the model of Operator_Ext_ConductingSheet, with
  the thickness taken from the body, or from the property for a 2D primitive.
  */
class Operator_Ext_LossyMetal : public Operator_Ext_LorentzMaterial
{
public:
	Operator_Ext_LossyMetal(Operator* op, double f_max);

	virtual Operator_Extension* Clone(Operator* op);

	virtual bool BuildExtension();

	virtual bool IsCylinderCoordsSave(bool closedAlpha, bool R0_included) const {UNUSED(closedAlpha); UNUSED(R0_included); return false;}
	virtual bool IsCylindricalMultiGridSave(bool child) const {UNUSED(child); return false;}

	virtual std::string GetExtensionName() const
	{
		return std::string("Lossy Metal Surface Impedance Extension");
	}

	virtual void ShowStat(std::ostream &ostr) const;

protected:
	//! Copy constructor
	Operator_Ext_LossyMetal(Operator* op, Operator_Ext_LossyMetal* op_ext);

	//! A piece of metal face inside the dual face of an edge
	struct Contribution
	{
		int fit;	   //!< index into m_Fits
		double weight; //!< (surface length inside the dual face) / (edge length), times 2 for a sheet
	};

	//! Surface contributions of the E edge (n,pos) at \a coord, empty if it lies inside the metal
	bool AddContributions(int n, const unsigned int* pos, const double* coord, CSPropLossyMetal* prop, std::vector<Contribution>& contrib);

	//! Lossy metal body at this point?
	bool IsLossy(const double* coord) const;
	//! Any metal at this point?
	bool IsMetal(const double* coord) const;
	//! Length (drawing units) of the lossy metal from \a coord along direction \a dir, sign \a s, at most \a max
	double InsideLength(const double* coord, int dir, int s, double max) const;
	//! Extent (drawing units) of the lossy metal from \a coord along \a dir, sign \a s, following the mesh up to the end of the body
	double Extent(const double* coord, const unsigned int* pos, int dir, int s) const;

	//! Index of the fit for this property and face thickness (m)
	int GetFit(CSPropLossyMetal* prop, double thickness);

	double m_f_max;
	int m_NumBranches;
	std::vector<double> m_w_k;

	struct FitKey
	{
		CSPropLossyMetal* prop;
		long t_index;
		bool operator<(const FitKey& o) const {return (prop<o.prop) || ((prop==o.prop) && (t_index<o.t_index));}
	};
	std::map<FitKey,int> m_FitIndex;
	std::vector<SurfaceImpedance::Fit> m_Fits;
	std::vector<double> m_FitThickness;

	unsigned int m_numSurface, m_numFinite, m_numPEC, m_numH, m_numNoThickness;
};

#endif // OPERATOR_EXT_LOSSYMETAL_H
