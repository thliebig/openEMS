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

#ifndef OPERATOR_EXT_LUMPEDRLC_H_
#define OPERATOR_EXT_LUMPEDRLC_H_

#include "vector"

#include "FDTD/operator.h"
#include "operator_extension.h"
#include "operator_ext_cylinder.h"

#include "engine_ext_lumpedRLC.h"

#define LUMPED_RLC_Z_FACT 20.0

class Operator_Ext_LumpedRLC : public Operator_Extension
{
	friend class Engine_Ext_LumpedRLC;
public:
	Operator_Ext_LumpedRLC(Operator* op);

	virtual ~Operator_Ext_LumpedRLC();

	virtual Operator_Extension* Clone(Operator* op);

	virtual bool BuildExtension();

	virtual Engine_Extension* CreateEngineExtention();

	virtual bool IsCylinderCoordsSave(bool closedAlpha, bool R0_included) const {UNUSED(closedAlpha); UNUSED(R0_included); return true;}
	virtual bool IsCylindricalMultiGridSave(bool child) const {UNUSED(child); return true;}
	virtual bool IsMPISave() const {return true;}

	virtual string GetExtensionName() const {return string("Series\\Parallel Lumped RLC load");}

	virtual void ShowStat(ostream &ostr) const;

	virtual bool IsLElumpedRLC(const CSPropLumpedElement* const p_prop);

protected:
	//! Copy constructor
	Operator_Ext_LumpedRLC(Operator* op, Operator_Ext_LumpedRLC* op_ext);

	// ADE update coefficients, array setup: coeff[mesh_pos_index]

	// Parallel circuit coefficients
	FDTD_FLOAT		*v_RLC_ilv;
	FDTD_FLOAT		*v_RLC_i2v;

	// Series circuit coefficients
	FDTD_FLOAT		*v_RLC_vv2;		// Coefficient for [n-2] time of Vd update in Vd equation
	FDTD_FLOAT		*v_RLC_vj1;		// Coefficient for [n-1] time of J update in Vd equation
	FDTD_FLOAT		*v_RLC_vj2;		// Coefficient for [n-2] time of J update in Vd equation
	FDTD_FLOAT		*v_RLC_vvd;		// Coefficient to multiply all Vd in the Vd update equation
	FDTD_FLOAT		*v_RLC_ib0;		// Inverse of beta_0
	FDTD_FLOAT		*v_RLC_b1;		// beta_1
	FDTD_FLOAT		*v_RLC_b2;		// beta_2

	// Additional containers
	int				*v_RLC_dir;
	unsigned int			**v_RLC_pos;

	// Vector length indicator
	unsigned int			RLC_count;




};



#endif /* OPERATOR_EXT_LUMPEDRLC_H_ */
