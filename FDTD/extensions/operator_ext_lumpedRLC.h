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

#ifndef OPERATOR_EXT_LUMPEDRLC_H_
#define OPERATOR_EXT_LUMPEDRLC_H_

#include <vector>

#include "FDTD/operator.h"
#include "operator_extension.h"
#include "operator_ext_cylinder.h"

#include "engine_ext_lumpedRLC.h"

#define LUMPED_RLC_Z_FACT 20.0

class Operator_Ext_LumpedRLC : public Operator_Extension
{
	friend class Engine_Ext_LumpedRLC;
	friend class CellModel;   // its edges in fused dispersive materials
	friend class GPU_Ext_LumpedRLC;
public:
	Operator_Ext_LumpedRLC(Operator* op);

	virtual ~Operator_Ext_LumpedRLC();

	virtual Operator_Extension* Clone(Operator* op);

	virtual bool BuildExtension();

	virtual Engine_Extension* CreateEngineExtention();

	virtual bool IsCylinderCoordsSave(bool closedAlpha, bool R0_included) const {UNUSED(closedAlpha); UNUSED(R0_included); return true;}
	virtual bool IsCylindricalMultiGridSave(bool child) const {UNUSED(child); return true;}

	virtual std::string GetExtensionName() const
	{
		return std::string("Series/Parallel Lumped RLC load");
	}

	virtual void ShowStat(std::ostream &ostr) const;

	virtual bool IsLElumpedRLC(const CSPropLumpedElement* const p_prop);

protected:
	//! Copy constructor
	Operator_Ext_LumpedRLC(Operator* op, Operator_Ext_LumpedRLC* op_ext);

	// ADE update coefficients, array setup: coeff[mesh_pos_index]

	// Parallel circuit coefficients
	FDTD_FLOAT		*v_RLC_ilv;
	FDTD_FLOAT		*v_RLC_i2v;

	// Series circuit coefficients
	// Vd[n] = vvd*(V + vv1*Vd[n-1] + vv2*Vd[n-2] + vj1*J[n-1] + vj2*J[n-2]), V the engine's update
	FDTD_FLOAT		*v_RLC_vvd;
	FDTD_FLOAT		*v_RLC_vv1;
	FDTD_FLOAT		*v_RLC_vv2;
	FDTD_FLOAT		*v_RLC_vj1;
	FDTD_FLOAT		*v_RLC_vj2;
	// J[n] = jv0*Vd[n] + jv1*Vd[n-1] + jv2*Vd[n-2] + jj1*J[n-1] + jj2*J[n-2]
	FDTD_FLOAT		*v_RLC_jv0;
	FDTD_FLOAT		*v_RLC_jv1;
	FDTD_FLOAT		*v_RLC_jv2;
	FDTD_FLOAT		*v_RLC_jj1;
	FDTD_FLOAT		*v_RLC_jj2;

	// Additional containers
	int				*v_RLC_dir;
	unsigned int			**v_RLC_pos;

	// Vector length indicator
	size_t					RLC_count;




};



#endif /* OPERATOR_EXT_LUMPEDRLC_H_ */
