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

#ifndef OPERATOR_EXT_ABSORBING_BC_H
#define OPERATOR_EXT_ABSORBING_BC_H

#include "FDTD/operator.h"
#include "operator_extension.h"

#include "CSPropAbsorbingBC.h"

class Operator_Ext_Absorbing_BC : public Operator_Extension
{
	friend class Engine_Ext_Absorbing_BC;
public:
	//Operator_Ext_Absorbing_BC(Operator* op, Operator_Ext_Absorbing_BC* op_ext);
	Operator_Ext_Absorbing_BC(Operator* op);
	Operator_Ext_Absorbing_BC(Operator* op, uint numThreads);
	~Operator_Ext_Absorbing_BC();

	virtual Operator_Extension* Clone(Operator* op);



	virtual bool BuildExtension();

	virtual Engine_Extension* CreateEngineExtention();

	//virtual bool IsMPISave() const {return true;}

	virtual string GetExtensionName() const {return string("Local absorbing boundary condition sheet");}

	virtual void ShowStat(ostream &ostr) const;

protected:

	Operator_Ext_Absorbing_BC(Operator* op, Operator_Ext_Absorbing_BC* op_ext);
	void Initialize();

	uint 		m_numPrims;		// Number of primitives that will apply absorbing B.C.
	uint		*m_numCells;	// Number of cells in each primitive

	// Storage for the modulated direcions. m_dir[primIdx][0] is the normal, the rest are the tangential components.
	uint		**m_dir;

	// Storage for sheet bounding box start and stop.
	uint		**m_sheetX0;
	uint		**m_sheetX1;

	// Storage for the normal direction of each sheet
	int			*m_normDir;

	// Coefficients, to be initialized on-demand.
	FDTD_FLOAT	***m_K1;
	FDTD_FLOAT	***m_K2;

	int			*m_boundaryTypes;

	uint		m_Nthreads;

};

#endif // OPERATOR_EXT_ABSORBING_BC_H
