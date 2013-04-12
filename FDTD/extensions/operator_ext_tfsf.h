/*
*	Copyright (C) 2012 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#ifndef OPERATOR_EXT_TFSF_H
#define OPERATOR_EXT_TFSF_H

#include "operator_extension.h"
#include "FDTD/operator.h"
#include "tools/constants.h"

class Excitation;

class Operator_Ext_TFSF : public Operator_Extension
{
	friend class Engine_Ext_TFSF;
public:
	Operator_Ext_TFSF(Operator* op);
	~Operator_Ext_TFSF();

	virtual Operator_Extension* Clone(Operator* op);

	virtual bool BuildExtension();

	virtual Engine_Extension* CreateEngineExtention();

	virtual bool IsCylinderCoordsSave(bool closedAlpha, bool R0_included) const {UNUSED(closedAlpha); UNUSED(R0_included); return false;}
	virtual bool IsCylindricalMultiGridSave(bool child) const {UNUSED(child); return false;}

	// FIXME, this extension is not save to use with MPI
	virtual bool IsMPISave() const {return false;}

	virtual string GetExtensionName() const {return string("Total-Field/Scattered-Field Extension");}

	virtual void ShowStat(ostream &ostr) const;

	virtual void Init();
	virtual void Reset();

protected:
	Excitation* m_Exc;

	bool m_IncLow[3];
	bool m_ActiveDir[3][2]; // m_ActiveDir[direction][low/high]
	unsigned int m_Start[3];
	unsigned int m_Stop[3];
	unsigned int m_numLines[3];

	double m_PropDir[3];
	double m_E_Amp[3];
	double m_H_Amp[3];

	double m_Frequency;
	double m_PhVel;

	unsigned int m_maxDelay;

	// array setup [direction][low/high][component][ <mesh_position> ]
	unsigned int* m_VoltDelay[3][2][2];
	FDTD_FLOAT* m_VoltDelayDelta[3][2][2];
	FDTD_FLOAT* m_VoltAmp[3][2][2];

	unsigned int* m_CurrDelay[3][2][2];
	FDTD_FLOAT* m_CurrDelayDelta[3][2][2];
	FDTD_FLOAT* m_CurrAmp[3][2][2];

};

#endif // OPERATOR_EXT_TFSF_H
