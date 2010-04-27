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

#ifndef ENGINE_EXT_MUR_ABC_H
#define ENGINE_EXT_MUR_ABC_H

#include "engine_extension.h"
#include "operator.h"

class Operator_Ext_Mur_ABC;

class Engine_Ext_Mur_ABC : public Engine_Extension
{
public:
	Engine_Ext_Mur_ABC(Operator_Ext_Mur_ABC* op_ext);
	virtual ~Engine_Ext_Mur_ABC();

	virtual void DoPreVoltageUpdates();
	virtual void DoPostVoltageUpdates();
	virtual void Apply2Voltages();

protected:
	Operator_Ext_Mur_ABC* m_Op_mur;

	int m_ny;
	int m_nyP,m_nyPP;
	unsigned int m_LineNr;
	int m_LineNr_Shift;
	unsigned int m_numLines[2];

	FDTD_FLOAT m_Mur_Coeff;
	FDTD_FLOAT** m_volt_nyP; //n+1 direction
	FDTD_FLOAT** m_volt_nyPP; //n+2 direction
};

#endif // ENGINE_EXT_MUR_ABC_H
