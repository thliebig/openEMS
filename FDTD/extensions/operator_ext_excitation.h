/*
*	Copyright (C) 2011 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#ifndef OPERATOR_EXT_EXCITATION_H
#define OPERATOR_EXT_EXCITATION_H

#include "operator_extension.h"
#include "FDTD/operator.h"

class Excitation;

class Operator_Ext_Excitation : public Operator_Extension
{
	friend class Engine_Ext_Excitation;
	friend class Engine_Ext_Mur_ABC;
	friend class Operator;
public:
	Operator_Ext_Excitation(Operator* op);
	~Operator_Ext_Excitation();

	virtual Operator_Extension* Clone(Operator* op);

	virtual bool BuildExtension();

	virtual Engine_Extension* CreateEngineExtention();

	virtual bool IsCylinderCoordsSave(bool closedAlpha, bool R0_included) const {UNUSED(closedAlpha); UNUSED(R0_included); return true;}
	virtual bool IsCylindricalMultiGridSave(bool child) const {UNUSED(child); return true;}
	virtual bool IsMPISave() const {return true;}

	virtual string GetExtensionName() const {return string("Excitation Extension");}

	virtual void ShowStat(ostream &ostr) const;

	virtual void Init();
	virtual void Reset();

	unsigned int GetVoltCount() const {return Volt_Count;}
	unsigned int GetVoltCount(int ny) const {return Volt_Count_Dir[ny];}

	unsigned int GetCurrCount() const {return Curr_Count;}
	unsigned int GetCurrCount(int ny) const {return Curr_Count_Dir[ny];}

protected:
	Operator_Ext_Excitation(Operator* op, Operator_Ext_Excitation* op_ext);

	Excitation* m_Exc;

	void setupVoltageExcitation( vector<unsigned int> const volt_vIndex[3], vector<FDTD_FLOAT> const& volt_vExcit,
								 vector<unsigned int> const& volt_vDelay, vector<unsigned int> const& volt_vDir );
	void setupCurrentExcitation( vector<unsigned int> const curr_vIndex[3], vector<FDTD_FLOAT> const& curr_vExcit,
								 vector<unsigned int> const& curr_vDelay, vector<unsigned int> const& curr_vDir );
	//E-Field/voltage Excitation
	unsigned int Volt_Count;
	unsigned int Volt_Count_Dir[3];
	unsigned int* Volt_index[3];
	unsigned short* Volt_dir;
	FDTD_FLOAT* Volt_amp; //represented as edge-voltages!!
	unsigned int* Volt_delay;

	//H-Field/current Excitation
	unsigned int Curr_Count;
	unsigned int Curr_Count_Dir[3];
	unsigned int* Curr_index[3];
	unsigned short* Curr_dir;
	FDTD_FLOAT* Curr_amp; //represented as edge-currents!!
	unsigned int* Curr_delay;
};

#endif // OPERATOR_EXT_EXCITATION_H
