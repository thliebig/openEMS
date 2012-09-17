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

#ifndef OPERATOR_EXT_UPML_H
#define OPERATOR_EXT_UPML_H

#include "FDTD/operator.h"
#include "operator_extension.h"

class FunctionParser;

//! Operator extension implemention an uniaxial perfectly matched layer (upml)
/*
  The priority for this extension should be the highest of all extensions since this operator will use the main engine to perform vital parts in the upml implementation.
  Therefore the voltages and currents as well as the operator are replaced during these update process.
  This extension is propably incompatible with the most other extensions operating in the same regions.
  */
class Operator_Ext_UPML : public Operator_Extension
{
	friend class Engine_Ext_UPML;
public:
	virtual ~Operator_Ext_UPML();

	//! Returns always true, Create_UPML method will take care of creating a valid pml for the cylindrical fdtd
	virtual bool IsCylinderCoordsSave(bool closedAlpha, bool R0_included) const { UNUSED(closedAlpha); UNUSED(R0_included); return true;}

	//! Returns always true if base grid, Create_UPML will create proper child pml extensions.
	virtual bool IsCylindricalMultiGridSave(bool child) const {if (child) return false; return true;}

	virtual bool IsMPISave() const {return true;}

	void SetBoundaryCondition(const int* BCs, const unsigned int size[6]);

	void SetRange(const unsigned int start[3], const unsigned int stop[3]);

	//! Set the grading function for the pml
	/*!
		Define the pml grading grading function.
		Predefined variables in this grading function are:
			D  = depth in the pml in meter
			dl = mesh delta inside the pml in meter
			W  = width (length) of the pml in meter
			N  = number of cells for the pml
			Z  = wave impedance at the current depth and position
		example: SetGradingFunction(" -log(1e-6)*log(2.5)/(2*dl*Z*(pow(2.5,W/dl)-1)) * pow(2.5, D/dl) ");

		An empty function string will be ignored.
	*/
	virtual bool SetGradingFunction(string func);

	virtual bool BuildExtension();

	virtual Engine_Extension* CreateEngineExtention();

	virtual string GetExtensionName() const {return string("Uniaxial PML Extension");}

	virtual void ShowStat(ostream &ostr) const;

	//! Create the UPML
	static bool Create_UPML(Operator* op, const int ui_BC[6], const unsigned int ui_size[6], const string gradFunc);

protected:
	Operator_Ext_UPML(Operator* op);
	int m_BC[6];
	unsigned int m_Size[6];

	unsigned int m_StartPos[3];
	unsigned int m_numLines[3];

	string m_GradFunc;
	FunctionParser* m_GradingFunction;

	void CalcGradingKappa(int ny, unsigned int pos[3], double Zm, double kappa_v[3], double kappa_i[3]);

	void DeleteOp();

	virtual FDTD_FLOAT& GetVV(int ny, unsigned int pos[3]) {return vv[ny][pos[0]][pos[1]][pos[2]];}
	virtual FDTD_FLOAT& GetVVFO(int ny, unsigned int pos[3]) {return vvfo[ny][pos[0]][pos[1]][pos[2]];}
	virtual FDTD_FLOAT& GetVVFN(int ny, unsigned int pos[3]) {return vvfn[ny][pos[0]][pos[1]][pos[2]];}
	virtual FDTD_FLOAT& GetII(int ny, unsigned int pos[3]) {return ii[ny][pos[0]][pos[1]][pos[2]];}
	virtual FDTD_FLOAT& GetIIFO(int ny, unsigned int pos[3]) {return iifo[ny][pos[0]][pos[1]][pos[2]];}
	virtual FDTD_FLOAT& GetIIFN(int ny, unsigned int pos[3]) {return iifn[ny][pos[0]][pos[1]][pos[2]];}

	FDTD_FLOAT**** vv;   //calc new voltage from old voltage
	FDTD_FLOAT**** vvfo; //calc new voltage from old voltage flux
	FDTD_FLOAT**** vvfn; //calc new voltage from new voltage flux
	FDTD_FLOAT**** ii;   //calc new current from old current
	FDTD_FLOAT**** iifo; //calc new current from old current flux
	FDTD_FLOAT**** iifn; //calc new current from new current flux
};

#endif // OPERATOR_EXT_UPML_H
