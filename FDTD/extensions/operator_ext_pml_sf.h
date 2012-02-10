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

#ifndef OPERATOR_EXT_PML_SF_H
#define OPERATOR_EXT_PML_SF_H

#include "operator_extension.h"
#include "FDTD/operator.h"

class FunctionParser;

//! Insert split field pml planes, edges and corner as necessary by the given boundary conditions
bool Build_Split_Field_PML(Operator* op, int BC[6], int size[6], string gradFunc);

//! This is the abstract operator extension for truncating the FDTD domain with a split field pml
class Operator_Ext_PML_SF : public Operator_Extension
{
	friend class Engine_Ext_PML_SF;
	friend class Engine_Ext_PML_SF_Plane;
public:
	~Operator_Ext_PML_SF();

	virtual void SetBoundaryCondition(int* BCs) {for (int n=0; n<6; ++n) m_BC[n]=BCs[n];}

	inline virtual FDTD_FLOAT& GetVV(unsigned int nP, unsigned int n, unsigned int x, unsigned int y, unsigned int z ) const { return vv[nP][n][x][y][z]; }
	inline virtual FDTD_FLOAT& GetVI(unsigned int nP, unsigned int n, unsigned int x, unsigned int y, unsigned int z ) const { return vi[nP][n][x][y][z]; }
	inline virtual FDTD_FLOAT& GetII(unsigned int nP, unsigned int n, unsigned int x, unsigned int y, unsigned int z ) const { return ii[nP][n][x][y][z]; }
	inline virtual FDTD_FLOAT& GetIV(unsigned int nP, unsigned int n, unsigned int x, unsigned int y, unsigned int z ) const { return iv[nP][n][x][y][z]; }

	virtual double GetEdgeArea(int ny, unsigned int pos[3], bool dualMesh = false) const {UNUSED(ny); UNUSED(pos); UNUSED(dualMesh); return 0.0;}
	virtual double GetEdgeLength(int ny, unsigned int pos[3], bool dualMesh = false) const {UNUSED(ny); UNUSED(pos); UNUSED(dualMesh); return 0.0;}

	//! This will resturn the pml parameter grading
	virtual double GetKappaGraded(double depth, double Zm) const {UNUSED(depth); UNUSED(Zm); return 0.0;}

	//! Set the grading function for the pml
	/*!
		Define the pml grading grading function.
		Predefined variables in this grading function are:
			D  = depth in the pml in meter
			dl = mesh delta inside the pml in meter
			W  = width (length) of the pml in meter
			N  = number of cells for the pml
			Z  = wave impedance at the current depth and position
		example: SetGradingFunction("-log(1e-6)*log(2.5)/(2*dl*pow(2.5,W/dl)-1) * pow(2.5, D/dl) / Z");

		An empty function string will be ignored.
	*/
	virtual bool SetGradingFunction(string func);

	virtual bool BuildExtension();

	virtual string GetExtensionName() const {return string("Split Field PML Extension");}

//	virtual void ShowStat(ostream &ostr) const;

protected:
	Operator_Ext_PML_SF(Operator* op);

	virtual void ApplyBC() {};

	virtual bool Calc_ECPos(int nP, int n, unsigned int* pos, double* inEC) const {UNUSED(n); UNUSED(nP); UNUSED(pos); UNUSED(inEC); return true;};

	unsigned int m_numLines[3];
	bool m_SetupDone;

	int m_BC[6];

	string m_GradFunc;
	FunctionParser* m_GradingFunction;

	void InitOP();
	void DeleteOP();
	//split field EC operator
	//the first array-index is the splitted field part
	FDTD_FLOAT**** vv[2]; //calc new voltage from old voltage
	FDTD_FLOAT**** vi[2]; //calc new voltage from old current
	FDTD_FLOAT**** ii[2]; //calc new current from old current
	FDTD_FLOAT**** iv[2]; //calc new current from old voltage
};

//! This is an operator extension for truncating the FDTD domain with a split field pml plane
class Operator_Ext_PML_SF_Plane : public Operator_Ext_PML_SF
{
	friend class Engine_Ext_PML_SF_Plane;
public:
	Operator_Ext_PML_SF_Plane(Operator* op);
	~Operator_Ext_PML_SF_Plane();

	//! Define the direction of this PML plane: ny=0,1,2 -> x,y,z and if at bottom_ny -> e.g. x=0 or x=end
	void SetDirection(int ny, bool top_ny);
	void SetPMLLength(int width);

	virtual double GetEdgeArea(int ny, unsigned int pos[3], bool dualMesh = false) const;
	virtual double GetEdgeLength(int ny, unsigned int pos[3], bool dualMesh = false) const;

	virtual double GetKappaGraded(double depth, double Zm) const;

	virtual bool Calc_ECPos(int nP, int n, unsigned int* pos, double* inEC) const;

	virtual Engine_Extension* CreateEngineExtention();

	virtual bool IsCylinderCoordsSave(bool closedAlpha, bool R0_included) const;

	virtual string GetExtensionName() const {return string("Split Field PML Plane Extension");}

	virtual void ShowStat(ostream &ostr) const;

protected:
	virtual void ApplyBC();

	int m_ny,m_nyP,m_nyPP;
	bool m_top;
	unsigned int m_LineNr;

	double m_pml_delta;
	double m_pml_width;
};

#endif // OPERATOR_EXT_PML_SF_H
