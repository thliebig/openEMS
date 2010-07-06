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

#ifndef OPERATOR_H
#define OPERATOR_H

#include "ContinuousStructure.h"
#include "tools/AdrOp.h"
#include "tools/constants.h"
#include "excitation.h"

class Operator_Extension;
class Engine;

//! Abstract base-class for the FDTD-operator
class Operator
{
	friend class Engine;
	friend class Operator_Ext_LorentzMaterial; //we need to find a way around this... friend class Operator_Extension only would be nice
public:
	//! Create a new operator
	static Operator* New();
	virtual ~Operator();

	virtual Engine* CreateEngine() const;

	virtual bool SetGeometryCSX(ContinuousStructure* geo);
	virtual ContinuousStructure* GetGeometryCSX() {return CSX;}

	virtual int CalcECOperator();

	inline virtual FDTD_FLOAT& GetVV( unsigned int n, unsigned int x, unsigned int y, unsigned int z ) const { return vv[n][x][y][z]; }
	inline virtual FDTD_FLOAT& GetVI( unsigned int n, unsigned int x, unsigned int y, unsigned int z ) const { return vi[n][x][y][z]; }

	inline virtual FDTD_FLOAT& GetII( unsigned int n, unsigned int x, unsigned int y, unsigned int z ) const { return ii[n][x][y][z]; }
	inline virtual FDTD_FLOAT& GetIV( unsigned int n, unsigned int x, unsigned int y, unsigned int z ) const { return iv[n][x][y][z]; }

	virtual void SetBoundaryCondition(int* BCs) {for (int n=0;n<6;++n) m_BC[n]=BCs[n];}
	virtual void ApplyElectricBC(bool* dirs); //applied by default to all boundaries
	virtual void ApplyMagneticBC(bool* dirs);

	double GetTimestep() const {return dT;};
	double GetNumberCells() const;

	//! Returns the number of lines as needed for post-processing etc. (for the engine, use GetOriginalNumLines)
	virtual unsigned int GetNumberOfLines(int ny) const {return numLines[ny];}
	//! Returns the number of lines as needed for the engine etc. (for post-processing etc, use GetOriginalNumLines)
	virtual unsigned int GetOriginalNumLines(int ny) const {return numLines[ny];}

	virtual void ShowStat() const;

	virtual void DumpOperator2File(string filename);
	virtual void DumpMaterial2File(string filename);
	virtual void DumpPEC2File( string filename );

	//! Get the name for the given direction: 0 -> x, 1 -> y, 2 -> z
	virtual string GetDirName(int ny) const;

	virtual double GetGridDelta() const {return gridDelta;}
	//! Get the mesh delta times the grid delta for a 3D position
	virtual double GetMeshDelta(int n, const int* pos, bool dualMesh=false) const;
	virtual double GetMeshDelta(int n, const unsigned int* pos, bool dualMesh=false) const;

	//! Get the disc line in n direction
	virtual double GetDiscLine(int n, int pos, bool dualMesh=false) const;
	virtual double GetDiscLine(int n, unsigned int pos, bool dualMesh=false) const;
	virtual bool SnapToMesh(double* coord, unsigned int* uicoord, bool lower=false, bool* inside=NULL);

	virtual void AddExtension(Operator_Extension* op_ext);
	virtual size_t GetNumberOfExtentions() const {return m_Op_exts.size();}
	virtual Operator_Extension* GetExtension(size_t index) const {return m_Op_exts.at(index);}

protected:
	//! use New() for creating a new Operator
	Operator();

	virtual void Init();
	virtual void Reset();
	virtual void InitOperator();
	virtual void InitExcitation();

	struct Grid_Path
	{
		vector<unsigned int> posPath[3];
		vector<unsigned short> dir;
	};
	struct Grid_Path FindPath(double start[], double stop[]);

	ContinuousStructure* CSX;

	int m_BC[6];

	//! Calculate the field excitations.
	virtual bool CalcFieldExcitation();

	unsigned int m_Nr_PEC[3]; //count PEC edges
	virtual bool CalcPEC();
	virtual void CalcPEC_Range(unsigned int startX, unsigned int stopX, unsigned int* counter);	//internal to CalcPEC
	virtual void CalcPEC_Curves();	//internal to CalcPEC

	//Calc timestep only internal use
	virtual double CalcTimestep();
	double dT; //FDTD timestep!

	double CalcTimestep_Var1();
	double CalcTimestep_Var3();

	//! Calc operator at certain pos
	virtual void Calc_ECOperatorPos(int n, unsigned int* pos);

	//EC elements, internal only!
	virtual void Init_EC();
	virtual bool Calc_EC();
	virtual bool Calc_ECPos(int n, unsigned int* pos, double* inEC);
	virtual bool Calc_EffMatPos(int n, unsigned int* pos, double* inMat);
	double* EC_C[3];
	double* EC_G[3];
	double* EC_L[3];
	double* EC_R[3];

	int m_MeshType;
	unsigned int numLines[3];
	double* discLines[3];
	double gridDelta;
	AdrOp* MainOp;
	AdrOp* DualOp;

	vector<Operator_Extension*> m_Op_exts;

	// engine/post-proc needs access
public:
	//EC operator
	FDTD_FLOAT**** vv; //calc new voltage from old voltage
	FDTD_FLOAT**** vi; //calc new voltage from old current
	FDTD_FLOAT**** ii; //calc new current from old current
	FDTD_FLOAT**** iv; //calc new current from old voltage

	Excitation* Exc;
};

#endif // OPERATOR_H
