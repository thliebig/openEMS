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
class TiXmlElement;

//! Abstract base-class for the FDTD-operator
class Operator
{
	friend class Engine;
	friend class Operator_Ext_LorentzMaterial; //we need to find a way around this... friend class Operator_Extension only would be nice
	friend class Operator_Ext_PML_SF_Plane;
public:	
	enum DebugFlags {None=0,debugMaterial=1,debugOperator=2,debugPEC=4};

	//! Create a new operator
	static Operator* New();
	virtual ~Operator();

	virtual Engine* CreateEngine() const;

	virtual bool SetGeometryCSX(ContinuousStructure* geo);
	virtual ContinuousStructure* GetGeometryCSX() {return CSX;}

	virtual int CalcECOperator( DebugFlags debugFlags = None );

	//! Calculate the FDTD equivalent circuit parameter for the given position and direction ny. \sa Calc_EffMat_Pos
	virtual bool Calc_ECPos(int ny, const unsigned int* pos, double* EC) const;

	//! Calculate the effective/averaged material properties at the given position and direction ny.
	virtual bool Calc_EffMatPos(int ny, const unsigned int* pos, double* EffMat) const;

	virtual bool SetupExcitation(TiXmlElement* Excite, unsigned int maxTS) {return Exc->setupExcitation(Excite,maxTS);};

	// the next four functions need to be reimplemented in a derived class
	inline virtual FDTD_FLOAT GetVV( unsigned int n, unsigned int x, unsigned int y, unsigned int z ) const { return vv[n][x][y][z]; }
	inline virtual FDTD_FLOAT GetVI( unsigned int n, unsigned int x, unsigned int y, unsigned int z ) const { return vi[n][x][y][z]; }
	inline virtual FDTD_FLOAT GetII( unsigned int n, unsigned int x, unsigned int y, unsigned int z ) const { return ii[n][x][y][z]; }
	inline virtual FDTD_FLOAT GetIV( unsigned int n, unsigned int x, unsigned int y, unsigned int z ) const { return iv[n][x][y][z]; }

	// convenient access functions
	inline virtual FDTD_FLOAT GetVV( unsigned int n, unsigned int pos[3] ) const { return GetVV(n,pos[0],pos[1],pos[2]); }
	inline virtual FDTD_FLOAT GetVI( unsigned int n, unsigned int pos[3] ) const { return GetVI(n,pos[0],pos[1],pos[2]); }
	inline virtual FDTD_FLOAT GetII( unsigned int n, unsigned int pos[3] ) const { return GetII(n,pos[0],pos[1],pos[2]); }
	inline virtual FDTD_FLOAT GetIV( unsigned int n, unsigned int pos[3] ) const { return GetIV(n,pos[0],pos[1],pos[2]); }

	// the next four functions need to be reimplemented in a derived class
	inline virtual void SetVV( unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value ) { vv[n][x][y][z] = value; }
	inline virtual void SetVI( unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value ) { vi[n][x][y][z] = value; }
	inline virtual void SetII( unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value ) { ii[n][x][y][z] = value; }
	inline virtual void SetIV( unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value ) { iv[n][x][y][z] = value; }

	virtual void SetBoundaryCondition(int* BCs) {for (int n=0;n<6;++n) m_BC[n]=BCs[n];}
	virtual void ApplyElectricBC(bool* dirs); //applied by default to all boundaries
	virtual void ApplyMagneticBC(bool* dirs);

	//! Set a forced timestep to use by the operator
	virtual void SetTimestep(double ts) {dT = ts;}
	double GetTimestep() const {return dT;};
	bool GetTimestepValid() const {return !m_InvaildTimestep;}
	virtual double GetNumberCells() const;

	//! Returns the number of lines as needed for post-processing etc. (for the engine, use GetOriginalNumLines())
	virtual unsigned int GetNumberOfLines(int ny) const {return numLines[ny];}
	//! Returns the number of lines as needed for the engine etc. (for post-processing etc, use GetOriginalNumLines())
	virtual unsigned int GetOriginalNumLines(int ny) const {return numLines[ny];}

	virtual void ShowStat() const;
	virtual void ShowExtStat() const;

	//! Get the name for the given direction: 0 -> x, 1 -> y, 2 -> z
	virtual string GetDirName(int ny) const;

	virtual double GetGridDelta() const {return gridDelta;}
	//! Get the mesh delta times the grid delta for a 3D position (unit is meter)
	virtual double GetMeshDelta(int n, const unsigned int* pos, bool dualMesh=false) const;

	//! Get the disc line in \a n direction (in drawing units)
	virtual double GetDiscLine(int n, unsigned int pos, bool dualMesh=false) const;

	//! Get the node width for a given direction \a n and a given mesh position \a pos
	virtual double GetNodeWidth(int ny, const unsigned int pos[3], bool dualMesh = false) const {return GetMeshDelta(ny,pos,!dualMesh);}

	//! Get the node area for a given direction \a n and a given mesh position \a pos
	virtual double GetNodeArea(int ny, const unsigned int pos[3], bool dualMesh = false) const;

	//! Get the length of an FDTD edge (unit is meter).
	virtual double GetEdgeLength(int ny, const unsigned int pos[3], bool dualMesh = false) const {return GetMeshDelta(ny,pos,dualMesh);}

	//! Get the area around an edge for a given direction \a n and a given mesh posisition \a pos
	/*!
		This will return the area around an edge with a given direction, measured at the middle of the edge.
		In a cartesian mesh this is equal to the NodeArea, may be different in other coordinate systems.
	*/
	virtual double GetEdgeArea(int ny, const unsigned int pos[3], bool dualMesh = false) const {return GetNodeArea(ny,pos,dualMesh);}

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

	// debug
	virtual void DumpOperator2File(string filename);
	virtual void DumpMaterial2File(string filename);
	virtual void DumpPEC2File( string filename );

	unsigned int m_Nr_PEC[3]; //count PEC edges
	virtual bool CalcPEC();
	virtual void CalcPEC_Range(unsigned int startX, unsigned int stopX, unsigned int* counter);	//internal to CalcPEC
	virtual void CalcPEC_Curves();	//internal to CalcPEC

	//Calc timestep only internal use
	virtual double CalcTimestep();
	double dT; //FDTD timestep!
	double opt_dT;
	bool m_InvaildTimestep;
	string m_Used_TS_Name;

	double CalcTimestep_Var1();
	double CalcTimestep_Var3();

	//! Calc operator at certain \a pos
	virtual void Calc_ECOperatorPos(int n, unsigned int* pos);

	//EC elements, internal only!
	virtual void Init_EC();
	virtual bool Calc_EC();
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

inline Operator::DebugFlags operator|( Operator::DebugFlags a, Operator::DebugFlags b ) { return static_cast<Operator::DebugFlags>(static_cast<int>(a) | static_cast<int>(b)); }
inline Operator::DebugFlags& operator|=( Operator::DebugFlags& a, const Operator::DebugFlags& b ) { return a = a | b; }

#endif // OPERATOR_H
