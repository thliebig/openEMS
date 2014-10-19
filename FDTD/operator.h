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

#include "tools/AdrOp.h"
#include "tools/constants.h"
#include "excitation.h"
#include "Common/operator_base.h"

class Operator_Extension;
class Operator_Ext_Excitation;
class Engine;
class TiXmlElement;

//! Basic FDTD-operator
class Operator : public Operator_Base
{
	friend class Engine;
	friend class Engine_Interface_FDTD;
	friend class Operator_Ext_LorentzMaterial; //we need to find a way around this... friend class Operator_Extension only would be nice
	friend class Operator_Ext_ConductingSheet; //we need to find a way around this... friend class Operator_Extension only would be nice
	friend class Operator_Ext_PML_SF_Plane;
	friend class Operator_Ext_Excitation;
	friend class Operator_Ext_UPML;
	friend class Operator_Ext_Cylinder;
public:
	enum DebugFlags {None=0,debugMaterial=1,debugOperator=2,debugPEC=4};

	enum MatAverageMethods {QuarterCell=0, CentralCell=1};

	//! Create a new operator
	static Operator* New();
	virtual ~Operator();

	virtual Engine* CreateEngine();
	virtual Engine* GetEngine() const {return m_Engine;}

	virtual bool SetGeometryCSX(ContinuousStructure* geo);

	virtual int CalcECOperator( DebugFlags debugFlags = None );

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

	virtual void ApplyElectricBC(bool* dirs); //applied by default to all boundaries
	virtual void ApplyMagneticBC(bool* dirs);

	virtual void SetBCSize(int dir, int size) {m_BC_Size[dir]=size;}
	virtual int GetBCSize(int dir) {return m_BC_Size[dir];}

	//! Set a forced timestep to use by the operator
	virtual void SetTimestep(double ts) {dT = ts;}
	virtual void SetTimestepFactor(double factor);
	bool GetTimestepValid() const {return !m_InvaildTimestep;}

	//! Choose a time step method (0=auto, 1=CFL, 3=Rennings)
	void SetTimeStepMethod(int var) {m_TimeStepVar=var;}

	//! Set the material averaging method /sa MatAverageMethods
	void SetMaterialAvgMethod(MatAverageMethods method);

	//! Set material averaging method to the advanced quarter cell material interpolation (default)
	void SetQuarterCellMaterialAvg() {m_MatAverageMethod=QuarterCell;}

	//! Set operator to assume a constant material inside a cell (material probing in the cell center)
	void SetCellConstantMaterial() {m_MatAverageMethod=CentralCell;}

	virtual double GetNumberCells() const;

	virtual unsigned int GetNumberOfNyquistTimesteps() const {return m_Exc->GetNyquistNum();}

	virtual unsigned int GetNumberOfLines(int ny, bool full=false) const {UNUSED(full);return numLines[ny];}

	virtual void ShowStat() const;
	virtual void ShowExtStat() const;

	virtual double GetGridDelta() const {return gridDelta;}

	//! Get the disc line in \a n direction (in drawing units)
	virtual double GetDiscLine(int n, unsigned int pos, bool dualMesh=false) const;

	//! Get the disc line delta in \a n direction (in drawing units)
	virtual double GetDiscDelta(int n, unsigned int pos, bool dualMesh=false) const;

	//! Get the coordinates for a given node index and component, according to the yee-algorithm. Returns true if inside the FDTD domain.
	virtual bool GetYeeCoords(int ny, unsigned int pos[3], double* coords, bool dualMesh) const;

	virtual bool GetNodeCoords(const unsigned int pos[3], double* coords, bool dualMesh=false, CoordinateSystem c_system=UNDEFINED_CS) const;

	//! Get the node width for a given direction \a n and a given mesh position \a pos
	virtual double GetNodeWidth(int ny, const unsigned int pos[3], bool dualMesh = false) const {return GetEdgeLength(ny,pos,!dualMesh);}
	//! Get the node width for a given direction \a n and a given mesh position \a pos
	virtual double GetNodeWidth(int ny, const int pos[3], bool dualMesh = false) const;

	//! Get the node area for a given direction \a n and a given mesh position \a pos
	virtual double GetNodeArea(int ny, const unsigned int pos[3], bool dualMesh = false) const;
	//! Get the node area for a given direction \a n and a given mesh position \a pos
	virtual double GetNodeArea(int ny, const int pos[3], bool dualMesh = false) const;

	//! Get the length of an FDTD edge (unit is meter).
	virtual double GetEdgeLength(int ny, const unsigned int pos[3], bool dualMesh = false) const;

	//! Get the volume of an FDTD cell
	virtual double GetCellVolume(const unsigned int pos[3], bool dualMesh = false) const;

	//! Get the area around an edge for a given direction \a n and a given mesh posisition \a pos
	/*!
		This will return the area around an edge with a given direction, measured at the middle of the edge.
		In a cartesian mesh this is equal to the NodeArea, may be different in other coordinate systems.
	*/
	virtual double GetEdgeArea(int ny, const unsigned int pos[3], bool dualMesh = false) const {return GetNodeArea(ny,pos,dualMesh);}

	virtual unsigned int SnapToMeshLine(int ny, double coord, bool &inside, bool dualMesh=false, bool fullMesh=false) const;

	//! Snap the given coodinates to mesh indices
	virtual bool SnapToMesh(const double* coord, unsigned int* uicoord, bool dualMesh=false, bool fullMesh=false, bool* inside=NULL) const;

	//! Snap a given box to the FDTD mesh
	virtual int SnapBox2Mesh(const double* start, const double* stop, unsigned int* uiStart, unsigned int* uiStop, bool dualMesh=false, bool fullMesh=false, int SnapMethod=0, bool* bStartIn=NULL, bool* bStopIn=NULL) const;

	//! Snap a given line to the operator mesh
	/*!
	  \param[in] start coorindate of the line
	  \param[in] stop coorindate of the line
	  \param[out] uiStart the snapped line-start coorindate index
	  \param[out] uiStop the snapped line-stop coorindate index
	  \param[in] dualMesh snap to main or dual mesh (default is main mesh)
	  \return returns a status, 0 = success, 1 = start outside, 2 = stop outside, 3 = both outside
	  */
	virtual int SnapLine2Mesh(const double* start, const double* stop, unsigned int* uiStart, unsigned int* uiStop, bool dualMesh=false, bool fullMesh=false) const;

	virtual void AddExtension(Operator_Extension* op_ext);
	virtual void DeleteExtension(Operator_Extension* op_ext);
	virtual size_t GetNumberOfExtentions() const {return m_Op_exts.size();}
	virtual Operator_Extension* GetExtension(size_t index) const {return m_Op_exts.at(index);}

	virtual void CleanupMaterialStorage();

	virtual double GetDiscMaterial(int type, int ny, const unsigned int pos[3]) const;

	//! Get the cell center coordinate usable for material averaging (Warning, may not be the yee cell center)
	virtual bool GetCellCenterMaterialAvgCoord(const int pos[3], double coord[3]) const;

	virtual void SetExcitationSignal(Excitation* exc);
	virtual Excitation* GetExcitationSignal() const {return m_Exc;}

	Operator_Ext_Excitation* GetExcitationExtension() const;

	virtual double CalcNumericPhaseVelocity(unsigned int start[3], unsigned int stop[3], double propDir[3], float freq) const;

	virtual vector<CSPrimitives*> GetPrimitivesBoundBox(int posX, int posY, int posZ, CSProperties::PropertyType type=CSProperties::ANY) const;

protected:
	//! use New() for creating a new Operator
	Operator();

	virtual void Init();
	void Delete();
	virtual void Reset();
	virtual void InitOperator();
	virtual void InitDataStorage();

	virtual bool SetupCSXGrid(CSRectGrid* grid);

	virtual Grid_Path FindPath(double start[], double stop[]);

	// debug
	virtual void DumpOperator2File(string filename);
	virtual void DumpMaterial2File(string filename);
	virtual void DumpPEC2File( string filename, unsigned int *range = NULL );

	unsigned int m_Nr_PEC[3]; //count PEC edges
	virtual bool CalcPEC();
	virtual void CalcPEC_Range(unsigned int startX, unsigned int stopX, unsigned int* counter);	//internal to CalcPEC
	virtual void CalcPEC_Curves();	//internal to CalcPEC

	//Calc timestep only internal use
	int m_TimeStepVar;
	double m_TimeStepFactor;
	virtual double CalcTimestep();
	double opt_dT;
	bool m_InvaildTimestep;
	string m_Used_TS_Name;

	double CalcTimestep_Var1();
	double CalcTimestep_Var3();

	//! Calculate the FDTD equivalent circuit parameter for the given position and direction ny. \sa Calc_EffMat_Pos
	virtual bool Calc_ECPos(int ny, const unsigned int* pos, double* EC, vector<CSPrimitives *> vPrims) const;

	//! Get the FDTD raw disc delta, needed by Calc_EffMatPos() \sa Calc_EffMatPos
	/*!
	  Get the raw disc delta for a given position and direction.
	  The result will be positive if a disc delta inside the simulation domain is requested.
	  The result will be the negative value of the first or last disc delta respectivly if the position is outside the field domain.
	  */
	virtual double GetRawDiscDelta(int ny, const int pos) const;

	//! Get the material at a given coordinate, direction and type from CSX (internal use only)
	virtual double GetMaterial(int ny, const double coords[3], int MatType, vector<CSPrimitives*> vPrims, bool markAsUsed=true) const;

	MatAverageMethods m_MatAverageMethod;

	//! Calculate the effective/averaged material properties at the given position and direction ny.
	virtual bool Calc_EffMatPos(int ny, const unsigned int* pos, double* EffMat, vector<CSPrimitives*> vPrims) const;

	virtual bool AverageMatCellCenter(int ny, const unsigned int* pos, double* EffMat, vector<CSPrimitives*> vPrims) const;
	virtual bool AverageMatQuarterCell(int ny, const unsigned int* pos, double* EffMat, vector<CSPrimitives*> vPrims) const;

	//! Calc operator at certain \a pos
	virtual void Calc_ECOperatorPos(int n, unsigned int* pos);

	//! Calculate and setup lumped elements
	virtual bool Calc_LumpedElements();

	//! Store the size of the applied boundary conditions
	int m_BC_Size[6];

	//store material properties for post-processing
	float**** m_epsR;
	float**** m_kappa;
	float**** m_mueR;
	float**** m_sigma;

	//EC elements, internal only!
	virtual void Init_EC();
	virtual bool Calc_EC();
	virtual void Calc_EC_Range(unsigned int xStart, unsigned int xStop);
	FDTD_FLOAT* EC_C[3];
	FDTD_FLOAT* EC_G[3];
	FDTD_FLOAT* EC_L[3];
	FDTD_FLOAT* EC_R[3];

	AdrOp* MainOp;

	vector<Operator_Extension*> m_Op_exts;

	Engine* m_Engine;

	// excitation classes
	Excitation* m_Exc; // excitation time signal class
//	Operator_Ext_Excitation* m_Op_Ext_Exc; // excitation extension

	// engine/post-proc needs access
public:
	//EC operator
	FDTD_FLOAT**** vv; //calc new voltage from old voltage
	FDTD_FLOAT**** vi; //calc new voltage from old current
	FDTD_FLOAT**** ii; //calc new current from old current
	FDTD_FLOAT**** iv; //calc new current from old voltage
};

inline Operator::DebugFlags operator|( Operator::DebugFlags a, Operator::DebugFlags b ) { return static_cast<Operator::DebugFlags>(static_cast<int>(a) | static_cast<int>(b)); }
inline Operator::DebugFlags& operator|=( Operator::DebugFlags& a, const Operator::DebugFlags& b ) { return a = a | b; }

#endif // OPERATOR_H
