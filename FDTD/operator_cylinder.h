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

#ifndef OPERATOR_CYLINDER_H
#define OPERATOR_CYLINDER_H

#define OPERATOR_CYLINDER_CLOSED_ALPHA_THRESHOLD 1e-6

#include "operator_multithread.h"

class Operator_Ext_Cylinder;

//! This class creates an operator for a cylindrical FDTD.
/*!
This class creates an operator for a cylindrical FDTD. No special engine is necessary,
all special cases e.g. a closed alpha mesh or an included r=0 case is treated by an operator/engine extension \sa operator_ext_cylinder.
*/
class Operator_Cylinder : public Operator_Multithread
{
	friend class Operator_CylinderMultiGrid;
	friend class Operator_Ext_Cylinder;
	friend class Operator_Ext_LorentzMaterial;
public:
	static Operator_Cylinder* New(unsigned int numThreads = 0);
	virtual ~Operator_Cylinder();

	virtual void ApplyMagneticBC(bool* dirs);

	virtual unsigned int GetNumberOfLines(int ny, bool full=false) const;

	//! Get the name for the given direction: 0 -> rho, 1 -> alpha, 2 -> z
	virtual string GetDirName(int ny) const;

	//! Get the coordinates for a given node index and component, according to the cylindrical yee-algorithm. Returns true if inside the FDTD domain.
	virtual bool GetYeeCoords(int ny, unsigned int pos[3], double* coords, bool dualMesh) const;

	//! Get the node width for a given direction \a n and a given mesh posisition \a pos
	virtual double GetNodeWidth(int ny, const unsigned int pos[3], bool dualMesh = false) const;
	//! Get the node width for a given direction \a n and a given mesh posisition \a pos
	virtual double GetNodeWidth(int ny, const int pos[3], bool dualMesh = false) const;

	//! Get the node area for a given direction \a n and a given mesh posisition \a pos
	virtual double GetNodeArea(int n, const unsigned int* pos, bool dualMesh=false) const;
	//! Get the node area for a given direction \a n and a given mesh posisition \a pos
	virtual double GetNodeArea(int ny, const int pos[3], bool dualMesh = false) const;

	//! Get the length of an FDTD edge, including radius corrected alpha-mesh width.
	virtual double GetEdgeLength(int ny, const unsigned int pos[3], bool dualMesh = false) const;

	//! Get the volume of an FDTD cell
	virtual double GetCellVolume(const unsigned int pos[3], bool dualMesh = false) const;

	//! Get the area around an edge for a given direction \a n and a given mesh posisition \a pos
	/*!
		This will return the area around an edge with a given direction, measured at the middle of the edge.
		In a cartesian mesh this is equal to the NodeArea, may be different in other coordinate systems.
	*/
	virtual double GetEdgeArea(int ny, const unsigned int pos[3], bool dualMesh = false) const;

	virtual double FitToAlphaRange(double a_coord, bool fullMesh=false) const;

	//! Map a negative or out of range index in alpha direction back into the closed alpha mesh
	virtual int MapAlphaIndex2Range(int pos) const;

	virtual bool GetCellCenterMaterialAvgCoord(const int pos[3], double coord[3]) const;

	virtual unsigned int SnapToMeshLine(int ny, double coord, bool &inside, bool dualMesh=false, bool fullMesh=false) const;

	//! Snap a given box to the FDTD mesh
	virtual int SnapBox2Mesh(const double* start, const double* stop, unsigned int* uiStart, unsigned int* uiStop, bool dualMesh=false, bool fullMesh=false, int SnapMethod=0, bool* bStartIn=NULL, bool* bStopIn=NULL) const;

	virtual int SnapLine2Mesh(const double* start, const double* stop, unsigned int* uiStart, unsigned int* uiStop, bool dualMesh=false, bool fullMesh=false) const;

	bool GetClosedAlpha() const {return CC_closedAlpha;}
	bool GetR0Included() const {return CC_R0_included;}

	virtual void AddExtension(Operator_Extension* op_ext);

	virtual Engine* CreateEngine();

protected:
	Operator_Cylinder();
	virtual void Init();

	virtual bool SetupCSXGrid(CSRectGrid* grid);

	virtual Grid_Path FindPath(double start[], double stop[]);

	virtual double GetRawDiscDelta(int ny, const int pos) const;

	virtual double GetMaterial(int ny, const double coords[3], int MatType, vector<CSPrimitives*> vPrims, bool markAsUsed=true) const;

	virtual int CalcECOperator( DebugFlags debugFlags = None );
	virtual double CalcTimestep();
	bool CC_closedAlpha;
	bool CC_R0_included;
	Operator_Ext_Cylinder* m_Cyl_Ext;

#ifdef MPI_SUPPORT
	bool CC_MPI_Alpha;
#endif
};

#endif // OPERATOR_CYLINDER_H
