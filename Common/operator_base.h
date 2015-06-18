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

#ifndef OPERATOR_BASE_H
#define OPERATOR_BASE_H

#include "ContinuousStructure.h"
#include "tools/global.h"
#include "Common/processing.h"
#include "string"

typedef struct
{
	std::vector<unsigned int> posPath[3];
	std::vector<unsigned short> dir;
} Grid_Path;

//! Abstract base-class for a common operator
class Operator_Base
{
public:
	virtual ~Operator_Base();

	virtual bool SetGeometryCSX(ContinuousStructure* geo);
	virtual ContinuousStructure* GetGeometryCSX() const {return CSX;}

	//! Get the timestep used by this operator
	virtual double GetTimestep() const {return dT;}

	//! Get the number of cells or nodes defined by this operator
	virtual double GetNumberCells() const =0;

	//! Get the number of timesteps satisfying the nyquist condition (may depend on the excitation)
	virtual unsigned int GetNumberOfNyquistTimesteps() const =0;

	//! Returns the number of lines as needed for post-processing etc.
	virtual unsigned int GetNumberOfLines(int ny, bool full=false) const =0;

	//! Get the name for the given direction: 0 -> x, 1 -> y, 2 -> z
	virtual std::string GetDirName(int ny) const;

	//! Get the grid drawing unit in m
	virtual double GetGridDelta() const =0;

	//! Get the disc line in \a n direction (in drawing units)
	virtual double GetDiscLine(int n, unsigned int pos, bool dualMesh=false) const =0;

	//! Get the disc line delta in \a n direction (in drawing units)
	virtual double GetDiscDelta(int n, unsigned int pos, bool dualMesh=false) const =0;

	//! Get the node width for a given direction \a n and a given mesh position \a pos
	virtual double GetNodeWidth(int ny, const unsigned int pos[3], bool dualMesh = false) const =0;

	//! Get the node area for a given direction \a n and a given mesh position \a pos
	virtual double GetNodeArea(int ny, const unsigned int pos[3], bool dualMesh = false) const =0;

	//! Get the length of an FDTD edge (unit is meter).
	virtual double GetEdgeLength(int ny, const unsigned int pos[3], bool dualMesh = false) const =0;

	//! Get the area around an edge for a given direction \a n and a given mesh posisition \a pos
	/*!
		This will return the area around an edge with a given direction, measured at the middle of the edge.
		In a cartesian mesh this is equal to the NodeArea, may be different in other coordinate systems.
	*/
	virtual double GetEdgeArea(int ny, const unsigned int pos[3], bool dualMesh = false) const =0;

	//! Get the volume of an FDTD cell
	virtual double GetCellVolume(const unsigned int pos[3], bool dualMesh = false) const =0;

	//! Snap the given coodinates to mesh indices, return box dimension
	virtual bool SnapToMesh(const double* coord, unsigned int* uicoord, bool dualMesh=false, bool fullMesh=false, bool* inside=NULL) const =0;

	//! Snap a given box to the operator mesh, uiStart will be always <= uiStop
	/*!
	  \param[in] start the box-start coorindate
	  \param[in] stop  the box-stopt coorindate
	  \param[out] uiStart the snapped box-start coorindate index
	  \param[out] uiStop the snapped box-stop coorindate index
	  \param[in] dualMesh snap to main or dual mesh (default is main mesh)
	  \param[in] SnapMethod Snapping method, 0=snap to closest line, 1/(2)=snap such that given box is inside (outside) the snapped lines
	  \return returns the box dimension or -1 if box is not inside the simulation domain
	  */
	virtual int SnapBox2Mesh(const double* start, const double* stop, unsigned int* uiStart, unsigned int* uiStop, bool dualMesh=false, bool fullMesh=false, int SnapMethod=0, bool* bStartIn=NULL, bool* bStopIn=NULL) const =0;

	//! Set the boundary conditions
	virtual void SetBoundaryCondition(int* BCs) {for (int n=0; n<6; ++n) m_BC[n]=BCs[n];}

	//! Set flags to store material data for post-processing
	virtual void SetMaterialStoreFlags(int type, bool val);

	//! Check storage flags and cleanup
	virtual void CleanupMaterialStorage() = 0;

	//! Get stored discrete material (if storage is enabled).
	virtual double GetDiscMaterial(int type, int ny, const unsigned int pos[3]) const = 0;

	//! Get the cell center coordinate usable for material averaging (Warning, may not be the yee cell center)
	bool GetCellCenterMaterialAvgCoord(const unsigned int pos[3], double coord[3]) const;

	//! Get the cell center coordinate usable for material averaging (Warning, may not be the yee cell center)
	virtual bool GetCellCenterMaterialAvgCoord(const int pos[3], double coord[3]) const = 0;

	virtual std::vector<CSPrimitives*> GetPrimitivesBoundBox(int posX, int posY, int posZ, CSProperties::PropertyType type=CSProperties::ANY) const = 0;

	//! Set the background material (default is vacuum)
	virtual void SetBackgroundMaterial(double epsR=0, double mueR=0, double kappa=0, double sigma=0, double density=0);

	//! Get background rel. electric permittivity
	double GetBackgroundEpsR() const {return m_BG_epsR;}
	//! Set background rel. electric permittivity
	void SetBackgroundEpsR(double val);

	//! Get background rel. magnetic permeability
	double GetBackgroundMueR() const {return m_BG_mueR;}
	//! Set background rel. magnetic permeability
	void SetBackgroundMueR(double val);

	//! Get background electric conductivity
	double GetBackgroundKappa() const {return m_BG_kappa;}
	//! Set background electric conductivity
	void SetBackgroundKappa(double val);

	//! Get background magnetic conductivity (artificial)
	double GetBackgroundSigma() const {return m_BG_sigma;}
	//! Set background magnetic conductivity (artificial)
	void SetBackgroundSigma(double val);

	//! Get background mass density
	double GetBackgroundDensity() const {return m_BG_density;}
	//! Set background mass density
	void SetBackgroundDensity(double val);

protected:
	Operator_Base();

	ContinuousStructure* CSX;

	virtual void Init();
	//! Cleanup data and reset
	void Delete();
	virtual void Reset();

	//! boundary conditions
	int m_BC[6];

	//! The operator timestep
	double dT;

	//! bool flag array to store material data for post-processing
	bool m_StoreMaterial[4];

	//! background materials
	double m_BG_epsR;
	double m_BG_mueR;
	double m_BG_kappa;
	double m_BG_sigma;
	double m_BG_density;

	CoordinateSystem m_MeshType;
	unsigned int numLines[3];
	double* discLines[3];
	double gridDelta;
};

#endif // OPERATOR_BASE_H
