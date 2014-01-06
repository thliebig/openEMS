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

#ifndef OPERATOR_CYLINDERMULTIGRID_H
#define OPERATOR_CYLINDERMULTIGRID_H

#define CYLIDINDERMULTIGRID_LIMIT 20

#include "operator_cylinder.h"

//! This is a cylindrical FDTD operator using a simple multi-grid approach.
/*!
  This cylindrical multi-grid operator itself is not calculating any real operator, instead it is hosting two separate "child" operators of type "Operator_Cylinder".
  This operator class (or the corresponding engine) will perform the interpolation and connection between these two child-operator/engines.
  One of the child operators itself may be another multi-grid operator to allow for a cascaded multi-grid approach.
  */
class Operator_CylinderMultiGrid : public Operator_Cylinder
{
	friend class Engine_CylinderMultiGrid;
public:
	static Operator_CylinderMultiGrid* New(vector<double> Split_Radii, unsigned int numThreads = 0, unsigned int level = 0);
	virtual ~Operator_CylinderMultiGrid();

	virtual double GetNumberCells() const;

	virtual Engine* CreateEngine();

	virtual bool SetGeometryCSX(ContinuousStructure* geo);

	//! Get the coordinates for a given node index and component, according to the cylindrical yee-algorithm. Returns true if inside the FDTD domain.
	virtual bool GetYeeCoords(int ny, unsigned int pos[3], double* coords, bool dualMesh) const;

	virtual unsigned int GetSplitPos() const {return m_Split_Pos;}

	virtual void SetBoundaryCondition(int* BCs);

	virtual void AddExtension(Operator_Extension* op_ext);

	//! Get the multi grid level of this operator, e.g. 0 is main grid --> no parent grid
	int GetMultiGridLevel() const {return m_MultiGridLevel;}

	Operator_Cylinder* GetInnerOperator() const {return m_InnerOp;}

	virtual void SetExcitationSignal(Excitation* exc);

	virtual void ShowStat() const;

	//! Get the cell center coordinate usable for material averaging (Warning, may not be the yee cell center)
	virtual bool GetCellCenterMaterialAvgCoord(const int pos[3], double coord[3]) const;

#ifdef MPI_SUPPORT
	virtual void SetTag(int tag);
	virtual void SetNeighborUp(int ny, int id);
	virtual void SetNeighborDown(int ny, int id);
#endif

protected:
	Operator_CylinderMultiGrid(vector<double> Split_Radii, unsigned int level);
	virtual void Init();
	void Delete();
	virtual void Reset();

	virtual bool SetupCSXGrid(CSRectGrid* grid);

	virtual int CalcECOperator( DebugFlags debugFlags = None );

	virtual void DumpPEC2File( string filename, unsigned int *range = NULL );

	//! The material data storage in the sub-grid area's will not be filled by the base-operator. Check and do this here!
	void FillMissingDataStorage();

	unsigned int m_MultiGridLevel;
	double m_Split_Rad;
	vector<double> m_Split_Radii;
	unsigned int m_Split_Pos;

	Operator_Cylinder* m_InnerOp;

	// sub-grid to base interpolation coefficients
	unsigned int* m_interpol_pos_v_2p[2];
	f4vector* f4_interpol_v_2p[2];
	unsigned int* m_interpol_pos_v_2pp[2];
	f4vector* f4_interpol_v_2pp[2];

	unsigned int* m_interpol_pos_i_2p[2];
	f4vector* f4_interpol_i_2p[2];
	unsigned int* m_interpol_pos_i_2pp[2];
	f4vector* f4_interpol_i_2pp[2];

	void SetupInterpolation();

	virtual void CalcStartStopLines(unsigned int &numThreads, vector<unsigned int> &start, vector<unsigned int> &stop) const;
};

#endif // OPERATOR_CYLINDERMULTIGRID_H
