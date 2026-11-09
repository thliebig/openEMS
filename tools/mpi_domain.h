/*
*	Copyright (C) 2026 Ivan Iudice
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

#ifndef MPI_DOMAIN_H
#define MPI_DOMAIN_H

#include <string>
#include <vector>

#include "tools/global.h"

class CSRectGrid;
class Operator_Base;

//! The MPI domain decomposition of a Cartesian mesh into blocks, one block per rank
/*!
  Every rank runs a complete openEMS simulation on its block of mesh lines,
  plus one ghost line towards each neighbour: the two blocks of a split line
  \a s share the lines \a s-1 to \a s+1. The lower rank owns the tangential
  voltages on line \a s, the upper rank the tangential currents between \a s
  and \a s+1; Engine_Ext_MPI copies them across after every half step, so the
  planes the other rank cannot compute are overwritten with the right values.
  With the ghost lines nothing a rank owns lies on the first or last line of
  its mesh, which the operator and the extensions treat as the mesh boundary.

  Without MPI_SUPPORT, or with a single rank, Enabled() is false and the
  reductions return their argument.
  */
class OPENEMS_EXPORT MPI_Domain
{
public:
	MPI_Domain();

	//! Initialize MPI (multithreaded, serialized), call first thing in main()
	void Init(int* argc, char*** argv);
	//! Finalize MPI, call last thing in main()
	void Finalize();
	//! Terminate all ranks
	void Abort(int code);

	int Rank() const {return m_Rank;}
	int Size() const {return m_Size;}
	//! More than one rank: decompose the mesh and exchange the fields
	bool Enabled() const {return m_Size>1;}

	//! The number of blocks per direction, 0 for a choice by MPI_Dims_create(); set before Decompose()
	void SetSplitRequest(const std::string& splits);

	//! A box that no split may cross (a voltage, current or mode probe), in drawing units; set before Decompose()
	void AddUnsplitBox(const double start[3], const double stop[3]);

	//! Replace the lines of \a grid by the block of this rank, given the boundary conditions and their sizes in cells (PML)
	/*!
	  No split is placed inside a boundary layer of \a bcSize cells or through
	  a box of AddUnsplitBox(): a split moves to the nearest allowed line. Returns
	  false (on all ranks) if the requested decomposition is impossible.
	  */
	bool Decompose(CSRectGrid* grid, const unsigned int bcSize[6]);

	//! The neighbour rank below (\a up false) or above in direction \a n, -1 for none
	int GetNeighbor(int n, bool up) const {return up ? m_NeighborUp[n] : m_NeighborDown[n];}
	//! Whether the face \a idx (0..5 as for the boundary conditions) is an inner face shared with a neighbour
	bool IsInnerFace(int idx) const {return GetNeighbor(idx/2, idx%2)>=0;}

	//! The number of cells of the whole (undecomposed) mesh
	double GetGlobalNumberCells() const {return m_GlobalCells;}
	//! The whole mesh (only its lines, see openEMS::SetupFDTD()), NULL without MPI
	void SetMeshOperator(const Operator_Base* op) {m_MeshOp=op;}
	const Operator_Base* GetMeshOperator() const {return Enabled() ? m_MeshOp : NULL;}
	//! The global line of the first local line in direction \a n (0 without MPI)
	unsigned int GetFirstLine(int n) const {return Enabled() ? m_FirstLine[n] : 0;}

	//! The number of lines of the whole mesh in direction \a n; \a local without MPI
	unsigned int GetGlobalNumberOfLines(int n, unsigned int local) const {return Enabled() ? m_GlobalLines[n] : local;}
	//! The local mesh position of the global node \a gpos, if this rank owns it (as for OwnsBox); always true without MPI
	bool GlobalToLocal(const unsigned int gpos[3], unsigned int lpos[3]) const;
	//! The local line of the global line \a g in direction \a n of a field dump, if this rank dumps it
	/*!
	  Primary mesh: the lower rank dumps the split line, it has the fields around it.
	  Dual mesh: the upper rank dumps the cell above the split line, the lower rank
	  has no fields on its last (ghost) line.
	  */
	bool OwnsDumpLine(int n, unsigned int g, bool dualMesh, unsigned int& local) const;

	//! Collect the blocks of the ranks into the array \a global of rank 0
	/*!
	  Every rank has the block \a local of \a lsize nodes at \a offset of the whole array
	  of \a gsize nodes, both in the N-I-J-K order of ArrayNIJK with \a ncomp components
	  of \a elem floats (2 for a complex value). A rank may have no nodes (\a local NULL).
	  */
	void GatherBlocks(const float* local, const unsigned int lsize[3], const unsigned int offset[3],
					  float* global, const unsigned int gsize[3], int ncomp, int elem) const;

	//! Whether this rank owns the box \a start..\a stop (drawing units), see OwnsBox
	enum BoxOwner {BOX_OWNED, BOX_OTHER, BOX_SPLIT};
	//! The owner of a probe box: the rank whose block contains the lower corner; BOX_SPLIT if the box crosses into another block
	BoxOwner OwnsBox(const double start[3], const double stop[3]) const;

	//! A file name with the rank of this process prepended, unchanged without MPI
	std::string RankFileName(const std::string& name) const;

	// collective reductions over all ranks, the argument unchanged without MPI
	double AllreduceMin(double val) const;
	double AllreduceMax(double val) const;
	double AllreduceSum(double val) const;
	int AllreduceMin(int val) const;
	bool AllreduceOr(bool val) const;
	void Barrier() const;

	//! Exchange \a count floats with the neighbours in direction \a n: send \a send to rank \a to, receive \a recv from rank \a from (-1: no partner)
	void SendRecv(const float* send, int to, float* recv, int from, int count, int tag) const;

protected:
	int m_Rank;
	int m_Size;
	bool m_Initialized;

	std::vector<std::vector<double> > m_UnsplitBoxes;
	const Operator_Base* m_MeshOp;
	int m_SplitRequest[3];   //!< requested blocks per direction, 0 for automatic
	int m_Splits[3];         //!< blocks per direction
	int m_Coord[3];          //!< the block position of this rank
	int m_NeighborUp[3];
	int m_NeighborDown[3];
	double m_GlobalCells;
	unsigned int m_GlobalLines[3];
	unsigned int m_SplitLo[3];   //!< the global line of the lower split of this block
	unsigned int m_SplitHi[3];   //!< ... and of the upper one (the last line at the end of the mesh)
	unsigned int m_FirstLine[3]; //!< the global line of the first local line (the ghost line, if any)

	//! The owned interval of this rank in drawing units, per direction: [lo, hi), closed at the upper end of the mesh
	double m_OwnLo[3];
	double m_OwnHi[3];
	bool m_OwnHiClosed[3];
};

OPENEMS_EXPORT extern MPI_Domain g_mpi;

#endif // MPI_DOMAIN_H
