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

//! The MPI domain decomposition of a Cartesian mesh into blocks, one block per rank
/*!
  Every rank runs a complete openEMS simulation on its block of mesh lines. A
  block that has a neighbour above in direction n gets one more line, the
  first line of that neighbour's block: the two blocks then share the planes
  of the lines \a s and \a s+1, \a s being the split line. The lower rank owns
  the tangential voltages on line \a s, the upper rank the tangential currents
  between \a s and \a s+1; Engine_Ext_MPI copies them across after every half
  step, so the planes the other rank cannot compute are overwritten with the
  right values.

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

	//! Replace the lines of \a grid by the block of this rank, given the boundary conditions and their sizes in cells (PML)
	/*!
	  No split is placed inside a boundary layer of \a bcSize cells. Returns
	  false (on all ranks) if the requested decomposition is impossible.
	  */
	bool Decompose(CSRectGrid* grid, const unsigned int bcSize[6]);

	//! The neighbour rank below (\a up false) or above in direction \a n, -1 for none
	int GetNeighbor(int n, bool up) const {return up ? m_NeighborUp[n] : m_NeighborDown[n];}
	//! Whether the face \a idx (0..5 as for the boundary conditions) is an inner face shared with a neighbour
	bool IsInnerFace(int idx) const {return GetNeighbor(idx/2, idx%2)>=0;}

	//! The number of cells of the whole (undecomposed) mesh
	double GetGlobalNumberCells() const {return m_GlobalCells;}

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

	int m_SplitRequest[3];   //!< requested blocks per direction, 0 for automatic
	int m_Splits[3];         //!< blocks per direction
	int m_Coord[3];          //!< the block position of this rank
	int m_NeighborUp[3];
	int m_NeighborDown[3];
	double m_GlobalCells;

	//! The owned interval of this rank in drawing units, per direction: [lo, hi), closed at the upper end of the mesh
	double m_OwnLo[3];
	double m_OwnHi[3];
	bool m_OwnHiClosed[3];
};

OPENEMS_EXPORT extern MPI_Domain g_mpi;

#endif // MPI_DOMAIN_H
