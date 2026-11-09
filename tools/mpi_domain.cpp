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

#include "mpi_domain.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <utility>

#include "CSRectGrid.h"
#include "tools/useful.h"

#ifdef MPI_SUPPORT
// the C API only: the C++ bindings are gone since MPI-3
#define OMPI_SKIP_MPICXX 1
#define MPICH_SKIP_MPICXX 1
#include <mpi.h>
#endif

using std::cerr;
using std::cout;
using std::endl;

MPI_Domain g_mpi;

MPI_Domain::MPI_Domain()
{
	m_Rank = 0;
	m_Size = 1;
	m_Initialized = false;
	m_GlobalCells = 0;
	m_MeshOp = NULL;
	for (int n=0; n<3; ++n)
	{
		m_SplitRequest[n] = 0;
		m_Splits[n] = 1;
		m_Coord[n] = 0;
		m_NeighborUp[n] = -1;
		m_NeighborDown[n] = -1;
		m_OwnLo[n] = 0;
		m_OwnHi[n] = 0;
		m_OwnHiClosed[n] = true;
	}
}

void MPI_Domain::Init(int* argc, char*** argv)
{
#ifdef MPI_SUPPORT
	// the halo exchange runs on a worker thread of the engine, the reductions on the main
	// thread, never both at once
	int provided = 0;
	MPI_Init_thread(argc, argv, MPI_THREAD_SERIALIZED, &provided);
	m_Initialized = true;
	MPI_Comm_rank(MPI_COMM_WORLD, &m_Rank);
	MPI_Comm_size(MPI_COMM_WORLD, &m_Size);
	if ((provided<MPI_THREAD_SERIALIZED) && (m_Size>1))
	{
		if (m_Rank==0)
			cerr << "MPI_Domain::Init: Error: the MPI library does not provide MPI_THREAD_SERIALIZED" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
#else
	UNUSED(argc);
	UNUSED(argv);
#endif
}

void MPI_Domain::Finalize()
{
#ifdef MPI_SUPPORT
	if (m_Initialized)
		MPI_Finalize();
	m_Initialized = false;
#endif
}

void MPI_Domain::Abort(int code)
{
#ifdef MPI_SUPPORT
	if (m_Initialized)
		MPI_Abort(MPI_COMM_WORLD, code);
#else
	UNUSED(code);
#endif
}

void MPI_Domain::SetSplitRequest(const std::string& splits)
{
	std::vector<double> v = SplitString2Double(splits, ",");
	if (v.size()!=3)
		throw std::invalid_argument("--mpi-split expects three numbers, e.g. --mpi-split=4,2,1 (0 for automatic)");
	for (int n=0; n<3; ++n)
		m_SplitRequest[n] = (int)v.at(n);
}

// the split lines of \a cells cells into \a parts blocks, none inside the first \a lo or the last \a hi cells
// and none inside one of the line intervals \a keep; a split moves to the nearest line allowed
static bool CalcSplitLines(unsigned int cells, int parts, unsigned int lo, unsigned int hi,
						   const std::vector<std::pair<unsigned int, unsigned int> >& keep, std::vector<unsigned int>& lines)
{
	lines.clear();
	lines.push_back(0);
	// a split line s must keep the boundary layers in one block, and the line s+1 inside the mesh
	int minLine = lo+1;
	int maxLine = (cells>hi+1) ? cells-hi-1 : 0;
	auto allowed = [&](int s, int c)
	{
		// at least two cells per block, also for the blocks still to come
		if ((s<minLine) || (s>maxLine) || (s<(int)lines.back()+2) || (s>(int)cells-2*(parts-c)))
			return false;
		for (size_t k=0; k<keep.size(); ++k)
			if ((s>=(int)keep[k].first) && (s<=(int)keep[k].second))
				return false;
		return true;
	};
	for (int c=1; c<parts; ++c)
	{
		int s0 = (int)std::lround((double)c*cells/parts);
		int s = -1;
		for (int d=0; d<=(int)cells && s<0; ++d)
		{
			if (allowed(s0+d, c))
				s = s0+d;
			else if (allowed(s0-d, c))
				s = s0-d;
		}
		if (s<0)
			return false;
		lines.push_back(s);
	}
	if (cells<lines.back()+2)
		return false;
	lines.push_back(cells);
	return true;
}

void MPI_Domain::AddUnsplitBox(const double start[3], const double stop[3])
{
	std::vector<double> box(6);
	for (int n=0; n<3; ++n)
	{
		box[2*n] = std::min(start[n], stop[n]);
		box[2*n+1] = std::max(start[n], stop[n]);
	}
	m_UnsplitBoxes.push_back(box);
}

bool MPI_Domain::Decompose(CSRectGrid* grid, const unsigned int bcSize[6])
{
	if (!Enabled())
		return true;

	unsigned int numLines[3];
	for (int n=0; n<3; ++n)
		numLines[n] = grid->GetQtyLines(n);
	m_GlobalCells = (double)(numLines[0]-1)*(numLines[1]-1)*(numLines[2]-1);

	// the blocks per direction: as requested, or the factorization with the least inner interface area
	int fixed[3];
	for (int n=0; n<3; ++n)
		fixed[n] = m_SplitRequest[n];
	int best[3] = {0,0,0};
	double bestArea = -1;
	for (int p0=1; p0<=m_Size; ++p0)
	{
		if ((m_Size%p0) || (fixed[0] && fixed[0]!=p0))
			continue;
		for (int p1=1; p1<=m_Size/p0; ++p1)
		{
			if (((m_Size/p0)%p1) || (fixed[1] && fixed[1]!=p1))
				continue;
			int p2 = m_Size/p0/p1;
			if (fixed[2] && fixed[2]!=p2)
				continue;
			int p[3] = {p0,p1,p2};
			bool fits = true;
			double area = 0;
			for (int n=0; n<3; ++n)
			{
				if ((unsigned int)p[n]*2 > numLines[n]-1)
					fits = false;
				area += (double)(p[n]-1)*(numLines[(n+1)%3]-1)*(numLines[(n+2)%3]-1);
			}
			if (fits && ((bestArea<0) || (area<bestArea)))
			{
				bestArea = area;
				for (int n=0; n<3; ++n)
					best[n] = p[n];
			}
		}
	}
	if (bestArea<0)
	{
		if (m_Rank==0)
			cerr << "MPI_Domain::Decompose: Error: no decomposition of the mesh into " << m_Size << " blocks"
				 << " matches --mpi-split and leaves two cells per block" << endl;
		return false;
	}

	std::vector<unsigned int> split[3];
	for (int n=0; n<3; ++n)
	{
		m_Splits[n] = best[n];
		// no split through a probe: the line range it snaps to, with two lines to spare for the dual mesh
		std::vector<std::pair<unsigned int, unsigned int> > keep;
		for (size_t k=0; k<m_UnsplitBoxes.size(); ++k)
		{
			bool inside;
			unsigned int i0 = grid->Snap2LineNumber(n, m_UnsplitBoxes[k][2*n], inside);
			unsigned int i1 = grid->Snap2LineNumber(n, m_UnsplitBoxes[k][2*n+1], inside);
			keep.push_back(std::make_pair(i0>2 ? i0-2 : 0, i1+2));
		}
		if (!CalcSplitLines(numLines[n]-1, best[n], bcSize[2*n], bcSize[2*n+1], keep, split[n]))
		{
			if (m_Rank==0)
				cerr << "MPI_Domain::Decompose: Error: cannot split direction " << n << " into " << best[n]
					 << " blocks outside the boundary layers (PML) and the voltage/current probes;"
					 << " use fewer blocks in this direction (--mpi-split)" << endl;
			return false;
		}
	}

	// block position and neighbours, z fastest
	m_Coord[2] = m_Rank % m_Splits[2];
	m_Coord[1] = (m_Rank / m_Splits[2]) % m_Splits[1];
	m_Coord[0] = m_Rank / m_Splits[2] / m_Splits[1];
	for (int n=0; n<3; ++n)
	{
		int c[3] = {m_Coord[0], m_Coord[1], m_Coord[2]};
		m_NeighborDown[n] = -1;
		m_NeighborUp[n] = -1;
		if (m_Coord[n]>0)
		{
			c[n] = m_Coord[n]-1;
			m_NeighborDown[n] = (c[0]*m_Splits[1] + c[1])*m_Splits[2] + c[2];
		}
		if (m_Coord[n]<m_Splits[n]-1)
		{
			c[n] = m_Coord[n]+1;
			m_NeighborUp[n] = (c[0]*m_Splits[1] + c[1])*m_Splits[2] + c[2];
		}
	}

	// replace the mesh lines by the block of this rank, plus one ghost line towards each neighbour
	for (int n=0; n<3; ++n)
	{
		unsigned int qty = 0;
		double* lines = grid->GetLines(n, NULL, qty, true);
		unsigned int first = split[n].at(m_Coord[n]);
		unsigned int last = split[n].at(m_Coord[n]+1);
		if (m_NeighborDown[n]>=0)
			--first;
		if (m_NeighborUp[n]>=0)
			++last;
		m_GlobalLines[n] = numLines[n];
		m_SplitLo[n] = split[n].at(m_Coord[n]);
		m_SplitHi[n] = split[n].at(m_Coord[n]+1);
		m_FirstLine[n] = first;
		grid->ClearLines(n);
		for (unsigned int i=first; i<=last; ++i)
			grid->AddDiscLine(n, lines[i]);

		// probes from the line above the lower split up to and including the upper split line
		m_OwnLo[n] = lines[split[n].at(m_Coord[n])];
		m_OwnHi[n] = lines[split[n].at(m_Coord[n]+1)];
		m_OwnHiClosed[n] = (m_NeighborUp[n]<0);
		delete[] lines;
	}

	if (m_Rank==0)
		cout << "MPI_Domain: decomposition into " << m_Splits[0] << "x" << m_Splits[1] << "x" << m_Splits[2]
			 << " blocks on " << m_Size << " ranks" << endl;
	cout << "MPI_Domain: rank " << m_Rank << " block (" << m_Coord[0] << "," << m_Coord[1] << "," << m_Coord[2]
		 << "): " << grid->GetQtyLines(0) << "x" << grid->GetQtyLines(1) << "x" << grid->GetQtyLines(2) << " lines" << endl;
	return true;
}

MPI_Domain::BoxOwner MPI_Domain::OwnsBox(const double start[3], const double stop[3]) const
{
	if (!Enabled())
		return BOX_OWNED;
	for (int n=0; n<3; ++n)
	{
		double lo = std::min(start[n], stop[n]);
		double hi = std::max(start[n], stop[n]);
		// the lower rank owns the split line itself: it has the right values on it
		bool firstBlock = (m_NeighborDown[n]<0);
		bool aboveLo = firstBlock ? true : (lo>m_OwnLo[n]);
		bool belowHi = m_OwnHiClosed[n] ? true : (lo<=m_OwnHi[n]);
		if (!(aboveLo && belowHi))
			return BOX_OTHER;
		if (!m_OwnHiClosed[n] && (hi>m_OwnHi[n]))
			return BOX_SPLIT;
	}
	return BOX_OWNED;
}

bool MPI_Domain::GlobalToLocal(const unsigned int gpos[3], unsigned int lpos[3]) const
{
	for (int n=0; n<3; ++n)
	{
		if (!Enabled())
		{
			lpos[n] = gpos[n];
			continue;
		}
		// the lower rank owns the split line itself, like the probe boxes
		if ((m_NeighborDown[n]>=0) && (gpos[n]<=m_SplitLo[n]))
			return false;
		if (gpos[n]>m_SplitHi[n])
			return false;
		lpos[n] = gpos[n] - m_FirstLine[n];
	}
	return true;
}

bool MPI_Domain::OwnsDumpLine(int n, unsigned int g, bool dualMesh, unsigned int& local) const
{
	if (!Enabled())
	{
		local = g;
		return true;
	}
	if (dualMesh)
	{
		if ((m_NeighborDown[n]>=0) && (g<m_SplitLo[n]))
			return false;
		if ((m_NeighborUp[n]>=0) && (g>=m_SplitHi[n]))
			return false;
	}
	else
	{
		if ((m_NeighborDown[n]>=0) && (g<=m_SplitLo[n]))
			return false;
		if ((m_NeighborUp[n]>=0) && (g>m_SplitHi[n]))
			return false;
	}
	local = g - m_FirstLine[n];
	return true;
}

void MPI_Domain::GatherBlocks(const float* local, const unsigned int lsize[3], const unsigned int offset[3],
							  float* global, const unsigned int gsize[3], int ncomp, int elem) const
{
#ifdef MPI_SUPPORT
	int meta[6] = {(int)lsize[0], (int)lsize[1], (int)lsize[2], (int)offset[0], (int)offset[1], (int)offset[2]};
	int count = ncomp*elem*(int)((size_t)lsize[0]*lsize[1]*lsize[2]);
	std::vector<int> metas, counts, displs;
	if (m_Rank==0)
	{
		metas.resize(6*m_Size);
		counts.resize(m_Size);
		displs.resize(m_Size);
	}
	MPI_Gather(meta, 6, MPI_INT, metas.data(), 6, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Gather(&count, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
	std::vector<float> buf;
	if (m_Rank==0)
	{
		int total = 0;
		for (int r=0; r<m_Size; ++r)
		{
			displs[r] = total;
			total += counts[r];
		}
		buf.resize(total);
	}
	MPI_Gatherv(local, count, MPI_FLOAT, buf.data(), counts.data(), displs.data(), MPI_FLOAT, 0, MPI_COMM_WORLD);
	if (m_Rank!=0)
		return;
	for (int r=0; r<m_Size; ++r)
	{
		const int* m = &metas[6*r];
		const float* src = buf.data() + displs[r];
		for (int c=0; c<ncomp; ++c)
			for (int i=0; i<m[0]; ++i)
				for (int j=0; j<m[1]; ++j)
				{
					size_t g = (((size_t)c*gsize[0] + m[3]+i)*gsize[1] + m[4]+j)*gsize[2] + m[5];
					std::copy(src, src + (size_t)m[2]*elem, global + g*elem);
					src += (size_t)m[2]*elem;
				}
	}
#else
	UNUSED(local);
	UNUSED(lsize);
	UNUSED(offset);
	UNUSED(global);
	UNUSED(gsize);
	UNUSED(ncomp);
	UNUSED(elem);
#endif
}

std::string MPI_Domain::RankFileName(const std::string& name) const
{
	if (!Enabled())
		return name;
	std::stringstream ss;
	ss << "rank" << m_Rank << "_" << name;
	return ss.str();
}

double MPI_Domain::AllreduceMin(double val) const
{
#ifdef MPI_SUPPORT
	if (Enabled())
	{
		double res = val;
		MPI_Allreduce(&val, &res, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
		return res;
	}
#endif
	return val;
}

double MPI_Domain::AllreduceMax(double val) const
{
#ifdef MPI_SUPPORT
	if (Enabled())
	{
		double res = val;
		MPI_Allreduce(&val, &res, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		return res;
	}
#endif
	return val;
}

double MPI_Domain::AllreduceSum(double val) const
{
#ifdef MPI_SUPPORT
	if (Enabled())
	{
		double res = val;
		MPI_Allreduce(&val, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		return res;
	}
#endif
	return val;
}

int MPI_Domain::AllreduceMin(int val) const
{
#ifdef MPI_SUPPORT
	if (Enabled())
	{
		int res = val;
		MPI_Allreduce(&val, &res, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
		return res;
	}
#endif
	return val;
}

bool MPI_Domain::AllreduceOr(bool val) const
{
#ifdef MPI_SUPPORT
	if (Enabled())
	{
		int v = val ? 1 : 0;
		int res = v;
		MPI_Allreduce(&v, &res, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
		return res!=0;
	}
#endif
	return val;
}

void MPI_Domain::Barrier() const
{
#ifdef MPI_SUPPORT
	if (Enabled())
		MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void MPI_Domain::SendRecv(const float* send, int to, float* recv, int from, int count, int tag) const
{
#ifdef MPI_SUPPORT
	MPI_Sendrecv(send, count, MPI_FLOAT, (to<0) ? MPI_PROC_NULL : to, tag,
				 recv, count, MPI_FLOAT, (from<0) ? MPI_PROC_NULL : from, tag,
				 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#else
	UNUSED(send);
	UNUSED(to);
	UNUSED(recv);
	UNUSED(from);
	UNUSED(count);
	UNUSED(tag);
#endif
}
