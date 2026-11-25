/*
*	Copyright (C) 2026 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#ifndef GPU_EXT_H
#define GPU_EXT_H

#include <vector>
#include <cstdint>

#include "gpu_backend.h"
#include "engine_plan.h"
#include "engine_subgrid.h"

class Engine_Extension;
class Engine_GPU;

//! The device extensions of every backend
/*!
  One class per engine extension, written once: it builds its point lists
  from the operator extension, uploads them with GPU_Backend::DeviceAlloc() and
  runs its per-point operations (kernels/ext_ops.h) with
  GPU_Backend::RunPoints(). The points are sorted along x, so that the box of
  a hook (a tile of the temporal blocking, or the whole domain) is a range of
  them, filtered along y and z on the device where the box cuts those
  (GPU_PointList::Run()).

  Returns NULL for an engine extension it has no device implementation of:
  Engine_GPU refuses the model then (every extension needs one).

  On a grid split over several backends (GPU_Backend_Split) it runs one
  device extension per part, see GPU_CreatePartExtension().
  */
GPU_Extension* GPU_CreateExtension(Engine_Extension* eng_ext, GPU_Backend* backend, const unsigned int numLines[3]);

//! The x lines of a grid a backend holds: all of them, or a part of a grid split along x over several backends (GPU_Backend_Split)
/*!
  A part owns the x lines [x0, x1) of the whole grid and holds them and, at
  each cut, the line on the other side of it (a halo line, which the neighbour
  updates); its local x line 0 is the line off of the whole grid.
  */
struct GPU_Slab
{
	unsigned int lines[3];   //!< of the whole grid
	unsigned int x0, x1;     //!< the x lines it owns
	unsigned int off;        //!< the x line of the whole grid at its local line 0
	unsigned int nx;         //!< its local x lines

	//! All lines of a grid of \a numLines
	static GPU_Slab Whole(const unsigned int numLines[3]);
	bool Owns(unsigned int x) const {return (x>=x0) && (x<x1);}
	bool Holds(unsigned int x) const {return (x>=off) && (x<off+nx);}
	//! Its local lines
	void Local(unsigned int l[3]) const {l[0] = nx; l[1] = lines[1]; l[2] = lines[2];}
	//! The local index of the flat ArrayNIJK index \a flat (or a node, component 0) of the whole grid, on an x line it holds
	uint64_t Local(uint64_t flat) const;
};

//! The device implementation of \a eng_ext on the part \a slab of a grid split over several backends, on \a backend (GPU_Backend_Split)
/*!
  Its points are those on the x lines the part owns. \a first is the device
  extension of \a eng_ext on the first part, NULL on the first part itself:
  the parts of the steady-state detection check the probes of all of them
  together. A part's extension may read its halo lines (GPU_Extension::HaloReads()).
  */
GPU_Extension* GPU_CreatePartExtension(Engine_Extension* eng_ext, GPU_Backend* backend, const GPU_Slab& slab, GPU_Extension* first);

//! The coupling of two GPU engines, a cylindrical multi-grid level and its sub-grid (Engine_SubGridCoupling), on the device: as SubGrid_Coupling_Basic
/*!
  The sub-grid's backend shares the work stream of the base grid's
  (GPU_Backend::NewSubGridBackend()). The coupling works on the device data of
  both, so the engine deletes it before either of them (Engine::ClearSubGrid()).
  */
class SubGrid_Coupling_GPU : public Engine_SubGridCoupling
{
public:
	SubGrid_Coupling_GPU(Engine_GPU* base, Engine_GPU* sub, const Operator_CylinderMultiGrid* op);

	//! Base grid voltages at r line split_pos-1 to the sub-grid, after the voltage updates of both grids
	virtual void SyncVoltages();
	//! Sub-grid currents to the base grid at r line split_pos-2, after the current updates of both grids
	virtual void SyncCurrents();
	//! Interpolate the sub-grid fields to the base grid inside the sub-grid (r lines below split_pos-1 / split_pos-2), at the end of a batch
	virtual void Interpolate();

protected:
	//! the interpolation of the voltages or of the currents on the device, [n][alpha] with n 0 and 1
	struct Weights {const uint32_t* p2p; const uint32_t* p2pp; const float* w2p; const float* w2pp;};
	void SetWeights(Weights& w, const std::vector<unsigned int> p2p[2], const std::vector<unsigned int> p2pp[2], const std::vector<float> w2p[2], const std::vector<float> w2pp[2]);
	void Interpolate(bool currents, unsigned int r0, unsigned int count);

	GPU_Backend* m_Base;
	GPU_Backend* m_Sub;
	unsigned int m_SplitPos;
	unsigned int m_BaseLines[3], m_SubLines[3];
	Weights m_Volt, m_Curr;
};

//! Points of a device extension sorted along x (stable: points at one x keep their order), and the points a box covers
class GPU_PointList
{
public:
	//! The points at the flat ArrayNIJK indices (or the nodes) \a flat of the whole grid of \a slab that it owns, sorted by x; Order() gives the original index of each sorted point
	/*!
	  Their y and z go to the device of \a b, for the boxes cut along y or z.
	  */
	void Init(GPU_Backend* b, const GPU_Slab& slab, const std::vector<uint64_t>& flat);
	size_t Size() const {return m_Order.size();}
	const std::vector<size_t>& Order() const {return m_Order;}
	//! The local flat indices of the sorted points
	const std::vector<uint64_t>& Index() const {return m_Index;}
	//! The local nodes (x*ny + y)*nz + z of the sorted points, added to \a nodes
	void AddNodes(std::vector<uint64_t>& nodes) const;
	//! Run the point operation \a op of type \a kind (GPU_Backend::RunPoints()) for the sorted points in \a box (local lines): their range along x, filtered along y and z if the box cuts those
	void Run(GPU_PointOp kind, const void* op, const Engine_Box& box, bool serial) const;

	//! The values \a v of the points in the sorted order
	template <typename T>
	std::vector<T> Sorted(const std::vector<T>& v) const
	{
		std::vector<T> s(m_Order.size());
		for (size_t n=0; n<m_Order.size(); ++n)
			s[n] = v[m_Order[n]];
		return s;
	}

protected:
	GPU_Backend* m_B = NULL;
	std::vector<size_t> m_Order;
	std::vector<uint64_t> m_Index;   //!< local flat index of the sorted points
	std::vector<unsigned int> m_X;   //!< local x of the sorted points
	unsigned int m_Lines[3] = {0, 0, 0};   //!< local
	const uint32_t* m_YZ = NULL;   //!< on the device: y and z of the sorted points
};

#endif // GPU_EXT_H
