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
  a hook (a tile of the temporal blocking, whose boxes span y and z, or the
  whole domain) is a range of them.

  Returns NULL for an engine extension it has no device implementation of:
  Engine_GPU refuses the model then (every extension needs one).
  */
GPU_Extension* GPU_CreateExtension(Engine_Extension* eng_ext, GPU_Backend* backend, const unsigned int numLines[3]);

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

//! Points of a device extension sorted along x (stable: points at one x keep their order), and the range a box covers
class GPU_PointList
{
public:
	//! Sort the points by \a x; Order() gives the original index of each sorted point
	void Init(const std::vector<unsigned int>& x);
	size_t Size() const {return m_Order.size();}
	const std::vector<size_t>& Order() const {return m_Order;}
	//! The sorted points in [box.start[0], box.stop[0]) along x
	void Range(const Engine_Box& box, uint64_t& begin, uint64_t& end) const;

	//! The values \a v of the points in the sorted order
	template <typename T>
	std::vector<T> Sorted(const std::vector<T>& v) const
	{
		std::vector<T> s(v.size());
		for (size_t n=0; n<m_Order.size(); ++n)
			s[n] = v[m_Order[n]];
		return s;
	}

protected:
	std::vector<size_t> m_Order;
	std::vector<unsigned int> m_X;   //!< x of the sorted points
};

#endif // GPU_EXT_H
