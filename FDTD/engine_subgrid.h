/*
*	Copyright (C) 2010-2026 Thorsten Liebig (Thorsten.Liebig@gmx.de)
*	Copyright (C) 2026 Sean Mollet (sean@malmoset.com)
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

#ifndef ENGINE_SUBGRID_H
#define ENGINE_SUBGRID_H

#include "tools/constants.h"
#include "tools/arraylib/array_nijk.h"
#include "tools/global.h"
#include "tools/array_ops.h"
#include "tools/aligned_allocator.h"
#include "multigrid_interpolation.h"

#include <vector>

class Engine;
class Engine_sse;
class Engine_GPU;
class GPU_MultiGridLink;
class Operator_CylinderMultiGrid;

//! The coupling of a cylindrical multi-grid level (base grid) and its sub-grid, see Engine::SetSubGrid()
/*!
  The sub-grid has the same r and z lines as the base grid up to the split
  radius, and every other alpha line. After the voltage half step of both grids
  the base voltages on the split radius go to the sub-grid, after the current
  half step the sub-grid currents next to it go to the base grid; at the end of
  a batch the base grid is filled inside the split radius from the sub-grid.
  One implementation per field layout of the engines.
  */
class Engine_SubGridCoupling
{
public:
	//! The coupling for the field layout of \a base and \a sub (Engine::GetType()); throws if the two differ or there is none for it
	static Engine_SubGridCoupling* New(Engine* base, Engine* sub, const Operator_CylinderMultiGrid* op);

	virtual ~Engine_SubGridCoupling() {}
	virtual void SyncVoltages() = 0;
	virtual void SyncCurrents() = 0;
	virtual void Interpolate() = 0;
};

//! On the basic field layout (Engine), also the host fallback of SubGrid_Coupling_GPU
class SubGrid_Coupling_Basic : public Engine_SubGridCoupling
{
public:
	SubGrid_Coupling_Basic(Engine* base, Engine* sub, const Operator_CylinderMultiGrid* op);

	virtual void SyncVoltages();
	virtual void SyncCurrents();
	virtual void Interpolate();

	//! base voltages on r line split_pos-1 to the sub-grid
	static void SyncVoltages(const ArrayLib::ArrayNIJK<FDTD_FLOAT>& base, ArrayLib::ArrayNIJK<FDTD_FLOAT>& sub, const MultiGridInterpolation& ip, const unsigned int numLines[3]);
	//! base grid fields on r line \a r from the sub-grid (\a currents: the current interpolation)
	static void Interpolate(ArrayLib::ArrayNIJK<FDTD_FLOAT>& base, const ArrayLib::ArrayNIJK<FDTD_FLOAT>& sub, const MultiGridInterpolation& ip, bool currents, unsigned int r, const unsigned int numLines[3]);

protected:
	Engine* m_Base;
	Engine* m_Sub;
	MultiGridInterpolation m_Interpol;
	unsigned int m_numLines[3];
};

//! Of two GPU engines (Engine_GPU): on the device with the GPU_MultiGridLink of the backend, else on the host mirrors
/*!
  The link needs both grids on the device. Without one (one of them in the host
  fallback, or the backend has no link), both grids and their sub-grids are
  forced into the host fallback, and the coupling of SubGrid_Coupling_Basic
  works on the host mirrors, which are then up to date after every half step.
  */
class SubGrid_Coupling_GPU : public SubGrid_Coupling_Basic
{
public:
	SubGrid_Coupling_GPU(Engine_GPU* base, Engine_GPU* sub, const Operator_CylinderMultiGrid* op);
	virtual ~SubGrid_Coupling_GPU();

	virtual void SyncVoltages();
	virtual void SyncCurrents();
	virtual void Interpolate();

	//! Couple on the host mirrors from now on (see Engine_GPU::ForceHostFallback())
	void DropLink();

protected:
	GPU_MultiGridLink* m_Link;   //!< refers to the backends of both grids
};

//! On the SSE field layout (Engine_sse)
class SubGrid_Coupling_SSE : public Engine_SubGridCoupling
{
public:
	SubGrid_Coupling_SSE(Engine_sse* base, Engine_sse* sub, const Operator_CylinderMultiGrid* op);

	virtual void SyncVoltages();
	virtual void SyncCurrents();
	virtual void Interpolate();

protected:
	void InterpolVoltChild2Base(unsigned int rPos);
	void InterpolCurrChild2Base(unsigned int rPos);

	Engine_sse* m_Base;
	Engine_sse* m_Sub;
	MultiGridInterpolation m_Interpol;   //!< the positions; the weights in all four lanes below
	typedef std::vector<f4vector, aligned_allocator<f4vector> > f4weights;
	f4weights m_w_v_2p[2], m_w_v_2pp[2], m_w_i_2p[2], m_w_i_2pp[2];
	unsigned int m_numLines[3];
	unsigned int m_numVectors;
};

#endif // ENGINE_SUBGRID_H
