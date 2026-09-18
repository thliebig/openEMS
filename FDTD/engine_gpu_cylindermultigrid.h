/*
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

#ifndef ENGINE_GPU_CYLINDERMULTIGRID_H
#define ENGINE_GPU_CYLINDERMULTIGRID_H

#include "engine_gpu.h"
#include "gpu_backend.h"

class Operator_CylinderMultiGrid;

//! GPU engine for a cylindrical multi-grid level, see Engine_CylinderMultiGrid
/*!
  Owns the engine of the sub-grid (created by the inner operator, so the field
  processing finds it), whose backend shares the device work stream. The CPU
  engine runs both grids in threads synchronized by barriers; the order they
  enforce is run here directly, in every half-step:
  this grid, the sub-grid (recursively with its own sub-grids), then the
  coupling (base voltages to the sub-grid, sub-grid currents to the base).
  At the end of a batch the sub-grid fields are interpolated to the base grid.

  The coupling runs on the device if the backend has a GPU_MultiGridLink,
  otherwise all levels run in the host fallback and the coupling works on the
  host mirrors.
  */
class Engine_GPU_CylinderMultiGrid : public Engine_GPU
{
public:
	static Engine_GPU_CylinderMultiGrid* New(const Operator_CylinderMultiGrid* op, const std::string& backend, GPU_Backend* parent=NULL);
	virtual ~Engine_GPU_CylinderMultiGrid();

	virtual void Init();
	virtual void Reset();

	virtual void VoltageHalfStep();
	virtual void CurrentHalfStep();
	virtual void NextTimestep();
	virtual void FinishBatch();

	virtual void ForceHostFallback();

protected:
	Engine_GPU_CylinderMultiGrid(const Operator_CylinderMultiGrid* op, const std::string& backend, GPU_Backend* parent);

	//! host fallback of the coupling, see GPU_MultiGridLink
	void HostSyncVoltages();
	void HostInterpolate(ArrayLib::ArrayNIJK<FDTD_FLOAT>& base, const ArrayLib::ArrayNIJK<FDTD_FLOAT>& sub, bool currents, unsigned int r);

	const Operator_CylinderMultiGrid* Op_CMG;
	Engine_GPU* m_InnerEngine;
	GPU_MultiGridInterpolation m_Interpol;
	GPU_MultiGridLink* m_Link;
};

#endif // ENGINE_GPU_CYLINDERMULTIGRID_H
