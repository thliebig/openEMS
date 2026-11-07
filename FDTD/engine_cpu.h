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

#ifndef ENGINE_CPU_H
#define ENGINE_CPU_H

#include <algorithm>
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "engine.h"
#include "engine_cpu_kernel.h"
#include "cell_model.h"

//! CPU engine on plain linear field arrays, the coefficients as runs of cell classes
/*!
  The cell model (CellModel) turns the operator into runs of cell classes
  along the fast axis, with the UPML and the dispersive materials fused into
  the update; the engine drops the engine extensions of what the model fuses
  and updates the runs with a kernel per class (engine_cpu_kernels.inc),
  compiled for the generic instruction set and for AVX2 with FMA, chosen at
  start (OPENEMS_CPU_ISA=generic|avx2 overrides).

  The fields are two ArrayNIJK in a layout chosen per model: N slowest, then
  the two slow axes, the fast axis last (stride 1). The fast axis is the one
  with the fewest runs (CellModel::RecommendFastAxis()), unless the box plans
  need another: a closed cylinder's alpha, the local absorbers' common normal
  (OPENEMS_CPU_FAST_AXIS=x|y|z overrides both); the threads split the slow
  axis with more lines (GetSplitAxis()). OPENEMS_CPU_INERT=0 updates the
  inert currents too (CellModel::Analyse()). The engine keeps the state of the runs (UPML flux, ADE
  states) and the fast geometry factor as float for the kernels.

  The extensions see the fields through Engine::GetVolt() and friends, which
  follow the layout, and run in the dispatch of the basic engine.

  Temporal blocking (BuildBlockedPlans()): the threads advance cache-sized
  tiles along the split axis up to 8 timesteps at a time (Engine_BlockedPlan,
  walked by Engine_ThreadPool), the extensions per step on the box of the
  step (Engine_Extension::SupportsBoxes()). It needs every extension to
  support boxes, no cylindrical multi-grid, and a model larger than the cache;
  else the engine runs the flat sweep. The results are the same either way.
  */
class Engine_CPU : public Engine
{
public:
	//! \a modelBuilt is called once the runs are built, before the fields are allocated: the operator frees the coefficients then
	static Engine_CPU* New(const Operator* op, unsigned int numThreads = 0, const std::function<void()>& modelBuilt = nullptr);
	virtual ~Engine_CPU();

	virtual void Init();
	virtual void Reset();

	//! The fast axis (stride 1 in memory)
	int GetFastAxis() const {return m_FastAxis;}

	//! In memory order and on several threads; the same sums in the same order whatever the threads, per line of the split axis
	/*!
	  The interface walked the cells x, y, z, across the layout when the fast
	  axis is x or y: about 9 % of a run with a Gaussian excitation, whose end
	  criteria takes the energy every Nyquist period during the first minute.
	  */
	virtual bool CalcFastEnergy(const unsigned int numNodes[3], double& E_energy, double& H_energy) const;

	//! The bytes left for the weight tables of the field dumps (Field_Gather_FDTD::BuildTable()), shared by all dumps
	/*!
	  As much as the fields and the state of the engine, 256 MiB at least: the
	  tables scale with the model (and so with the machine that runs it), and
	  the same model gets the same tables on any machine.
	  */
	size_t GatherBudget() const {return m_GatherBudget;}
	void UseGatherBudget(size_t bytes) {m_GatherBudget -= std::min(bytes, m_GatherBudget);}
	//! The threads for work at setup: --numThreads, or DefaultThreads()
	unsigned int SetupThreads() const {return (m_numThreads>0) ? m_numThreads : DefaultThreads();}

protected:
	Engine_CPU(const Operator* op);

	using Engine::UpdateVoltages;   // the x ranges of the basic kernels: never called here, the box updates are overridden
	using Engine::UpdateCurrents;
	virtual void UpdateVoltages(const Engine_Box& box);
	virtual void UpdateCurrents(const Engine_Box& box);
	//! Any box: the kernels walk the runs of each line clipped to it
	virtual bool CanUpdateBox(const Engine_Box& box, bool currents) const {(void)box; (void)currents; return true;}
	virtual int GetSplitAxis() const {return m_SplitAxis;}

	//! Trapezoid tiles along the split axis sized for the L3, if every extension supports boxes and there is no sub-grid
	/*!
	  The tile from the estimate of the cache and, as candidates, tiles 4x and 2x,
	  half and a quarter as wide; the thread pool searches them from the estimate
	  and tries the flat sweep at the start of the run (OPENEMS_CPU_TRIAL seconds
	  each, default 0.5, 0 none). OPENEMS_CPU_TBLOCK=0 keeps the flat sweep,
	  =k:W[,k:W...] sets the timesteps per block and the tile widths in lines,
	  without trials.
	  */
	virtual bool BuildBlockedPlans(unsigned int numThreads, Engine_BlockedChoice& choice);
	//! All physical cores, no ramp
	virtual unsigned int DefaultThreads() const;

	virtual std::string GetFieldLayout() const {return m_Layout;}

	//! the axes and the layout: the model's fast axis, the constraints of the box plans, OPENEMS_CPU_FAST_AXIS=x|y|z
	void ChooseLayout();
	//! drop the engine extensions of the operator extensions the model fuses
	void DropFusedExtensions();

	//! the kernel data and the kernels of the instruction set of this CPU, after CellModel::Build()
	void InitKernels();
	Engine_CPU_Kernel m_Kernel;
	Engine_CPU_Update m_UpdateVoltages;
	Engine_CPU_Update m_UpdateCurrents;
	std::string m_ISA;

	std::unique_ptr<CellModel> m_Model;
	std::function<void()> m_ModelBuilt;   //!< see New(), called once
	//! [field][component] the state of the UPML and DISP runs (CellModel_Run::u)
	std::vector<float> m_State[2][3];
	//! [field][component] the fast geometry factor as float, for the kernels
	std::vector<float> m_Gf[2][3];

	int m_FastAxis;
	int m_SplitAxis;     //!< the slowest axis
	int m_MidAxis;       //!< the other slow axis
	std::string m_Layout;
	size_t m_GatherBudget;   //!< see GatherBudget()
};

#endif // ENGINE_CPU_H
