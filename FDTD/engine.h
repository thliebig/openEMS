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

#ifndef ENGINE_H
#define ENGINE_H

#include <fstream>
#include "operator.h"
#include "engine_type.h"
#include "engine_plan.h"

#include "tools/arraylib/array_nijk.h"

class Engine_Extension;
class Engine_ThreadPool;
class Engine_SubGridCoupling;

class Engine : public Engine_Types
{
public:
	//! A basic engine on \a numThreads threads, 0 for finding the best number while running (see NextInterval())
	static Engine* New(const Operator* op, unsigned int numThreads = 0);
	virtual ~Engine();

	virtual void Init();
	virtual void Reset();

	//!Iterate a number of timesteps
	virtual bool IterateTS(unsigned int iterTS);

	virtual unsigned int GetNumberOfTimesteps() {return numTS;}

	//! The speed of the last interval, for the thread count of the CPU engines
	virtual void NextInterval(float curr_speed);

	//! Number of threads of a CPU engine, 0 (default) for finding the best number while running; set before the first IterateTS()
	void SetNumThreads(unsigned int numThreads) {m_numThreads = numThreads;}

	//this access functions muss be overloaded by any new engine using a different storage model
	inline virtual FDTD_FLOAT GetVolt(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const
	{
		ArrayLib::ArrayNIJK<FDTD_FLOAT>& volt = *volt_ptr;
		return volt(n,x,y,z);
	}

	inline virtual FDTD_FLOAT GetVolt(unsigned int n, const unsigned int pos[3]) const
	{
		ArrayLib::ArrayNIJK<FDTD_FLOAT>& volt = *volt_ptr;
		return volt(n, pos[0], pos[1], pos[2]);
	}

	inline virtual FDTD_FLOAT GetCurr(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const
	{
		ArrayLib::ArrayNIJK<FDTD_FLOAT>& curr = *curr_ptr;
		return curr(n,x,y,z);
	}

	inline virtual FDTD_FLOAT GetCurr(unsigned int n, const unsigned int pos[3]) const
	{
		ArrayLib::ArrayNIJK<FDTD_FLOAT>& curr = *curr_ptr;
		return curr(n, pos[0], pos[1], pos[2]);
	}

	inline virtual void SetVolt(unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value)
	{
		ArrayLib::ArrayNIJK<FDTD_FLOAT>& volt = *volt_ptr;
		volt(n,x,y,z) = value;
	}

	inline virtual void SetVolt(unsigned int n, const unsigned int pos[3], FDTD_FLOAT value)
	{
		ArrayLib::ArrayNIJK<FDTD_FLOAT>& volt = *volt_ptr;
		volt(n, pos[0], pos[1], pos[2]) = value;
	}

	inline virtual void SetCurr(unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value)
	{
		ArrayLib::ArrayNIJK<FDTD_FLOAT>& curr = *curr_ptr;
		curr(n,x,y,z) = value;
	}

	inline virtual void SetCurr(unsigned int n, const unsigned int pos[3], FDTD_FLOAT value )
	{
		ArrayLib::ArrayNIJK<FDTD_FLOAT>& curr = *curr_ptr;
		curr(n, pos[0], pos[1], pos[2]) = value;
	}

	//! Execute Pre-Voltage extension updates
	virtual void DoPreVoltageUpdates();
	//! Main FDTD engine voltage updates
	virtual void UpdateVoltages(unsigned int startX, unsigned int numX);
	//! The voltage update of \a box (Engine_Plan); the default takes boxes spanning all of y and z, through UpdateVoltages(startX, numX)
	virtual void UpdateVoltages(const Engine_Box& box);
	//! Execute Post-Voltage extension updates
	virtual void DoPostVoltageUpdates();
	//! Apply extension voltage changes
	virtual void Apply2Voltages();

	//! Execute Pre-Current extension updates
	virtual void DoPreCurrentUpdates();
	//! Main FDTD engine current updates
	virtual void UpdateCurrents(unsigned int startX, unsigned int numX);
	//! The current update of \a box, see UpdateVoltages(const Engine_Box&)
	virtual void UpdateCurrents(const Engine_Box& box);

	//! The axis the flat plan splits over the threads: one of the slowest in memory
	virtual int GetSplitAxis() const {return 0;}
	//! Execute Post-Current extension updates
	virtual void DoPostCurrentUpdates();
	//! Apply extension current changes
	virtual void Apply2Current();

	inline size_t GetExtensionCount() {return m_Eng_exts.size();}
	inline Engine_Extension* GetExtension(size_t nr) {return m_Eng_exts.at(nr);}
	virtual void SortExtensionByPriority();

	EngineType GetType() const {return m_type;}

	//! The sub-grid of a cylindrical multi-grid level: its engine (of the same type) and the coupling of the two, both taken over; before the first IterateTS(), throws after
	void SetSubGrid(Engine* sub, Engine_SubGridCoupling* coupling);
	//! The engine of the sub-grid, NULL without one; this engine steps it right after itself in every half step (the threads of Engine_ThreadPool, or Engine_GPU)
	Engine* GetSubGridEngine() const {return m_SubGrid;}
	//! Couple this grid and its sub-grid, once both did their voltage half step
	void SyncSubGridVoltages();
	//! Couple this grid and its sub-grid, once both did their current half step
	void SyncSubGridCurrents();
	//! Fill this grid from its sub-grid, at the end of IterateTS() (the sub-grid filled from its own first)
	void InterpolateFromSubGrid();

protected:
	EngineType m_type;

	Engine(const Operator* op);
	const Operator* Op;

	unsigned int numLines[3];

	ArrayLib::ArrayNIJK<FDTD_FLOAT>* volt_ptr;
	ArrayLib::ArrayNIJK<FDTD_FLOAT>* curr_ptr;
	unsigned int numTS;

	//! The memory layout of the fields (ArrayNIJK), allocated by Init()
	virtual std::string GetFieldLayout() const {return "NIJK";}

	virtual void InitExtensions();
	virtual void ClearExtensions();
	std::vector<Engine_Extension*> m_Eng_exts;

	friend class Engine_ThreadPool;
	friend class Engine_ThreadPool_Worker;
	friend class SubGrid_Coupling_Basic;

	Engine* m_SubGrid;                        //!< see SetSubGrid()
	Engine_SubGridCoupling* m_SubGridCoupling;
	//! Delete the coupling and the sub-grid engine
	void ClearSubGrid();

	unsigned int m_numThreads;        //!< see SetNumThreads()
	Engine_ThreadPool* m_ThreadPool;  //!< the worker threads of a CPU engine, created by the first IterateTS()
};

#endif // ENGINE_H
