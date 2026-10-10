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

#include <complex>
#include <fstream>
#include <vector>
#include "operator.h"

#include "tools/arraylib/array_nijk.h"

namespace NS_Engine_Multithread
{
class thread; // evil hack to access numTS from multithreading context
}

class Engine_Extension;
struct Engine_GatherEntry;

class Engine
{
public:
	enum EngineType
	{
		BASIC, SSE, GPU, UNKNOWN
	};

	static Engine* New(const Operator* op);
	virtual ~Engine();

	virtual void Init();
	virtual void Reset();

	//!Iterate a number of timesteps
	virtual bool IterateTS(unsigned int iterTS);

	virtual unsigned int GetNumberOfTimesteps() {return numTS;}

	virtual void NextInterval(float curr_speed) {UNUSED(curr_speed);};

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
	//! Execute Post-Voltage extension updates
	virtual void DoPostVoltageUpdates();
	//! Apply extension voltage changes
	virtual void Apply2Voltages();

	//! Execute Pre-Current extension updates
	virtual void DoPreCurrentUpdates();
	//! Main FDTD engine current updates
	virtual void UpdateCurrents(unsigned int startX, unsigned int numX);
	//! Execute Post-Current extension updates
	virtual void DoPostCurrentUpdates();
	//! Apply extension current changes
	virtual void Apply2Current();

	//! Calc (roughly) the total energy of the first numNodes[n] nodes in each direction, see Engine_Interface_Base::CalcFastEnergy()
	/*!
	  To be overloaded by an engine with a different storage model, like the field access functions.
	  */
	virtual double CalcFastEnergy(const unsigned int numNodes[3]) const;

	//! Called before several threads read the fields concurrently (GetVolt()/GetCurr()), see Engine_Interface_Base::PrepareFieldAccess()
	virtual void PrepareFieldAccess() {}

	//! Copy the current fields into snapshot \a slot (0 or 1), see Engine_Interface_Base::TakeFieldSnapshot(); false if the engine has no snapshots
	virtual bool SnapshotFields(unsigned int slot, const FDTD_FLOAT* &volt, const FDTD_FLOAT* &curr) {UNUSED(slot); UNUSED(volt); UNUSED(curr); return false;}
	//! Wait until snapshot \a slot can be read, may be called from another thread
	virtual void WaitSnapshot(unsigned int slot) const {UNUSED(slot);}

	//! The fields in the basic engine layout for a direct evaluation of the dumped nodes (see Field_Gather_FDTD), valid after PrepareFieldAccess()
	/*!
	  False for an engine with another storage model, or that does not offer it.
	  */
	virtual bool GetBasicFieldArrays(const ArrayLib::ArrayNIJK<FDTD_FLOAT>* &volt, const ArrayLib::ArrayNIJK<FDTD_FLOAT>* &curr) const {UNUSED(volt); UNUSED(curr); return false;}
	//! Whether the snapshots evaluate the dumped nodes themselves (see AddSnapshotGather()), instead of copying the fields
	virtual bool CanSnapshotGather() const {return false;}
	//! Evaluate the dump entries \a entries of the E (\a h_field false) or H field at every snapshot
	/*!
	  The snapshot then holds their values from \a offset on, in the voltages (E) or currents (H)
	  returned by SnapshotFields(). Only before the first snapshot, returns false otherwise.
	  */
	virtual bool AddSnapshotGather(bool h_field, const std::vector<Engine_GatherEntry>& entries, size_t& offset) {UNUSED(h_field); UNUSED(entries); UNUSED(offset); return false;}

	//! Sums over timesteps of weighted dump entries, kept by the engine, see Engine_Interface_Base::CreateFieldDFT(); -1 if not available
	virtual int AddFieldDFT(bool h_field, const std::vector<Engine_GatherEntry>& entries, unsigned int count) {UNUSED(h_field); UNUSED(entries); UNUSED(count); return -1;}
	virtual void AccumulateFieldDFT(int id, const std::vector<std::complex<float>>& weights) {UNUSED(id); UNUSED(weights);}
	virtual bool ReadFieldDFT(int id, std::vector<std::complex<float>>& sums) {UNUSED(id); UNUSED(sums); return false;}

	inline size_t GetExtensionCount() {return m_Eng_exts.size();}
	inline Engine_Extension* GetExtension(size_t nr) {return m_Eng_exts.at(nr);}
	virtual void SortExtensionByPriority();

	EngineType GetType() const {return m_type;}

protected:
	EngineType m_type;

	Engine(const Operator* op);
	const Operator* Op;

	unsigned int numLines[3];

	ArrayLib::ArrayNIJK<FDTD_FLOAT>* volt_ptr;
	ArrayLib::ArrayNIJK<FDTD_FLOAT>* curr_ptr;
	unsigned int numTS;

	virtual void InitExtensions();
	virtual void ClearExtensions();
	std::vector<Engine_Extension*> m_Eng_exts;

	friend class NS_Engine_Multithread::thread; // evil hack to access numTS from multithreading context
};

#endif // ENGINE_H
