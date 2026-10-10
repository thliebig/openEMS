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

#ifndef ASYNC_DUMPS_H
#define ASYNC_DUMPS_H

#include <condition_variable>
#include <deque>
#include <functional>
#include <map>
#include <mutex>
#include <thread>

#include "engine_interface_base.h"

//! Background work of the field dumps, one thread in the order of the tasks.
//! - the HDF5 writes of the TD dumps: during the run it is the only user of HDF5 (the
//!   other dumps write in InitProcess()/PostProcess(), after FinishAsync()), so the writes
//!   and the simulation overlap,
//! - with field snapshots of the engine (see Engine_Interface_Base::TakeFieldSnapshot()),
//!   the dumped fields from the snapshot (TD dumps) and the FD sums (FD dumps), while the
//!   engine continues.
class AsyncDumps
{
public:
	static AsyncDumps& Get()
	{
		static AsyncDumps* dumps = new AsyncDumps();   // never deleted: no joining at exit
		return *dumps;
	}

	//! Snapshot of the fields of timestep \a ts, shared by all dumps of that timestep; false if the engine has none
	bool Snapshot(Engine_Interface_Base* eng_if, unsigned int ts, int &slot, const float* &volt, const float* &curr)
	{
		std::unique_lock<std::mutex> lock(m_Mutex);
		const void* id = eng_if->GetEngineID();
		if ((m_Current>=0) && (m_Slots[m_Current].id==id) && (m_Slots[m_Current].ts==ts))
		{
			slot = m_Current;
			volt = m_Slots[slot].volt;
			curr = m_Slots[slot].curr;
			return true;
		}
		// the two slots alternate, wait for the tasks reading the other one
		const int next = (m_Current<0) ? 0 : 1-m_Current;
		m_Changed.wait(lock, [&] {return m_Slots[next].refs==0;});
		m_Slots[next].id = NULL;
		lock.unlock();
		const bool ok = eng_if->TakeFieldSnapshot(next, volt, curr);   // on the main thread, which drives the engine
		lock.lock();
		if (!ok)
			return false;
		m_Slots[next].id = id;
		m_Slots[next].ts = ts;
		m_Slots[next].volt = volt;
		m_Slots[next].curr = curr;
		m_Current = next;
		slot = next;
		return true;
	}

	//! Run \a fn in the background for \a owner, reading snapshot \a slot (-1: none)
	void Push(const void* owner, int slot, const std::function<void()> &fn)
	{
		std::unique_lock<std::mutex> lock(m_Mutex);
		if (!m_Started)
		{
			std::thread(&AsyncDumps::Work, this).detach();
			m_Started = true;
		}
		// bounded: tasks without a snapshot hold their field data
		m_Changed.wait(lock, [&] {return m_Tasks.size()<64;});
		if (slot>=0)
			++m_Slots[slot].refs;
		++m_Pending[owner];
		m_Tasks.push_back(Task{owner, slot, fn});
		m_Changed.notify_all();
	}

	//! Wait for all tasks of \a owner
	void Wait(const void* owner)
	{
		std::unique_lock<std::mutex> lock(m_Mutex);
		m_Changed.wait(lock, [&] {return m_Pending[owner]==0;});
		m_Pending.erase(owner);
		// no tasks left (e.g. at the end of a run): forget the snapshots, the engine may be
		// deleted and a new one get the same address
		if (m_Tasks.empty() && (m_Slots[0].refs==0) && (m_Slots[1].refs==0))
		{
			m_Slots[0] = Slot();
			m_Slots[1] = Slot();
			m_Current = -1;
		}
	}

protected:
	struct Task
	{
		const void* owner;
		int slot;
		std::function<void()> fn;
	};
	struct Slot
	{
		const void* id = NULL;
		unsigned int ts = 0;
		const float* volt = NULL;
		const float* curr = NULL;
		int refs = 0;
	};

	void Work()
	{
		std::unique_lock<std::mutex> lock(m_Mutex);
		while (true)
		{
			m_Changed.wait(lock, [&] {return !m_Tasks.empty();});
			Task task = m_Tasks.front();
			m_Tasks.pop_front();
			m_Changed.notify_all();
			lock.unlock();
			task.fn();
			lock.lock();
			if (task.slot>=0)
				--m_Slots[task.slot].refs;
			--m_Pending[task.owner];
			m_Changed.notify_all();
		}
	}

	std::mutex m_Mutex;
	std::condition_variable m_Changed;
	std::deque<Task> m_Tasks;
	std::map<const void*, int> m_Pending;
	Slot m_Slots[2];
	int m_Current = -1;
	bool m_Started = false;
};

#endif // ASYNC_DUMPS_H
