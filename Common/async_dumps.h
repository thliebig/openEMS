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
	bool Snapshot(Engine_Interface_Base* eng_if, unsigned int ts, const float* &volt, const float* &curr)
	{
		std::unique_lock<std::mutex> lock(m_Mutex);
		const void* id = eng_if->GetEngineID();
		if ((m_Snap.id==id) && (m_Snap.ts==ts))
		{
			volt = m_Snap.volt;
			curr = m_Snap.curr;
			return true;
		}
		// one snapshot: wait for the tasks still reading the last one
		m_Changed.wait(lock, [&] {return m_Snap.refs==0;});
		m_Snap.id = NULL;
		lock.unlock();
		const bool ok = eng_if->TakeFieldSnapshot(volt, curr);   // on the main thread, which drives the engine
		lock.lock();
		if (!ok)
			return false;
		m_Snap.id = id;
		m_Snap.ts = ts;
		m_Snap.volt = volt;
		m_Snap.curr = curr;
		return true;
	}

	//! Run \a fn in the background for \a owner, reading the snapshot if \a snapshot
	void Push(const void* owner, bool snapshot, const std::function<void()> &fn)
	{
		std::unique_lock<std::mutex> lock(m_Mutex);
		if (!m_Started)
		{
			std::thread(&AsyncDumps::Work, this).detach();
			m_Started = true;
		}
		// one task per dump: tasks without a snapshot hold their field data, so a dump whose
		// previous write is not done yet makes the engine wait instead of queueing more copies
		m_Changed.wait(lock, [&] {return m_Pending[owner]==0;});
		if (snapshot)
			++m_Snap.refs;
		++m_Pending[owner];
		m_Tasks.push_back(Task{owner, snapshot, fn});
		m_Changed.notify_all();
	}

	//! Wait for all tasks of \a owner
	void Wait(const void* owner)
	{
		std::unique_lock<std::mutex> lock(m_Mutex);
		m_Changed.wait(lock, [&] {return m_Pending[owner]==0;});
		m_Pending.erase(owner);
		// no tasks left (e.g. at the end of a run): forget the snapshot, the engine may be
		// deleted and a new one get the same address
		if (m_Tasks.empty() && (m_Snap.refs==0))
			m_Snap = Snap();
	}

protected:
	struct Task
	{
		const void* owner;
		bool snapshot;
		std::function<void()> fn;
	};
	struct Snap
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
			if (task.snapshot)
				--m_Snap.refs;
			--m_Pending[task.owner];
			m_Changed.notify_all();
		}
	}

	std::mutex m_Mutex;
	std::condition_variable m_Changed;
	std::deque<Task> m_Tasks;
	std::map<const void*, int> m_Pending;
	Snap m_Snap;
	bool m_Started = false;
};

#endif // ASYNC_DUMPS_H
