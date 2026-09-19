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

#include "processfields_td.h"
#include "Common/operator_base.h"
#include "tools/vtk_file_writer.h"
#include "tools/hdf5_file_writer.h"
#include <condition_variable>
#include <deque>
#include <functional>
#include <iomanip>
#include <map>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>

using namespace std;

namespace
{
//! Background writer of the HDF5 TD dumps. During the run it is the only user of HDF5
//! (the other dumps write in InitProcess()/PostProcess(), after FinishAsync()), so the
//! writes and the simulation overlap. With field snapshots of the engine (see
//! Engine_Interface_Base::TakeFieldSnapshot()) it also calculates the dumped fields
//! from the snapshot, while the engine continues.
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
}

ProcessFieldsTD::ProcessFieldsTD(Engine_Interface_Base* eng_if) : ProcessFields(eng_if)
{
	pad_length = 8;
	m_AsyncFailed = false;
	m_AsyncUsed = false;
	m_Snapshots = false;
}

ProcessFieldsTD::~ProcessFieldsTD()
{
	FinishAsync();
}

void ProcessFieldsTD::FinishAsync()
{
	if (!m_AsyncUsed)
		return;
	AsyncDumps::Get().Wait(this);
	m_AsyncUsed = false;
	if (m_AsyncFailed)
	{
		SetEnable(false);
		cerr << "ProcessFieldsTD::Process: can't dump to file... disabled! " << endl;
		m_AsyncFailed = false;
	}
}

// on the background writer: calculate the field from the snapshot \a src (if \a field is NULL) and write it
bool ProcessFieldsTD::WriteHDF5(unsigned int ts, float time, ArrayLib::ArrayNIJK<float>* field, const float* src)
{
	ArrayLib::ArrayNIJK<float> calc_field;
	if (!field)
	{
		calc_field.Init("TD_field", numLines);
		if (!CalcField(calc_field, src))
			return false;
		field = &calc_field;
	}
	stringstream ss;
	ss << std::setw( pad_length ) << std::setfill( '0' ) << ts;
	bool success = m_HDF5_Dump_File->WriteVectorField<float>(ss.str(), *field, g_settings.GetLegacyHDF5Dumps());
	float time_attr[1] = {time};
	success &= m_HDF5_Dump_File->WriteAttribute("/FieldData/TD/"+ss.str(),"time",time_attr,1);
	return success;
}

void ProcessFieldsTD::InitProcess()
{
	m_Snapshots = false;
	if (Enabled==false) return;

	ProcessFields::InitProcess();

	// before the first snapshot, e.g. to evaluate the dumped nodes on the device
	if (m_fileType==HDF5_FILETYPE)
	{
		GetGather();
		m_Snapshots = m_Gather && m_Eng_Interface->PrepareSnapshotGather(m_Gather);
	}

	if (m_Vtk_Dump_File)
		m_Vtk_Dump_File->SetHeader(string("openEMS TD Field Dump -- Interpolation: ")+m_Eng_Interface->GetInterpolationTypeString());

	if (m_HDF5_Dump_File)
	{
		m_HDF5_Dump_File->SetCurrentGroup("/FieldData/TD");
		// a dataset every few timesteps: keep the file open, see PostProcess()
		m_HDF5_Dump_File->SetKeepOpen(true);
	}
}

void ProcessFieldsTD::PostProcess()
{
	FinishAsync();
	ProcessFields::PostProcess();
	if (m_HDF5_Dump_File)
		m_HDF5_Dump_File->Close();
}

int ProcessFieldsTD::Process()
{
	if (Enabled==false) return -1;
	if (m_AsyncFailed)
	{
		FinishAsync();   // reports and disables
		return -1;
	}
	if (CheckTimestep()==false) return GetNextInterval();

	string filename = m_filename;

	// HDF5: written in the background, see AsyncDumps; with a snapshot of the fields also calculated there
	if (m_fileType==HDF5_FILETYPE)
	{
		const unsigned int ts = m_Eng_Interface->GetNumberOfTimesteps();
		const float time = (float)m_Eng_Interface->GetTime(m_dualTime);
		int slot = -1;
		const float *volt = NULL, *curr = NULL;
		const float* src = NULL;
		if (m_Snapshots && AsyncDumps::Get().Snapshot(m_Eng_Interface, ts, slot, volt, curr))
			src = (m_DumpType==H_FIELD_DUMP) ? curr : volt;
		ArrayLib::ArrayNIJK<float>* field = NULL;
		if (!src)
		{
			field = new ArrayLib::ArrayNIJK<float>("TD_field", numLines);
			if (!CalcField(*field))
			{
				delete field;
				SetEnable(false);
				cerr << "ProcessFieldsTD::Process: can't dump to file... disabled! " << endl;
				return GetNextInterval();
			}
		}
		m_AsyncUsed = true;
		AsyncDumps::Get().Push(this, slot, [this, ts, time, field, src, slot]()
		{
			if (src)
				m_Eng_Interface->WaitFieldSnapshot(slot);   // the copy may still run on the device
			if (!WriteHDF5(ts, time, field, src))
				m_AsyncFailed = true;
			delete field;
		});
		return GetNextInterval();
	}

	ArrayLib::ArrayNIJK<float> field("TD_field", numLines);
	bool success = CalcField(field);

	if (m_fileType==VTK_FILETYPE)
	{
		m_Vtk_Dump_File->SetTimestep(m_Eng_Interface->GetNumberOfTimesteps());
		m_Vtk_Dump_File->ClearAllFields();
		m_Vtk_Dump_File->AddVectorField(GetFieldNameByType(m_DumpType),field);
		success &= m_Vtk_Dump_File->Write();
	}
	else
	{
		success = false;
		cerr << "ProcessFieldsTD::Process: unknown File-Type" << endl;
	}

	if (success==false)
	{
		SetEnable(false);
		cerr << "ProcessFieldsTD::Process: can't dump to file... disabled! " << endl;
	}

	return GetNextInterval();
}
