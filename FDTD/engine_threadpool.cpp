/*
*	Copyright (C) 2010 Sebastian Held (sebastian.held@gmx.de)
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

#include <iostream>
#include <vector>

#include "engine_threadpool.h"
#include "engine.h"
#include "extensions/engine_extension.h"
#include "tools/denormal.h"
#include "tools/useful.h"

using std::cout;
using std::endl;

//! One worker thread of Engine_ThreadPool
class Engine_ThreadPool_Worker
{
public:
	Engine_ThreadPool_Worker(Engine_ThreadPool* pool, unsigned int threadID) : m_pool(pool), m_threadID(threadID) {}
	void operator()();

protected:
	Engine_ThreadPool* m_pool;
	unsigned int m_threadID;
};

void Engine_ThreadPool_Worker::operator()()
{
	// speed up the calculation of denormal floating point values (flush-to-zero)
	Denormal::Disable();

	const size_t numGrids = m_pool->m_Grids.size();
	while (!m_pool->m_stopThreads)
	{
		// wait for start
		m_pool->m_startBarrier->wait();

		if (m_pool->m_stopThreads)
			return;

		for (unsigned int iter=0; iter<m_pool->m_iterTS; ++iter)
		{
			for (size_t g=0; g<numGrids; ++g)
				m_pool->VoltageHalfStep(g, m_threadID);
			m_pool->SyncSubGrids(m_threadID, true);

			for (size_t g=0; g<numGrids; ++g)
				m_pool->CurrentHalfStep(g, m_threadID);
			m_pool->SyncSubGrids(m_threadID, false);

			if (m_threadID == 0)
				for (size_t g=0; g<numGrids; ++g)
					++m_pool->m_Grids[g]->numTS; // only the first thread increments numTS
		}

		m_pool->m_stopBarrier->wait();
	}
}

void Engine_ThreadPool::VoltageHalfStep(size_t g, unsigned int threadID)
{
	Engine* eng = m_Grids[g];
	// pre voltage stuff...
	DoPreVoltageUpdates(eng, threadID);

	//voltage updates
	if (threadID<m_Plans[g].GetNumThreads())
		for (const Engine_Box& box : m_Plans[g].volt[threadID])
			eng->UpdateVoltages(box);

	m_IterateBarrier->wait();

	//post voltage stuff...
	DoPostVoltageUpdates(eng, threadID);
	Apply2Voltages(eng, threadID);
}

void Engine_ThreadPool::CurrentHalfStep(size_t g, unsigned int threadID)
{
	Engine* eng = m_Grids[g];
	//pre current stuff
	DoPreCurrentUpdates(eng, threadID);

	//current updates
	if (threadID<m_Plans[g].GetNumThreads())
		for (const Engine_Box& box : m_Plans[g].curr[threadID])
			eng->UpdateCurrents(box);

	m_IterateBarrier->wait();

	//post current stuff
	DoPostCurrentUpdates(eng, threadID);
	Apply2Current(eng, threadID);
}

void Engine_ThreadPool::SyncSubGrids(unsigned int threadID, bool voltages)
{
	if (m_Grids.size()<2)
		return;
	m_IterateBarrier->wait();   // all grids done
	if (threadID==0)
		for (size_t g=m_Grids.size()-1; g-->0; )
		{
			if (voltages)
				m_Grids[g]->SyncSubGridVoltages();
			else
				m_Grids[g]->SyncSubGridCurrents();
		}
	m_IterateBarrier->wait();   // coupling done
}

// ----------------------------------------------------------------------------

Engine_ThreadPool::Engine_ThreadPool(Engine* eng, unsigned int numThreads)
{
	m_Eng = eng;
	m_thread_group = NULL;
	m_startBarrier = NULL;
	m_stopBarrier = NULL;
	m_IterateBarrier = NULL;
	m_iterTS = 0;
	m_max_numThreads = AvailableThreads();
	m_last_speed = 0;
	m_opt_speed = false;
	m_stopThreads = false;

	for (Engine* g=eng; g; g=g->GetSubGridEngine())
		m_Grids.push_back(g);

	if (numThreads == 0)
	{
		m_opt_speed = true;
		numThreads = 1;
	}
	else if (numThreads > m_max_numThreads)
		numThreads = m_max_numThreads;
	m_numThreads = 0;

	ChangeNumThreads(numThreads);
}

Engine_ThreadPool::~Engine_ThreadPool()
{
	StopThreads();
	delete m_IterateBarrier;
	delete m_startBarrier;
	delete m_stopBarrier;
}

void Engine_ThreadPool::StopThreads()
{
	if (m_thread_group==NULL)
		return;
	m_thread_group->interrupt_all();
	m_thread_group->join_all(); // wait for termination
	delete m_thread_group;
	m_thread_group = NULL;
}

void Engine_ThreadPool::ChangeNumThreads(unsigned int numThreads)
{
	StopThreads();

	m_numThreads = numThreads;

	if (g_settings.GetVerboseLevel()>0)
		cout << "Multithreaded engine using " << m_numThreads << " threads. Utilization: (";

	// the plan of each grid; the utilization reported is the one of the engine itself
	m_Plans.clear();
	for (size_t g=0; g<m_Grids.size(); ++g)
	{
		Engine* eng = m_Grids[g];
		m_Plans.push_back(Engine_Plan::Flat(eng->numLines, m_numThreads, eng->GetSplitAxis()));
		if (g==0)
			m_numThreads = m_Plans[0].GetNumThreads();   // the engine itself may use fewer threads (few lines)
	}
	if (g_settings.GetVerboseLevel()>0)
	{
		const int a = m_Eng->GetSplitAxis();
		for (unsigned int n=0; n<m_numThreads; n++)
		{
			const Engine_Box& box = m_Plans[0].volt[n].front();
			cout << box.stop[a]-box.start[a] << ((n == m_numThreads-1) ? ")" : ";");
		}
	}
	if (g_settings.GetVerboseLevel()>0)
		cout << endl;

	delete m_IterateBarrier;
	// make sure all threads are waiting
	m_IterateBarrier = new boost::barrier(m_numThreads); // numThread workers

	delete m_startBarrier;
	m_startBarrier = new boost::barrier(m_numThreads+1); // numThread workers + 1 controller

	delete m_stopBarrier;
	m_stopBarrier = new boost::barrier(m_numThreads+1); // numThread workers + 1 controller

	m_thread_group = new boost::thread_group();
	for (unsigned int n=0; n<m_numThreads; n++)
	{
		boost::thread *t = new boost::thread( Engine_ThreadPool_Worker(this,n) );
		m_thread_group->add_thread( t );
	}

	for (size_t g=0; g<m_Grids.size(); ++g)
		for (size_t n=0; n<m_Grids[g]->m_Eng_exts.size(); ++n)
			m_Grids[g]->m_Eng_exts.at(n)->SetNumberOfThreads(m_numThreads);
}

void Engine_ThreadPool::IterateTS(unsigned int iterTS)
{
	m_iterTS = iterTS;
	m_startBarrier->wait(); // start the threads
	m_stopBarrier->wait(); // wait for the threads to finish <iterTS> time steps

	// the grids from their sub-grids, the innermost first
	for (size_t g=m_Grids.size()-1; g-->0; )
		m_Grids[g]->InterpolateFromSubGrid();
}

void Engine_ThreadPool::NextInterval(float curr_speed)
{
	if (!m_opt_speed) return;
	if (curr_speed<m_last_speed)
	{
		ChangeNumThreads(m_numThreads-1);
		cout << "Multithreaded Engine: Best performance found using " << m_numThreads << " threads." << std::endl;
		m_opt_speed = false;
	}
	else if (m_numThreads<m_max_numThreads)
	{
		m_last_speed = curr_speed;
		ChangeNumThreads(m_numThreads+1);
	}
}

void Engine_ThreadPool::DoPreVoltageUpdates(Engine* eng, int threadID)
{
	//execute extensions in reverse order -> highest priority gets access to the voltages last
	for (int n=eng->m_Eng_exts.size()-1; n>=0; --n)
	{
		eng->m_Eng_exts.at(n)->DoPreVoltageUpdates(threadID);
		m_IterateBarrier->wait();
	}
}

void Engine_ThreadPool::DoPostVoltageUpdates(Engine* eng, int threadID)
{
	//execute extensions in normal order -> highest priority gets access to the voltages first
	for (size_t n=0; n<eng->m_Eng_exts.size(); ++n)
	{
		eng->m_Eng_exts.at(n)->DoPostVoltageUpdates(threadID);
		m_IterateBarrier->wait();
	}
}

void Engine_ThreadPool::Apply2Voltages(Engine* eng, int threadID)
{
	//execute extensions in normal order -> highest priority gets access to the voltages first
	for (size_t n=0; n<eng->m_Eng_exts.size(); ++n)
	{
		eng->m_Eng_exts.at(n)->Apply2Voltages(threadID);
		m_IterateBarrier->wait();
	}
}

void Engine_ThreadPool::DoPreCurrentUpdates(Engine* eng, int threadID)
{
	//execute extensions in reverse order -> highest priority gets access to the currents last
	for (int n=eng->m_Eng_exts.size()-1; n>=0; --n)
	{
		eng->m_Eng_exts.at(n)->DoPreCurrentUpdates(threadID);
		m_IterateBarrier->wait();
	}
}

void Engine_ThreadPool::DoPostCurrentUpdates(Engine* eng, int threadID)
{
	//execute extensions in normal order -> highest priority gets access to the currents first
	for (size_t n=0; n<eng->m_Eng_exts.size(); ++n)
	{
		eng->m_Eng_exts.at(n)->DoPostCurrentUpdates(threadID);
		m_IterateBarrier->wait();
	}
}

void Engine_ThreadPool::Apply2Current(Engine* eng, int threadID)
{
	//execute extensions in normal order -> highest priority gets access to the currents first
	for (size_t n=0; n<eng->m_Eng_exts.size(); ++n)
	{
		eng->m_Eng_exts.at(n)->Apply2Current(threadID);
		m_IterateBarrier->wait();
	}
}
