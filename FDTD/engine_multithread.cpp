/*
*	Copyright (C) 2010 Sebastian Held (sebastian.held@gmx.de)
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

//#define ENABLE_DEBUG_TIME

#ifdef ENABLE_DEBUG_TIME
#define DEBUG_TIME(x) x;
#else
#define DEBUG_TIME(x) ;
#endif



#include "engine_multithread.h"
#include "extensions/engine_extension.h"
#include "tools/array_ops.h"

#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/date_time/gregorian/gregorian.hpp"
#include <iomanip>

#ifndef SSE_CORRECT_DENORMALS
#include <xmmintrin.h>
#endif

//! \brief construct an Engine_Multithread instance
//! it's the responsibility of the caller to free the returned pointer
Engine_Multithread* Engine_Multithread::New(const Operator_Multithread* op, unsigned int numThreads)
{
	cout << "Create FDTD engine (compressed SSE + multi-threading)" << endl;
	Engine_Multithread* e = new Engine_Multithread(op);
	e->setNumThreads( numThreads );
	e->Init();
	return e;
}

Engine_Multithread::Engine_Multithread(const Operator_Multithread* op) : ENGINE_MULTITHREAD_BASE(op)
{
	m_Op_MT = op;
	m_type = SSE;
	m_IterateBarrier = 0;
	m_startBarrier = 0;
	m_stopBarrier = 0;

#ifdef ENABLE_DEBUG_TIME
	m_MPI_Barrier = 0;
#endif
}

Engine_Multithread::~Engine_Multithread()
{
#ifdef ENABLE_DEBUG_TIME
	NS_Engine_Multithread::DBG().cout() << "Engine_Multithread::~Engine_Multithread()" << endl;
	std::map<boost::thread::id, std::vector<double> >::iterator it;
	for (it=m_timer_list.begin(); it!=m_timer_list.end(); it++)
	{
		NS_Engine_Multithread::DBG().cout() << "*** DEBUG Thread: " << it->first << std::endl;
		std::vector<double>::iterator it2;
		for (it2=it->second.begin(); it2<it->second.end();)
		{
			NS_Engine_Multithread::DBG().cout() << "after voltage update, before barrier1: " << fixed << setprecision(6) << *(it2++) << std::endl;
			NS_Engine_Multithread::DBG().cout() << "after barrier1, before barrier2: "       << fixed << setprecision(6) << *(it2++) << std::endl;
			NS_Engine_Multithread::DBG().cout() << "after barrier2, before current update: " << fixed << setprecision(6) << *(it2++) << std::endl;
			NS_Engine_Multithread::DBG().cout() << "after current update, before barrier3: " << fixed << setprecision(6) << *(it2++) << std::endl;
			NS_Engine_Multithread::DBG().cout() << "after barrier3: "                        << fixed << setprecision(6) << *(it2++) << std::endl;
		}
	}
#endif

	Reset();
}

void Engine_Multithread::setNumThreads( unsigned int numThreads )
{
	m_numThreads = numThreads;
}

void Engine_Multithread::Init()
{
	m_stopThreads = true;
	ENGINE_MULTITHREAD_BASE::Init();

	// initialize threads
	m_stopThreads = false;
	if (m_numThreads == 0)
		m_numThreads = boost::thread::hardware_concurrency();

	vector<unsigned int> m_Start_Lines;
	vector<unsigned int> m_Stop_Lines;
	m_Op_MT->CalcStartStopLines( m_numThreads, m_Start_Lines, m_Stop_Lines );

	if (g_settings.GetVerboseLevel()>0)
		cout << "Multithreaded engine using " << m_numThreads << " threads. Utilization: (";
	m_IterateBarrier = new boost::barrier(m_numThreads); // numThread workers

	m_startBarrier = new boost::barrier(m_numThreads+1); // numThread workers + 1 controller
	m_stopBarrier = new boost::barrier(m_numThreads+1); // numThread workers + 1 controller
#ifdef MPI_SUPPORT
	m_MPI_Barrier = 0;
#endif

	for (unsigned int n=0; n<m_numThreads; n++)
	{
		unsigned int start = m_Start_Lines.at(n);
		unsigned int stop = m_Stop_Lines.at(n);
		unsigned int stop_h = stop;
		if (n == m_numThreads-1)
		{
			// last thread
			stop_h = stop-1;
			if (g_settings.GetVerboseLevel()>0)
				cout << stop-start+1 << ")" << endl;
		}
		else
			if (g_settings.GetVerboseLevel()>0)
				cout << stop-start+1 << ";";
//		NS_Engine_Multithread::DBG().cout() << "###DEBUG## Thread " << n << ": start=" << start << " stop=" << stop  << " stop_h=" << stop_h << std::endl;
		boost::thread *t = new boost::thread( NS_Engine_Multithread::thread(this,start,stop,stop_h,n) );
		m_thread_group.add_thread( t );
	}

	for (size_t n=0; n<m_Eng_exts.size(); ++n)
		m_Eng_exts.at(n)->SetNumberOfThreads(m_numThreads);
}

void Engine_Multithread::Reset()
{
	if (!m_stopThreads) // prevent multiple invocations
	{
		ClearExtensions(); //prevent extensions from interfering with thread reset...

		// stop the threads
		//NS_Engine_Multithread::DBG().cout() << "stopping all threads" << endl;
		m_iterTS = 1;
		m_startBarrier->wait(); // start the threads
		m_stopThreads = true;
		m_stopBarrier->wait(); // wait for the threads to finish
		m_thread_group.join_all(); // wait for termination
		delete m_IterateBarrier;
		m_IterateBarrier = 0;
		delete m_startBarrier;
		m_startBarrier = 0;
		delete m_stopBarrier;
		m_stopBarrier = 0;
	}

	ENGINE_MULTITHREAD_BASE::Reset();
}

bool Engine_Multithread::IterateTS(unsigned int iterTS)
{
	m_iterTS = iterTS;

	//cout << "bool Engine_Multithread::IterateTS(): starting threads ...";
	m_startBarrier->wait(); // start the threads

	//cout << "... threads started";

	m_stopBarrier->wait(); // wait for the threads to finish <iterTS> time steps
	return true;
}

void Engine_Multithread::DoPreVoltageUpdates(int threadID)
{
	//execute extensions in reverse order -> highest priority gets access to the voltages last
	for (int n=m_Eng_exts.size()-1; n>=0; --n)
	{
		m_Eng_exts.at(n)->DoPreVoltageUpdates(threadID);
		m_IterateBarrier->wait();
	}

}

void Engine_Multithread::DoPostVoltageUpdates(int threadID)
{
	//execute extensions in normal order -> highest priority gets access to the voltages first
	for (size_t n=0; n<m_Eng_exts.size(); ++n)
	{
		m_Eng_exts.at(n)->DoPostVoltageUpdates(threadID);
		m_IterateBarrier->wait();
	}
}

void Engine_Multithread::Apply2Voltages(int threadID)
{
	//execute extensions in normal order -> highest priority gets access to the voltages first
	for (size_t n=0; n<m_Eng_exts.size(); ++n)
	{
		m_Eng_exts.at(n)->Apply2Voltages(threadID);
		m_IterateBarrier->wait();
	}
}

void Engine_Multithread::DoPreCurrentUpdates(int threadID)
{
	//execute extensions in reverse order -> highest priority gets access to the currents last
	for (int n=m_Eng_exts.size()-1; n>=0; --n)
	{
		m_Eng_exts.at(n)->DoPreCurrentUpdates(threadID);
		m_IterateBarrier->wait();
	}
}

void Engine_Multithread::DoPostCurrentUpdates(int threadID)
{
	//execute extensions in normal order -> highest priority gets access to the currents first
	for (size_t n=0; n<m_Eng_exts.size(); ++n)
	{
		m_Eng_exts.at(n)->DoPostCurrentUpdates(threadID);
		m_IterateBarrier->wait();
	}
}

void Engine_Multithread::Apply2Current(int threadID)
{
	//execute extensions in normal order -> highest priority gets access to the currents first
	for (size_t n=0; n<m_Eng_exts.size(); ++n)
	{
		m_Eng_exts.at(n)->Apply2Current(threadID);
		m_IterateBarrier->wait();
	}
}

//
// *************************************************************************************************************************
//
namespace NS_Engine_Multithread
{

thread::thread( Engine_Multithread* ptr, unsigned int start, unsigned int stop, unsigned int stop_h, unsigned int threadID )
{
	m_enginePtr = ptr;
	m_start = start;
	m_stop = stop;
	m_stop_h = stop_h;
	m_threadID = threadID;
}

void thread::operator()()
{
	//std::cout << "thread::operator() Parameters: " << m_start << " " << m_stop << std::endl;
	//DBG().cout() << "Thread " << m_threadID << " (" << boost::this_thread::get_id() << ") started." << endl;

	// speed up the calculation of denormal floating point values (flush-to-zero)
#ifndef SSE_CORRECT_DENORMALS
	unsigned int oldMXCSR = _mm_getcsr(); //read the old MXCSR setting
	unsigned int newMXCSR = oldMXCSR | 0x8040; // set DAZ and FZ bits
	_mm_setcsr( newMXCSR ); //write the new MXCSR setting to the MXCSR
#endif

	while (!m_enginePtr->m_stopThreads)
	{
		// wait for start
		//DBG().cout() << "Thread " << m_threadID << " (" << boost::this_thread::get_id() << ") waiting..." << endl;
		m_enginePtr->m_startBarrier->wait();
		//cout << "Thread " << boost::this_thread::get_id() << " waiting... started." << endl;

		DEBUG_TIME( Timer timer1 );

		for (unsigned int iter=0; iter<m_enginePtr->m_iterTS; ++iter)
		{
			// pre voltage stuff...
			m_enginePtr->DoPreVoltageUpdates(m_threadID);

			//voltage updates
			m_enginePtr->UpdateVoltages(m_start,m_stop-m_start+1);

			// record time
			DEBUG_TIME( m_enginePtr->m_timer_list[boost::this_thread::get_id()].push_back( timer1.elapsed() ); )

			//cout << "Thread " << boost::this_thread::get_id() << " m_barrier1 waiting..." << endl;
			m_enginePtr->m_IterateBarrier->wait();

			// record time
			DEBUG_TIME( m_enginePtr->m_timer_list[boost::this_thread::get_id()].push_back( timer1.elapsed() ); )

			//post voltage stuff...
			m_enginePtr->DoPostVoltageUpdates(m_threadID);
			m_enginePtr->Apply2Voltages(m_threadID);

#ifdef MPI_SUPPORT
			if (m_threadID==0)
			{
				if (m_enginePtr->m_MPI_Barrier)
					m_enginePtr->m_MPI_Barrier->wait();
				m_enginePtr->SendReceiveVoltages();
			}
			m_enginePtr->m_IterateBarrier->wait();
#endif

			// record time
			DEBUG_TIME( m_enginePtr->m_timer_list[boost::this_thread::get_id()].push_back( timer1.elapsed() ); )

			//pre current stuff
			m_enginePtr->DoPreCurrentUpdates(m_threadID);

			//current updates
			m_enginePtr->UpdateCurrents(m_start,m_stop_h-m_start+1);

			// record time
			DEBUG_TIME( m_enginePtr->m_timer_list[boost::this_thread::get_id()].push_back( timer1.elapsed() ); )
			m_enginePtr->m_IterateBarrier->wait();

			// record time
			DEBUG_TIME( m_enginePtr->m_timer_list[boost::this_thread::get_id()].push_back( timer1.elapsed() ); )

			//post current stuff
			m_enginePtr->DoPostCurrentUpdates(m_threadID);
			m_enginePtr->Apply2Current(m_threadID);

#ifdef MPI_SUPPORT
			if (m_threadID==0)
			{
				if (m_enginePtr->m_MPI_Barrier)
					m_enginePtr->m_MPI_Barrier->wait();
				m_enginePtr->SendReceiveCurrents();
			}
			m_enginePtr->m_IterateBarrier->wait();
#endif

			if (m_threadID == 0)
				++m_enginePtr->numTS; // only the first thread increments numTS
		}

		m_enginePtr->m_stopBarrier->wait();
	}

	//DBG().cout() << "Thread " << m_threadID << " (" << boost::this_thread::get_id() << ") finished." << endl;
}

} // namespace

