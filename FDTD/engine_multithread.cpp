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
#include "tools/array_ops.h"

#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/date_time/gregorian/gregorian.hpp"
#include <iomanip>

//! \brief construct an Engine_Multithread instance
//! it's the responsibility of the caller to free the returned pointer
Engine_Multithread* Engine_Multithread::createEngine(const Operator* op, unsigned int numThreads)
{
	Engine_Multithread* e = new Engine_Multithread(op);
	e->setNumThreads( numThreads );
	e->Init();
	return e;
}

Engine_Multithread::Engine_Multithread(const Operator* op) : Engine(op)
{
}

Engine_Multithread::~Engine_Multithread()
{
#ifdef ENABLE_DEBUG_TIME
	NS_Engine_Multithread::DBG().cout() << "Engine_Multithread::~Engine_Multithread()" << endl;
	std::map<boost::thread::id, std::vector<double> >::iterator it;
	for (it=m_timer_list.begin(); it!=m_timer_list.end(); it++) {
		NS_Engine_Multithread::DBG().cout() << "*** DEBUG Thread: " << it->first << std::endl;
		std::vector<double>::iterator it2;
		for (it2=it->second.begin(); it2<it->second.end();) {
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
	Engine::Init(); // gets cleaned up by Engine::~Engine()

	// initialize threads
	m_stopThreads = false;
	if (m_numThreads == 0)
		m_numThreads = boost::thread::hardware_concurrency();
	cout << "using " << m_numThreads << " threads" << std::endl;
	m_barrier1 = new boost::barrier(m_numThreads+1); // numThread workers + 1 excitation thread
	m_barrier2 = new boost::barrier(m_numThreads+1); // numThread workers + 1 excitation thread
	m_barrier3 = new boost::barrier(m_numThreads); // numThread workers
	m_startBarrier = new boost::barrier(m_numThreads+1); // numThread workers + 1 controller
	m_stopBarrier = new boost::barrier(m_numThreads+1); // numThread workers + 1 controller

	unsigned int linesPerThread = round((float)Op->numLines[0] / (float)m_numThreads);
	for (unsigned int n=0; n<m_numThreads; n++) {
		unsigned int start = n * linesPerThread;
		unsigned int stop = (n+1) * linesPerThread - 1;
		unsigned int stop_h = stop;
		if (n == m_numThreads-1) {
			// last thread
			stop = Op->numLines[0]-1;
			stop_h = stop-1;
		}
		//NS_Engine_Multithread::DBG().cout() << "###DEBUG## Thread " << n << ": start=" << start << " stop=" << stop  << " stop_h=" << stop_h << std::endl;
		boost::thread *t = new boost::thread( NS_Engine_Multithread::thread(this,start,stop,stop_h,n) );
		m_thread_group.add_thread( t );
	}
	boost::thread *t = new boost::thread( NS_Engine_Multithread::thread_e_excitation(this) );
	m_thread_group.add_thread( t );
}

void Engine_Multithread::Reset()
{
	if (!m_stopThreads) {
		// prevent multiple invocations

		// stop the threads
		//NS_Engine_Multithread::DBG().cout() << "stopping all threads" << endl;
		m_iterTS = 1;
		m_startBarrier->wait(); // start the threads
		m_stopThreads = true;
		m_stopBarrier->wait(); // wait for the threads to finish
		m_thread_group.join_all(); // wait for termination
		delete m_barrier1; m_barrier1 = 0;
		delete m_barrier2; m_barrier2 = 0;
		delete m_barrier3; m_barrier3 = 0;
		delete m_startBarrier; m_startBarrier = 0;
		delete m_stopBarrier; m_stopBarrier = 0;
	}

	Engine::Reset();
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

//
// *************************************************************************************************************************
//
namespace NS_Engine_Multithread {

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

	unsigned int pos[3];
	bool shift[3];

	while (!m_enginePtr->m_stopThreads) {
		// wait for start
		//DBG().cout() << "Thread " << m_threadID << " (" << boost::this_thread::get_id() << ") waiting..." << endl;
		m_enginePtr->m_startBarrier->wait();
		//cout << "Thread " << boost::this_thread::get_id() << " waiting... started." << endl;

		DEBUG_TIME( Timer timer1 );

		for (unsigned int iter=0;iter<m_enginePtr->m_iterTS;++iter)
		{
			//voltage updates
			for (pos[0]=m_start;pos[0]<=m_stop;++pos[0])
			{
				shift[0]=pos[0];
				for (pos[1]=0;pos[1]<m_enginePtr->Op->numLines[1];++pos[1])
				{
					shift[1]=pos[1];
					for (pos[2]=0;pos[2]<m_enginePtr->Op->numLines[2];++pos[2])
					{
						shift[2]=pos[2];
						//do the updates here
						//for x
						m_enginePtr->volt[0][pos[0]][pos[1]][pos[2]] *= m_enginePtr->Op->vv[0][pos[0]][pos[1]][pos[2]];
						m_enginePtr->volt[0][pos[0]][pos[1]][pos[2]] += m_enginePtr->Op->vi[0][pos[0]][pos[1]][pos[2]] * ( m_enginePtr->curr[2][pos[0]][pos[1]][pos[2]] - m_enginePtr->curr[2][pos[0]][pos[1]-shift[1]][pos[2]] - m_enginePtr->curr[1][pos[0]][pos[1]][pos[2]] + m_enginePtr->curr[1][pos[0]][pos[1]][pos[2]-shift[2]]);

						//for y
						m_enginePtr->volt[1][pos[0]][pos[1]][pos[2]] *= m_enginePtr->Op->vv[1][pos[0]][pos[1]][pos[2]];
						m_enginePtr->volt[1][pos[0]][pos[1]][pos[2]] += m_enginePtr->Op->vi[1][pos[0]][pos[1]][pos[2]] * ( m_enginePtr->curr[0][pos[0]][pos[1]][pos[2]] - m_enginePtr->curr[0][pos[0]][pos[1]][pos[2]-shift[2]] - m_enginePtr->curr[2][pos[0]][pos[1]][pos[2]] + m_enginePtr->curr[2][pos[0]-shift[0]][pos[1]][pos[2]]);

						//for x
						m_enginePtr->volt[2][pos[0]][pos[1]][pos[2]] *= m_enginePtr->Op->vv[2][pos[0]][pos[1]][pos[2]];
						m_enginePtr->volt[2][pos[0]][pos[1]][pos[2]] += m_enginePtr->Op->vi[2][pos[0]][pos[1]][pos[2]] * ( m_enginePtr->curr[1][pos[0]][pos[1]][pos[2]] - m_enginePtr->curr[1][pos[0]-shift[0]][pos[1]][pos[2]] - m_enginePtr->curr[0][pos[0]][pos[1]][pos[2]] + m_enginePtr->curr[0][pos[0]][pos[1]-shift[1]][pos[2]]);
					}
				}
			}

			// record time
			DEBUG_TIME( m_enginePtr->m_timer_list[boost::this_thread::get_id()].push_back( timer1.elapsed() ); )

			//cout << "Thread " << boost::this_thread::get_id() << " m_barrier1 waiting..." << endl;
			m_enginePtr->m_barrier1->wait();

			// record time
			DEBUG_TIME( m_enginePtr->m_timer_list[boost::this_thread::get_id()].push_back( timer1.elapsed() ); )

			// e-field excitation (thread thread_e_excitation)

			m_enginePtr->m_barrier2->wait();
			// e_excitation finished

			// record time
			DEBUG_TIME( m_enginePtr->m_timer_list[boost::this_thread::get_id()].push_back( timer1.elapsed() ); )

			//current updates
			for (pos[0]=m_start;pos[0]<=m_stop_h;++pos[0])
			{
				for (pos[1]=0;pos[1]<m_enginePtr->Op->numLines[1]-1;++pos[1])
				{
					for (pos[2]=0;pos[2]<m_enginePtr->Op->numLines[2]-1;++pos[2])
					{
						//do the updates here
						//for x
						m_enginePtr->curr[0][pos[0]][pos[1]][pos[2]] *= m_enginePtr->Op->ii[0][pos[0]][pos[1]][pos[2]];
						m_enginePtr->curr[0][pos[0]][pos[1]][pos[2]] += m_enginePtr->Op->iv[0][pos[0]][pos[1]][pos[2]] * ( m_enginePtr->volt[2][pos[0]][pos[1]][pos[2]] - m_enginePtr->volt[2][pos[0]][pos[1]+1][pos[2]] - m_enginePtr->volt[1][pos[0]][pos[1]][pos[2]] + m_enginePtr->volt[1][pos[0]][pos[1]][pos[2]+1]);

						//for y
						m_enginePtr->curr[1][pos[0]][pos[1]][pos[2]] *= m_enginePtr->Op->ii[1][pos[0]][pos[1]][pos[2]];
						m_enginePtr->curr[1][pos[0]][pos[1]][pos[2]] += m_enginePtr->Op->iv[1][pos[0]][pos[1]][pos[2]] * ( m_enginePtr->volt[0][pos[0]][pos[1]][pos[2]] - m_enginePtr->volt[0][pos[0]][pos[1]][pos[2]+1] - m_enginePtr->volt[2][pos[0]][pos[1]][pos[2]] + m_enginePtr->volt[2][pos[0]+1][pos[1]][pos[2]]);

						//for x
						m_enginePtr->curr[2][pos[0]][pos[1]][pos[2]] *= m_enginePtr->Op->ii[2][pos[0]][pos[1]][pos[2]];
						m_enginePtr->curr[2][pos[0]][pos[1]][pos[2]] += m_enginePtr->Op->iv[2][pos[0]][pos[1]][pos[2]] * ( m_enginePtr->volt[1][pos[0]][pos[1]][pos[2]] - m_enginePtr->volt[1][pos[0]+1][pos[1]][pos[2]] - m_enginePtr->volt[0][pos[0]][pos[1]][pos[2]] + m_enginePtr->volt[0][pos[0]][pos[1]+1][pos[2]]);
					}
				}
			}

			// record time
			DEBUG_TIME( m_enginePtr->m_timer_list[boost::this_thread::get_id()].push_back( timer1.elapsed() ); )

			m_enginePtr->m_barrier3->wait();

			// record time
			DEBUG_TIME( m_enginePtr->m_timer_list[boost::this_thread::get_id()].push_back( timer1.elapsed() ); )

			//soft current excitation here (H-field excite)

			if (m_threadID == 0)
				++m_enginePtr->numTS; // only the first thread increments numTS
		}

		m_enginePtr->m_stopBarrier->wait();
	}

	//DBG().cout() << "Thread " << m_threadID << " (" << boost::this_thread::get_id() << ") finished." << endl;
}

} // namespace

//
// *************************************************************************************************************************
//
namespace NS_Engine_Multithread {

thread_e_excitation::thread_e_excitation( Engine_Multithread* ptr )
{
	m_enginePtr = ptr;
}

void thread_e_excitation::operator()()
{
	//std::cout << "thread_e_excitation::operator()" << std::endl;
	//DBG().cout() << "Thread e_excitation (" << boost::this_thread::get_id() << ") started." << endl;

	int exc_pos;
	const unsigned int E_Exc_Count = m_enginePtr->Op->E_Exc_Count;

	while (!m_enginePtr->m_stopThreads)
	{
		// waiting on NS_Engine_Multithread::thread
		m_enginePtr->m_barrier1->wait();

		// soft voltage excitation here (E-field excite)
		for (unsigned int n=0;n<E_Exc_Count;++n)
		{
			exc_pos = (int)m_enginePtr->numTS - (int)m_enginePtr->Op->E_Exc_delay[n];
			exc_pos*= (exc_pos>0 && exc_pos<=(int)m_enginePtr->Op->ExciteLength);
			m_enginePtr->volt[m_enginePtr->Op->E_Exc_dir[n]][m_enginePtr->Op->E_Exc_index[0][n]][m_enginePtr->Op->E_Exc_index[1][n]][m_enginePtr->Op->E_Exc_index[2][n]] += m_enginePtr->Op->E_Exc_amp[n]*m_enginePtr->Op->ExciteSignal[exc_pos];
		}

		// continue NS_Engine_Multithread::thread
		m_enginePtr->m_barrier2->wait();
	}

	//DBG().cout() << "Thread e_excitation (" << boost::this_thread::get_id() << ") finished." << endl;
}

} // namespace
