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

#include "engine_multithread.h"
#include "tools/array_ops.h"

#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/date_time/gregorian/gregorian.hpp"
#include <iomanip>


//! \brief construct an Engine_Multithread instance
//! it's the responsibility of the caller to free the returned pointer
Engine_Multithread* Engine_Multithread::createEngine(Operator* op)
{
	Engine_Multithread* e = new Engine_Multithread(op);
	e->Init();
	return e;
}

Engine_Multithread::Engine_Multithread(Operator* op) : Engine(op)
{
}

Engine_Multithread::~Engine_Multithread()
{
	//DEBUG
	cout << "Engine_Multithread::~Engine_Multithread()" << endl;
	std::map<boost::thread::id, std::vector<double> >::iterator it;
	for (it=m_timer_list.begin(); it!=m_timer_list.end(); it++) {
		std::cout << "*** DEBUG Thread: " << it->first << std::endl;
		std::vector<double>::iterator it2;
		for (it2=it->second.begin(); it2<it->second.end();) {
			std::cout << "after voltage update, before barrier1: " << fixed << setprecision(6) << *(it2++) << std::endl;
			std::cout << "after barrier1, before barrier2: "       << fixed << setprecision(6) << *(it2++) << std::endl;
			std::cout << "after barrier2, before current update: " << fixed << setprecision(6) << *(it2++) << std::endl;
			std::cout << "after current update, before barrier3: " << fixed << setprecision(6) << *(it2++) << std::endl;
			std::cout << "after barrier3: "                        << fixed << setprecision(6) << *(it2++) << std::endl;
		}
	}
	//DEBUG

}

void Engine_Multithread::Init()
{
	Engine::Init();

	numTS = 0;
	m_numTS_times_threads = 0;

	// initialize threads
	m_numThreads = boost::thread::hardware_concurrency();
	std::cout << "using " << m_numThreads << " threads" << std::endl;
	m_barrier1 = new boost::barrier(m_numThreads+1); // numThread workers + 1 excitation thread
	m_barrier2 = new boost::barrier(m_numThreads+1); // numThread workers + 1 excitation thread
	m_barrier3 = new boost::barrier(m_numThreads); // numThread workers
	m_startBarrier = new boost::barrier(m_numThreads+1); // numThread workers + 1 controller
	m_stopBarrier = new boost::barrier(m_numThreads+1); // numThread workers + 1 controller

	for (int n=0; n<m_numThreads; n++) {
		unsigned int linesPerThread = (Op->numLines[0]+m_numThreads-1) / m_numThreads;
		unsigned int start = n * linesPerThread;
		unsigned int stop = min( (n+1) * linesPerThread - 1, Op->numLines[0]-1 );
		//std::cout << "### " << Op->numLines[0] << " " << linesPerThread << " " << start << " " << stop << std::endl;
		boost::thread *t = new boost::thread( thread(this,start,stop) );
		m_thread_group.add_thread( t );
	}
	boost::thread *t = new boost::thread( thread_e_excitation(this) );
	m_thread_group.add_thread( t );
}

void Engine_Multithread::Reset()
{

	Engine::Reset();
}

bool Engine_Multithread::IterateTS(unsigned int iterTS)
{
	m_iterTS = iterTS;
	//cout << "bool Engine_Multithread::IterateTS(): starting threads ...";
	m_startBarrier->wait(); // start the threads
	//cout << "... threads started";

	m_stopBarrier->wait(); // wait for the threads to finish <iterTS> time steps
	numTS = m_numTS_times_threads;
	return true;
}







thread::thread( Engine_Multithread* ptr, unsigned int start, unsigned int stop ) : m_enginePtr(ptr), m_start(start), m_stop(stop), m_stopThread(false)
{
	Op = m_enginePtr->Op;
	volt = m_enginePtr->volt;
	curr = m_enginePtr->curr;
}

void thread::operator()()
{
	//std::cout << "thread::operator() Parameters: " << m_start << " " << m_stop << std::endl;

	unsigned int pos[3];
	bool shift[3];

	while (!m_stopThread) {
		// wait for start
		//cout << "Thread " << boost::this_thread::get_id() << " waiting..." << endl;
		m_enginePtr->m_startBarrier->wait();
		//cout << "Thread " << boost::this_thread::get_id() << " waiting... started." << endl;

		Timer timer1;

		for (unsigned int iter=0;iter<m_enginePtr->m_iterTS;++iter)
		{
			//voltage updates
			for (pos[0]=m_start;pos[0]<=m_stop;++pos[0])
			{
				shift[0]=pos[0];
				for (pos[1]=0;pos[1]<Op->numLines[1];++pos[1])
				{
					shift[1]=pos[1];
					for (pos[2]=0;pos[2]<Op->numLines[2];++pos[2])
					{
						shift[2]=pos[2];
						//do the updates here
						//for x
						volt[0][pos[0]][pos[1]][pos[2]] *= Op->vv[0][pos[0]][pos[1]][pos[2]];
						volt[0][pos[0]][pos[1]][pos[2]] += Op->vi[0][pos[0]][pos[1]][pos[2]] * ( curr[2][pos[0]][pos[1]][pos[2]] - curr[2][pos[0]][pos[1]-shift[1]][pos[2]] - curr[1][pos[0]][pos[1]][pos[2]] + curr[1][pos[0]][pos[1]][pos[2]-shift[2]]);

						//for y
						volt[1][pos[0]][pos[1]][pos[2]] *= Op->vv[1][pos[0]][pos[1]][pos[2]];
						volt[1][pos[0]][pos[1]][pos[2]] += Op->vi[1][pos[0]][pos[1]][pos[2]] * ( curr[0][pos[0]][pos[1]][pos[2]] - curr[0][pos[0]][pos[1]][pos[2]-shift[2]] - curr[2][pos[0]][pos[1]][pos[2]] + curr[2][pos[0]-shift[0]][pos[1]][pos[2]]);

						//for x
						volt[2][pos[0]][pos[1]][pos[2]] *= Op->vv[2][pos[0]][pos[1]][pos[2]];
						volt[2][pos[0]][pos[1]][pos[2]] += Op->vi[2][pos[0]][pos[1]][pos[2]] * ( curr[1][pos[0]][pos[1]][pos[2]] - curr[1][pos[0]-shift[0]][pos[1]][pos[2]] - curr[0][pos[0]][pos[1]][pos[2]] + curr[0][pos[0]][pos[1]-shift[1]][pos[2]]);
					}
				}
			}

			// record time
			m_enginePtr->m_timer_list[boost::this_thread::get_id()].push_back( timer1.elapsed() );

			//cout << "Thread " << boost::this_thread::get_id() << " m_barrier1 waiting..." << endl;
			m_enginePtr->m_barrier1->wait();

			// record time
			m_enginePtr->m_timer_list[boost::this_thread::get_id()].push_back( timer1.elapsed() );

			// e-field excitation (thread thread_e_excitation)

			m_enginePtr->m_barrier2->wait();
			// e_excitation finished

			// record time
			m_enginePtr->m_timer_list[boost::this_thread::get_id()].push_back( timer1.elapsed() );

			//current updates
			for (pos[0]=m_start;pos[0]<=m_stop-1;++pos[0])
			{
				for (pos[1]=0;pos[1]<Op->numLines[1]-1;++pos[1])
				{
					for (pos[2]=0;pos[2]<Op->numLines[2]-1;++pos[2])
					{
						//do the updates here
						//for x
						curr[0][pos[0]][pos[1]][pos[2]] *= Op->ii[0][pos[0]][pos[1]][pos[2]];
						curr[0][pos[0]][pos[1]][pos[2]] += Op->iv[0][pos[0]][pos[1]][pos[2]] * ( volt[2][pos[0]][pos[1]][pos[2]] - volt[2][pos[0]][pos[1]+1][pos[2]] - volt[1][pos[0]][pos[1]][pos[2]] + volt[1][pos[0]][pos[1]][pos[2]+1]);

						//for y
						curr[1][pos[0]][pos[1]][pos[2]] *= Op->ii[1][pos[0]][pos[1]][pos[2]];
						curr[1][pos[0]][pos[1]][pos[2]] += Op->iv[1][pos[0]][pos[1]][pos[2]] * ( volt[0][pos[0]][pos[1]][pos[2]] - volt[0][pos[0]][pos[1]][pos[2]+1] - volt[2][pos[0]][pos[1]][pos[2]] + volt[2][pos[0]+1][pos[1]][pos[2]]);

						//for x
						curr[2][pos[0]][pos[1]][pos[2]] *= Op->ii[2][pos[0]][pos[1]][pos[2]];
						curr[2][pos[0]][pos[1]][pos[2]] += Op->iv[2][pos[0]][pos[1]][pos[2]] * ( volt[1][pos[0]][pos[1]][pos[2]] - volt[1][pos[0]+1][pos[1]][pos[2]] - volt[0][pos[0]][pos[1]][pos[2]] + volt[0][pos[0]][pos[1]+1][pos[2]]);
					}
				}
			}

			// record time
			m_enginePtr->m_timer_list[boost::this_thread::get_id()].push_back( timer1.elapsed() );

			m_enginePtr->m_barrier3->wait();

			// record time
			m_enginePtr->m_timer_list[boost::this_thread::get_id()].push_back( timer1.elapsed() );

			//soft current excitation here (H-field excite)

			++m_enginePtr->m_numTS_times_threads;
		}

		m_enginePtr->m_stopBarrier->wait();
	}
}


thread_e_excitation::thread_e_excitation( Engine_Multithread* ptr ) : m_enginePtr(ptr), m_stopThread(false)
{
	Op = m_enginePtr->Op;
	volt = m_enginePtr->volt;
	curr = m_enginePtr->curr;
}

void thread_e_excitation::operator()()
{
	//std::cout << "thread_e_excitation::operator()" << std::endl;

	while (!m_stopThread) {

		// waiting on thread
		m_enginePtr->m_barrier1->wait();

		int exc_pos;
		unsigned int numTS = m_enginePtr->m_numTS_times_threads / m_enginePtr->m_numThreads;

		//soft voltage excitation here (E-field excite)
		for (unsigned int n=0;n<Op->E_Exc_Count;++n)
		{
			exc_pos = (int)numTS - (int)Op->E_Exc_delay[n];
			exc_pos*= (exc_pos>0 && exc_pos<(int)Op->ExciteLength);
	//			if (n==0) cerr << numTS << " => " << Op->ExciteSignal[exc_pos] << endl;
			volt[Op->E_Exc_dir[n]][Op->E_Exc_index[0][n]][Op->E_Exc_index[1][n]][Op->E_Exc_index[2][n]] += Op->E_Exc_amp[n]*Op->ExciteSignal[exc_pos];
		}

		// continueing thread
		m_enginePtr->m_barrier2->wait();
	}
}
