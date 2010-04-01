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

#ifndef ENGINE_MULTITHREAD_H
#define ENGINE_MULTITHREAD_H

#include "operator.h"
#include "engine.h"

#include <boost/thread.hpp>
#include <boost/fusion/include/list.hpp>
#include <boost/fusion/container/list/list_fwd.hpp>
#include <boost/fusion/include/list_fwd.hpp>

class Engine_Multithread;

namespace NS_Engine_Multithread {

	class DBG { // debug
	public:
		DBG() {}
		~DBG() { std::cout << os.str();}
		std::ostringstream& cout() {return os;}
	protected:
		std::ostringstream os;
	};

	class Timer { //debug
	public:
		Timer() {gettimeofday(&t1,NULL);}
		double elapsed() {gettimeofday(&t2,NULL); return (t2.tv_sec-t1.tv_sec) + (t2.tv_usec-t1.tv_usec)*1e-6;}
	protected:
		timeval t1,t2;
	};

	class thread {
	public:
		thread( Engine_Multithread* ptr, unsigned int start, unsigned int stop, unsigned int stop_h, unsigned int threadID );
		void operator()();

	protected:
		unsigned int m_start, m_stop, m_stop_h, m_threadID;
		Engine_Multithread *m_enginePtr;
	};

	class thread_e_excitation {
	public:
		thread_e_excitation( Engine_Multithread* ptr);
		void operator()();

	protected:
		Engine_Multithread *m_enginePtr;
	};
} // namespace


class Engine_Multithread : public Engine
{
	friend class NS_Engine_Multithread::thread;
	friend class NS_Engine_Multithread::thread_e_excitation;
public:
	static Engine_Multithread* createEngine(const Operator* op, unsigned int numThreads = 0);
	virtual ~Engine_Multithread();

	virtual void setNumThreads( unsigned int numThreads );
	virtual void Init();
	virtual void Reset();

	//!Iterate a number of timesteps
	virtual bool IterateTS(unsigned int iterTS);

protected:
	Engine_Multithread(const Operator* op);
	boost::thread_group m_thread_group;
	boost::barrier *m_barrier1, *m_barrier2, *m_barrier3, *m_startBarrier, *m_stopBarrier;
	volatile unsigned int m_iterTS;
	unsigned int m_numThreads; //!< number of worker threads
	volatile bool m_stopThreads;

#ifdef ENABLE_DEBUG_TIME
	std::map<boost::thread::id, std::vector<double> > m_timer_list;
#endif
};

#endif // ENGINE_MULTITHREAD_H
