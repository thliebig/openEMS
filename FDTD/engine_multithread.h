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

//debug
class Timer {
public:
	Timer() {gettimeofday(&t1,NULL);}
	double elapsed() {gettimeofday(&t2,NULL); return (t2.tv_sec-t1.tv_sec) + (t2.tv_usec-t1.tv_usec)*1e-6;}
protected:
	timeval t1,t2;
};

class Engine_Multithread : public Engine
{
	friend class Processing;
	friend class ProcessVoltage;
	friend class ProcessCurrent;
	friend class ProcessFields;
	friend class ProcessFieldsTD;
	friend class thread;
	friend class thread_e_excitation;
public:
	static Engine_Multithread* createEngine(Operator* op);
	virtual ~Engine_Multithread();

	virtual void Init();
	virtual void Reset();

	//!Iterate a number of timesteps
	virtual bool IterateTS(unsigned int iterTS);

	virtual unsigned int GetNumberOfTimesteps() {return m_numTS_times_threads / m_numThreads;}

protected:
	Engine_Multithread(Operator* op);
	boost::thread_group m_thread_group;
	boost::barrier *m_barrier1, *m_barrier2, *m_barrier3, *m_startBarrier, *m_stopBarrier;
	volatile unsigned int m_iterTS;
	volatile unsigned int m_numTS_times_threads; //!< numTS times the number of worker threads
	unsigned int m_numThreads; //!< number of worker threads

#ifdef ENABLE_DEBUG_TIME
	std::map<boost::thread::id, std::vector<double> > m_timer_list;
#endif
};



class thread {
public:
	thread( Engine_Multithread* ptr, unsigned int start, unsigned int stop );
	void operator()();

protected:
	unsigned int m_start, m_stop;
	volatile bool m_stopThread;
	Engine_Multithread *m_enginePtr;
	Operator *Op;
	FDTD_FLOAT**** volt;
	FDTD_FLOAT**** curr;
};

class thread_e_excitation {
public:
	thread_e_excitation( Engine_Multithread* ptr);
	void operator()();

protected:
	volatile bool m_stopThread;
	Engine_Multithread *m_enginePtr;
	Operator *Op;
	FDTD_FLOAT**** volt;
	FDTD_FLOAT**** curr;
};

#endif // ENGINE_MULTITHREAD_H
