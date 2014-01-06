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

#ifndef OPERATOR_MULTITHREAD_H
#define OPERATOR_MULTITHREAD_H

#include "operator_sse_compressed.h"

#include <boost/thread.hpp>

#ifdef MPI_SUPPORT
	#define OPERATOR_MULTITHREAD_BASE Operator_MPI
	#include "operator_mpi.h"
#else
	#define OPERATOR_MULTITHREAD_BASE Operator_SSE_Compressed
#endif

class Operator_Multithread : public OPERATOR_MULTITHREAD_BASE
{
	friend class Engine_Multithread;
	friend class Operator_Thread;
public:
	//! Create a new operator
	static Operator_Multithread* New(unsigned int numThreads = 0);
	virtual ~Operator_Multithread();

	virtual void setNumThreads( unsigned int numThreads );

	virtual Engine* CreateEngine();

protected:
	Operator_Multithread();
	virtual void Init();
	void Delete();
	virtual void Reset();

	virtual bool Calc_EC(); //this method is using multi-threading

	unsigned int (*m_Nr_PEC_thread)[3]; //count PEC edges per thread
	virtual bool CalcPEC(); //this method is using multi-threading

	virtual int CalcECOperator( DebugFlags debugFlags = None );

	//Calc_EC barrier
	boost::barrier* m_CalcEC_Start;
	boost::barrier* m_CalcEC_Stop;
	//CalcPEC barrier
	boost::barrier* m_CalcPEC_Start;
	boost::barrier* m_CalcPEC_Stop;

	boost::thread_group m_thread_group;
	unsigned int m_numThreads; // number of worker threads

	//! Calculate the start/stop lines for the multithreading operator and engine.
	/*!
		It depends on the number of threads and number of lines to simulate.
		This method is also used by the multithreading engine!
		This method may also reduce the usable number of thread in case of too few lines or otherwise bad utilization.
	*/
	virtual void CalcStartStopLines(unsigned int &numThreads, vector<unsigned int> &start, vector<unsigned int> &stop) const;
};

class Operator_Thread
{
public:
	Operator_Thread( Operator_Multithread* ptr, unsigned int start, unsigned int stop, unsigned int threadID );
	void operator()();

protected:
	unsigned int m_start, m_stop, m_threadID;
	Operator_Multithread *m_OpPtr;
};


#endif // OPERATOR_MULTITHREAD_H
