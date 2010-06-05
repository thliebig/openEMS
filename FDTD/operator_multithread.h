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

class Operator_Multithread : public Operator_SSE_Compressed
{
	friend class Operator_Thread;
public:
	//! Create a new operator
	static Operator_Multithread* New(unsigned int numThreads = 0);
	virtual ~Operator_Multithread();

	virtual int CalcECOperator();

	virtual void setNumThreads( unsigned int numThreads );

	virtual Engine* CreateEngine() const;

protected:
	Operator_Multithread();
	virtual void Init();
	virtual void Reset();

	virtual bool Calc_EC(); //this method is using multi-threading

	boost::barrier* m_CalcEC_Start;
	boost::barrier* m_CalcEC_Stop;
	boost::thread_group m_thread_group;
	unsigned int m_numThreads; // number of worker threads
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
