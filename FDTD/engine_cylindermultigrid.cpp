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

#include "engine_cylindermultigrid.h"
#include "operator_cylindermultigrid.h"
#include "extensions/engine_ext_cylindermultigrid.h"

Engine_CylinderMultiGrid* Engine_CylinderMultiGrid::New(const Operator_CylinderMultiGrid* op, unsigned int numThreads)
{
	cout << "Create FDTD engine (cylindrical multi grid mesh using sse compression + multithreading)" << endl;
	Engine_CylinderMultiGrid* e = new Engine_CylinderMultiGrid(op);
	e->setNumThreads( numThreads );
	e->Init();
	return e;
}

Engine_CylinderMultiGrid::Engine_CylinderMultiGrid(const Operator_CylinderMultiGrid* op) : Engine_Cylinder(op)
{
	Op_CMG = op;

	m_WaitOnBase = new boost::barrier(2);
	m_WaitOnChild = new boost::barrier(2);
	m_WaitOnSync = new boost::barrier(2);

	m_Eng_Ext_MG = new Engine_Ext_CylinderMultiGrid(NULL,true);
	m_Eng_Ext_MG->SetBarrier(m_WaitOnBase, m_WaitOnChild, m_WaitOnSync);
	m_Eng_Ext_MG->SetEngine(this);

	Engine* eng = op->GetInnerOperator()->CreateEngine();
	m_InnerEngine = dynamic_cast<Engine_Multithread*>(eng);

	Engine_Ext_CylinderMultiGrid* m_InnerEng_Ext_MG = new Engine_Ext_CylinderMultiGrid(NULL,false);
	m_InnerEng_Ext_MG->SetBarrier(m_WaitOnBase, m_WaitOnChild, m_WaitOnSync);

	// if already has a base extension, switch places ... seems to be faster...
	for (size_t n=0; n<m_InnerEngine->m_Eng_exts.size(); ++n)
	{
		Engine_Ext_CylinderMultiGrid* eng_mg = dynamic_cast<Engine_Ext_CylinderMultiGrid*>(m_InnerEngine->m_Eng_exts.at(n));
		if (eng_mg)
		{
			m_InnerEngine->m_Eng_exts.at(n) = m_InnerEng_Ext_MG;
			m_InnerEng_Ext_MG = eng_mg;
			break;
		}
	}
	m_InnerEngine->m_Eng_exts.push_back(m_InnerEng_Ext_MG);
}

Engine_CylinderMultiGrid::~Engine_CylinderMultiGrid()
{
#ifdef MPI_SUPPORT
	delete m_InnerEngine->m_MPI_Barrier;
	m_InnerEngine->m_MPI_Barrier = NULL;
#endif

	m_Thread_NumTS = 0;
	m_startBarrier->wait();

	m_IteratorThread_Group.join_all();

	delete m_InnerEngine;
	m_InnerEngine = NULL;

	delete m_WaitOnBase;
	m_WaitOnBase = NULL;
	delete m_WaitOnChild;
	m_WaitOnChild = NULL;
	delete m_WaitOnSync;
	m_WaitOnSync = NULL;

	delete m_startBarrier;
	m_startBarrier = NULL;
	delete m_stopBarrier;
	m_stopBarrier = NULL;
}

void Engine_CylinderMultiGrid::Init()
{
	Engine_Multithread::Init();

	m_Eng_exts.push_back(m_Eng_Ext_MG);

	m_startBarrier = new boost::barrier(3); //both engines + organizer
	m_stopBarrier = new boost::barrier(3); //both engines + organizer

	boost::thread *t = NULL;

	t = new boost::thread( Engine_CylinderMultiGrid_Thread(this,m_startBarrier,m_stopBarrier,&m_Thread_NumTS, true) );
	m_IteratorThread_Group.add_thread( t );

	t = new boost::thread( Engine_CylinderMultiGrid_Thread(m_InnerEngine,m_startBarrier,m_stopBarrier,&m_Thread_NumTS, false) );
	m_IteratorThread_Group.add_thread( t );

	m_InnerEngine->SortExtensionByPriority();
	SortExtensionByPriority();

#ifdef MPI_SUPPORT
	//assign an MPI barrier to inner Engine
	m_InnerEngine->m_MPI_Barrier  = new boost::barrier(2);
#endif
}

bool Engine_CylinderMultiGrid::IterateTS(unsigned int iterTS)
{
	m_Thread_NumTS = iterTS;

	m_startBarrier->wait(); //start base and child iterations

	m_stopBarrier->wait();  //tell base and child to wait for another start event...

	//interpolate child data to base mesh...
	for (unsigned int n=0; n<Op_CMG->m_Split_Pos-1; ++n)
		InterpolVoltChild2Base(n);
	for (unsigned int n=0; n<Op_CMG->m_Split_Pos-2; ++n)
		InterpolCurrChild2Base(n);

	return true;
}

void Engine_CylinderMultiGrid::InterpolVoltChild2Base(unsigned int rPos)
{
	//interpolate voltages from child engine to the base engine...
	unsigned int pos[3];
	pos[0] = rPos;
	for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
	{
		for (pos[2]=0; pos[2]<numVectors; ++pos[2])
		{
			//r - direction
			f4_volt[0][pos[0]][pos[1]][pos[2]].v  = Op_CMG->f4_interpol_v_2p[0][pos[1]].v * m_InnerEngine->f4_volt[0][pos[0]][Op_CMG->m_interpol_pos_v_2p[0][pos[1]]][pos[2]].v
													+ Op_CMG->f4_interpol_v_2pp[0][pos[1]].v * m_InnerEngine->f4_volt[0][pos[0]][Op_CMG->m_interpol_pos_v_2pp[0][pos[1]]][pos[2]].v;

			//z - direction
			f4_volt[2][pos[0]][pos[1]][pos[2]].v  = Op_CMG->f4_interpol_v_2p[0][pos[1]].v * m_InnerEngine->f4_volt[2][pos[0]][Op_CMG->m_interpol_pos_v_2p[0][pos[1]]][pos[2]].v
													+ Op_CMG->f4_interpol_v_2pp[0][pos[1]].v * m_InnerEngine->f4_volt[2][pos[0]][Op_CMG->m_interpol_pos_v_2pp[0][pos[1]]][pos[2]].v;

			//alpha - direction
			f4_volt[1][pos[0]][pos[1]][pos[2]].v  = Op_CMG->f4_interpol_v_2p[1][pos[1]].v * m_InnerEngine->f4_volt[1][pos[0]][Op_CMG->m_interpol_pos_v_2p[1][pos[1]]][pos[2]].v
													+ Op_CMG->f4_interpol_v_2pp[1][pos[1]].v * m_InnerEngine->f4_volt[1][pos[0]][Op_CMG->m_interpol_pos_v_2pp[1][pos[1]]][pos[2]].v;
		}
	}
}

void Engine_CylinderMultiGrid::InterpolCurrChild2Base(unsigned int rPos)
{
	//interpolate voltages from child engine to the base engine...
	unsigned int pos[3];
	pos[0] = rPos;
	for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
	{
		for (pos[2]=0; pos[2]<numVectors; ++pos[2])
		{
			//r - direction
			f4_curr[0][pos[0]][pos[1]][pos[2]].v  = Op_CMG->f4_interpol_i_2p[0][pos[1]].v * m_InnerEngine->f4_curr[0][pos[0]][Op_CMG->m_interpol_pos_i_2p[0][pos[1]]][pos[2]].v
													+ Op_CMG->f4_interpol_i_2pp[0][pos[1]].v * m_InnerEngine->f4_curr[0][pos[0]][Op_CMG->m_interpol_pos_i_2pp[0][pos[1]]][pos[2]].v;

			//z - direction
			f4_curr[2][pos[0]][pos[1]][pos[2]].v  = Op_CMG->f4_interpol_i_2p[0][pos[1]].v * m_InnerEngine->f4_curr[2][pos[0]][Op_CMG->m_interpol_pos_i_2p[0][pos[1]]][pos[2]].v
													+ Op_CMG->f4_interpol_i_2pp[0][pos[1]].v * m_InnerEngine->f4_curr[2][pos[0]][Op_CMG->m_interpol_pos_i_2pp[0][pos[1]]][pos[2]].v;

			//alpha - direction
			f4_curr[1][pos[0]][pos[1]][pos[2]].v  = Op_CMG->f4_interpol_i_2p[1][pos[1]].v * m_InnerEngine->f4_curr[1][pos[0]][Op_CMG->m_interpol_pos_i_2p[1][pos[1]]][pos[2]].v
													+ Op_CMG->f4_interpol_i_2pp[1][pos[1]].v * m_InnerEngine->f4_curr[1][pos[0]][Op_CMG->m_interpol_pos_i_2pp[1][pos[1]]][pos[2]].v;
		}
	}
}

#ifdef MPI_SUPPORT
	void Engine_CylinderMultiGrid::SendReceiveVoltages()
	{
		//do the local voltage sync, child is waiting...
		Engine_Multithread::SendReceiveVoltages();

		//run inner voltage sync
		m_InnerEngine->m_MPI_Barrier->wait();
	}

	void Engine_CylinderMultiGrid::SendReceiveCurrents()
	{
		//do the local current sync, child is waiting...
		Engine_Multithread::SendReceiveCurrents();

		//run inner voltage sync
		m_InnerEngine->m_MPI_Barrier->wait();
	}

#endif

/****************************************************************************************/
Engine_CylinderMultiGrid_Thread::Engine_CylinderMultiGrid_Thread( Engine_Multithread* engine, boost::barrier *start, boost::barrier *stop, volatile unsigned int* numTS, bool isBase)
{
	m_startBarrier = start;
	m_stopBarrier = stop;
	m_Eng=engine;
	m_isBase=isBase;
	m_numTS = numTS;
}

void Engine_CylinderMultiGrid_Thread::operator()()
{
	m_startBarrier->wait(); //wait for Base engine to start the iterations...

	while (*m_numTS>0)	//m_numTS==0 request to terminate this thread...
	{
		if (m_isBase)
			m_Eng->Engine_Multithread::IterateTS(*m_numTS);
		else
			m_Eng->IterateTS(*m_numTS);
		m_stopBarrier->wait();		//sync all workers after iterations are performed
		m_startBarrier->wait();		//wait for Base engine to start the iterations again ...
	}
}
