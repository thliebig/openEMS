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

#include "FDTD/engine.h"
#include "engine_ext_cylindermultigrid.h"
#include "FDTD/engine_cylindermultigrid.h"

Engine_Ext_CylinderMultiGrid::Engine_Ext_CylinderMultiGrid(Operator_Extension* op_ext, bool isBase) : Engine_Extension(op_ext)
{
	m_IsBase = isBase;
	m_Eng_MG = NULL;

	// the multi-grid should be applies last?
	m_Priority = ENG_EXT_PRIO_CYLINDERMULTIGRID;
}

Engine_Ext_CylinderMultiGrid::~Engine_Ext_CylinderMultiGrid()
{
}

void Engine_Ext_CylinderMultiGrid::SetBarrier(boost::barrier* waitBase, boost::barrier* waitChild, boost::barrier* waitSync)
{
	m_WaitOnBase = waitBase;
	m_WaitOnChild = waitChild;
	m_WaitOnSync = waitSync;
}

void Engine_Ext_CylinderMultiGrid::SetEngine(Engine* eng)
{
	m_Eng_MG = dynamic_cast<Engine_CylinderMultiGrid*>(eng);
	if (m_Eng_MG==NULL)
	{
		cerr << "Engine_Ext_CylinderMultiGrid::SetEngine(): Error" << endl;
		exit(0);
	}
}

void Engine_Ext_CylinderMultiGrid::DoPreVoltageUpdates()
{
	//cerr << "Engine_Ext_CylinderMultiGrid::DoPreVoltageUpdates() for " << m_IsBase << endl;
	if (!m_IsBase)
	{
		//cerr << "child: volt wait on base " << endl;
		m_WaitOnBase->wait(); //wait on base to finisch current sync and/or to finisch voltage updates, than start child voltage updates
	}
}

void Engine_Ext_CylinderMultiGrid::DoPostVoltageUpdates()
{

}

void Engine_Ext_CylinderMultiGrid::Apply2Voltages()
{
	if (m_IsBase)
	{
		m_WaitOnBase->wait();	//base voltage updates are done, tell child to start its voltage updates
		m_WaitOnChild->wait();	//wait for child to finisch its updates
		SyncVoltages();			//child is finisch, run sync and go to current updates next
		m_WaitOnSync->wait();	//sync is done... move on and tell child to move on...
	}
	else
	{
		m_WaitOnChild->wait();	//child is finished voltage updates, will tell base to run sync
		m_WaitOnSync->wait();	//wait for base to finisch sync before going to wait for current updates
	}
}

void Engine_Ext_CylinderMultiGrid::SyncVoltages()
{
	if (m_Eng_MG==NULL)
	{
		cerr << "Engine_Ext_CylinderMultiGrid::SyncVoltages: Error engine is NULL" << endl;
		return;
	}

	unsigned int* numLines = m_Eng_MG->numLines;

	Engine_Multithread* m_InnerEng = m_Eng_MG->m_InnerEngine;

	//interpolate voltages from base engine to child engine...
	unsigned int pos[3];
	pos[0] = m_Eng_MG->Op_CMG->GetSplitPos()-1;
	unsigned int pos1_half = 0;
	f4vector v_null;
	v_null.f[0] = 0;
	v_null.f[1] = 0;
	v_null.f[2] = 0;
	v_null.f[3] = 0;
	for (pos[1]=0; pos[1]<numLines[1]-1; pos[1]+=2)
	{
		pos1_half = pos[1]/2;
		for (pos[2]=0; pos[2]<m_Eng_MG->numVectors; ++pos[2])
		{
			//r - direczion
			m_InnerEng->f4_volt[0][pos[0]][pos1_half][pos[2]].v = v_null.v;

			//z - direction
			m_InnerEng->f4_volt[2][pos[0]][pos1_half][pos[2]].v = m_Eng_MG->f4_volt[2][pos[0]][pos[1]][pos[2]].v;

			//alpha - direction
			m_InnerEng->f4_volt[1][pos[0]][pos1_half][pos[2]].v  = m_Eng_MG->f4_volt[1][pos[0]][pos[1]][pos[2]].v;
			m_InnerEng->f4_volt[1][pos[0]][pos1_half][pos[2]].v += m_Eng_MG->f4_volt[1][pos[0]][pos[1]+1][pos[2]].v;
		}
	}
}

void Engine_Ext_CylinderMultiGrid::DoPreCurrentUpdates()
{
	//cerr << "Engine_Ext_CylinderMultiGrid::DoPreCurrentUpdates() for " << m_IsBase << endl;
	if (!m_IsBase)
	{
		//cerr << "child: curr wait on base " << endl;
		m_WaitOnBase->wait(); //wait on base to finisch voltage sync and current updates, than start child current updates
	}
}

void Engine_Ext_CylinderMultiGrid::DoPostCurrentUpdates()
{

}

void Engine_Ext_CylinderMultiGrid::Apply2Current()
{
	if (m_IsBase)
	{
		//cerr << "Base: curr wait on base done, wait on sync" << endl;
		m_WaitOnBase->wait();	//base current updates are done, tell child to start its current updates
		m_WaitOnChild->wait();	//wait for child to finisch its updates
		SyncCurrents();			//child is finisch, run sync and go to voltage updates next
		m_WaitOnSync->wait();	//sync is done... move on and tell child to move on...
	}
	else
	{
		m_WaitOnChild->wait();	//child is finished current updates, will tell base to run sync...
		m_WaitOnSync->wait();	//wait for base to finisch sync before going to wait for next voltage updates
		//cerr << "Child: curr done, wait on sync" << endl;
	}
}

void Engine_Ext_CylinderMultiGrid::SyncCurrents()
{
	if (m_Eng_MG==NULL)
	{
		cerr << "Engine_Ext_CylinderMultiGrid::SyncCurrents: Error engine is NULL" << endl;
		return;
	}

	m_Eng_MG->InterpolCurrChild2Base(m_Eng_MG->Op_CMG->GetSplitPos()-2);
	return;
}
