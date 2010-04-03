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

#include "processing.h"

Processing::Processing(Operator* op, Engine* eng)
{
	Op=op;
	Eng=eng;
	Enabled = true;
	m_PS_pos = 0;
	ProcessInterval=0;
}

Processing::~Processing()
{
}

void Processing::Reset()
{
	m_PS_pos=0;
}

bool Processing::CheckTimestep()
{
	if (m_ProcessSteps.size()>m_PS_pos)
	{
		if (m_ProcessSteps.at(m_PS_pos)==Eng->GetNumberOfTimesteps())
		{
			++m_PS_pos;
			return true;
		}
	}
	if (ProcessInterval)
	{
		if (Eng->GetNumberOfTimesteps()%ProcessInterval==0) return true;
	}
	return false;
}

int Processing::GetNextInterval() const
{
	if (Enabled==false) return -1;
	unsigned int next=-1;
	if (m_ProcessSteps.size()>m_PS_pos)
	{
		next = m_ProcessSteps.at(m_PS_pos)-Eng->GetNumberOfTimesteps();
	}
	if (ProcessInterval==0) return next;
	unsigned int next_Interval = ProcessInterval - Eng->GetNumberOfTimesteps()%ProcessInterval;
	if (next_Interval<next)
		next = next_Interval;
	return next;
}

void Processing::AddStep(unsigned int step)
{
	if (m_ProcessSteps.size()==0)
		m_ProcessSteps.push_back(step);
	else if (find(m_ProcessSteps.begin(), m_ProcessSteps.end(),step)==m_ProcessSteps.end())
		m_ProcessSteps.push_back(step);
}

void Processing::AddSteps(vector<unsigned int> steps)
{
	for (size_t n=0;n<steps.size();++n)
	{
		AddStep(steps.at(n));
	}
}

void Processing::DefineStartStopCoord(double* dstart, double* dstop)
{
	if (Op->SnapToMesh(dstart,start)==false) cerr << "Processing::DefineStartStopCoord: Warning: Snapping problem, check start value!!" << endl;
	if (Op->SnapToMesh(dstop,stop)==false) cerr << "Processing::DefineStartStopCoord: Warning: Snapping problem, check stop value!!" << endl;
}

double Processing::CalcLineIntegral(unsigned int* start, unsigned int* stop, int field)
{
	double result=0;
	FDTD_FLOAT**** array;
	if (field==0)
		array=Eng->GetVoltages();
	else if (field==1)
		array=Eng->GetCurrents();
	else return 0.0;

	for (int n=0;n<3;++n)
	{
		if (start[n]<stop[n])
		{
			unsigned int pos[3]={start[0],start[1],start[2]};
			for (;pos[n]<stop[n];++pos[n])
			{
				result+=array[n][pos[0]][pos[1]][pos[2]];
			}
		}
		else
		{
			unsigned int pos[3]={stop[0],stop[1],stop[2]};
			for (;pos[n]<start[n];++pos[n])
			{
				result-=array[n][pos[0]][pos[1]][pos[2]];
			}
		}
	}
	return result;
}


void ProcessingArray::AddProcessing(Processing* proc)
{
	ProcessArray.push_back(proc);
}

void ProcessingArray::Reset()
{
	for (size_t i=0;i<ProcessArray.size();++i)
	{
		ProcessArray.at(i)->Reset();
	}
}

void ProcessingArray::DeleteAll()
{
	for (size_t i=0;i<ProcessArray.size();++i)
	{
		delete ProcessArray.at(i);
	}
	ProcessArray.clear();
}

int ProcessingArray::Process()
{
	int nextProcess=maxInterval;
	//this could be done nicely in parallel??
	for (size_t i=0;i<ProcessArray.size();++i)
	{
		int step = ProcessArray.at(i)->Process();
		if ((step>0) && (step<nextProcess))
			nextProcess=step;
	}
	return nextProcess;
}
