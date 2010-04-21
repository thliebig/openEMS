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
	SetPrecision(12);
	ProcessInterval=0;
}

Processing::~Processing()
{
	file.close();
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
	if (Op->SnapToMesh(dstart,start)==false) cerr << "Processing::DefineStartStopCoord: Warning: Snapped line outside field domain!!" << endl;
	if (Op->SnapToMesh(dstop,stop)==false) cerr << "Processing::DefineStartStopCoord: Warning: Snapped line outside field domain!!" << endl;
}

double Processing::CalcLineIntegral(unsigned int* start, unsigned int* stop, int field) const
{
	switch (field) {
	case 0:
		return CalcLineIntegral_V( start, stop );
	case 1:
		return CalcLineIntegral_I( start, stop );
	}
	return 0;
}

double Processing::CalcLineIntegral_I(unsigned int* start, unsigned int* stop) const
{
	double result=0;
	for (int n=0;n<3;++n)
	{
		if (start[n]<stop[n])
		{
			unsigned int pos[3]={start[0],start[1],start[2]};
			for (;pos[n]<stop[n];++pos[n])
				result += Eng->GetCurr(n,pos[0],pos[1],pos[2]);
		}
		else
		{
			unsigned int pos[3]={stop[0],stop[1],stop[2]};
			for (;pos[n]<start[n];++pos[n])
				result -= Eng->GetCurr(n,pos[0],pos[1],pos[2]);
		}
	}
	return result;
}
double Processing::CalcLineIntegral_V(unsigned int* start, unsigned int* stop) const
{
	double result=0;
	for (int n=0;n<3;++n)
	{
		if (start[n]<stop[n])
		{
			unsigned int pos[3]={start[0],start[1],start[2]};
			for (;pos[n]<stop[n];++pos[n])
				result += Eng->GetVolt(n,pos[0],pos[1],pos[2]);
		}
		else
		{
			unsigned int pos[3]={stop[0],stop[1],stop[2]};
			for (;pos[n]<start[n];++pos[n])
				result -= Eng->GetVolt(n,pos[0],pos[1],pos[2]);
		}
	}
	return result;
}

void Processing::OpenFile( string outfile )
{
	if (file.is_open())
		file.close();

	file.open( outfile.c_str() );
	if (!file.is_open())
		cerr << "Can't open file: " << outfile << endl;

	m_filename = outfile;
}

void Processing::DumpBox2File( string vtkfilenameprefix, bool dualMesh ) const
{
	string vtkfilename = vtkfilenameprefix + m_filename + ".vtk";

	ofstream file( vtkfilename.c_str() );
	if (!file.is_open())
	{
		cerr << "Processing::DumpBoxes2File(): Can't open file: " << vtkfilename << endl;
		return;
	}

	// normalize coordinates
	double s1[3], s2[3];
	for (int i=0; i<3; i++) {
		s1[i] = min(Op->GetDiscLine(i,start[i],dualMesh),Op->GetDiscLine(i,stop[i],dualMesh));
		s2[i] = max(Op->GetDiscLine(i,start[i],dualMesh),Op->GetDiscLine(i,stop[i],dualMesh));
	}

	// fix degenerate box/plane -> line (paraview display problem)
	if (((s1[0] == s2[0]) && (s1[1] == s2[1])) || ((s1[0] == s2[0]) && (s1[2] == s2[2])) || ((s1[2] == s2[2]) && (s1[1] == s2[1]))) {
		// line are not displayed correctly -> enlarge
		for (int i=0; i<3; i++) {
			double delta = min( Op->GetMeshDelta( i, start,dualMesh ), Op->GetMeshDelta( i, stop,dualMesh ) ) / Op->GetGridDelta() / 4.0;
			s1[i] -= delta;
			s2[i] += delta;
		}
	}

	file << "# vtk DataFile Version 2.0" << endl;
	file << "" << endl;
	file << "ASCII" << endl;
	file << "DATASET POLYDATA" << endl;

	file << "POINTS 8 float" << endl;
	file << s1[0] << " " << s1[1] << " " << s1[2] << endl;
	file << s2[0] << " " << s1[1] << " " << s1[2] << endl;
	file << s2[0] << " " << s2[1] << " " << s1[2] << endl;
	file << s1[0] << " " << s2[1] << " " << s1[2] << endl;
	file << s1[0] << " " << s1[1] << " " << s2[2] << endl;
	file << s2[0] << " " << s1[1] << " " << s2[2] << endl;
	file << s2[0] << " " << s2[1] << " " << s2[2] << endl;
	file << s1[0] << " " << s2[1] << " " << s2[2] << endl;

	file << "POLYGONS 6 30" << endl;
	file << "4 0 1 2 3" << endl;
	file << "4 4 5 6 7" << endl;
	file << "4 7 6 2 3" << endl;
	file << "4 4 5 1 0" << endl;
	file << "4 0 4 7 3" << endl;
	file << "4 5 6 2 1" << endl;

	file.close();
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

void ProcessingArray::DumpBoxes2File( string vtkfilenameprefix ) const
{
	for (size_t i=0;i<ProcessArray.size();++i)
		ProcessArray.at(i)->DumpBox2File( vtkfilenameprefix );
}
