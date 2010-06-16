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

#include "tools/global.h"
#include "processcurrent.h"
#include <iomanip>

ProcessCurrent::ProcessCurrent(Operator* op, Engine* eng) : Processing(op, eng)
{
}

ProcessCurrent::~ProcessCurrent()
{
}

void ProcessCurrent::DefineStartStopCoord(double* dstart, double* dstop)
{
	if (Op->SnapToMesh(dstart,start,true,m_start_inside)==false)
		cerr << "ProcessCurrent::DefineStartStopCoord: Warning: Snapped line outside field domain!!" << endl;
	if (Op->SnapToMesh(dstop,stop,true,m_stop_inside)==false)
		cerr << "ProcessCurrent::DefineStartStopCoord: Warning: Snapped line outside field domain!!" << endl;

	if (g_settings.showProbeDiscretization()) {
		cerr << m_filename << ": snapped coords: (" << Op->GetDiscLine( 0, start[0], true ) << ","
				<< Op->GetDiscLine( 1, start[1], true ) << "," << Op->GetDiscLine( 2, start[2], true ) << ") -> ("
				<< Op->GetDiscLine( 0, stop[0], true ) << ","<< Op->GetDiscLine( 1, stop[1], true ) << ","
				<< Op->GetDiscLine( 2, stop[2], true ) << ")";
		cerr << "   [" << start[0] << "," << start[1] << "," << start[2] << "] -> ["
				<< stop[0] << "," << stop[1] << "," << stop[2] << "]" << endl;
	}
}

int ProcessCurrent::Process()
{
	if (Enabled==false) return -1;
	if (CheckTimestep()==false) return GetNextInterval();
	FDTD_FLOAT current=0;

	for (int n=0;n<3;++n)
	{
		if (start[n]>stop[n])
		{
			unsigned int help=start[n];
			start[n]=stop[n];
			stop[n]=help;
			bool b_help=m_start_inside[n];
			m_start_inside[n] = m_stop_inside[n];
			m_stop_inside[n] = b_help;
		}
		if (m_stop_inside[n]==false) // integrate up to the wall, Operator::SnapToMesh would give numLines[n]-2
			stop[n]=Op->GetNumberOfLines(n)-1;
	}
	//x-current
	if (m_start_inside[1] && m_start_inside[2])
		for (unsigned int i=start[0]+1;i<=stop[0];++i)
			current+=Eng->GetCurr(0,i,start[1],start[2]);
	//y-current
	if (m_stop_inside[0] && m_start_inside[2])
		for (unsigned int i=start[1]+1;i<=stop[1];++i)
			current+=Eng->GetCurr(1,stop[0],i,start[2]);
	//z-current
	if (m_stop_inside[0] && m_stop_inside[1])
		for (unsigned int i=start[2]+1;i<=stop[2];++i)
			current+=Eng->GetCurr(2,stop[0],stop[1],i);
	//x-current
	if (m_stop_inside[1] && m_stop_inside[2])
		for (unsigned int i=start[0]+1;i<=stop[0];++i)
			current-=Eng->GetCurr(0,i,stop[1],stop[2]);
	//y-current
	if (m_start_inside[0] && m_stop_inside[2])
		for (unsigned int i=start[1]+1;i<=stop[1];++i)
			current-=Eng->GetCurr(1,start[0],i,stop[2]);
	//z-current
	if (m_start_inside[0] && m_start_inside[1])
		for (unsigned int i=start[2]+1;i<=stop[2];++i)
			current-=Eng->GetCurr(2,start[0],start[1],i);

//	cerr << "ts: " << Eng->numTS << " i: " << current << endl;
	current*=m_weight;
	v_current.push_back(current);
	//current is sampled half a timestep later then the voltages
	file  << setprecision(m_precision) << (0.5 + (double)Eng->GetNumberOfTimesteps())*Op->GetTimestep() << "\t" << current << endl;
	return GetNextInterval();
}

void ProcessCurrent::DumpBox2File( string vtkfilenameprefix, bool /*dualMesh*/ ) const
{
	Processing::DumpBox2File( vtkfilenameprefix, true );
}
