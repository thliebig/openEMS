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

#include <iomanip>
#include "tools/global.h"
#include "process_hfield.h"

ProcessHField::ProcessHField(Operator_Base* op, Engine* eng) : ProcessEField(op, eng)
{
}

ProcessHField::~ProcessHField()
{
}

void ProcessHField::InitProcess()
{
	OpenFile(m_Name);
	for (int n=0; n<3; n++)
		FD_Values[n].assign(m_FD_Samples.size(),double_complex(0.0,0.0));

	file << "% time-domain H-field probe by openEMS " GIT_VERSION << endl;
	file << "% coords: (" << Op->GetDiscLine(0,start[0],true)*Op->GetGridDelta() << ","
	<< Op->GetDiscLine(1,start[1],true)*Op->GetGridDelta() << ","
	<< Op->GetDiscLine(2,start[2],true)*Op->GetGridDelta() << ") m -> [" << start[0] << "," << start[1] << "," << start[2] << "]" << endl;
	file << "% t/s\tEx/(A/m)\tEy/(A/m)\tEz/(A/m)" << endl;
}

void ProcessHField::DefineStartStopCoord(double* dstart, double* dstop)
{
	if (Op->SnapToMesh(dstart,start,true,m_start_inside)==false)
		cerr << "ProcessHField::DefineStartStopCoord: Warning: Snapped line outside field domain!!" << endl;
	if (Op->SnapToMesh(dstop,stop,true,m_stop_inside)==false)
		cerr << "ProcessHField::DefineStartStopCoord: Warning: Snapped line outside field domain!!" << endl;

	if (g_settings.showProbeDiscretization())
	{
		cerr << m_Name << ": snapped coords: (" << Op->GetDiscLine( 0, start[0], true ) << ","
		     << Op->GetDiscLine( 1, start[1], true ) << "," << Op->GetDiscLine( 2, start[2], true ) << ") -> ("
		     << Op->GetDiscLine( 0, stop[0], true ) << ","<< Op->GetDiscLine( 1, stop[1], true ) << ","
		     << Op->GetDiscLine( 2, stop[2], true ) << ")";
		cerr << "   [" << start[0] << "," << start[1] << "," << start[2] << "] -> ["
		     << stop[0] << "," << stop[1] << "," << stop[2] << "]" << endl;
	}
}


int ProcessHField::Process()
{
	if (!Enabled)
		return -1;
	if (CheckTimestep()==false)
		return GetNextInterval();


	if (ProcessInterval)
	{
		// time-domain processing
		if (Eng->GetNumberOfTimesteps()%ProcessInterval==0)
		{
			file << setprecision(m_precision) << ((double)Eng->GetNumberOfTimesteps()+0.5)*Op->GetTimestep();
			for (int n=0; n<3; n++)
			{
				FDTD_FLOAT field = Eng->GetCurr(n,start) / Op->GetMeshDelta(n,start,true);
				field *= m_weight;
//				TD_Values.push_back(voltage);
				file << "\t" << field;
			}
			file << endl;
		}
	}

	if (m_FD_Interval)
	{
		// frequency-domain processing
		if (Eng->GetNumberOfTimesteps()%m_FD_Interval==0)
		{
			double T = (double)Eng->GetNumberOfTimesteps() * Op->GetTimestep();
			for (int pol=0; pol<3; pol++)
			{
				FDTD_FLOAT field = Eng->GetCurr(pol,start) / Op->GetMeshDelta(pol,start,true);
				field *= m_weight;
				for (size_t n=0; n<m_FD_Samples.size(); ++n)
				{
					FD_Values[pol].at(n) += (double)field * std::exp( -2.0 * _I * M_PI * m_FD_Samples.at(n) * T );
				}
				++m_FD_SampleCount;
			}
			if (m_Flush)
				FlushData();
			m_Flush = false;
		}
	}

	return GetNextInterval();
}
