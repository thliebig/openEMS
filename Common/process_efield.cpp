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
#include "time.h"
#include "process_efield.h"

ProcessEField::ProcessEField(Operator_Base* op, Engine* eng) : Processing(op)
{
	Eng = eng;
}

ProcessEField::~ProcessEField()
{
	FlushData();
}

void ProcessEField::InitProcess()
{
	OpenFile(m_Name);
	for (int n=0; n<3; n++)
		FD_Values[n].assign(m_FD_Samples.size(),double_complex(0.0,0.0));

	file << "% time-domain E-field probe by openEMS " GIT_VERSION << endl;
	file << "% coords: (" << Op->GetDiscLine(0,start[0])*Op->GetGridDelta() << ","
						  << Op->GetDiscLine(1,start[1])*Op->GetGridDelta() << ","
						  << Op->GetDiscLine(2,start[2])*Op->GetGridDelta() << ") m -> [" << start[0] << "," << start[1] << "," << start[2] << "]" << endl;
	file << "% t/s\tEx/(V/m)\tEy/(V/m)\tEz/(V/m)" << endl;
}

void ProcessEField::FlushData()
{
	if (m_FD_Samples.size())
		Dump_FD_Data(FD_Values,1.0/(double)m_FD_SampleCount,m_filename + "_FD");
}

void ProcessEField::Dump_FD_Data(vector<double_complex> value[3], double factor, string filename)
{
	if (value[0].size()==0)
		return;
	if (value[0].size()!=m_FD_Samples.size())
	{
		cerr << "Processing::Dump_FD_Data: Error: Complex value and frequency vector have different size! This should never happend!!!" << endl;
		return;
	}
	ofstream file;
	file.open( filename.c_str() );
	if (!file.is_open())
		cerr << "Can't open file: " << filename << endl;

	time_t rawTime;
	time(&rawTime);
	file << "%dump by openEMS @" << ctime(&rawTime) << "%frequency\treal_x\timag_x\treal_y\timag_y\treal_z\timag_z\n";
	for (size_t n=0;n<value[0].size();++n)
	{
		file << m_FD_Samples.at(n)
			 << "\t" << 2.0 * std::real(value[0].at(n))*factor << "\t" << 2.0 * std::imag(value[0].at(n))*factor
			 << "\t" << 2.0 * std::real(value[1].at(n))*factor << "\t" << 2.0 * std::imag(value[1].at(n))*factor
			 << "\t" << 2.0 * std::real(value[2].at(n))*factor << "\t" << 2.0 * std::imag(value[2].at(n))*factor << "\n";
	}
	file.close();
}

int ProcessEField::Process()
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
			file << setprecision(m_precision) << (double)Eng->GetNumberOfTimesteps()*Op->GetTimestep();
			for (int n=0; n<3; n++)
			{
				FDTD_FLOAT field = Eng->GetVolt(n,start) / Op->GetMeshDelta(n,start);
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
				FDTD_FLOAT field = Eng->GetVolt(pol,start) / Op->GetMeshDelta(pol,start);
				field *= m_weight;
				for (size_t n=0;n<m_FD_Samples.size();++n)
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
