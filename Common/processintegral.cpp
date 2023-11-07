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

#include "processintegral.h"
#include "Common/operator_base.h"
#include "time.h"
#include <iomanip>

using namespace std;

ProcessIntegral::ProcessIntegral(Engine_Interface_Base* eng_if)  : Processing(eng_if)
{
	m_Results=NULL;
	m_FD_Results=NULL;
	m_normDir = -1;
}

ProcessIntegral::~ProcessIntegral()
{
	delete[] m_Results;
	delete[] m_FD_Results;
	m_Results = NULL;
	m_FD_Results = NULL;
}


void ProcessIntegral::InitProcess()
{
	delete[] m_Results; m_Results = NULL;
	delete[] m_FD_Results; m_FD_Results = NULL;

	if (!Enabled)
		return;

	m_Results = new double[GetNumberOfIntegrals()];
	m_FD_Results = new vector<double_complex>[GetNumberOfIntegrals()];

	m_filename = m_Name;
	OpenFile(m_filename);

	//write header
	time_t rawTime;
	time(&rawTime);
	file << "% time-domain " << GetProcessingName() << " by openEMS " << GIT_VERSION << " @" << ctime(&rawTime);
	file << "% start-coordinates: ("
	<< Op->GetDiscLine(0,start[0])*Op->GetGridDelta() << ","
	<< Op->GetDiscLine(1,start[1])*Op->GetGridDelta() << ","
	<< Op->GetDiscLine(2,start[2])*Op->GetGridDelta() << ") m -> [" << start[0] << "," << start[1] << "," << start[2] << "]" << endl;
	file << "% stop-coordinates: ("
	<< Op->GetDiscLine(0,stop[0])*Op->GetGridDelta() << ","
	<< Op->GetDiscLine(1,stop[1])*Op->GetGridDelta() << ","
	<< Op->GetDiscLine(2,stop[2])*Op->GetGridDelta() << ") m -> [" << stop[0] << "," << stop[1] << "," << stop[2] << "]" << endl;
	file << "% t/s";
	for (int n=0;n<GetNumberOfIntegrals();++n)
	{
		file << "\t" << GetIntegralName(n);
	}
	file << endl;

	for (int i=0;i<GetNumberOfIntegrals();++i)
	{
		for (size_t n=0; n<m_FD_Samples.size(); ++n)
		{
			m_FD_Results[i].push_back(0);
		}
	}
}

void ProcessIntegral::FlushData()
{
	if (!Enabled)
		return;
	if (m_FD_Samples.size())
		Dump_FD_Data(1.0,m_filename + "_FD");
}


void ProcessIntegral::Dump_FD_Data(double factor, string filename)
{
	if (m_FD_Samples.size()==0)
		return;
	ofstream file;
	file.open( filename.c_str() );
	if (!file.is_open())
		cerr << "ProcessIntegral::Dump_FD_Data: Error: Can't open file: " << filename << endl;

	//write header
	time_t rawTime;
	time(&rawTime);
	file << "% frequency-domain " << GetProcessingName() << " by openEMS " << GIT_VERSION << " @" << ctime(&rawTime);
	file << "% start-coordinates: ("
	<< Op->GetDiscLine(0,start[0])*Op->GetGridDelta() << ","
	<< Op->GetDiscLine(1,start[1])*Op->GetGridDelta() << ","
	<< Op->GetDiscLine(2,start[2])*Op->GetGridDelta() << ") m -> [" << start[0] << "," << start[1] << "," << start[2] << "]" << endl;
	file << "% stop-coordinates: ("
	<< Op->GetDiscLine(0,stop[0])*Op->GetGridDelta() << ","
	<< Op->GetDiscLine(1,stop[1])*Op->GetGridDelta() << ","
	<< Op->GetDiscLine(2,stop[2])*Op->GetGridDelta() << ") m -> [" << stop[0] << "," << stop[1] << "," << stop[2] << "]" << endl;
	file << "% f/Hz";
	for (int n=0;n<GetNumberOfIntegrals();++n)
	{
		file << "\t" << GetIntegralName(n) << "\t";
	}
	file << endl << "%";
	for (int i = 0; i < GetNumberOfIntegrals();++i)
		file << "\treal\timag";
	file << endl;

	for (size_t n=0; n<m_FD_Samples.size(); ++n)
	{
		file << m_FD_Samples.at(n) ;
		for (int i = 0; i < GetNumberOfIntegrals();++i)
			file << "\t" << std::real(m_FD_Results[i].at(n))*factor << "\t" << std::imag(m_FD_Results[i].at(n))*factor;
		file << "\n";
	}

	file.close();
}

int ProcessIntegral::Process()
{
	if (Enabled==false) return -1;
	if (CheckTimestep()==false) return GetNextInterval();

	CalcMultipleIntegrals();
	int NrInt = GetNumberOfIntegrals();
	double time = m_Eng_Interface->GetTime(m_dualTime);

	if (ProcessInterval)
	{
		if (m_Eng_Interface->GetNumberOfTimesteps()%ProcessInterval==0)
		{
			file << setprecision(m_precision) << time;
			for (int n=0; n<NrInt; ++n)
				file << "\t" << m_Results[n] * m_weight;
			file << endl;
		}
	}

	if (m_FD_Interval)
	{
		if (m_Eng_Interface->GetNumberOfTimesteps()%m_FD_Interval==0)
		{
			for (size_t n=0; n<m_FD_Samples.size(); ++n)
			{
				for (int i=0; i<NrInt; ++i)
					m_FD_Results[i].at(n) += (double)m_Results[i] * m_weight * std::exp( -2.0 * _I * M_PI * m_FD_Samples.at(n) * time ) * 2.0 * Op->GetTimestep() * (double)m_FD_Interval;
			}
			++m_FD_SampleCount;
			if (m_Flush)
				FlushData();
			m_Flush = false;
		}
	}

	return GetNextInterval();
}

double* ProcessIntegral::CalcMultipleIntegrals()
{
	m_Results[0] = CalcIntegral();
	return m_Results;
}
