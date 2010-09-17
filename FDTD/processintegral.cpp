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
#include <iomanip>

ProcessIntegral::ProcessIntegral(Operator* op, Engine* eng)  : Processing(op, eng)
{
	m_TimeShift = 0.0;
	m_Results=NULL;
}

ProcessIntegral::~ProcessIntegral()
{
	delete[] m_Results;
	m_Results = NULL;
	ProcessIntegral::FlushData();
}


void ProcessIntegral::InitProcess()
{
	m_filename = m_Name;
	OpenFile(m_filename);
	FD_Values.clear();
	for (size_t n=0;n<m_FD_Samples.size();++n)
		FD_Values.push_back(0);
}

void ProcessIntegral::FlushData()
{
	if (m_FD_Samples.size())
		Dump_FD_Data(FD_Values,1.0/(double)m_FD_SampleCount,m_filename + "_FD");
}

int ProcessIntegral::Process()
{
	if (Enabled==false) return -1;
	if (CheckTimestep()==false) return GetNextInterval();

	CalcMultipleIntegrals();
	int NrInt = GetNumberOfIntegrals();
	double integral = m_Results[0] * m_weight;

	double time = (double)Eng->GetNumberOfTimesteps()*Op->GetTimestep() + m_TimeShift;

	if (ProcessInterval)
	{
		if (Eng->GetNumberOfTimesteps()%ProcessInterval==0)
		{
			TD_Values.push_back(integral);
			file << setprecision(m_precision) << time;
			for (int n=0;n<NrInt;++n)
				file << "\t" << m_Results[n] * m_weight;
			file << endl;
		}
	}

	if (m_FD_Interval)
	{
		if (Eng->GetNumberOfTimesteps()%m_FD_Interval==0)
		{
			double T = time;
			for (size_t n=0;n<m_FD_Samples.size();++n)
			{
				FD_Values.at(n) += (double)integral * std::exp( -2.0 * II * M_PI * m_FD_Samples.at(n) * T );
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
	if (m_Results==NULL)
		m_Results = new double[1];
	m_Results[0] = CalcIntegral();
	return m_Results;
}
