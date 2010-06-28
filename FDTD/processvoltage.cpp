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

#include "processvoltage.h"
#include <complex.h>
#include <iomanip>

ProcessVoltage::ProcessVoltage(Operator* op, Engine* eng) : Processing(op, eng)
{
}

ProcessVoltage::~ProcessVoltage()
{
	ProcessVoltage::FlushData();
}

void ProcessVoltage::InitProcess()
{
	m_filename = m_Name;
	OpenFile(m_filename);
	FD_voltages.clear();
	for (size_t n=0;n<m_FD_Samples.size();++n)
		FD_voltages.push_back(0);
}

int ProcessVoltage::Process()
{
	if (Enabled==false) return -1;
	if (CheckTimestep()==false) return GetNextInterval();

	FDTD_FLOAT voltage=CalcLineIntegral(start,stop,0);
	voltage*=m_weight;

	if (ProcessInterval)
	{
		if (Eng->GetNumberOfTimesteps()%ProcessInterval==0)
		{
			voltages.push_back(voltage);
			file << setprecision(m_precision) << (double)Eng->GetNumberOfTimesteps()*Op->GetTimestep() << "\t" << voltage << endl;
		}
	}

	if (m_FD_Interval)
	{
		if (Eng->GetNumberOfTimesteps()%m_FD_Interval==0)
		{
			double T = (double)Eng->GetNumberOfTimesteps() * Op->GetTimestep();
			for (size_t n=0;n<m_FD_Samples.size();++n)
			{
				FD_voltages.at(n) += voltage * cexp( -2.0 * 1.0i * M_PI * m_FD_Samples.at(n) * T );
			}
			++m_FD_SampleCount;
			if (m_Flush)
					FlushData();
			m_Flush = false;
		}
	}

	return GetNextInterval();
}

void ProcessVoltage::FlushData()
{
	Dump_FD_Data(FD_voltages,1.0/(double)m_FD_SampleCount,m_filename + "_FD");
}
