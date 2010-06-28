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

ProcessIntegral::ProcessIntegral(Operator* op, Engine* eng)  : Processing(op, eng)
{
}

ProcessIntegral::~ProcessIntegral()
{
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

