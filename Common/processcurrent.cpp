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
#include "FDTD/engine_interface_fdtd.h"
#include <iomanip>

ProcessCurrent::ProcessCurrent(Engine_Interface_Base* eng_if) : ProcessIntegral(eng_if)
{
	m_SnapMethod=1;
}

ProcessCurrent::~ProcessCurrent()
{
}

string ProcessCurrent::GetIntegralName(int row) const
{
	if (row==0)
		return "current";
	return "unknown";
}

void ProcessCurrent::DefineStartStopCoord(double* dstart, double* dstop)
{
	ProcessIntegral::DefineStartStopCoord(dstart, dstop);

	if ((m_normDir<0) || (m_normDir>2))
	{
		if (m_Dimension!=2)
		{
			cerr << "ProcessCurrent::DefineStartStopCoord(): Warning Current Integration Box \"" << m_filename << "\" is not a surface (found dimension: " << m_Dimension << ") nor has it set a valid normal direction! --> disabled" << endl;
			SetEnable(false);
			return;
		}

		for (int n=0; n<3; ++n)
		{
			if (stop[n] == start[n])
				m_normDir = n;
		}
	}
	else
	{
		//expand dimension to 2 if possible
		m_Dimension = 0;
		for (int n=0; n<3; ++n)
		{
			if (n==m_normDir)
				continue;
			if (dstart[n]==dstop[n])
			{
				if ((Op->GetDiscLine( n, start[n], true ) > dstart[n]) && (start[n]>0))
					--start[n];
				if ((Op->GetDiscLine( n, start[n], true ) < dstart[n]) && (stop[n]<Op->GetNumberOfLines(n)-1))
					++stop[n];
			}
			if (stop[n] > start[n])
				++m_Dimension;
		}
	}

	if (start[m_normDir]!=stop[m_normDir])
	{
		cerr << "ProcessCurrent::DefineStartStopCoord(): Warning Current Integration Box \"" << m_filename << "\" has an expansion in normal direction! --> disabled" << endl;
		SetEnable(false);
		return;
	}

	if ((m_normDir<0) || (m_normDir>2))
	{
		cerr << "ProcessCurrent::DefineStartStopCoord(): Warning Current Integration Box \"" << m_filename << "\" has an invalid normal direction --> disabled" << endl;
		SetEnable(false);
		return;
	}
}

double ProcessCurrent::CalcIntegral()
{
	FDTD_FLOAT current=0;

	Engine_Interface_FDTD* EI_FDTD = dynamic_cast<Engine_Interface_FDTD*>(m_Eng_Interface);

	if (EI_FDTD)
	{
		const Engine* Eng = EI_FDTD->GetFDTDEngine();


		switch (m_normDir)
		{
		case 0:
			//y-current
			if (m_stop_inside[0] && m_start_inside[2])
				for (unsigned int i=start[1]+1; i<=stop[1]; ++i)
					current+=Eng->GetCurr(1,stop[0],i,start[2]);
			//z-current
			if (m_stop_inside[0] && m_stop_inside[1])
				for (unsigned int i=start[2]+1; i<=stop[2]; ++i)
					current+=Eng->GetCurr(2,stop[0],stop[1],i);
			//y-current
			if (m_start_inside[0] && m_stop_inside[2])
				for (unsigned int i=start[1]+1; i<=stop[1]; ++i)
					current-=Eng->GetCurr(1,start[0],i,stop[2]);
			//z-current
			if (m_start_inside[0] && m_start_inside[1])
				for (unsigned int i=start[2]+1; i<=stop[2]; ++i)
					current-=Eng->GetCurr(2,start[0],start[1],i);
			break;
		case 1:
			//z-current
			if (m_start_inside[0] && m_start_inside[1])
				for (unsigned int i=start[2]+1; i<=stop[2]; ++i)
					current+=Eng->GetCurr(2,start[0],start[1],i);
			//x-current
			if (m_stop_inside[1] && m_stop_inside[2])
				for (unsigned int i=start[0]+1; i<=stop[0]; ++i)
					current+=Eng->GetCurr(0,i,stop[1],stop[2]);
			//z-current
			if (m_stop_inside[0] && m_stop_inside[1])
				for (unsigned int i=start[2]+1; i<=stop[2]; ++i)
					current-=Eng->GetCurr(2,stop[0],stop[1],i);
			//x-current
			if (m_start_inside[1] && m_start_inside[2])
				for (unsigned int i=start[0]+1; i<=stop[0]; ++i)
					current-=Eng->GetCurr(0,i,start[1],start[2]);
			break;
		case 2:
			//x-current
			if (m_start_inside[1] && m_start_inside[2])
				for (unsigned int i=start[0]+1; i<=stop[0]; ++i)
					current+=Eng->GetCurr(0,i,start[1],start[2]);
			//y-current
			if (m_stop_inside[0] && m_start_inside[2])
				for (unsigned int i=start[1]+1; i<=stop[1]; ++i)
					current+=Eng->GetCurr(1,stop[0],i,start[2]);
			//x-current
			if (m_stop_inside[1] && m_stop_inside[2])
				for (unsigned int i=start[0]+1; i<=stop[0]; ++i)
					current-=Eng->GetCurr(0,i,stop[1],stop[2]);
			//y-current
			if (m_start_inside[0] && m_stop_inside[2])
				for (unsigned int i=start[1]+1; i<=stop[1]; ++i)
					current-=Eng->GetCurr(1,start[0],i,stop[2]);
			break;
		default:
			//this cannot happen...
			return 0.0;
			break;
		}
	}

	return current;
}
