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

#include "useful.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <climits>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

unsigned int CalcNyquistNum(double fmax, double dT)
{
	if (fmax==0) return UINT_MAX;
	if (dT==0) return 1;
	double T0 = 1/fmax;
	return floor(T0/2/dT);
}

std::vector<unsigned int> AssignJobs2Threads(unsigned int jobs, unsigned int nrThreads, bool RemoveEmpty)
{
	std::vector<unsigned int> jpt; //jobs per thread

	unsigned int ui_jpt = jobs/nrThreads;
	for (unsigned int n=0;n<nrThreads;++n)
	{
		jpt.push_back(ui_jpt);
		jobs-=ui_jpt;
	}

	for (unsigned int n=0;n<nrThreads;++n)
	{
		if (jobs>0)
		{
			++jpt.at(n);
			--jobs;
		}
	}

	if (jobs>0)
		std::cerr << "AssignJobs2Threads: Error, " << jobs << " remain to be assigned, this should not have happend..." << std::endl;

	if (RemoveEmpty)
	{
		while (jpt.back()==0)
			jpt.pop_back();
	}

	return jpt;
}
