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

#ifndef USEFUL_H
#define USEFUL_H

#include <vector>

//! Calc the nyquist number of timesteps for a given frequency and timestep
unsigned int CalcNyquistNum(double fmax, double dT);

//! Calc the highest frequency allowed for a given nyquist number of timesteps and timestep
double CalcNyquistFrequency(unsigned int nyquist, double dT);

//! Calculate an optimal job distribution to a given number of threads. Will return a vector with the jobs for each thread.
std::vector<unsigned int> AssignJobs2Threads(unsigned int jobs, unsigned int nrThreads, bool RemoveEmpty=false);

#endif // USEFUL_H
