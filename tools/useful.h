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
#include <string>

//! Calc the nyquist number of timesteps for a given frequency and timestep
unsigned int CalcNyquistNum(double fmax, double dT);

//! Calc the highest frequency allowed for a given nyquist number of timesteps and timestep
double CalcNyquistFrequency(unsigned int nyquist, double dT);

//! Calculate an optimal job distribution to a given number of threads. Will return a vector with the jobs for each thread.
std::vector<unsigned int> AssignJobs2Threads(unsigned int jobs, unsigned int nrThreads, bool RemoveEmpty=false);

std::vector<float> SplitString2Float(std::string str, std::string delimiter=",");
std::vector<double> SplitString2Double(std::string str, std::string delimiter=",");

bool CrossProd(const double* v1, const double* v2, double* out);
double ScalarProd(const double* v1, const double* v2);

double Determinant(const double* mat);
double* Invert(const double* in, double* out);

int LinePlaneIntersection(const double *p0, const double* p1, const double* p2, const double* l_start, const double* l_stop, double* is_point, double &dist);

#ifndef __GNUC__
int gettimeofday(struct timeval* tp, struct timezone* tzp);
#endif // _WIN32

#endif // USEFUL_H
