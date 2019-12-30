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
#include <iomanip>
#include <sstream>
#include <boost/algorithm/string.hpp>
#include <iostream>

unsigned int CalcNyquistNum(double fmax, double dT)
{
	if (fmax==0) return UINT_MAX;
	if (dT==0) return 1;
	double T0 = 1/fmax;
	return floor(T0/2/dT);
}

double CalcNyquistFrequency(unsigned int nyquist, double dT)
{
	if (nyquist==0) return 0;
	if (dT==0) return 0;
	return floor(1/(double)nyquist/2/dT);
}

std::vector<unsigned int> AssignJobs2Threads(unsigned int jobs, unsigned int nrThreads, bool RemoveEmpty)
{
	std::vector<unsigned int> jpt; //jobs per thread

	unsigned int ui_jpt = jobs/nrThreads;
	for (unsigned int n=0; n<nrThreads; ++n)
	{
		jpt.push_back(ui_jpt);
		jobs-=ui_jpt;
	}

	for (unsigned int n=0; n<nrThreads; ++n)
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

std::vector<float> SplitString2Float(std::string str, std::string delimiter)
{
	std::vector<float> v_f;
	std::vector<std::string> results;
	boost::split(results, str, boost::is_any_of(delimiter));

	for (size_t n=0;n<results.size();++n)
	{
		std::istringstream is(results.at(n));
		float num;
		if (is >> num)
			v_f.push_back(num);
	}
	return v_f;
}

std::vector<double> SplitString2Double(std::string str, std::string delimiter)
{
	std::vector<double> v_f;
	std::vector<std::string> results;
	boost::split(results, str, boost::is_any_of(delimiter));

	for (size_t n=0;n<results.size();++n)
	{
		std::istringstream is(results.at(n));
		double num;
		if (is >> num)
			v_f.push_back(num);
	}
	return v_f;
}

bool CrossProd(const double *v1, const double *v2, double* out)
{
	int nP,nPP;
	for (int n=0;n<3;++n)
	{
		nP = (n+1)%3;
		nPP = (n+2)%3;
		out[n] = v1[nP]*v2[nPP] - v1[nPP]*v2[nP];
	}
	return ((out[0]+out[1]+out[2])>0);
}

double ScalarProd(const double *v1, const double *v2)
{
	double out=0;
	for (int n=0;n<3;++n)
		out+=v1[n]*v2[n];
	return out;
}

double Determinant(const double *mat)
{
	return mat[0]*mat[4]*mat[8]+mat[1]*mat[5]*mat[6]+mat[2]*mat[3]*mat[7]-mat[2]*mat[4]*mat[6]-mat[1]*mat[3]*mat[8]-mat[0]*mat[5]*mat[7];
}

double* Invert(const double* in, double* out)
{
	double det = Determinant(in);
	out[0] = (in[4]*in[8]-in[5]*in[7])/det;
	out[1] = (in[2]*in[7]-in[1]*in[8])/det;
	out[2] = (in[1]*in[5]-in[2]*in[4])/det;
	out[3] = (in[5]*in[6]-in[3]*in[8])/det;
	out[4] = (in[0]*in[8]-in[2]*in[6])/det;
	out[5] = (in[2]*in[3]-in[0]*in[5])/det;
	out[6] = (in[3]*in[7]-in[4]*in[6])/det;
	out[7] = (in[1]*in[6]-in[0]*in[7])/det;
	out[8] = (in[0]*in[4]-in[1]*in[3])/det;
	return out;
}

int LinePlaneIntersection(const double *p0, const double *p1, const double *p2, const double *l_start, const double *l_stop, double* is_point, double &dist)
{
	dist = 0;
	double mat[9];
	for (int n=0;n<3;++n)
	{
		is_point[n] = 0;
		mat[3*n] = l_start[n]-l_stop[n];
		mat[3*n+1] = p1[n]-p0[n];
		mat[3*n+2] = p2[n]-p0[n];
	}
	double det = Determinant(mat);
	if (fabs(det)<1e-50)
		return -1;

	double inv_mat[9];
	Invert(mat, inv_mat);

	double t=0,u=0,v=0;
	for (int n=0;n<3;++n)
	{
		t+=inv_mat[n]*(l_start[n]-p0[n]);
		u+=inv_mat[3+n]*(l_start[n]-p0[n]);
		v+=inv_mat[6+n]*(l_start[n]-p0[n]);
	}
	dist = t;

	for (int n=0;n<3;++n)
		is_point[n] = l_start[n]*(1-dist) + l_stop[n]*dist;

	if ((u<0) || (u>1) || (v<0) || (v>1))
		return 1;
	if ((t<0) || (t>1))
		return 2;

	return 0;
}

#ifndef __GNUC__
#include <chrono>
#include <Winsock2.h> // for struct timeval

int gettimeofday(struct timeval* tp, struct timezone* tzp) {
  namespace sc = std::chrono;
  sc::system_clock::duration d = sc::system_clock::now().time_since_epoch();
  sc::seconds s = sc::duration_cast<sc::seconds>(d);
  tp->tv_sec = s.count();
  tp->tv_usec = sc::duration_cast<sc::microseconds>(d - s).count();

  return 0;
}

#endif // _WIN32