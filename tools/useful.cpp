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
#include <iomanip>
#include <sstream>
#include <fstream>
#include <thread>
#ifdef __linux__
#include <sched.h>
#endif
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

namespace
{
#ifdef __linux__
//! true if the comma-separated cgroup v1 controller list contains "name"
bool HasController(const std::string& list, const std::string& name)
{
	for (size_t pos=0; pos<=list.size(); )
	{
		size_t comma = list.find(',', pos);
		size_t len = (comma==std::string::npos) ? std::string::npos : comma-pos;
		if (list.compare(pos, len, name)==0)
			return true;
		if (comma==std::string::npos)
			break;
		pos = comma+1;
	}
	return false;
}

//! This process' cgroup path for a given v1 controller (e.g. "cpu"), or empty if
//! the controller isn't mounted / cgroup v1 isn't in use
std::string OwnCgroupV1Path(const std::string& controller)
{
	std::ifstream self("/proc/self/cgroup");
	std::string line;
	while (std::getline(self, line))
	{
		size_t c1 = line.find(':');
		size_t c2 = (c1==std::string::npos) ? std::string::npos : line.find(':', c1+1);
		if (c2==std::string::npos)
			continue;
		if (HasController(line.substr(c1+1, c2-c1-1), controller))
			return line.substr(c2+1);
	}
	return "";
}

//! CPUs of a cgroup v1 CPU quota (rounded up) for this process' own cgroup, 0 if none
unsigned int CgroupV1Quota()
{
	std::string rel = OwnCgroupV1Path("cpu");
	if (rel.empty())
		return 0;
	// the "cpu" controller may be mounted combined with "cpuacct"
	for (const char* base : {"/sys/fs/cgroup/cpu", "/sys/fs/cgroup/cpu,cpuacct"})
	{
		std::string dir = std::string(base) + rel;
		std::ifstream q1(dir + "/cpu.cfs_quota_us"), p1(dir + "/cpu.cfs_period_us");
		double q = -1, p = 0;
		if ((q1 >> q) && (p1 >> p) && (q>0) && (p>0))
			return (unsigned int)std::ceil(q/p);
	}
	return 0;
}

//! CPUs allowed by one cgroup v2 "cpu.max" file, 0 if unlimited or missing
unsigned int CgroupV2QuotaAt(const std::string& dir)
{
	std::ifstream f(dir + "/cpu.max");
	std::string quota;
	double period = 0;
	if (!(f >> quota >> period))
		return 0;
	if ((quota=="max") || (period<=0))
		return 0;
	return (unsigned int)std::ceil(std::atof(quota.c_str())/period);
}

//! CPUs of a cgroup v2 CPU quota (rounded up), 0 if there is none. cgroup v2 quotas
//! are hierarchical (a parent's quota always constrains its children), so this
//! walks from this process' own cgroup up to the mount root and returns the
//! tightest limit found.
unsigned int CgroupV2Quota()
{
	std::ifstream self("/proc/self/cgroup");
	std::string line, rel;
	while (std::getline(self, line))
	{
		if (line.compare(0, 3, "0::")==0) // unified hierarchy: "0::<path>"
		{
			rel = line.substr(3);
			break;
		}
	}
	if (rel.empty())
		return 0;

	const std::string root = "/sys/fs/cgroup";
	unsigned int quota = 0;
	std::string dir = root + rel;
	while (true)
	{
		unsigned int q = CgroupV2QuotaAt(dir);
		if ((q>0) && ((quota==0) || (q<quota)))
			quota = q;
		if (dir==root)
			break;
		size_t slash = dir.find_last_of('/');
		dir = (slash>root.size()) ? dir.substr(0, slash) : root;
	}
	return quota;
}

//! CPUs of this process' cgroup CPU quota (rounded up), 0 if there is none
unsigned int CgroupCPUQuota()
{
	unsigned int q = CgroupV2Quota();
	return (q>0) ? q : CgroupV1Quota();
}
#endif
}

unsigned int AvailableThreads()
{
	static unsigned int threads = 0;
	if (threads)
		return threads;

	unsigned int n = std::thread::hardware_concurrency();
#ifdef __linux__
	unsigned int visible = n;
	cpu_set_t set;
	if (sched_getaffinity(0, sizeof(set), &set)==0)
	{
		const unsigned int affinity = CPU_COUNT(&set);
		if ((affinity>0) && ((n==0) || (affinity<n)))
			n = affinity;
	}
	const unsigned int quota = CgroupCPUQuota();
	if ((quota>0) && ((n==0) || (quota<n)))
		n = quota;
	if ((n>0) && (n!=visible))
		std::cerr << "Note: limiting the default thread count to " << n << " of "
				   << visible << " visible CPUs (CPU affinity and/or a cgroup CPU quota)"
				   << std::endl;
#endif
	threads = (n>0) ? n : 1;
	return threads;
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
		std::cerr << "AssignJobs2Threads: Error, " << jobs << " remain to be assigned, this should not have happened..." << std::endl;

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

#if defined(_WIN32) && !defined(__GNUC__)
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

#endif // defined(_WIN32) && !defined(__GNUC__)
