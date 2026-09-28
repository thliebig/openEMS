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
#include <exception>
#include <functional>
#include <thread>

//! Calc the nyquist number of timesteps for a given frequency and timestep
unsigned int CalcNyquistNum(double fmax, double dT);

//! Calc the highest frequency allowed for a given nyquist number of timesteps and timestep
double CalcNyquistFrequency(unsigned int nyquist, double dT);

//! Number of threads this process can usefully run in parallel: the visible CPUs,
//! limited by the CPU affinity and a cgroup CPU quota (e.g. a container or systemd
//! unit with a CPU limit). Falls back to hardware_concurrency() where none of that
//! applies or cannot be determined (non-Linux).
unsigned int AvailableThreads();

//! Calculate an optimal job distribution to a given number of threads. Will return a vector with the jobs for each thread.
std::vector<unsigned int> AssignJobs2Threads(unsigned int jobs, unsigned int nrThreads, bool RemoveEmpty=false);

//! Split [0, count) into \a workers contiguous ranges and run \a body(start, stop)
//! (stop exclusive) on each in its own thread. All threads are joined before the
//! first worker exception is rethrown.
inline void ParallelRanges(unsigned int count, unsigned int workers, const std::function<void(unsigned int, unsigned int)>& body)
{
	if (workers > count) workers = count;
	if (workers <= 1)
	{
		if (count) body(0, count);
		return;
	}
	std::vector<std::thread> threads;
	std::vector<std::exception_ptr> errors(workers);
	auto run = [&](unsigned int worker) {
		try { body(count*worker/workers, count*(worker+1)/workers); }
		catch (...) { errors[worker] = std::current_exception(); }
	};
	try {
		for (unsigned int i=0; i<workers; ++i) threads.emplace_back(run, i);
	} catch (...) {
		for (auto& thread : threads) thread.join();
		throw;
	}
	for (auto& thread : threads) thread.join();
	for (auto& error : errors) if (error) std::rethrow_exception(error);
}

std::vector<float> SplitString2Float(std::string str, std::string delimiter=",");
std::vector<double> SplitString2Double(std::string str, std::string delimiter=",");

bool CrossProd(const double* v1, const double* v2, double* out);
double ScalarProd(const double* v1, const double* v2);

double Determinant(const double* mat);
double* Invert(const double* in, double* out);

int LinePlaneIntersection(const double *p0, const double* p1, const double* p2, const double* l_start, const double* l_stop, double* is_point, double &dist);

#if defined(_WIN32) && !defined(__GNUC__)
int gettimeofday(struct timeval* tp, struct timezone* tzp);
#endif // defined(_WIN32) && !defined(__GNUC__)

#endif // USEFUL_H
