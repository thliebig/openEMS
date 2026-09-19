/*
*	Copyright (C) 2026 Sean Mollet (sean@malmoset.com)
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

#include <thread>

#include "engine_gpu.h"
#include "operator_gpu.h"
#include "tools/useful.h"

using std::cout;
using std::cerr;
using std::endl;

Operator_GPU* Operator_GPU::New(const std::string& backend)
{
	cout << "Create FDTD operator (GPU)" << endl;
	Operator_GPU* op = new Operator_GPU(backend);
	op->Init();
	return op;
}

Operator_GPU::Operator_GPU(const std::string& backend) : Operator()
{
	m_Backend = backend;
}

Operator_GPU::~Operator_GPU()
{
}

Engine* Operator_GPU::CreateEngine()
{
	m_Engine = Engine_GPU::New(this, m_Backend);
	return m_Engine;
}

void Operator_GPU::ThreadRanges(std::vector<unsigned int>& start, std::vector<unsigned int>& stop) const
{
	std::vector<unsigned int> jobs = AssignJobs2Threads(numLines[0], AvailableCPUs(), true);
	unsigned int line = 0;
	for (size_t n=0; n<jobs.size(); ++n)
	{
		start.push_back(line);
		line += jobs.at(n);
		stop.push_back(line-1);
	}
}

// see Operator::Calc_EC() and Operator_Multithread
bool Operator_GPU::Calc_EC()
{
	if (CSX==NULL)
	{
		cerr << "Operator_GPU::Calc_EC: CSX not given or invalid!!!" << endl;
		return false;
	}
	MainOp->SetPos(0,0,0);

	std::vector<unsigned int> start, stop;
	ThreadRanges(start, stop);
	std::vector<std::thread> threads;
	for (size_t n=0; n<start.size(); ++n)
		threads.push_back(std::thread(&Operator_GPU::Calc_EC_Range, this, start[n], stop[n]));
	for (size_t n=0; n<threads.size(); ++n)
		threads[n].join();
	return true;
}

// see Operator::CalcPEC() and Operator_Multithread
bool Operator_GPU::CalcPEC()
{
	std::vector<unsigned int> start, stop;
	ThreadRanges(start, stop);
	std::vector<unsigned int> counter(3*start.size(), 0);
	std::vector<std::thread> threads;
	for (size_t n=0; n<start.size(); ++n)
		threads.push_back(std::thread(&Operator_GPU::CalcPEC_Range, this, start[n], stop[n], &counter[3*n]));
	for (size_t n=0; n<threads.size(); ++n)
		threads[n].join();

	for (int n=0; n<3; ++n)
	{
		m_Nr_PEC[n] = 0;
		for (size_t t=0; t<start.size(); ++t)
			m_Nr_PEC[n] += counter[3*t+n];
	}
	CalcPEC_Curves();
	return true;
}
