/*
*	Copyright (C) 2026 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

/*
  Tests for Engine_CPU_NodeSet, the lookup of the CPU engine from a mesh node
  to its index in an extension's cell list (the coefficients of a dispersive
  cell): for nodes added in ascending order, Find() must give each its index
  and -1 for every other node. The sets cover the word boundaries of the bit
  set (nodes 0, 63, 64, the last), empty and full words, and a sparse random
  set over many words.

  Run through ctest, see test_engine_plan.cpp. Exits non-zero and prints
  "FAIL: ..." per failed check.
*/

#include "FDTD/engine_cpu.h"

#include <iostream>
#include <random>
#include <string>
#include <vector>

static int g_failures = 0;

static void Check(const std::string& name, size_t nodes, const std::vector<size_t>& members)
{
	Engine_CPU_NodeSet set;
	set.Init(nodes);
	for (size_t lin : members)
		if (!set.Add(lin))
		{
			std::cerr << "FAIL: " << name << ": node " << lin << " added twice" << std::endl;
			++g_failures;
			return;
		}
	if (!members.empty() && set.Add(members.front()))
	{
		std::cerr << "FAIL: " << name << ": a node in the set added again" << std::endl;
		++g_failures;
	}
	set.Finish();
	std::vector<int64_t> expect(nodes, -1);
	for (size_t i=0; i<members.size(); ++i)
		expect[members[i]] = (int64_t)i;
	for (size_t lin=0; lin<nodes; ++lin)
		if (set.Find(lin)!=expect[lin])
		{
			std::cerr << "FAIL: " << name << ": Find(" << lin << ") = " << set.Find(lin) << ", expected " << expect[lin] << std::endl;
			++g_failures;
			return;
		}
}

int main()
{
	Check("edges", 200, {0, 63, 64, 127, 128, 199});
	std::vector<size_t> full;
	for (size_t lin=64; lin<192; ++lin)
		full.push_back(lin);
	Check("full words", 256, full);
	Check("empty", 100, {});
	Check("one node", 1, {0});
	std::mt19937 rng(1);
	std::vector<size_t> sparse;
	for (size_t lin=0; lin<100000; ++lin)
		if (rng()%7==0)
			sparse.push_back(lin);
	Check("sparse", 100000, sparse);
	if (g_failures==0)
		std::cout << "test_engine_cpu: all checks passed" << std::endl;
	return g_failures ? 1 : 0;
}
