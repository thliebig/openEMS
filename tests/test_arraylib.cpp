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
  Tests for the memory layout of ArrayLib::ArrayNIJK: in every layout each
  element has its own place in memory and the last axis of the layout string
  is the fastest; the default stays N-I-J-K, which code walking data()
  linearly (the GPU copies) relies on; a string that is not a permutation of
  NIJK is refused.

  Build with -DOPENEMS_BUILD_TESTS=ON and run it through ctest:

    cmake -DOPENEMS_BUILD_TESTS=ON <the usual options> ..
    make && ctest --output-on-failure

  Exits non-zero and prints "FAIL: ..." per failed check.
*/

#include "tools/arraylib/array_nijk.h"

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

static int fails = 0;
#define CHECK(cond, msg) do { if (!(cond)) { std::cout << "FAIL: " << msg << "\n"; ++fails; } } while (0)

using Array = ArrayLib::ArrayNIJK<float>;

// odd extents, so a swapped pair of axes cannot go unnoticed
static const unsigned int ext[3] = {5, 7, 11};

// a distinct value per element
static float Value(unsigned int n, unsigned int i, unsigned int j, unsigned int k)
{
	return (float)(((n*ext[0] + i)*ext[1] + j)*ext[2] + k);
}

// fill through operator() and check every element reads back, through both
// access paths, and that the elements cover the memory exactly once
static void CheckAccess(Array& a, const std::string& what)
{
	for (unsigned int n=0; n<3; ++n)
		for (unsigned int i=0; i<ext[0]; ++i)
			for (unsigned int j=0; j<ext[1]; ++j)
				for (unsigned int k=0; k<ext[2]; ++k)
					a(n, i, j, k) = Value(n, i, j, k);
	int bad = 0;
	for (unsigned int n=0; n<3; ++n)
		for (unsigned int i=0; i<ext[0]; ++i)
			for (unsigned int j=0; j<ext[1]; ++j)
				for (unsigned int k=0; k<ext[2]; ++k)
					bad += (a(n, i, j, k)!=Value(n, i, j, k)) + (a({n, i, j, k})!=Value(n, i, j, k));
	CHECK(bad==0, what << ": " << bad << " elements read back wrong");

	std::vector<float> mem(a.data(), a.data() + a.size());
	std::sort(mem.begin(), mem.end());
	bool bijective = true;
	for (size_t m=0; m<mem.size(); ++m)
		bijective = bijective && mem[m]==(float)m;
	CHECK(bijective, what << ": the elements do not cover the memory exactly once");
}

int main()
{
	Array a("default", {ext[0], ext[1], ext[2]});
	CheckAccess(a, "default layout");
	// N-I-J-K, K fastest: the memory in the order of Value()
	bool linear = true;
	for (size_t m=0; m<a.size(); ++m)
		linear = linear && a.data()[m]==(float)m;
	CHECK(linear, "default layout is not N-I-J-K");
	CHECK(a.layout()=="NIJK", "default layout() is " << a.layout());

	// every permutation of the four axes
	std::string layout = "IJKN";
	std::sort(layout.begin(), layout.end());
	int count = 0;
	do
	{
		Array b("permuted", {ext[0], ext[1], ext[2]}, layout);
		CheckAccess(b, layout);
		CHECK(b.layout()==layout, layout << ": layout() is " << b.layout());
		CHECK(b.stride(std::string("NIJK").find(layout[3]))==1, layout << ": the last axis is not the fastest");
		++count;
	} while (std::next_permutation(layout.begin(), layout.end()));
	CHECK(count==24, "tested " << count << " layouts, not 24");

	for (const char* bad : {"", "NIJ", "NIJKK", "NIIK", "XIJK", "nijk"})
	{
		bool refused = false;
		try {Array c("bad", {2, 2, 2}, bad);}
		catch (const std::invalid_argument&) {refused = true;}
		CHECK(refused, "layout \"" << bad << "\" accepted");
	}

	if (fails==0)
		std::cout << "all ArrayLib checks passed\n";
	return fails ? 1 : 0;
}
