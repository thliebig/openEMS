// Host unit test for FDTD/metal_predicates.h.
//
// Compile (Boost is a build dependency of openEMS):
//   clang++ -std=c++11 -O2 -I<repo-root> \
//     TESTSUITE/enginetests/metal_predicates_test.cpp -o /tmp/metal_predicates_test \
//     -I/opt/homebrew/include
//   /tmp/metal_predicates_test
//
// It checks that whenever mp_orient2d_sign() returns a non-zero sign it equals
// the exact double-precision determinant sign computed with 100-bit arithmetic,
// and reports how often the predicate correctly declines (returns 0).

#include "FDTD/metal_predicates.h"

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <random>

using boost::multiprecision::cpp_bin_float_100;

static mp_df split(double d)
{
	mp_df r;
	r.hi = static_cast<float>(d);
	r.lo = static_cast<float>(d - static_cast<double>(r.hi));
	return r;
}

static int exact_sign(double x1, double y1, double x2, double y2, double px, double py)
{
	cpp_bin_float_100 a(x1), b(y1), c(x2), d(y2), e(px), f(py);
	cpp_bin_float_100 det = (a - e) * (d - f) - (b - f) * (c - e);
	if (det > 0) return 1;
	if (det < 0) return -1;
	return 0;
}

int main()
{
	std::mt19937_64 rng(0x9E3779B97F4A7C15ull);
	auto rand_double = [&]() {
		double u = static_cast<double>(rng() >> 11) / static_cast<double>(1ull << 53);
		return u * 2000.0 - 1000.0;
	};

	long random_total = 0, random_certain = 0, edge_certain = 0, edge_total = 0;
	long uncertain = 0, mismatches = 0;

	for (long i = 0; i < 4000000; ++i)
	{
		double x1 = rand_double(), y1 = rand_double();
		double x2 = rand_double(), y2 = rand_double();
		double px, py;
		bool edge = false;
		if (i % 2 == 0)
		{
			// Nearly collinear: point on the segment plus a tiny offset.
			double t = (rng() >> 11) / static_cast<double>(1ull << 53);
			px = x1 + t * (x2 - x1);
			py = y1 + t * (y2 - y1);
			// Perturb by a small but resolvable amount so the exact sign is
			// usually non-zero at double precision. A few cases remain ties.
			double scale = 1.0 + std::fabs(px) + std::fabs(py);
			int ulps = static_cast<int>(rng() % 9) - 4;
			px += ulps * 1e-10 * scale;
			py += (static_cast<int>(rng() % 9) - 4) * 1e-10 * scale;
			edge = true;
		}
		else
		{
			px = rand_double();
			py = rand_double();
		}

		int ref = exact_sign(x1, y1, x2, y2, px, py);
		int got = mp_orient2d_sign(split(x1), split(y1), split(x2), split(y2),
		                           split(px), split(py));
		if (edge) ++edge_total; else ++random_total;
		if (got == 0) ++uncertain;
		else
		{
			if (edge) ++edge_certain; else ++random_certain;
			if (got != ref)
			{
				++mismatches;
				if (mismatches <= 10)
					printf("MISMATCH ref=%d got=%d  (%g,%g) (%g,%g) (%g,%g)\n",
					       ref, got, x1, y1, x2, y2, px, py);
			}
		}
	}

	printf("random: %ld cases, %ld certified (%.3f%%)\n",
	       random_total, random_certain, 100.0 * random_certain / random_total);
	printf("near-edge: %ld cases, %ld certified (%.3f%%), uncertain(total)=%ld\n",
	       edge_total, edge_certain, 100.0 * edge_certain / edge_total, uncertain);
	printf("mismatches=%ld\n", mismatches);
	if (mismatches)
	{
		printf("FAIL: certified sign disagreed with exact arithmetic\n");
		return 1;
	}
	// General cases must almost always be certified; a low rate means the
	// filter is too conservative and would keep falling back to the CPU.
	if (random_certain < random_total * 99 / 100)
	{
		printf("FAIL: too few certified general cases\n");
		return 1;
	}
	printf("PASS\n");
	return 0;
}
