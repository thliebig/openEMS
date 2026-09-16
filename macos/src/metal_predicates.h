/*
 * Copyright (C) 2026 openEMS contributors
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * Higher-precision geometric predicates for the Metal PEC mapping.
 *
 * Apple GPUs have no FP64 in Metal, so coordinates are carried as "double
 * float" pairs (hi + lo, ~48 bits). The polygon winding predicate mirrors
 * CSXCAD's CSPrimPolygon::IsInside. It returns 0 when the result is within a
 * conservative bound of the decision boundary, in which case the caller must
 * fall back to the CPU so CSXCAD's exact double result is preserved.
 *
 * This header is compiled both as C++ (for the host unit test) and as part of
 * the Metal shader source, so it must stay valid in both languages.
 */

#ifndef OPENEMS_METAL_PREDICATES_H
#define OPENEMS_METAL_PREDICATES_H

#ifdef __METAL_VERSION__
  #define MP_FMA(a, b, c) fma(a, b, c)
  #define MP_FABS(a) fabs(a)
#else
  #include <cmath>
  #define MP_FMA(a, b, c) std::fma(a, b, c)
  #define MP_FABS(a) std::fabs(a)
#endif

// Double-float (two-float) value hi + lo.
typedef struct { float hi, lo; } mp_df;

static inline mp_df mp_quick_two_sum(float a, float b)
{
	mp_df r;
	r.hi = a + b;
	r.lo = b - (r.hi - a);
	return r;
}

static inline mp_df mp_two_sum(float a, float b)
{
	mp_df r;
	r.hi = a + b;
	float bv = r.hi - a;
	r.lo = (a - (r.hi - bv)) + (b - bv);
	return r;
}

static inline mp_df mp_df_add(mp_df a, mp_df b)
{
	mp_df s = mp_two_sum(a.hi, b.hi);
	mp_df t = mp_two_sum(a.lo, b.lo);
	s.lo += t.hi;
	s = mp_quick_two_sum(s.hi, s.lo);
	s.lo += t.lo;
	return mp_quick_two_sum(s.hi, s.lo);
}

static inline mp_df mp_df_sub(mp_df a, mp_df b)
{
	mp_df nb;
	nb.hi = -b.hi;
	nb.lo = -b.lo;
	return mp_df_add(a, nb);
}

static inline mp_df mp_df_mul(mp_df a, mp_df b)
{
	float p = a.hi * b.hi;
	float e = MP_FMA(a.hi, b.hi, -p) + a.hi * b.lo + a.lo * b.hi;
	return mp_quick_two_sum(p, e);
}

static inline int mp_df_ge(mp_df a, mp_df b)
{
	if (a.hi > b.hi) return 1;
	if (a.hi < b.hi) return 0;
	return a.lo >= b.lo;
}

static inline int mp_df_eq(mp_df a, mp_df b)
{
	return a.hi == b.hi && a.lo == b.lo;
}

// Sign of orient2d(prev, cur, point) = (x1-px)*(y2-py) - (y1-py)*(x2-px).
// Returns +1 / -1, or 0 when the sign cannot be certified against CSXCAD's
// double-precision decision. The caller treats 0 as "query the CPU".
static inline int mp_orient2d_sign(mp_df x1, mp_df y1, mp_df x2, mp_df y2, mp_df px, mp_df py)
{
	mp_df A = mp_df_sub(x1, px);
	mp_df B = mp_df_sub(y2, py);
	mp_df C = mp_df_sub(y1, py);
	mp_df D = mp_df_sub(x2, px);
	mp_df det = mp_df_sub(mp_df_mul(A, B), mp_df_mul(C, D));

	// Error budget: double-float arithmetic (~2^-46) plus the hi/lo input
	// representation error (~2^-48) scale like scale^2. 1e-11 leaves a large
	// margin over both while still resolving ~2^13 finer than the old FP32
	// filter, so genuine near-edge cases stay on the GPU.
	float scale = MP_FABS(x1.hi) + MP_FABS(y1.hi) + MP_FABS(x2.hi) +
	              MP_FABS(y2.hi) + MP_FABS(px.hi) + MP_FABS(py.hi);
	float bound = scale * scale * 1.0e-11f;

	if (MP_FABS(det.hi) > bound)
		return det.hi > 0.0f ? 1 : -1;
	if (det.hi == 0.0f && MP_FABS(det.lo) > bound)
		return det.lo > 0.0f ? 1 : -1;
	return 0;
}

#endif
