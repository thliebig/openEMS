/*
*	Copyright (C) 2026 Roel Budé (roel.bude@proton.me)
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

#ifndef NF2FF_KERNELS_H
#define NF2FF_KERNELS_H

// The far-field row sum of nf2ff_calc (nf2ff_kernels.inc), compiled per instruction set.
// Only plain types here: whatever a header defines inline is compiled into the AVX2
// kernels too, and the linker may keep that copy for every caller.

#include <cstddef>

// the q axis of the currents is padded to a multiple of this, the lanes of every kernel
#define NF2FF_STRIDE_ALIGN 8

// S_re/S_im[c] += (a_re + j a_im) * sum_q b[q]*F[c*fstride+q] for c<numComp; F and b are
// split complex arrays of n values, n a multiple of NF2FF_STRIDE_ALIGN
typedef void (*nf2ff_row_func)(const float* F_re, const float* F_im, size_t fstride, unsigned int numComp,
                               const float* b_re, const float* b_im, unsigned int n,
                               double a_re, double a_im, double* S_re, double* S_im);

namespace nf2ff_generic
{
void Row(const float* F_re, const float* F_im, size_t fstride, unsigned int numComp,
         const float* b_re, const float* b_im, unsigned int n,
         double a_re, double a_im, double* S_re, double* S_im);
}

#ifdef NF2FF_AVX2
namespace nf2ff_avx2
{
void Row(const float* F_re, const float* F_im, size_t fstride, unsigned int numComp,
         const float* b_re, const float* b_im, unsigned int n,
         double a_re, double a_im, double* S_re, double* S_im);
}
#endif

#endif // NF2FF_KERNELS_H
