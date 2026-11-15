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

// The far-field row sum of nf2ff_calc (declared in nf2ff_kernels.h), included once per
// instruction set with NF2FF_KERNEL_NAMESPACE and NF2FF_LANES defined. Plain arithmetic
// only, see nf2ff_kernels.h.

namespace NF2FF_KERNEL_NAMESPACE
{

// S[c] += a * sum_q b[q]*f[c*fstride+q] for c<NC, split complex arrays, n a multiple of NF2FF_LANES
template <unsigned int NC>
static inline void CDot(const float* b_re, const float* b_im, const float* f_re, const float* f_im, size_t fstride, unsigned int n,
                        double a_re, double a_im, double* S_re, double* S_im)
{
	float s_re[NC][NF2FF_LANES] = {};
	float s_im[NC][NF2FF_LANES] = {};
	for (unsigned int q=0; q<n; q+=NF2FF_LANES)
		for (unsigned int c=0; c<NC; ++c)
			for (unsigned int l=0; l<NF2FF_LANES; ++l)
			{
				float fr = f_re[c*fstride+q+l];
				float fi = f_im[c*fstride+q+l];
				s_re[c][l] += b_re[q+l]*fr - b_im[q+l]*fi;
				s_im[c][l] += b_re[q+l]*fi + b_im[q+l]*fr;
			}
	for (unsigned int c=0; c<NC; ++c)
	{
		double re=0, im=0;
		for (unsigned int l=0; l<NF2FF_LANES; ++l)
		{
			re += s_re[c][l];
			im += s_im[c][l];
		}
		S_re[c] += a_re*re - a_im*im;
		S_im[c] += a_re*im + a_im*re;
	}
}

// all current components in one pass, so the phase factors are loaded once
void Row(const float* F_re, const float* F_im, size_t fstride, unsigned int numComp,
         const float* b_re, const float* b_im, unsigned int n,
         double a_re, double a_im, double* S_re, double* S_im)
{
	switch (numComp)
	{
	case 4:
		CDot<4>(b_re, b_im, F_re, F_im, fstride, n, a_re, a_im, S_re, S_im);
		break;
	case 6:
		CDot<6>(b_re, b_im, F_re, F_im, fstride, n, a_re, a_im, S_re, S_im);
		break;
	default:
		for (unsigned int c=0; c<numComp; ++c)
			CDot<1>(b_re, b_im, F_re+c*fstride, F_im+c*fstride, 0, n, a_re, a_im, S_re+c, S_im+c);
	}
}

}
