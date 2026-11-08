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

// The nf2ff kernels for AVX2 and FMA: compiled with -mavx2 -mfma (CMakeLists.txt), run only
// on a CPU that has both (nf2ff_calc::AddSinglePlane()).

#include "nf2ff_kernels.h"   // and nothing else, see there

#define NF2FF_KERNEL_NAMESPACE nf2ff_avx2
#define NF2FF_LANES 8
#include "nf2ff_kernels.inc"
