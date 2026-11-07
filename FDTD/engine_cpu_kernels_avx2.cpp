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

// The kernels of Engine_CPU for AVX2 and FMA: compiled with -mavx2 -mfma (CMakeLists.txt), run only on a CPU that has both (Engine_CPU::Init()).

#include "engine_cpu_kernel.h"

#define ENGINE_CPU_NAMESPACE Engine_CPU_AVX2
#include "engine_cpu_kernels.inc"
