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

#ifndef ENGINE_CPU_KERNEL_H
#define ENGINE_CPU_KERNEL_H

// The data and the update kernels of Engine_CPU. The kernels are compiled
// once per instruction set (engine_cpu_kernels*.cpp) and chosen at run time,
// so this header holds plain data only: an inline function here, compiled
// with AVX2 in one file and without in another, could be merged by the linker
// into the AVX2 version for all callers.

#include <cstddef>
#include <cstdint>

//! The update class and the coefficients a run of cells shares
struct Engine_CPU_Set
{
	enum Class : uint8_t {ZERO, LOSSLESS, LOSSY, UPML};

	double m;              //!< the material part of vi (iv): vi = m * g, see Engine_CPU
	float vv;              //!< vv (ii), 1 for LOSSLESS
	float uvv, ufo, ufn;   //!< UPML: Operator_Ext_UPML::vv, vvfo, vvfn (ii, iifo, iifn); vv and m are the main update on the flux
	Class cls;
};

//! A run of cells along the fast axis with one coefficient set
struct Engine_CPU_Run
{
	uint32_t start;   //!< first line along the fast axis
	uint32_t stop;    //!< one past the last
	uint32_t set;     //!< index in Engine_CPU_Runs::sets
	uint32_t u;       //!< UPML: the index of the first cell's flux in Engine_CPU_UPML
};

//! What the kernels of one engine read and write, filled by Engine_CPU::Init()
struct Engine_CPU_Kernel
{
	int F, S, M;                 //!< the fast, the split and the other slow axis
	unsigned int numLines[3];
	float* field[2];             //!< the voltages and the currents (ArrayNIJK data)
	size_t sN;                   //!< the stride of a component
	size_t s[3];                 //!< the stride of each mesh axis
	struct Component
	{
		const Engine_CPU_Run* runs;
		const uint64_t* first;   //!< per line of the fast axis its first run, see Engine_CPU_Runs
		const Engine_CPU_Set* sets;
		float* flux;             //!< the UPML flux of the component
		const float* gf;         //!< the geometry factor along the fast axis
		const double* gS;        //!< the geometry factors along the slow axes
		const double* gM;
	} comp[2][3];                //!< [voltages, currents][component]
};

//! A voltage or current update of the box [start, stop) (mesh lines)
typedef void (*Engine_CPU_Update)(const Engine_CPU_Kernel& k, const unsigned int start[3], const unsigned int stop[3]);

namespace Engine_CPU_Generic
{
	void UpdateVoltages(const Engine_CPU_Kernel& k, const unsigned int start[3], const unsigned int stop[3]);
	void UpdateCurrents(const Engine_CPU_Kernel& k, const unsigned int start[3], const unsigned int stop[3]);
}

#ifdef OPENEMS_CPU_AVX2
namespace Engine_CPU_AVX2
{
	void UpdateVoltages(const Engine_CPU_Kernel& k, const unsigned int start[3], const unsigned int stop[3]);
	void UpdateCurrents(const Engine_CPU_Kernel& k, const unsigned int start[3], const unsigned int stop[3]);
}
#endif

#endif // ENGINE_CPU_KERNEL_H
