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

#include <iostream>

#include "gpu_backend.h"
#include "gpu_backend_reference.h"
#ifdef OPENEMS_WITH_METAL
#include "metal/gpu_backend_metal.h"
#endif

GPU_Backend* GPU_Backend::New(const std::string& name)
{
	if (name=="reference")
		return new GPU_Backend_Reference();
	if (name!="auto")
		std::cerr << "GPU_Backend::New: Warning: unknown GPU backend \"" << name << "\", using the best available" << std::endl;

#ifdef OPENEMS_WITH_METAL
	GPU_Backend* metal = GPU_Backend_Metal::New();
	if (metal)
		return metal;
	std::cerr << "GPU_Backend::New: Warning: no Metal device found" << std::endl;
#endif

	std::cerr << "GPU_Backend::New: Warning: no GPU backend available, using the reference backend on the CPU" << std::endl;
	return new GPU_Backend_Reference();
}
