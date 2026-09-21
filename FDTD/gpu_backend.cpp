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

#include <algorithm>
#include <iostream>

#include "gpu_backend.h"
#include "gpu_backend_reference.h"
#include "engine.h"
#include "extensions/engine_ext_upml.h"
#include "extensions/engine_ext_steadystate.h"
#ifdef OPENEMS_WITH_METAL
#include "metal/gpu_backend_metal.h"
#endif
#ifdef OPENEMS_WITH_HIP
#include "hip/gpu_backend_hip.h"
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

#ifdef OPENEMS_WITH_HIP
	GPU_Backend* hip = GPU_Backend_HIP::New();
	if (hip)
		return hip;
	std::cerr << "GPU_Backend::New: Warning: no HIP device found" << std::endl;
#endif

	std::cerr << "GPU_Backend::New: Warning: no GPU backend available, using the reference backend on the CPU" << std::endl;
	return new GPU_Backend_Reference();
}

bool GPU_UPMLFusionBox(Engine* eng, const std::vector<GPU_UPMLRegion>& regions, const unsigned int numLines[3],
                       unsigned int start[3], unsigned int stop[3])
{
	// the extension order: UPML (and steady-state) first, so that nothing else
	// touches the fields between the UPML hooks and the main update
	size_t num_upml = 0;
	size_t last_upml = 0;
	for (size_t n=0; n<eng->GetExtensionCount(); ++n)
		if (dynamic_cast<Engine_Ext_UPML*>(eng->GetExtension(n)))
		{
			++num_upml;
			last_upml = n;
		}
	if ((num_upml==0) || (num_upml!=regions.size()))
		return false;
	for (size_t n=0; n<last_upml; ++n)
	{
		Engine_Extension* ext = eng->GetExtension(n);
		if (!dynamic_cast<Engine_Ext_UPML*>(ext) && !dynamic_cast<Engine_Ext_SteadyState*>(ext))
			return false;
	}

	// mark the region nodes, the others must form a box
	const size_t nx = numLines[0], ny = numLines[1], nz = numLines[2];
	std::vector<unsigned char> mask(nx*ny*nz, 0);
	for (size_t r=0; r<regions.size(); ++r)
	{
		const GPU_UPMLRegion& R = regions[r];
		for (int n=0; n<3; ++n)
			if (R.start[n]+R.size[n]>numLines[n])
				return false;
		for (size_t x=R.start[0]; x<R.start[0]+R.size[0]; ++x)
			for (size_t y=R.start[1]; y<R.start[1]+R.size[1]; ++y)
				for (size_t z=R.start[2]; z<R.start[2]+R.size[2]; ++z)
				{
					unsigned char& m = mask[(x*ny + y)*nz + z];
					if (m)
						return false;   // overlapping regions
					m = 1;
				}
	}
	unsigned int lo[3] = {numLines[0], numLines[1], numLines[2]};
	unsigned int hi[3] = {0, 0, 0};
	size_t free_nodes = 0;
	for (unsigned int x=0; x<nx; ++x)
		for (unsigned int y=0; y<ny; ++y)
			for (unsigned int z=0; z<nz; ++z)
				if (!mask[((size_t)x*ny + y)*nz + z])
				{
					++free_nodes;
					const unsigned int pos[3] = {x, y, z};
					for (int n=0; n<3; ++n)
					{
						lo[n] = std::min(lo[n], pos[n]);
						hi[n] = std::max(hi[n], pos[n]+1);
					}
				}
	if (free_nodes==0)
		return false;
	// the nodes left over must fill their bounding box, else they are no single box
	if (free_nodes!=(size_t)(hi[0]-lo[0])*(hi[1]-lo[1])*(hi[2]-lo[2]))
		return false;

	for (int n=0; n<3; ++n)
	{
		start[n] = lo[n];
		stop[n] = hi[n];
	}
	return true;
}
