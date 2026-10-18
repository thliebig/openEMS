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
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "gpu_backend.h"
#include "gpu_backend_reference.h"
#include "gpu_plugin.h"
#ifdef OPENEMS_WITH_HIP
#include "hip/gpu_backend_hip.h"
#endif
#include "engine.h"
#include "operator.h"
#include "extensions/engine_ext_upml.h"
#include "extensions/engine_ext_steadystate.h"
#ifdef OPENEMS_WITH_METAL
#include "metal/gpu_backend_metal.h"
#endif

namespace
{
struct GPU_Device
{
	const GPU_Plugin* plugin;
	int index;              //!< device number within the plugin
	GPU_DeviceInfo info;
};

//! Discrete GPUs first, then the one with more memory
bool BetterDevice(const GPU_Device& a, const GPU_Device& b)
{
	if (a.info.integrated!=b.info.integrated)
		return !a.info.integrated;
	return a.info.memory>b.info.memory;
}

std::string DeviceLine(size_t n, const GPU_Device& dev)
{
	std::ostringstream line;
	line << n << ": " << dev.info.name << " (" << dev.plugin->name << ", "
	     << std::fixed << std::setprecision(1) << dev.info.memory/1073741824.0 << " GiB"
	     << (dev.info.integrated ? ", integrated" : "") << ")";
	return line.str();
}
}

GPU_Backend* GPU_Backend::New(const std::string& name)
{
	if (name=="reference")
		return new GPU_Backend_Reference();
	if (name!="auto")
		std::cerr << "GPU_Backend::New: Warning: unknown GPU backend \"" << name << "\", using the best available" << std::endl;

	// every backend there is, built in or a plugin; the device numbers of --gpu-device
	// follow this order
	std::vector<std::string> messages;
	std::vector<const GPU_Plugin*> plugins;
#ifdef OPENEMS_WITH_HIP
	plugins.push_back(GPU_Backend_HIP::Plugin());   // OPENEMS_GPU_PLUGINS=OFF
#endif
	const std::vector<const GPU_Plugin*> loaded = GPU_LoadPlugins(messages);
	plugins.insert(plugins.end(), loaded.begin(), loaded.end());
#ifdef OPENEMS_WITH_METAL
	plugins.push_back(GPU_Backend_Metal::Plugin());
#endif

	std::vector<GPU_Device> devices;
	for (const GPU_Plugin* plugin : plugins)
	{
		const int count = plugin->DeviceCount();
		if (count==0)
			messages.push_back(std::string("no ") + plugin->name + " device");
		for (int n=0; n<count; ++n)
		{
			GPU_Device dev;
			dev.plugin = plugin;
			dev.index = n;
			if (plugin->DeviceInfo(n, &dev.info))
				devices.push_back(dev);
		}
	}

	const int chosen = g_settings.GetGPUDevice();
	if (g_settings.GetVerboseLevel()>0)
	{
		std::cout << "GPU devices:" << std::endl;
		for (size_t n=0; n<devices.size(); ++n)
			std::cout << "  " << DeviceLine(n, devices.at(n)) << std::endl;
	}
	if ((chosen<-1) || (chosen>=(int)devices.size()))
	{
		std::ostringstream err;
		err << "--gpu-device=" << chosen << ", but there " << (devices.size()==1 ? "is " : "are ")
		    << devices.size() << " GPU device" << (devices.size()==1 ? "" : "s");
		for (size_t n=0; n<devices.size(); ++n)
			err << (n==0 ? ": " : ", ") << DeviceLine(n, devices.at(n));
		throw std::runtime_error(err.str());
	}

	// the chosen device, else the best one that can run
	std::vector<size_t> order;
	if (chosen>=0)
		order.push_back(chosen);
	else
	{
		for (size_t n=0; n<devices.size(); ++n)
			order.push_back(n);
		std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {return BetterDevice(devices.at(a), devices.at(b));});
	}
	for (size_t n : order)
	{
		const GPU_Device& dev = devices.at(n);
		if (GPU_Backend* backend = dev.plugin->NewBackend(dev.index))
			return backend;
		messages.push_back("cannot run on device " + DeviceLine(n, dev));
	}

	for (const std::string& msg : messages)
		std::cerr << "GPU_Backend::New: Warning: " << msg << std::endl;
	// not the reference backend: it is far slower than the CPU engines, the caller has to say so
	return NULL;
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

static FDTD_FLOAT GPU_Coefficient(const Operator* op, int c, int n, unsigned int x, unsigned int y, unsigned int z)
{
	switch (c)
	{
	case 0: return op->GetVV(n, x, y, z);
	case 1: return op->GetVI(n, x, y, z);
	case 2: return op->GetII(n, x, y, z);
	default: return op->GetIV(n, x, y, z);
	}
}

void GPU_NodeCoefficients(const Operator* op, const unsigned int numLines[3], size_t node, float values[12])
{
	const unsigned int x = node/((size_t)numLines[1]*numLines[2]);
	const unsigned int y = (node/numLines[2]) % numLines[1];
	const unsigned int z = node % numLines[2];
	for (int c=0; c<4; ++c)
		for (int n=0; n<3; ++n)
			values[3*c+n] = GPU_Coefficient(op, c, n, x, y, z);
}

void GPU_CoefficientArray(const Operator* op, const unsigned int numLines[3], int c, std::vector<float>& values)
{
	values.resize(3*(size_t)numLines[0]*numLines[1]*numLines[2]);
	size_t idx = 0;
	for (int n=0; n<3; ++n)
		for (unsigned int x=0; x<numLines[0]; ++x)
			for (unsigned int y=0; y<numLines[1]; ++y)
				for (unsigned int z=0; z<numLines[2]; ++z, ++idx)
					values[idx] = GPU_Coefficient(op, c, n, x, y, z);
}
