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
