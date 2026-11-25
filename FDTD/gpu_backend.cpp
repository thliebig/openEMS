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
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "gpu_backend.h"
#include "gpu_cell_model.h"
#include "gpu_ext.h"
#include "operator.h"
#include "gpu_backend_reference.h"
#include "gpu_split.h"
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

bool GPU_Backend::FitsOperator(const Operator* op, const GPU_CellModel& model)
{
	for (int n=0; op && (n<3); ++n)
		if (model.numLines[n]!=op->GetNumberOfLines(n, true))
			return false;
	return true;
}

void GPU_Backend::ReadLine(bool currents, unsigned int n, unsigned int x, unsigned int y, const unsigned int numLines[3], FDTD_FLOAT* line)
{
	const size_t key = ((size_t)n*numLines[0] + x)*numLines[1] + y;
	DeviceDownload(line, DeviceField(currents) + key*numLines[2], numLines[2]*sizeof(FDTD_FLOAT));
}

GPU_Extension* GPU_Backend::CreateExtension(Engine_Extension* eng_ext, const unsigned int numLines[3])
{
	return GPU_CreatePartExtension(eng_ext, this, GPU_Slab::Whole(numLines), NULL);
}

GPU_Backend* GPU_Backend::New(const std::string& name)
{
	// --gpu-devices, else OPENEMS_GPU_DEVICES: the grid split over several devices (GPU_Backend_Split).
	// The environment variable gives way to --gpu-device; the two options must agree.
	std::string option;
	bool from_env = false;
	std::vector<unsigned int> split = GPU_NumberList(g_settings.GetGPUDevices(), "--gpu-devices", "OPENEMS_GPU_DEVICES", option, from_env);
	const int device = g_settings.GetGPUDevice();
	if ((device>=0) && !split.empty())
	{
		if (from_env)
			split.clear();
		else if ((split.size()>1) || ((int)split.front()!=device))
			throw std::runtime_error(option + " and --gpu-device=" + std::to_string(device) + ": use one of them");
	}

	if (name=="reference")
	{
		if (split.size()<2)
			return new GPU_Backend_Reference();
		// for tests: only the number of parts counts
		std::vector<GPU_Backend*> parts;
		for (size_t n=0; n<split.size(); ++n)
			parts.push_back(new GPU_Backend_Reference());
		return new GPU_Backend_Split(parts);
	}
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

	int chosen = g_settings.GetGPUDevice();
	if (g_settings.GetVerboseLevel()>0)
	{
		std::cout << "GPU devices:" << std::endl;
		for (size_t n=0; n<devices.size(); ++n)
			std::cout << "  " << DeviceLine(n, devices.at(n)) << std::endl;
	}
	auto no_such_device = [&](const std::string& what, int device)
	{
		std::ostringstream err;
		err << what << device << ", but there " << (devices.size()==1 ? "is " : "are ")
		    << devices.size() << " GPU device" << (devices.size()==1 ? "" : "s");
		for (size_t n=0; n<devices.size(); ++n)
			err << (n==0 ? ": " : ", ") << DeviceLine(n, devices.at(n));
		throw std::runtime_error(err.str());
	};
	for (unsigned int dev : split)
		if (dev>=devices.size())
			no_such_device(option + ": device ", dev);
	if (split.size()==1)
		chosen = split.front();   // one device: as --gpu-device
	else if (split.size()>1)
	{
		// a backend per part, on its device; a device may hold several parts
		std::vector<GPU_Backend*> parts;
		for (unsigned int dev : split)
		{
			GPU_Backend* backend = devices.at(dev).plugin->NewBackend(devices.at(dev).index);
			if (backend==NULL)
			{
				for (GPU_Backend* b : parts)
					delete b;
				for (const std::string& msg : messages)
					std::cerr << "GPU_Backend::New: Warning: " << msg << std::endl;
				throw std::runtime_error(option + ": cannot run on device " + DeviceLine(dev, devices.at(dev)));
			}
			parts.push_back(backend);
		}
		return new GPU_Backend_Split(parts);
	}
	if ((chosen<-1) || (chosen>=(int)devices.size()))
		no_such_device("--gpu-device=", chosen);

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
