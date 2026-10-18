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

#ifndef GPU_PLUGIN_H
#define GPU_PLUGIN_H

#include <string>
#include <vector>

class GPU_Backend;

//! What a GPU plugin provides: a shared library with a GPU backend, loaded at runtime
/*!
  A plugin is built from the same sources as libopenEMS and reads its classes
  directly, so it only fits the build it belongs to: libopenEMS refuses a plugin
  whose build_id differs from its own (OPENEMS_GPU_BUILD_ID from CMake: the version
  and a random id per build directory, or only the version of a release). The
  plugin exports GPU_PLUGIN_ENTRY, which returns this.
  */
struct GPU_Plugin
{
	const char* build_id;           //!< OPENEMS_GPU_BUILD_ID of the build the plugin belongs to
	const char* name;               //!< short name, e.g. "hip-amd"
	GPU_Backend* (*NewBackend)();   //!< a backend on a usable device, or NULL if there is none
};

//! Name of the function a plugin exports, of type GPU_PluginEntry
#define GPU_PLUGIN_ENTRY "openEMS_gpu_plugin"
typedef const GPU_Plugin* (*GPU_PluginEntry)();

#ifdef _WIN32
#define GPU_PLUGIN_EXPORT extern "C" __declspec(dllexport)
#else
#define GPU_PLUGIN_EXPORT extern "C" __attribute__((visibility("default")))
#endif

//! A backend of the first GPU plugin with a usable device, or NULL
/*!
  The plugins are the files openEMS_gpu_*.so (.dll on Windows) in the directories
  of OPENEMS_PLUGIN_PATH if set, else in openEMS/plugins next to libopenEMS. They
  are loaded once per process and never unloaded. \a messages says, for every
  directory and plugin, why it gave no backend.
  */
GPU_Backend* GPU_NewPluginBackend(std::vector<std::string>& messages);

#endif // GPU_PLUGIN_H
