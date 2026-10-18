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

#include "gpu_plugin.h"
#include "gpu_build_id.h"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <filesystem>

#ifdef _WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#endif

namespace fs = std::filesystem;

namespace
{
#ifdef _WIN32
const char PATH_SEPARATOR = ';';
const char* PLUGIN_SUFFIX = ".dll";
#else
const char PATH_SEPARATOR = ':';
const char* PLUGIN_SUFFIX = ".so";
#endif
const char* PLUGIN_PREFIX = "openEMS_gpu_";

//! Directory of the loaded libopenEMS itself
fs::path CoreDirectory()
{
#ifdef _WIN32
	HMODULE self = NULL;
	if (!GetModuleHandleExW(GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS | GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
	                        reinterpret_cast<LPCWSTR>(&CoreDirectory), &self))
		return fs::path();
	wchar_t path[MAX_PATH];
	if (GetModuleFileNameW(self, path, MAX_PATH)==0)
		return fs::path();
	return fs::path(path).parent_path();
#else
	Dl_info info;
	if (!dladdr(reinterpret_cast<void*>(&CoreDirectory), &info) || !info.dli_fname)
		return fs::path();
	return fs::absolute(info.dli_fname).parent_path();
#endif
}

std::vector<fs::path> PluginDirectories()
{
	std::vector<fs::path> dirs;
	if (const char* env = getenv("OPENEMS_PLUGIN_PATH"))
	{
		std::string list = env;
		size_t start = 0;
		while (start<=list.size())
		{
			size_t end = list.find(PATH_SEPARATOR, start);
			if (end==std::string::npos)
				end = list.size();
			if (end>start)
				dirs.push_back(list.substr(start, end-start));
			start = end+1;
		}
		return dirs;
	}
	const fs::path core = CoreDirectory();
	if (!core.empty())
		dirs.push_back(core / "openEMS" / "plugins");
	return dirs;
}

#ifdef _WIN32
std::string LastError()
{
	const DWORD err = GetLastError();
	char* buf = NULL;
	FormatMessageA(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
	               NULL, err, 0, reinterpret_cast<LPSTR>(&buf), 0, NULL);
	std::string msg = buf ? buf : "";
	LocalFree(buf);
	while (!msg.empty() && ((msg.back()=='\n') || (msg.back()=='\r')))
		msg.pop_back();
	if (err==ERROR_MOD_NOT_FOUND)
		msg += " (a library the plugin needs is missing, e.g. the GPU driver or runtime)";
	return msg;
}
#endif

//! Load \a file and return its plugin description, or NULL with the reason in \a error
const GPU_Plugin* LoadPlugin(const fs::path& file, std::string& error)
{
	GPU_PluginEntry entry = NULL;
#ifdef _WIN32
	// an absolute path, else LOAD_LIBRARY_SEARCH_DLL_LOAD_DIR is refused; the plugin's
	// own directory is searched for its dependencies, and no error box pops up
	const UINT oldMode = SetErrorMode(SEM_FAILCRITICALERRORS);
	HMODULE lib = LoadLibraryExW(fs::absolute(file).wstring().c_str(), NULL,
	                             LOAD_LIBRARY_SEARCH_DLL_LOAD_DIR | LOAD_LIBRARY_SEARCH_DEFAULT_DIRS);
	SetErrorMode(oldMode);
	if (!lib)
	{
		error = LastError();
		return NULL;
	}
	entry = reinterpret_cast<GPU_PluginEntry>(reinterpret_cast<void*>(GetProcAddress(lib, GPU_PLUGIN_ENTRY)));
#else
	// RTLD_NOW: a missing symbol is a load error here, not a crash in the middle of a run
	void* lib = dlopen(file.c_str(), RTLD_NOW | RTLD_LOCAL);
	if (!lib)
	{
		error = dlerror();
		return NULL;
	}
	entry = reinterpret_cast<GPU_PluginEntry>(dlsym(lib, GPU_PLUGIN_ENTRY));
#endif
	// a plugin is never unloaded: the backends it creates run its code
	if (!entry)
	{
		error = std::string("not an openEMS GPU plugin (no ") + GPU_PLUGIN_ENTRY + ")";
		return NULL;
	}
	const GPU_Plugin* plugin = entry();
	if (!plugin || !plugin->build_id || !plugin->name || !plugin->DeviceCount || !plugin->DeviceInfo || !plugin->NewBackend)
	{
		error = "invalid plugin description";
		return NULL;
	}
	if (strcmp(plugin->build_id, OPENEMS_GPU_BUILD_ID)!=0)
	{
		error = std::string("belongs to another openEMS build (") + plugin->build_id + "), this is " + OPENEMS_GPU_BUILD_ID;
		return NULL;
	}
	return plugin;
}

struct PluginFile
{
	fs::path file;
	const GPU_Plugin* plugin;   //!< NULL if it failed to load
	std::string error;
};

//! All plugin files, loaded on the first call; \a dir_messages says which directories had none
const std::vector<PluginFile>& Plugins(std::vector<std::string>& dir_messages)
{
	static std::vector<std::string> s_dir_messages;
	static const std::vector<PluginFile> plugins = []
	{
		std::vector<PluginFile> found;
		for (const fs::path& dir : PluginDirectories())
		{
			std::vector<fs::path> files;
			std::error_code ec;
			for (fs::directory_iterator it(dir, ec), end; !ec && (it!=end); it.increment(ec))
			{
				const std::string name = it->path().filename().string();
				if ((name.compare(0, strlen(PLUGIN_PREFIX), PLUGIN_PREFIX)==0) && (it->path().extension()==PLUGIN_SUFFIX))
					files.push_back(it->path());
			}
			if (files.empty())
				s_dir_messages.push_back("no GPU plugin in " + dir.string());
			std::sort(files.begin(), files.end());
			for (const fs::path& file : files)
			{
				PluginFile p;
				p.file = file;
				p.plugin = LoadPlugin(file, p.error);
				found.push_back(p);
			}
		}
		if (PluginDirectories().empty())
			s_dir_messages.push_back("no directory to look for GPU plugins in");
		return found;
	}();
	dir_messages.insert(dir_messages.end(), s_dir_messages.begin(), s_dir_messages.end());
	return plugins;
}
}

std::vector<const GPU_Plugin*> GPU_LoadPlugins(std::vector<std::string>& messages)
{
	std::vector<const GPU_Plugin*> loaded;
	for (const PluginFile& p : Plugins(messages))
	{
		if (p.plugin)
			loaded.push_back(p.plugin);
		else
			messages.push_back(p.file.string() + ": " + p.error);
	}
	return loaded;
}
