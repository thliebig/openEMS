// SPDX-License-Identifier: GPL-3.0-or-later
// Test-only accounting for the single-process fixture, Linux/shared CUDA runtime.
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <dlfcn.h>
#include <mutex>
#include <unordered_map>

namespace
{
struct Allocations
{
	std::mutex mutex;
	std::unordered_map<void*, size_t> live;
	size_t current = 0, peak = 0, calls = 0;
};

Allocations& allocations()
{
	// Keep the accounting alive through CUDA/library teardown.
	static Allocations* value = new Allocations;
	return *value;
}

void* runtime_symbol(const char* name)
{
	static void* handle = []()
	{
		const char* library = std::getenv("OPENEMS_CUDART_LIBRARY");
		return dlopen(library ? library : "libcudart.so", RTLD_NOW | RTLD_LOCAL);
	}();
	void* symbol = handle ? dlsym(handle, name) : nullptr;
	if (!symbol)
	{
		std::fprintf(stderr, "CUDA allocation trace: cannot resolve %s; set OPENEMS_CUDART_LIBRARY to the shared CUDA runtime\n", name);
		std::abort();
	}
	return symbol;
}
}

extern "C" int cudaMalloc(void** ptr, size_t bytes)
{
	static auto real = reinterpret_cast<int(*)(void**, size_t)>(runtime_symbol("cudaMalloc"));
	const int result = real(ptr, bytes);
	if (!result)
	{
		Allocations& a = allocations();
		std::lock_guard<std::mutex> lock(a.mutex);
		a.live[*ptr] = bytes;
		a.current += bytes;
		a.peak = std::max(a.peak, a.current);
		++a.calls;
	}
	return result;
}

extern "C" int cudaFree(void* ptr)
{
	static auto real = reinterpret_cast<int(*)(void*)>(runtime_symbol("cudaFree"));
	const int result = real(ptr);
	if (!result)
	{
		Allocations& a = allocations();
		std::lock_guard<std::mutex> lock(a.mutex);
		auto it = a.live.find(ptr);
		if (it!=a.live.end())
		{
			a.current -= it->second;
			a.live.erase(it);
		}
	}
	return result;
}

__attribute__((destructor)) static void report()
{
	const Allocations& a = allocations();
	std::fprintf(stderr, "CUDA_ALLOCATION_AUDIT calls=%zu peak_bytes=%zu remaining_bytes=%zu\n", a.calls, a.peak, a.current);
}
