/*
*	Copyright (C) 2026 Thorsten Liebig (Thorsten.Liebig@gmx.de)
*	Copyright (C) 2026 Alexandru Iercosan-Lucaci (alex.lucaci@alacrity.ro)
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
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

#include "sycl_internal.h"
#include "gpu_build_id.h"
#include "FDTD/operator.h"

/***************************** helpers *****************************/

SYCL_Context::SYCL_Context(const sycl::device& dev) : device(dev)
{
	// the asynchronous errors (a failed kernel) of both queues, thrown by ThrowAsyncError()
	auto handler = [this](sycl::exception_list errors) {
		for (const std::exception_ptr& e : errors)
		{
			try
			{
				std::rethrow_exception(e);
			}
			catch (const std::exception& ex)
			{
				std::lock_guard<std::mutex> lock(m_ErrorMutex);
				if (m_AsyncError.empty())
					m_AsyncError = ex.what();
			}
		}
	};
	context = sycl::context(device, handler);
	queue = sycl::queue(context, device, handler, sycl::property::queue::in_order());
	copy_queue = sycl::queue(context, device, handler, sycl::property::queue::in_order());
	name = device.get_info<sycl::info::device::name>();
	cache_bytes = device.get_info<sycl::info::device::global_mem_cache_size>();
	max_alloc = device.get_info<sycl::info::device::max_mem_alloc_size>();
	max_group = device.get_info<sycl::info::device::max_work_group_size>();
}

SYCL_Context::~SYCL_Context()
{
	try
	{
		queue.wait();
		copy_queue.wait();
		// an asynchronous error still pending goes through the handler, which writes
		// members destroyed before the queues: deliver it here
		queue.throw_asynchronous();
		copy_queue.throw_asynchronous();
	}
	catch (...)
	{
	}
}

void SYCL_Context::ThrowAsyncError()
{
	std::lock_guard<std::mutex> lock(m_ErrorMutex);
	if (!m_AsyncError.empty())
		throw std::runtime_error("GPU_Backend_SYCL: device error: " + m_AsyncError);
}

GPU_Backend_SYCL::Impl::Impl()
{
	dim.nx = dim.ny = dim.nz = 0;
}

GPU_Backend_SYCL::Impl::~Impl()
{
	if (!ctx)   // the context failed (New())
		return;
	// nothing may throw here
	try
	{
		Queue().wait();
		ctx->copy_queue.wait();
	}
	catch (...)
	{
	}
	FreeSnapshots();
	for (void* ptr : m_Allocations)
		sycl::free(ptr, ctx->context);
#ifdef SYCL_EXT_ONEAPI_COPY_OPTIMIZE
	for (void* ptr : m_Pinned)
	{
		try
		{
			sycl::ext::oneapi::experimental::release_from_device_copy(ptr, ctx->context);
		}
		catch (...)
		{
		}
	}
#endif
}

void GPU_Backend_SYCL::Impl::Free(void* ptr)
{
	if (ptr==NULL)
		return;
	// unlike hipFree, sycl::free does not wait for the work that uses the buffer
	try
	{
		Queue().wait();
	}
	catch (...)
	{
	}
	m_Allocations.erase(std::remove(m_Allocations.begin(), m_Allocations.end(), ptr), m_Allocations.end());
	sycl::free(ptr, ctx->context);
}

void GPU_Backend_SYCL::Impl::FreeSnapshots()
{
	try
	{
		ctx->copy_queue.wait();
	}
	catch (...)
	{
	}
	for (size_t r=0; r<snap_lines.size(); ++r)
		Free(snap_lines[r]);
	snap_lines.clear();
	snap_regions[0].clear();
	snap_regions[1].clear();
	Free(snap_dev);
	if (snap_host)
		sycl::free(snap_host, ctx->context);
	snap_dev = snap_host = NULL;
	snap_nv = snap_n = 0;
	snap_taken = false;
}

void GPU_Backend_SYCL::Impl::Flush()
{
	SYCL_Check("queue synchronization", [&] {Queue().wait_and_throw();});
	ctx->ThrowAsyncError();
}

void GPU_Backend_SYCL::Impl::PinHostMemory(void* ptr, size_t bytes)
{
#ifdef SYCL_EXT_ONEAPI_COPY_OPTIMIZE
	if (m_Pinned.count(ptr))
		return;
	// not fatal: unpinned memory only makes the transfers slower
	try
	{
		sycl::ext::oneapi::experimental::prepare_for_device_copy(ptr, bytes, ctx->context);
		m_Pinned.insert(ptr);
	}
	catch (...)
	{
	}
#else
	(void)ptr;
	(void)bytes;
#endif
}

/***************************** devices *****************************/

namespace
{
//! The Level Zero GPUs, found once; empty without a runtime or driver. With ONEAPI_DEVICE_SELECTOR
//! set, every device it selects instead, e.g. opencl:cpu to test the backend on a machine without
//! an Intel GPU (far slower than the CPU engines) or hip:gpu / cuda:gpu with a runtime that has them.
const std::vector<sycl::device>& Devices()
{
	// never freed: the SYCL runtime may be torn down before the static destructors run
	static const std::vector<sycl::device>* devices = []()
	{
		// the kernels are compiled for the device when the GPU engine starts (SPIR-V), seconds on
		// a GPU: keep them in SYCL's cache on disk, unless the user decided otherwise. Before the
		// first call into the runtime, which reads its configuration once.
		if (getenv("SYCL_CACHE_PERSISTENT")==NULL)
		{
#ifdef _WIN32
			_putenv_s("SYCL_CACHE_PERSISTENT", "1");
#else
			setenv("SYCL_CACHE_PERSISTENT", "1", 0);
#endif
		}
		std::vector<sycl::device>* list = new std::vector<sycl::device>();
		try
		{
			const bool selected = (getenv("ONEAPI_DEVICE_SELECTOR")!=NULL);
			for (const sycl::device& dev : sycl::device::get_devices())
				if (selected || (dev.is_gpu() && (dev.get_backend()==sycl::backend::ext_oneapi_level_zero)))
					list->push_back(dev);
		}
		catch (...)
		{
			list->clear();
		}
		return list;
	}();
	return *devices;
}

GPU_Backend* NewSYCLBackend(int n)
{
	return GPU_Backend_SYCL::New(n);
}
}

int GPU_Backend_SYCL::DeviceCount()
{
	return (int)Devices().size();
}

bool GPU_Backend_SYCL::DeviceInfo(int n, GPU_DeviceInfo* info)
{
	if ((n<0) || (n>=DeviceCount()))
		return false;
	try
	{
		const sycl::device& dev = Devices().at(n);
		snprintf(info->name, sizeof(info->name), "%s", dev.get_info<sycl::info::device::name>().c_str());
		info->memory = dev.get_info<sycl::info::device::global_mem_size>();
		// deprecated in SYCL 2020, without a replacement that says the same: no aspect tells
		// that the device memory is the host memory
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
		info->integrated = dev.get_info<sycl::info::device::host_unified_memory>();
#pragma clang diagnostic pop
		return true;
	}
	catch (...)
	{
		return false;
	}
}

const GPU_Plugin* GPU_Backend_SYCL::Plugin()
{
	static const GPU_Plugin plugin = {OPENEMS_GPU_BUILD_ID, OPENEMS_GPU_PLUGIN_NAME, DeviceCount, DeviceInfo, NewSYCLBackend};
	return &plugin;
}

GPU_Backend_SYCL* GPU_Backend_SYCL::New(int device)
{
	// called through the plugin entry: nothing may throw from here
	if ((device<0) || (device>=DeviceCount()))
		return NULL;
	Impl* impl = NULL;
	std::string name = "device " + std::to_string(device);
	try
	{
		const sycl::device& dev = Devices().at(device);
		name = dev.get_info<sycl::info::device::name>();
		impl = new Impl();
		impl->ctx = std::make_shared<SYCL_Context>(dev);
		// A device the build has no code for (no ahead-of-time image, and the SPIR-V does
		// not compile for it) fails when the first kernel is built for it: run the smallest
		// kernel there is, so that the CPU is left to fall back to.
		int* probe = impl->TryAlloc<int>(1);
		if (probe==NULL)
			throw std::runtime_error("no device memory");
		impl->Queue().single_task([=]() {*probe = 1;});
		impl->Flush();
	}
	catch (const std::exception& e)
	{
		std::cerr << "GPU_Backend_SYCL: cannot run on " << name << ": " << e.what() << std::endl;
		delete impl;
		return NULL;
	}
	catch (...)
	{
		std::cerr << "GPU_Backend_SYCL: cannot run on " << name << std::endl;
		delete impl;
		return NULL;
	}
	return new GPU_Backend_SYCL(impl);
}

/***************************** GPU_Backend_SYCL *****************************/

GPU_Backend* GPU_Backend_SYCL::NewSubGridBackend()
{
	Impl* impl = new Impl();
	impl->ctx = d->ctx;
	return new GPU_Backend_SYCL(impl);
}

GPU_Backend_SYCL::GPU_Backend_SYCL(Impl* impl)
{
	d = impl;
}

GPU_Backend_SYCL::~GPU_Backend_SYCL()
{
	delete d;
}

std::string GPU_Backend_SYCL::GetName() const
{
	return "SYCL (" + d->ctx->name + ")";
}

size_t GPU_Backend_SYCL::CacheBytes(std::string& source) const
{
	source = d->ctx->cache_bytes ? "the global memory cache of the SYCL device" : "";
	return d->ctx->cache_bytes;
}

bool GPU_Backend_SYCL::InitFields(const Operator* op)
{
	unsigned int numLines[3];
	for (int n=0; n<3; ++n)
		numLines[n] = op->GetNumberOfLines(n, true);
	d->dim.nx = numLines[0];
	d->dim.ny = numLines[1];
	d->dim.nz = numLines[2];
	d->numCells = (size_t)numLines[0]*numLines[1]*numLines[2];
	// one allocation per field, up to the device's limit: 4 GiB on many Intel GPUs, the
	// whole memory on a B580, whose driver addresses larger buffers (tested past 4 GiB)
	if (3*d->numCells*sizeof(float) > d->ctx->max_alloc)
	{
		std::cerr << "GPU_Backend_SYCL: a field of " << 3*d->numCells*sizeof(float)/1048576 << " MiB is more than the largest allocation of "
		          << d->ctx->name << " (" << d->ctx->max_alloc/1048576 << " MiB)" << std::endl;
		return false;
	}
	d->volt = d->Alloc<float>(3*d->numCells);
	d->curr = d->Alloc<float>(3*d->numCells);
	return true;
}

static void CheckSize(size_t host_count, size_t device_count)
{
	if (host_count!=device_count)
		throw std::runtime_error("GPU_Backend_SYCL: host and device field size mismatch");
}

void GPU_Backend_SYCL::Download(bool currents, ArrayLib::ArrayNIJK<FDTD_FLOAT>& field)
{
	CheckSize(field.size(), 3*d->numCells);
	d->PinHostMemory(field.data(), field.size()*sizeof(FDTD_FLOAT));
	SYCL_Check("download", [&] {d->Queue().memcpy(field.data(), currents ? d->curr : d->volt, field.size()*sizeof(float));});
	d->Flush();
}

void GPU_Backend_SYCL::Synchronize()
{
	d->Flush();
}

void* GPU_Backend_SYCL::Impl::UploadRegion(const Engine_GatherRegion& region, SYCL_GatherRegion& dev)
{
	std::vector<unsigned int> lines;
	for (int n=0; n<3; ++n)
		for (int k=0; k<3; ++k)
			lines.insert(lines.end(), region.lines[n][k].begin(), region.lines[n][k].end());
	unsigned int* buf = TryAlloc<unsigned int>(lines.size());
	if (!buf)
		return NULL;
	SYCL_Check("upload", [&] {Queue().memcpy(buf, lines.data(), lines.size()*sizeof(unsigned int));});
	Flush();
	size_t pos = 0, start = 0;
	for (int n=0; n<3; ++n)
	{
		for (int k=0; k<3; ++k)
		{
			dev.lines[n][k] = buf + pos;
			dev.size[n][k] = region.lines[n][k].size();
			pos += dev.size[n][k];
		}
		dev.start[n] = start;
		start += region.Size(n);
	}
	dev.total = start;
	return buf;
}

bool GPU_Backend_SYCL::SetSnapshotGather(const std::vector<Engine_GatherRegion>& volt_regions, const std::vector<Engine_GatherRegion>& curr_regions)
{
	d->FreeSnapshots();
	size_t nv = 0, n = 0;
	for (size_t r=0; r<volt_regions.size(); ++r)
		nv += volt_regions[r].Size();
	n = nv;
	for (size_t r=0; r<curr_regions.size(); ++r)
		n += curr_regions[r].Size();
	if (n==0)
		return false;
	// not fatal: without snapshots the dumps read the fields from the host mirror
	d->snap_dev = d->TryAlloc<float>(n);
	try
	{
		d->snap_host = sycl::malloc_host<float>(n, d->ctx->context);
	}
	catch (const sycl::exception&)
	{
		d->snap_host = NULL;
	}
	bool ok = d->snap_dev && d->snap_host;
	const std::vector<Engine_GatherRegion>* regions[2] = {&volt_regions, &curr_regions};
	for (int f=0; f<2 && ok; ++f)
		for (size_t r=0; r<regions[f]->size() && ok; ++r)
		{
			SYCL_GatherRegion dev;
			void* buf = d->UploadRegion(regions[f]->at(r), dev);
			ok = (buf!=NULL);
			if (ok)
			{
				d->snap_lines.push_back(buf);
				d->snap_regions[f].push_back(dev);
			}
		}
	if (!ok)
	{
		d->FreeSnapshots();
		std::cerr << "GPU_Backend_SYCL: not enough memory for field snapshots, the field dumps wait for the device" << std::endl;
		return false;
	}
	d->snap_nv = nv;
	d->snap_n = n;
	return true;
}

bool GPU_Backend_SYCL::SnapshotFields(const FDTD_FLOAT* &volt, const FDTD_FLOAT* &curr)
{
	if (!d->snap_n)
		return false;
	// the raw values on the work queue, once the last snapshot was downloaded from snap_dev
	std::vector<sycl::event> after;
	if (d->snap_taken)
		after.push_back(d->snap_done);
	size_t offset = 0;
	sycl::event evaluated;
	const SYCL_GridDim N = d->dim;
	for (int f=0; f<2; ++f)
		for (size_t r=0; r<d->snap_regions[f].size(); ++r)
		{
			const SYCL_GatherRegion R = d->snap_regions[f][r];
			if (R.total==0)
				continue;
			const float* field = f ? d->curr : d->volt;
			float* out = d->snap_dev + offset;
			evaluated = SYCL_Launch(d, "gather_region", R.total, 1, 1, [=](size_t e, size_t, size_t) {
				out[e] = field[region_field_index(R, N, e)];
			}, after);
			after.clear();   // the queue is in order
			offset += R.total;
		}
	// the download on the copy queue, while the work queue continues with the next timesteps
	d->snap_done = SYCL_Check("snapshot download", [&] {return d->ctx->copy_queue.memcpy(d->snap_host, d->snap_dev, d->snap_n*sizeof(float), evaluated);});
	d->snap_taken = true;
	volt = d->snap_host;
	curr = d->snap_host + d->snap_nv;
	return true;
}

void GPU_Backend_SYCL::WaitSnapshot()
{
	// called on the dump thread: report errors instead of throwing
	try
	{
		d->snap_done.wait_and_throw();
	}
	catch (const std::exception& e)
	{
		std::cerr << "GPU_Backend_SYCL: field snapshot failed: " << e.what() << std::endl;
	}
}

//! Weights per launch of the frequency domain sums, passed by value with the kernel
#define DFT_WEIGHTS 32
struct SYCL_DFTWeights
{
	SYCL_Float2 w[DFT_WEIGHTS];
};

int GPU_Backend_SYCL::AddFieldDFT(bool h_field, const Engine_GatherRegion& region, unsigned int count)
{
	const size_t n = region.Size();
	if ((n==0) || (count==0))
		return -1;
	Impl::FieldDFT dft;
	dft.count = count;
	dft.currents = h_field;
	void* lines = d->UploadRegion(region, dft.region);
	dft.sums = lines ? d->TryAlloc<SYCL_Float2>(n*count) : NULL;
	if (!dft.sums)
	{
		// not fatal: without device sums the dump sums the fields on the host
		d->Free(lines);
		std::cerr << "GPU_Backend_SYCL: not enough memory for a frequency domain dump on the device" << std::endl;
		return -1;
	}
	SYCL_Check("memset", [&] {d->Queue().memset(dft.sums, 0, n*count*sizeof(SYCL_Float2));});
	d->field_dfts.push_back(dft);
	return d->field_dfts.size()-1;
}

// frequency domain dumps: sums[k*total + e] += raw value(e) * w[k] for the weights of one timestep,
// a float times a complex float, then a complex sum; the host interpolates the sums
void GPU_Backend_SYCL::AccumulateFieldDFT(int id, const std::vector<std::complex<float>>& weights)
{
	if ((id<0) || ((size_t)id>=d->field_dfts.size()))
		return;
	const Impl::FieldDFT& dft = d->field_dfts[id];
	const unsigned int count = std::min((size_t)dft.count, weights.size());
	const SYCL_GridDim N = d->dim;
	const SYCL_GatherRegion R = dft.region;
	const float* f = dft.currents ? d->curr : d->volt;
	for (unsigned int k0=0; k0<count; k0+=DFT_WEIGHTS)
	{
		SYCL_DFTWeights W;
		const unsigned int nw = std::min(count-k0, (unsigned int)DFT_WEIGHTS);
		for (unsigned int k=0; k<nw; ++k)
			W.w[k] = SYCL_Float2{weights[k0+k].real(), weights[k0+k].imag()};
		SYCL_Float2* sums = dft.sums + k0*R.total;
		SYCL_Launch(d, "field_dft", R.total, 1, 1, [=](size_t e, size_t, size_t) {
			const float v = f[region_field_index(R, N, e)];
			for (unsigned int k=0; k<nw; ++k)
			{
				SYCL_Float2 s = sums[k*R.total + e];
				s.x += v*W.w[k].x;
				s.y += v*W.w[k].y;
				sums[k*R.total + e] = s;
			}
		});
	}
}

bool GPU_Backend_SYCL::ReadFieldDFT(int id, unsigned int k, std::vector<std::complex<float>>& sums)
{
	if ((id<0) || ((size_t)id>=d->field_dfts.size()) || (k>=d->field_dfts[id].count))
		return false;
	const Impl::FieldDFT& dft = d->field_dfts[id];
	sums.resize(dft.region.total);
	SYCL_Check("download", [&] {d->Queue().memcpy(sums.data(), dft.sums + k*dft.region.total, sums.size()*sizeof(SYCL_Float2));});
	d->Flush();
	return true;
}

// Field energy: the squared voltages and currents of the nodes below L, summed along x
// in float, one work-item per (y,z) line. The host sums the lines in double.
bool GPU_Backend_SYCL::CalcFastEnergy(const unsigned int numNodes[3], double& E_energy, double& H_energy)
{
	E_energy = H_energy = 0;
	const SYCL_GridDim L = {numNodes[0], numNodes[1], numNodes[2]};
	const SYCL_GridDim N = d->dim;
	if ((L.nx>N.nx) || (L.ny>N.ny) || (L.nz>N.nz))
		return false;
	const size_t count = (size_t)L.ny*L.nz;
	if (L.nx==0 || count==0)
		return true;
	if (d->energy_count<count)
	{
		d->energy = d->Alloc<SYCL_Float2>(count);
		d->energy_count = count;
		d->energy_host.resize(count);
	}
	const float* volt = d->volt;
	const float* curr = d->curr;
	SYCL_Float2* partial = d->energy;
	SYCL_Launch(d, "energy_lines", L.nz, L.ny, 1, [=](size_t z, size_t y, size_t) {
		const size_t sn = (size_t)N.nx*N.ny*N.nz;
		float E = 0, H = 0;
		for (unsigned int x=0; x<L.nx; ++x)
		{
			const size_t i = nijk(N, 0, x, y, z);
			E += volt[i]*volt[i];
			E += volt[sn+i]*volt[sn+i];
			E += volt[2*sn+i]*volt[2*sn+i];
			H += curr[i]*curr[i];
			H += curr[sn+i]*curr[sn+i];
			H += curr[2*sn+i]*curr[2*sn+i];
		}
		partial[y*L.nz + z] = SYCL_Float2{E, H};
	});
	SYCL_Check("energy download", [&] {d->Queue().memcpy(d->energy_host.data(), d->energy, count*sizeof(SYCL_Float2));});
	d->Flush();
	for (size_t k=0; k<count; ++k)
	{
		E_energy += d->energy_host[k].x;
		H_energy += d->energy_host[k].y;
	}
	return true;
}
