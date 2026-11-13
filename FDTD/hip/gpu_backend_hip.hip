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

#include "hip/hip_runtime.h"
#include <algorithm>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <map>
#include <set>
#include <iostream>
#include <limits>
#include <stdexcept>

#include "hip_internal.h"
#include "gpu_build_id.h"
#include "FDTD/operator.h"

// the raw values of a dump region (see HIP_GatherRegion) for a snapshot: only these go over
// PCIe, not the whole fields; the host interpolates them
__global__ void gather_region(const float* __restrict__ f, HIP_GatherRegion R, HIP_GridDim N, float* __restrict__ out)
{
	const size_t e = (size_t)blockIdx.x*blockDim.x + threadIdx.x;
	if (e<R.total)
		out[e] = f[region_field_index(R, N, e)];
}

// frequency domain dumps: sums[k*total + e] += raw value(e) * w[k] for the weights of one timestep,
// a float times a complex float, then a complex sum; the host interpolates the sums
#define DFT_WEIGHTS 32   // weights per launch, passed by value in the kernel parameters
struct HIP_DFTWeights { float2 w[DFT_WEIGHTS]; };

__global__ void field_dft(const float* __restrict__ f, HIP_GatherRegion R, HIP_GridDim N,
                          HIP_DFTWeights W, unsigned int nw, float2* __restrict__ sums)
{
	const size_t e = (size_t)blockIdx.x*blockDim.x + threadIdx.x;
	if (e>=R.total)
		return;
	const float v = f[region_field_index(R, N, e)];
	for (unsigned int k=0; k<nw; ++k)
	{
		float2 s = sums[k*R.total + e];
		s.x += v*W.w[k].x;
		s.y += v*W.w[k].y;
		sums[k*R.total + e] = s;
	}
}

// Field energy: the squared voltages and currents of the nodes below L, summed along x
// in float, one thread per (y,z) line. The host sums the lines in double.
__global__ void energy_lines(const float* __restrict__ volt, const float* __restrict__ curr, float2* partial,
                             HIP_GridDim N, HIP_GridDim L)
{
	const unsigned int z = blockIdx.x*blockDim.x + threadIdx.x;
	const unsigned int y = blockIdx.y*blockDim.y + threadIdx.y;
	if (z>=L.nz || y>=L.ny)
		return;
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
	partial[y*L.nz + z] = make_float2(E, H);
}

/***************************** helpers *****************************/

namespace
{
struct KernelTimeTable
{
	hipEvent_t start = NULL, stop = NULL;
	std::map<std::string, std::pair<double, size_t>> times;   // total ms, launches
	~KernelTimeTable()
	{
		if (times.empty())
			return;
		double total = 0;
		for (const auto& t : times)
			total += t.second.first;
		std::vector<std::pair<double, std::string>> sorted;
		for (const auto& t : times)
			sorted.push_back(std::make_pair(t.second.first, t.first));
		std::sort(sorted.rbegin(), sorted.rend());
		fprintf(stderr, "HIP kernel times (%.1f ms total):\n", total);
		for (const auto& t : sorted)
		{
			const size_t n = times[t.second].second;
			fprintf(stderr, "  %6.1f %%  %10.1f ms  %8zu x %8.1f us  %s\n", 100*t.first/total, t.first, n, 1000*t.first/n, t.second.c_str());
		}
	}
};
KernelTimeTable& KernelTimes() {static KernelTimeTable table; return table;}
}

bool HIP_KernelTimes::Enabled()
{
	static const bool enabled = getenv("OPENEMS_HIP_KERNEL_TIMES") && (atoi(getenv("OPENEMS_HIP_KERNEL_TIMES"))!=0);
	return enabled;
}

void HIP_KernelTimes::Begin(hipStream_t stream)
{
	KernelTimeTable& t = KernelTimes();
	if (!t.start)
	{
		(void)hipEventCreate(&t.start);
		(void)hipEventCreate(&t.stop);
	}
	(void)hipEventRecord(t.start, stream);
}

void HIP_KernelTimes::End(hipStream_t stream, const char* name)
{
	KernelTimeTable& t = KernelTimes();
	(void)hipEventRecord(t.stop, stream);
	(void)hipEventSynchronize(t.stop);
	float ms = 0;
	(void)hipEventElapsedTime(&ms, t.start, t.stop);
	std::pair<double, size_t>& entry = t.times[name];
	entry.first += ms;
	++entry.second;
}

void HIP_Check(hipError_t err, const char* what)
{
	if (err!=hipSuccess)
		throw std::runtime_error(std::string("GPU_Backend_HIP: ") + what + " failed: " + hipGetErrorString(err));
}

HIP_Context::HIP_Context()
{
	device = 0;
	stream = 0;
	copy_stream = 0;
}

HIP_Context::~HIP_Context()
{
	if (stream)
	{
		(void)hipStreamSynchronize(stream);
		(void)hipStreamDestroy(stream);
	}
	if (copy_stream)
	{
		(void)hipStreamSynchronize(copy_stream);
		(void)hipStreamDestroy(copy_stream);
	}
}

GPU_Backend_HIP::Impl::Impl()
{
	dim.nx = dim.ny = dim.nz = 0;
	numCells = 0;
	volt = curr = NULL;
	energy = NULL;
	energy_count = 0;
	snap_nv = snap_n = 0;
	snap_dev = NULL;
	snap_host = NULL;
	snap_taken = false;
	HIP_Check(hipEventCreateWithFlags(&snap_evaluated, hipEventDisableTiming), "hipEventCreate");
	HIP_Check(hipEventCreateWithFlags(&snap_done, hipEventDisableTiming), "hipEventCreate");
}

GPU_Backend_HIP::Impl::~Impl()
{
	if (ctx)
	{
		(void)hipStreamSynchronize(ctx->stream);
		(void)hipStreamSynchronize(ctx->copy_stream);
	}
	FreeSnapshots();
	(void)hipEventDestroy(snap_evaluated);
	(void)hipEventDestroy(snap_done);
	for (size_t n=0; n<m_Allocations.size(); ++n)
		(void)hipFree(m_Allocations.at(n));
	for (std::set<void*>::iterator it=m_Pinned.begin(); it!=m_Pinned.end(); ++it)
		(void)hipHostUnregister(*it);
}

void GPU_Backend_HIP::Impl::FreeSnapshots()
{
	for (size_t r=0; r<snap_lines.size(); ++r)
		Free(snap_lines[r]);
	snap_lines.clear();
	snap_regions[0].clear();
	snap_regions[1].clear();
	(void)hipFree(snap_dev);
	(void)hipHostFree(snap_host);
	snap_dev = snap_host = NULL;
	snap_nv = snap_n = 0;
	snap_taken = false;
}

void GPU_Backend_HIP::Impl::Flush()
{
	HIP_Check(hipStreamSynchronize(Stream()), "stream synchronization");
}

void GPU_Backend_HIP::Impl::CheckLaunch(const char* kernel)
{
	HIP_Check(hipGetLastError(), kernel);
}

void GPU_Backend_HIP::Impl::PinHostMemory(void* ptr, size_t bytes)
{
	if (m_Pinned.count(ptr))
		return;
	// not fatal: unpinned memory only makes the transfers slower
	if (hipHostRegister(ptr, bytes, hipHostRegisterDefault)==hipSuccess)
		m_Pinned.insert(ptr);
	else
		(void)hipGetLastError();
}

dim3 GPU_Backend_HIP::Impl::Block(size_t ni, size_t nj)
{
	if (nj<=1)
		return dim3(256, 1, 1);
	if (ni<=16)
		return dim3(16, 16, 1);
	return dim3(32, 8, 1);
}

dim3 GPU_Backend_HIP::Impl::Grid(dim3 block, size_t ni, size_t nj, size_t nk)
{
	return dim3((ni+block.x-1)/block.x, (nj+block.y-1)/block.y, nk);
}

/***************************** GPU_Backend_HIP *****************************/

//! Smallest kernel there is, to see whether this build has code for this device (see New())
__global__ void probe_device()
{
}

int GPU_Backend_HIP::DeviceCount()
{
	int count = 0;
	if (hipGetDeviceCount(&count)!=hipSuccess)
	{
		(void)hipGetLastError();
		return 0;
	}
	return count;
}

bool GPU_Backend_HIP::DeviceInfo(int n, GPU_DeviceInfo* info)
{
	hipDeviceProp_t prop;
	if (hipGetDeviceProperties(&prop, n)!=hipSuccess)
	{
		(void)hipGetLastError();
		return false;
	}
	snprintf(info->name, sizeof(info->name), "%s", prop.name);
	info->memory = prop.totalGlobalMem;
	info->integrated = prop.integrated!=0;
	return true;
}

namespace
{
GPU_Backend* NewHIPBackend(int n)
{
	return GPU_Backend_HIP::New(n);
}
}

const GPU_Plugin* GPU_Backend_HIP::Plugin()
{
	static const GPU_Plugin plugin = {OPENEMS_GPU_BUILD_ID, OPENEMS_GPU_PLUGIN_NAME, DeviceCount, DeviceInfo, NewHIPBackend};
	return &plugin;
}

#ifdef OPENEMS_HIP_ARCHITECTURES
namespace
{
//! Whether this build has code for the AMD architecture \a arch (gcnArchName, e.g. "gfx1100:sramecc-:xnack-")
bool BuiltForArch(std::string arch)
{
	arch = arch.substr(0, arch.find(':'));
	// a generic target covers a family of architectures (see the AMDGPU docs of LLVM)
	static const char* const generic[][2] = {
		{"gfx9-generic", "gfx90"}, {"gfx9-4-generic", "gfx94"}, {"gfx9-4-generic", "gfx95"},
		{"gfx10-1-generic", "gfx101"}, {"gfx10-3-generic", "gfx103"}, {"gfx11-generic", "gfx11"}, {"gfx12-generic", "gfx12"}};
	const std::string list = OPENEMS_HIP_ARCHITECTURES;
	size_t start = 0;
	while (start<list.size())
	{
		size_t end = list.find(',', start);
		if (end==std::string::npos)
			end = list.size();
		const std::string target = list.substr(start, std::min(list.find(':', start), end)-start);
		start = end+1;
		if (target==arch)
			return true;
		for (const auto& g : generic)
			if ((target==g[0]) && (arch.compare(0, strlen(g[1]), g[1])==0) && (arch!="gfx908") && (arch!="gfx90a"))
				return true;
	}
	return false;
}
}
#endif

#if defined(__HIP_PLATFORM_AMD__) && defined(__linux__)
namespace
{
// the largest cache of the highest level of the GPU at the PCI address of \a prop from the topology of the
// amdgpu driver (/sys/class/kfd): the Infinity Cache, which HIP's l2CacheSize leaves out; 0 if not found
size_t KFD_LastLevelCache(const hipDeviceProp_t& prop)
{
	const unsigned int location = ((unsigned int)prop.pciBusID << 8) | ((unsigned int)prop.pciDeviceID << 3);
	for (unsigned int node=0; node<256; ++node)
	{
		const std::string dir = "/sys/class/kfd/kfd/topology/nodes/" + std::to_string(node);
		std::ifstream props(dir + "/properties");
		if (!props)
			break;
		long long loc = -1, domain = -1;
		std::string key;
		long long value;
		while (props >> key >> value)
			if (key=="location_id")
				loc = value;
			else if (key=="domain")
				domain = value;
		if ((loc!=(long long)location) || (domain!=(long long)prop.pciDomainID))
			continue;
		unsigned int best_level = 0;
		size_t best = 0;
		for (unsigned int c=0; c<4096; ++c)
		{
			std::ifstream cache(dir + "/caches/" + std::to_string(c) + "/properties");
			if (!cache)
				break;
			unsigned int level = 0;
			size_t kib = 0;
			while (cache >> key >> value)
				if (key=="level")
					level = (unsigned int)value;
				else if (key=="size")
					kib = (size_t)value;
			if ((level>best_level) || ((level==best_level) && (kib*1024>best)))
			{
				best_level = level;
				best = kib*1024;
			}
		}
		return best;
	}
	return 0;
}
}
#endif

GPU_Backend_HIP* GPU_Backend_HIP::New(int device)
{
	if ((device<0) || (device>=DeviceCount()) || (hipSetDevice(device)!=hipSuccess))
	{
		(void)hipGetLastError();
		return NULL;
	}
	Impl* impl = new Impl();
	impl->ctx = std::make_shared<HIP_Context>();
	impl->ctx->device = device;
	hipDeviceProp_t prop;
	HIP_Check(hipGetDeviceProperties(&prop, impl->ctx->device), "hipGetDeviceProperties");
	impl->ctx->name = prop.name;
	// the last-level cache: NVIDIA's L2 is; on AMD the Infinity Cache behind it (RDNA 2 on, CDNA 3) is not in
	// the properties, only in the driver's topology on Linux: else unknown
#if defined(__HIP_PLATFORM_AMD__)
#if defined(__linux__)
	impl->ctx->cache_bytes = KFD_LastLevelCache(prop);
	impl->ctx->cache_source = "the amdgpu topology";
#endif
#else
	impl->ctx->cache_bytes = (prop.l2CacheSize>0) ? (size_t)prop.l2CacheSize : 0;
	impl->ctx->cache_source = "the L2 of the device properties";
#endif
#ifdef OPENEMS_HIP_ARCHITECTURES
	// Some AMD runtimes (ROCm 7.2 with assertions) abort in the query below on a device
	// the build has no code for, instead of returning an error: compare with the list first.
	if (!BuiltForArch(prop.gcnArchName))
	{
		std::cerr << "GPU_Backend_HIP: no kernels for " << impl->ctx->name << " (" << prop.gcnArchName
		          << ") in this build, which is for " << OPENEMS_HIP_ARCHITECTURES << std::endl;
		delete impl;
		return NULL;
	}
#endif
	HIP_Check(hipStreamCreateWithFlags(&impl->ctx->stream, hipStreamNonBlocking), "hipStreamCreate");
	HIP_Check(hipStreamCreateWithFlags(&impl->ctx->copy_stream, hipStreamNonBlocking), "hipStreamCreate");

	// A GPU whose architecture this build has no code for is found and described like any
	// other, and launching on it takes the process down with it rather than returning an
	// error. Ask about a kernel instead of running one: this reports the missing code and
	// leaves the CPU to fall back to.
	hipFuncAttributes attr;
	hipError_t err = hipFuncGetAttributes(&attr, reinterpret_cast<const void*>(probe_device));
	if (err!=hipSuccess)
	{
		// the architecture is the name to add to the build, but only the AMD platform
		// fills it in: the NVIDIA one leaves gcnArchName of the properties untouched
		std::string arch;
#ifdef __HIP_PLATFORM_AMD__
		arch = std::string(" (") + prop.gcnArchName + ")";
#endif
		std::cerr << "GPU_Backend_HIP: no kernels for " << impl->ctx->name << arch
		          << " in this build: " << hipGetErrorString(err) << std::endl;
		(void)hipGetLastError();
		delete impl;
		return NULL;
	}
	return new GPU_Backend_HIP(impl);
}

GPU_Backend* GPU_Backend_HIP::NewSubGridBackend()
{
	Impl* impl = new Impl();
	impl->ctx = d->ctx;
	return new GPU_Backend_HIP(impl);
}

GPU_Backend_HIP::GPU_Backend_HIP(Impl* impl)
{
	d = impl;
}

GPU_Backend_HIP::~GPU_Backend_HIP()
{
	delete d;
}

std::string GPU_Backend_HIP::GetName() const
{
	return "HIP (" + d->ctx->name + ")";
}

size_t GPU_Backend_HIP::CacheBytes(std::string& source) const
{
	source = d->ctx->cache_bytes ? d->ctx->cache_source : std::string();
	return d->ctx->cache_bytes;
}

bool GPU_Backend_HIP::InitFields(const Operator* op)
{
	unsigned int numLines[3];
	for (int n=0; n<3; ++n)
		numLines[n] = op->GetNumberOfLines(n, true);
	d->dim.nx = numLines[0];
	d->dim.ny = numLines[1];
	d->dim.nz = numLines[2];
	d->numCells = (size_t)numLines[0]*numLines[1]*numLines[2];

	// the main updates launch one block per x line along z, which the device limits (64k on NVIDIA): more
	// lines in launches of that many. From the properties: hipDeviceGetAttribute calls the CUDA driver API
	// on NVIDIA, which the plugin does not link (MSVC keeps that call, GCC folds it away)
	hipDeviceProp_t prop;
	HIP_Check(hipGetDeviceProperties(&prop, d->ctx->device), "hipGetDeviceProperties");
	if (prop.maxGridSize[2]>0)
		d->max_grid_z = (unsigned int)prop.maxGridSize[2];

	d->volt = d->Alloc<float>(3*d->numCells);
	d->curr = d->Alloc<float>(3*d->numCells);
	return true;
}

static void CheckSize(size_t host_count, size_t device_count)
{
	if (host_count!=device_count)
		throw std::runtime_error("GPU_Backend_HIP: host and device field size mismatch");
}

// the host mirror is page-locked on first use, for the full transfer speed
void GPU_Backend_HIP::Download(bool currents, ArrayLib::ArrayNIJK<FDTD_FLOAT>& field)
{
	CheckSize(field.size(), 3*d->numCells);
	d->PinHostMemory(field.data(), field.size()*sizeof(FDTD_FLOAT));
	HIP_Check(hipMemcpyAsync(field.data(), currents ? d->curr : d->volt, field.size()*sizeof(float), hipMemcpyDeviceToHost, d->Stream()), "download");
	d->Flush();
}

void GPU_Backend_HIP::Synchronize()
{
	d->Flush();
}

void* GPU_Backend_HIP::Impl::UploadRegion(const Engine_GatherRegion& region, HIP_GatherRegion& dev)
{
	std::vector<unsigned int> lines;
	for (int n=0; n<3; ++n)
		for (int k=0; k<3; ++k)
			lines.insert(lines.end(), region.lines[n][k].begin(), region.lines[n][k].end());
	unsigned int* buf = TryAlloc<unsigned int>(lines.size());
	if (!buf)
		return NULL;
	HIP_Check(hipMemcpy(buf, lines.data(), lines.size()*sizeof(unsigned int), hipMemcpyHostToDevice), "hipMemcpy");
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

bool GPU_Backend_HIP::SetSnapshotGather(const std::vector<Engine_GatherRegion>& volt_regions, const std::vector<Engine_GatherRegion>& curr_regions)
{
	HIP_Check(hipStreamSynchronize(d->ctx->copy_stream), "stream synchronization");
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
	bool ok = (hipMalloc(&d->snap_dev, n*sizeof(float))==hipSuccess)
	          && (hipHostMalloc(&d->snap_host, n*sizeof(float))==hipSuccess);
	const std::vector<Engine_GatherRegion>* regions[2] = {&volt_regions, &curr_regions};
	for (int f=0; f<2 && ok; ++f)
		for (size_t r=0; r<regions[f]->size() && ok; ++r)
		{
			HIP_GatherRegion dev;
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
		(void)hipGetLastError();
		d->FreeSnapshots();
		std::cerr << "GPU_Backend_HIP: not enough memory for field snapshots, the field dumps wait for the device" << std::endl;
		return false;
	}
	d->snap_nv = nv;
	d->snap_n = n;
	return true;
}

bool GPU_Backend_HIP::SnapshotFields(const FDTD_FLOAT* &volt, const FDTD_FLOAT* &curr)
{
	if (!d->snap_n)
		return false;
	// the raw values on the work stream, once the last snapshot was downloaded from snap_dev
	if (d->snap_taken)
		HIP_Check(hipStreamWaitEvent(d->Stream(), d->snap_done, 0), "hipStreamWaitEvent");
	size_t offset = 0;
	for (int f=0; f<2; ++f)
		for (size_t r=0; r<d->snap_regions[f].size(); ++r)
		{
			const HIP_GatherRegion& R = d->snap_regions[f][r];
			HIP_Launch(d, "gather_region", gather_region, R.total, 1, 1, (const float*)(f ? d->curr : d->volt), R, d->dim, d->snap_dev + offset);
			offset += R.total;
		}
	HIP_Check(hipEventRecord(d->snap_evaluated, d->Stream()), "hipEventRecord");
	// the download on the copy stream, while the work stream continues with the next timesteps
	hipStream_t copy = d->ctx->copy_stream;
	HIP_Check(hipStreamWaitEvent(copy, d->snap_evaluated, 0), "hipStreamWaitEvent");
	HIP_Check(hipMemcpyAsync(d->snap_host, d->snap_dev, d->snap_n*sizeof(float), hipMemcpyDeviceToHost, copy), "snapshot download");
	HIP_Check(hipEventRecord(d->snap_done, copy), "hipEventRecord");
	d->snap_taken = true;
	volt = d->snap_host;
	curr = d->snap_host + d->snap_nv;
	return true;
}

int GPU_Backend_HIP::AddFieldDFT(bool h_field, const Engine_GatherRegion& region, unsigned int count)
{
	const size_t n = region.Size();
	if ((n==0) || (count==0))
		return -1;
	Impl::FieldDFT dft;
	dft.count = count;
	dft.currents = h_field;
	void* lines = d->UploadRegion(region, dft.region);
	dft.sums = lines ? d->TryAlloc<float2>(n*count) : NULL;
	if (!dft.sums)
	{
		// not fatal: without device sums the dump sums the fields on the host
		d->Free(lines);
		std::cerr << "GPU_Backend_HIP: not enough memory for a frequency domain dump on the device" << std::endl;
		return -1;
	}
	HIP_Check(hipMemsetAsync(dft.sums, 0, n*count*sizeof(float2), d->Stream()), "hipMemset");
	d->field_dfts.push_back(dft);
	return d->field_dfts.size()-1;
}

void GPU_Backend_HIP::AccumulateFieldDFT(int id, const std::vector<std::complex<float>>& weights)
{
	if ((id<0) || ((size_t)id>=d->field_dfts.size()))
		return;
	const Impl::FieldDFT& dft = d->field_dfts[id];
	const unsigned int count = std::min((size_t)dft.count, weights.size());
	for (unsigned int k0=0; k0<count; k0+=DFT_WEIGHTS)
	{
		HIP_DFTWeights W;
		const unsigned int nw = std::min(count-k0, (unsigned int)DFT_WEIGHTS);
		for (unsigned int k=0; k<nw; ++k)
			W.w[k] = make_float2(weights[k0+k].real(), weights[k0+k].imag());
		HIP_Launch(d, "field_dft", field_dft, dft.region.total, 1, 1, (const float*)(dft.currents ? d->curr : d->volt),
		            dft.region, d->dim, W, nw, dft.sums + k0*dft.region.total);
	}
}

bool GPU_Backend_HIP::ReadFieldDFT(int id, unsigned int k, std::vector<std::complex<float>>& sums)
{
	if ((id<0) || ((size_t)id>=d->field_dfts.size()) || (k>=d->field_dfts[id].count))
		return false;
	const Impl::FieldDFT& dft = d->field_dfts[id];
	sums.resize(dft.region.total);
	HIP_Check(hipMemcpyAsync(sums.data(), dft.sums + k*dft.region.total, sums.size()*sizeof(float2), hipMemcpyDeviceToHost, d->Stream()), "download");
	d->Flush();
	return true;
}

void GPU_Backend_HIP::WaitSnapshot()
{
	// called on the dump thread: report errors instead of throwing
	const hipError_t err = hipEventSynchronize(d->snap_done);
	if (err!=hipSuccess)
		std::cerr << "GPU_Backend_HIP: field snapshot failed: " << hipGetErrorString(err) << std::endl;
}

bool GPU_Backend_HIP::CalcFastEnergy(const unsigned int numNodes[3], double& E_energy, double& H_energy)
{
	E_energy = H_energy = 0;
	HIP_GridDim L = {numNodes[0], numNodes[1], numNodes[2]};
	if ((L.nx>d->dim.nx) || (L.ny>d->dim.ny) || (L.nz>d->dim.nz))
		return false;
	const size_t count = (size_t)L.ny*L.nz;
	if (L.nx==0 || count==0)
		return true;
	if (d->energy_count<count)
	{
		d->energy = d->Alloc<float2>(count);
		d->energy_count = count;
		d->energy_host.resize(count);   // small, not pinned
	}
	HIP_Launch(d, "energy_lines", energy_lines, L.nz, L.ny, 1, (const float*)d->volt, (const float*)d->curr, d->energy, d->dim, L);
	HIP_Check(hipMemcpyAsync(d->energy_host.data(), d->energy, count*sizeof(float2), hipMemcpyDeviceToHost, d->Stream()), "energy download");
	d->Flush();
	for (size_t k=0; k<count; ++k)
	{
		E_energy += d->energy_host[k].x;
		H_energy += d->energy_host[k].y;
	}
	return true;
}
