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
#include <limits>
#include <stdexcept>

#include "hip_internal.h"
#include "FDTD/operator.h"

// Main FDTD updates, operation by operation the same as Engine::UpdateVoltages/UpdateCurrents,
// over the nodes in [B, E). One thread per mesh node: z along x-threads (contiguous), y along
// y-threads, x along the grid z-dimension. The coefficients are read with main_coeff().
__global__ void update_voltages(float* __restrict__ volt, const float* __restrict__ curr,
                                const void* index, const float* ca, const float* cb, unsigned int mode,
                                HIP_GridDim N, HIP_GridDim B, HIP_GridDim E)
{
	const unsigned int z = B.nz + blockIdx.x*blockDim.x + threadIdx.x;
	const unsigned int y = B.ny + blockIdx.y*blockDim.y + threadIdx.y;
	const unsigned int x = B.nx + blockIdx.z;
	if (z>=E.nz || y>=E.ny || x>=E.nx)
		return;
	const unsigned int sn = N.nx*N.ny*N.nz;
	const unsigned int i  = nijk(N, 0, x, y, z);
	const unsigned int xm = (x>0) ? N.ny*N.nz : 0;   // shift to the previous line, none on the first
	const unsigned int ym = (y>0) ? N.nz : 0;
	const unsigned int zm = (z>0) ? 1 : 0;
	const HIP_MainCoeff C = main_coeff(index, ca, cb, mode, sn, i, 0);
	const float* vv = C.a;
	const float* vi = C.b;
	const unsigned int s = C.s;
	float v;

	//for x
	v  = volt[i] * vv[0];
	v += vi[0] * (curr[2*sn+i] - curr[2*sn+i-ym] - curr[sn+i] + curr[sn+i-zm]);
	volt[i] = v;

	//for y
	v  = volt[sn+i] * vv[s];
	v += vi[s] * (curr[i] - curr[i-zm] - curr[2*sn+i] + curr[2*sn+i-xm]);
	volt[sn+i] = v;

	//for z
	v  = volt[2*sn+i] * vv[2*s];
	v += vi[2*s] * (curr[sn+i] - curr[sn+i-xm] - curr[i] + curr[i-ym]);
	volt[2*sn+i] = v;
}

// the currents on the last mesh line of each direction are not updated
__global__ void update_currents(float* __restrict__ curr, const float* __restrict__ volt,
                                const void* index, const float* ca, const float* cb, unsigned int mode,
                                HIP_GridDim N, HIP_GridDim B, HIP_GridDim E)
{
	const unsigned int z = B.nz + blockIdx.x*blockDim.x + threadIdx.x;
	const unsigned int y = B.ny + blockIdx.y*blockDim.y + threadIdx.y;
	const unsigned int x = B.nx + blockIdx.z;
	if (z>=E.nz || y>=E.ny || x>=E.nx || z+1>=N.nz || y+1>=N.ny || x+1>=N.nx)
		return;
	const unsigned int sn = N.nx*N.ny*N.nz;
	const unsigned int i  = nijk(N, 0, x, y, z);
	const unsigned int xp = N.ny*N.nz;
	const unsigned int yp = N.nz;
	const HIP_MainCoeff C = main_coeff(index, ca, cb, mode, sn, i, 6);
	const float* ii = C.a;
	const float* iv = C.b;
	const unsigned int s = C.s;
	float c;

	//for x
	c  = curr[i] * ii[0];
	c += iv[0] * (volt[2*sn+i] - volt[2*sn+i+yp] - volt[sn+i] + volt[sn+i+1]);
	curr[i] = c;

	//for y
	c  = curr[sn+i] * ii[s];
	c += iv[s] * (volt[i] - volt[i+1] - volt[2*sn+i] + volt[2*sn+i+xp]);
	curr[sn+i] = c;

	//for z
	c  = curr[2*sn+i] * ii[2*s];
	c += iv[2*s] * (volt[sn+i] - volt[sn+i+xp] - volt[i] + volt[i+yp]);
	curr[2*sn+i] = c;
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
	const unsigned int sn = N.nx*N.ny*N.nz;
	float E = 0, H = 0;
	for (unsigned int x=0; x<L.nx; ++x)
	{
		const unsigned int i = nijk(N, 0, x, y, z);
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

void HIP_Check(hipError_t err, const char* what)
{
	if (err!=hipSuccess)
		throw std::runtime_error(std::string("GPU_Backend_HIP: ") + what + " failed: " + hipGetErrorString(err));
}

HIP_Context::HIP_Context()
{
	device = 0;
	stream = 0;
}

HIP_Context::~HIP_Context()
{
	if (stream)
	{
		hipStreamSynchronize(stream);
		hipStreamDestroy(stream);
	}
}

GPU_Backend_HIP::Impl::Impl()
{
	dim.nx = dim.ny = dim.nz = 0;
	numCells = 0;
	volt = curr = vv = vi = ii = iv = NULL;
	coeff_mode = 0;
	index = NULL;
	coeff = NULL;
	energy = NULL;
	energy_count = 0;
}

GPU_Backend_HIP::Impl::~Impl()
{
	if (ctx)
		hipStreamSynchronize(ctx->stream);
	for (size_t n=0; n<m_Allocations.size(); ++n)
		hipFree(m_Allocations.at(n));
	for (std::set<void*>::iterator it=m_Pinned.begin(); it!=m_Pinned.end(); ++it)
		hipHostUnregister(*it);
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
		hipGetLastError();
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

GPU_Backend_HIP* GPU_Backend_HIP::New()
{
	int count = 0;
	if ((hipGetDeviceCount(&count)!=hipSuccess) || (count==0))
	{
		hipGetLastError();
		return NULL;
	}
	Impl* impl = new Impl();
	impl->ctx = std::make_shared<HIP_Context>();
	HIP_Check(hipGetDevice(&impl->ctx->device), "hipGetDevice");
	hipDeviceProp prop;
	HIP_Check(hipGetDeviceProperties(&prop, impl->ctx->device), "hipGetDeviceProperties");
	impl->ctx->name = prop.name;
	HIP_Check(hipStreamCreateWithFlags(&impl->ctx->stream, hipStreamNonBlocking), "hipStreamCreate");
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
	return "CUDA (" + d->ctx->name + ")";
}

bool GPU_Backend_HIP::Init(const Operator* op)
{
	unsigned int numLines[3];
	for (int n=0; n<3; ++n)
		numLines[n] = op->GetNumberOfLines(n, true);
	d->dim.nx = numLines[0];
	d->dim.ny = numLines[1];
	d->dim.nz = numLines[2];
	d->numCells = (size_t)numLines[0]*numLines[1]*numLines[2];
	d->main_start.nx = d->main_start.ny = d->main_start.nz = 0;
	d->main_stop = d->dim;
	d->upml_fused = -1;

	// the kernels index with 32 bit
	if (3*d->numCells > std::numeric_limits<unsigned int>::max())
	{
		std::cerr << "GPU_Backend_HIP::Init: Error: the mesh is too large for 32 bit indexing" << std::endl;
		return false;
	}

	const size_t count = 3*d->numCells;
	d->volt = d->Alloc<float>(count);
	d->curr = d->Alloc<float>(count);

	// upload the final operator coefficients, including all changes by operator extensions
	std::vector<float> vv(count), vi(count), ii(count), iv(count);
	unsigned int pos[3];
	size_t idx = 0;
	for (int n=0; n<3; ++n)
		for (pos[0]=0; pos[0]<numLines[0]; ++pos[0])
			for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
				for (pos[2]=0; pos[2]<numLines[2]; ++pos[2], ++idx)
				{
					vv[idx] = op->GetVV(n, pos[0], pos[1], pos[2]);
					vi[idx] = op->GetVI(n, pos[0], pos[1], pos[2]);
					ii[idx] = op->GetII(n, pos[0], pos[1], pos[2]);
					iv[idx] = op->GetIV(n, pos[0], pos[1], pos[2]);
				}

	// Most nodes share one of a few coefficient sets (same material and mesh spacing): store
	// every distinct set of the 12 coefficients of a node once and a set index per node
	const size_t sn = d->numCells;
	const std::vector<float>* src[4] = {&vv, &vi, &ii, &iv};
	GPU_CoeffSets sets;
	if (GPU_FindSets(sn, 12, [&](size_t i, float* values)
	    {
		    for (int c=0; c<4; ++c)
			    for (int n=0; n<3; ++n)
				    values[3*c+n] = (*src[c])[n*sn+i];
	    }, sets))
	{
		d->coeff_mode = sets.mode;
		d->index = d->AllocIndex(sets);
		d->coeff = d->Alloc<float>(sets.table.size(), sets.table.data());
		std::cout << "GPU_Backend_HIP: " << sets.count << " distinct coefficient sets, compressed update coefficients ("
		          << (sets.mode==1 ? 16 : 32) << " bit index)" << std::endl;
		return true;
	}

	d->coeff_mode = 0;
	d->vv = d->Alloc<float>(count, vv.data());
	d->vi = d->Alloc<float>(count, vi.data());
	d->ii = d->Alloc<float>(count, ii.data());
	d->iv = d->Alloc<float>(count, iv.data());
	return true;
}

void GPU_Backend_HIP::UpdateVoltages()
{
	const bool c = d->coeff_mode>0;
	const HIP_GridDim& B = d->main_start;
	const HIP_GridDim& E = d->main_stop;
	HIP_Launch(d, "update_voltages", update_voltages, E.nz-B.nz, E.ny-B.ny, E.nx-B.nx,
	            d->volt, (const float*)d->curr, (const void*)d->index, (const float*)(c ? d->coeff : d->vv), (const float*)(c ? d->coeff : d->vi), d->coeff_mode, d->dim, B, E);
}

void GPU_Backend_HIP::UpdateCurrents()
{
	const bool c = d->coeff_mode>0;
	const HIP_GridDim& B = d->main_start;
	const HIP_GridDim& E = d->main_stop;
	HIP_Launch(d, "update_currents", update_currents, E.nz-B.nz, E.ny-B.ny, E.nx-B.nx,
	            d->curr, (const float*)d->volt, (const void*)d->index, (const float*)(c ? d->coeff : d->ii), (const float*)(c ? d->coeff : d->iv), d->coeff_mode, d->dim, B, E);
}

static void CheckSize(size_t host_count, size_t device_count)
{
	if (host_count!=device_count)
		throw std::runtime_error("GPU_Backend_HIP: host and device field size mismatch");
}

void GPU_Backend_HIP::DownloadVoltages(ArrayLib::ArrayNIJK<FDTD_FLOAT>& volt)
{
	CheckSize(volt.size(), 3*d->numCells);
	d->PinHostMemory(volt.data(), volt.size()*sizeof(FDTD_FLOAT));
	HIP_Check(hipMemcpyAsync(volt.data(), d->volt, volt.size()*sizeof(float), hipMemcpyDeviceToHost, d->Stream()), "download");
	d->Flush();
}

void GPU_Backend_HIP::DownloadCurrents(ArrayLib::ArrayNIJK<FDTD_FLOAT>& curr)
{
	CheckSize(curr.size(), 3*d->numCells);
	d->PinHostMemory(curr.data(), curr.size()*sizeof(FDTD_FLOAT));
	HIP_Check(hipMemcpyAsync(curr.data(), d->curr, curr.size()*sizeof(float), hipMemcpyDeviceToHost, d->Stream()), "download");
	d->Flush();
}

void GPU_Backend_HIP::UploadVoltages(const ArrayLib::ArrayNIJK<FDTD_FLOAT>& volt)
{
	CheckSize(volt.size(), 3*d->numCells);
	d->PinHostMemory(volt.data(), volt.size()*sizeof(FDTD_FLOAT));
	HIP_Check(hipMemcpyAsync(d->volt, volt.data(), volt.size()*sizeof(float), hipMemcpyHostToDevice, d->Stream()), "upload");
	d->Flush();
}

void GPU_Backend_HIP::UploadCurrents(const ArrayLib::ArrayNIJK<FDTD_FLOAT>& curr)
{
	CheckSize(curr.size(), 3*d->numCells);
	d->PinHostMemory(curr.data(), curr.size()*sizeof(FDTD_FLOAT));
	HIP_Check(hipMemcpyAsync(d->curr, curr.data(), curr.size()*sizeof(float), hipMemcpyHostToDevice, d->Stream()), "upload");
	d->Flush();
}

bool GPU_Backend_HIP::DownloadRange(bool currents, size_t offset, size_t count, FDTD_FLOAT* dst)
{
	if (offset+count > 3*d->numCells)
		return false;
	const float* src = (currents ? d->curr : d->volt) + offset;
	HIP_Check(hipMemcpyAsync(dst, src, count*sizeof(float), hipMemcpyDeviceToHost, d->Stream()), "range download");
	d->Flush();
	return true;
}

void GPU_Backend_HIP::Synchronize()
{
	d->Flush();
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

// all HIP extensions, see hip_internal.h
static const HIP_ExtensionFactory HIP_EXTENSIONS[] = {
	HIP_CreateExt_Excitation,
	HIP_CreateExt_UPML,
	HIP_CreateExt_Mur_ABC,
	HIP_CreateExt_LorentzMaterial,
	HIP_CreateExt_LumpedRLC,
	HIP_CreateExt_TFSF,
	HIP_CreateExt_Absorbing_BC,
	HIP_CreateExt_SteadyState,
	HIP_CreateExt_Cylinder,
};

GPU_Extension* GPU_Backend_HIP::CreateExtension(Engine_Extension* eng_ext, Engine* eng)
{
	for (size_t n=0; n<sizeof(HIP_EXTENSIONS)/sizeof(HIP_EXTENSIONS[0]); ++n)
	{
		if (HIP_EXTENSIONS[n]==NULL)
			continue;
		GPU_Extension* gpu_ext = HIP_EXTENSIONS[n](d, eng_ext, eng);
		if (gpu_ext)
			return gpu_ext;
	}
	return NULL;
}
