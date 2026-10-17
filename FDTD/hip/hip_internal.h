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

// HIP only: shared by the HIP backend and the HIP extensions

#ifndef HIP_INTERNAL_H
#define HIP_INTERNAL_H

#include <hip/hip_runtime.h>

#include <algorithm>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "gpu_backend_hip.h"
#include "FDTD/gpu_coeff_sets.h"

//! Grid size (number of mesh lines) as passed to the kernels
struct HIP_GridDim
{
	unsigned int nx, ny, nz;
};

//! The kernels index the fields with size_t. The fused step (update_fused) is also built with
//! 32 bit indices, chosen in Init() (see GPU_Backend_HIP::Impl::index32); the helpers it uses take
//! the index type I as a template parameter, size_t by default.
template<typename T> struct HIP_NoDeduce { typedef T type; };

//! Fields and coefficients use the ArrayNIJK layout of the host:
//! index(n,x,y,z) = n*nx*ny*nz + x*ny*nz + y*nz + z
template<typename I = size_t>
__device__ __forceinline__ I nijk(const HIP_GridDim& N, unsigned int n, unsigned int x, unsigned int y, unsigned int z)
{
	return (((I)n*N.nx + x)*N.ny + y)*N.nz + z;
}

//! Main update coefficients of node i, see GPU_Backend_HIP::Impl::coeff_mode:
//! a[n*s] and b[n*s] are vv/vi (voltages, set_offset 0) or ii/iv (currents, set_offset 6)
//! of direction n. Mode 0: ca/cb are the full arrays, else ca is the table of sets.
//! The raw values of a dump on the device, see Engine_GatherRegion: value e of component n is
//! the edge (n, lines[n][0][a], lines[n][1][b], lines[n][2][c]) with e = start[n] + (a*size[n][1] + b)*size[n][2] + c
struct HIP_GatherRegion
{
	const unsigned int* lines[3][3];
	unsigned int size[3][3];
	size_t start[3];
	size_t total;
};

__device__ __forceinline__ size_t region_field_index(const HIP_GatherRegion& R, const HIP_GridDim& N, size_t e)
{
	const int n = (e>=R.start[2]) ? 2 : ((e>=R.start[1]) ? 1 : 0);
	size_t l = e - R.start[n];
	const unsigned int c = l % R.size[n][2];
	l /= R.size[n][2];
	const unsigned int b = l % R.size[n][1];
	const unsigned int a = l / R.size[n][1];
	return nijk(N, n, R.lines[n][0][a], R.lines[n][1][b], R.lines[n][2][c]);
}

template<typename I> struct HIP_MainCoeffT { const float* a; const float* b; I s; };
typedef HIP_MainCoeffT<size_t> HIP_MainCoeff;

template<typename I = size_t>
__device__ __forceinline__ HIP_MainCoeffT<I> main_coeff(const void* index, const float* ca, const float* cb, unsigned int mode,
                                                        typename HIP_NoDeduce<I>::type sn, typename HIP_NoDeduce<I>::type i,
                                                        unsigned int set_offset)
{
	HIP_MainCoeffT<I> c;
	if (mode==0)
	{
		c.a = ca + i; c.b = cb + i; c.s = sn;
		return c;
	}
	const unsigned int set = (mode==1) ? (unsigned int)((const unsigned short*)index)[i] : ((const unsigned int*)index)[i];
	c.a = ca + (size_t)12*set + set_offset; c.b = c.a + 3; c.s = 1;   // 64 bit: more than 358M sets wrap in 32
	return c;
}

#include "hip_upml.hiph"

//! Throw a std::runtime_error if \a err is an error
void HIP_Check(hipError_t err, const char* what);

//! Device and work stream, shared by the backends of all grids of a simulation (see NewSubGridBackend())
struct HIP_Context
{
	int device;
	hipStream_t stream;
	hipStream_t copy_stream;   //!< downloads of field snapshots, overlapping the work (see GPU_Backend_HIP::SnapshotFields())
	std::string name;
	unsigned int warp_size;    //!< threads that run in lockstep: 32 on NVIDIA, 64 on the CDNA GPUs of AMD

	HIP_Context();
	~HIP_Context();
};

//! State of one grid
struct GPU_Backend_HIP::Impl
{
	std::shared_ptr<HIP_Context> ctx;

	HIP_GridDim dim;
	size_t numCells;    //!< nx*ny*nz, the field buffers hold 3*numCells values

	float *volt, *curr;
	float2* energy;                      //!< per-line energy sums, see GPU_Backend_HIP::CalcFastEnergy()
	size_t energy_count;
	std::vector<float2> energy_host;
	//! Main update coefficients: 0: full arrays vv, vi, ii, iv; 1/2: a 16/32 bit set index per node and the sets (see main_coeff())
	unsigned int coeff_mode;
	float *vv, *vi, *ii, *iv;
	void* index;
	float* coeff;

	//! The main updates cover the nodes in [main_start, main_stop), the fused UPML kernels the others (see hip_ext_upml.cu)
	HIP_GridDim main_start, main_stop;
	//! UPML extensions of this grid and whether they run fused with the main updates (-1: not decided yet), see hip_ext_upml.cu
	std::vector<GPU_Extension*> upml;
	int upml_fused;
	//! UPML regions along z updated by the main kernels (nodes below main_start.nz / above main_stop.nz)
	HIP_ZSlabs zslabs;
	bool zslab_lo, zslab_hi;

	//! Fused step (see update_fused): voltages and currents of the main nodes in one pass, from the
	//! current fields (volt, curr) into the next buffers (volt_next, curr_next), swapped afterwards.
	int fused_step;                          //!< -1: not decided yet, 0: off, 1: on
	int fused_blockers;                      //!< extensions or grids that do not allow it
	float *volt_next, *curr_next;
	std::vector<size_t> volt_modified;       //!< flat indices of the voltages the extensions register as changed between the half-steps (excitation, lumped RLC)
	HIP_FusedRegions fregions;              //!< the UPML regions in the kernel (see HIP_Ext_UPML::CanFuse())
	size_t* fixup;                           //!< main nodes whose currents are recomputed after the voltage extensions
	unsigned int fixup_count;
	bool index32;                            //!< the fused step indexes with 32 bit (see Init())
	//! Decide once whether the fused step is used, and prepare it
	bool DecideFusedStep();

	//! Field snapshots (see GPU_Backend_HIP::SnapshotFields()): the raw values the dumps read, copied on
	//! the device (gather_region) and downloaded on the copy stream into a page-locked host buffer
	//! while the next timesteps run
	std::vector<HIP_GatherRegion> snap_regions[2];   //!< of the voltages, of the currents
	std::vector<void*> snap_lines;    //!< device line tables of the regions
	unsigned int snap_nv, snap_n;     //!< number of raw voltages, of all raw values
	float* snap_dev;                  //!< the raw values
	float* snap_host;                 //!< the downloaded raw values
	hipEvent_t snap_evaluated;       //!< snap_dev is written
	hipEvent_t snap_done;            //!< the download is done
	bool snap_taken;                  //!< a download from snap_dev was started
	void FreeSnapshots();

	//! Sums of weighted raw values (frequency domain dumps, see GPU_Backend_HIP::AddFieldDFT())
	struct FieldDFT
	{
		HIP_GatherRegion region;
		unsigned int count;   //!< number of sums per raw value
		bool currents;        //!< raw values of the currents, else of the voltages
		float2* sums;         //!< sum k of raw value e at [k*region.total + e]
	};
	//! Device copy of the line tables of \a region (one buffer, see TryAlloc()); NULL if out of memory
	void* UploadRegion(const Engine_GatherRegion& region, HIP_GatherRegion& dev);
	std::vector<FieldDFT> field_dfts;

	Impl();
	~Impl();

	hipStream_t Stream() const {return ctx->stream;}

	//! Device buffer of \a count elements, initialized with \a host or zero; freed with this Impl
	template <typename T>
	T* Alloc(size_t count, const T* host=NULL)
	{
		T* ptr = NULL;
		HIP_Check(hipMalloc(&ptr, std::max(count, (size_t)1)*sizeof(T)), "hipMalloc");
		m_Allocations.push_back(ptr);
		// on the work stream: it does not synchronize with the legacy default stream
		if (host && count)
			HIP_Check(hipMemcpyAsync(ptr, host, count*sizeof(T), hipMemcpyHostToDevice, Stream()), "hipMemcpy");
		else
			HIP_Check(hipMemsetAsync(ptr, 0, std::max(count, (size_t)1)*sizeof(T), Stream()), "hipMemset");
		Flush();
		return ptr;
	}

	//! Uninitialized device buffer of \a count elements for an optional optimization, freed with
	//! this Impl or by Free(); NULL if the device memory is not sufficient
	template <typename T>
	T* TryAlloc(size_t count)
	{
		T* ptr = NULL;
		if (hipMalloc(&ptr, std::max(count, (size_t)1)*sizeof(T))!=hipSuccess)
		{
			(void)hipGetLastError();   // not sticky, but reported by the next CheckLaunch()
			return NULL;
		}
		m_Allocations.push_back(ptr);
		return ptr;
	}

	//! Free a buffer of Alloc() or TryAlloc() before this Impl
	void Free(void* ptr)
	{
		if (ptr==NULL)
			return;
		m_Allocations.erase(std::remove(m_Allocations.begin(), m_Allocations.end(), ptr), m_Allocations.end());
		hipFree(ptr);
	}

	//! Device buffer of the set indices of \a sets (16 or 32 bit)
	void* AllocIndex(const GPU_CoeffSets& sets)
	{
		if (sets.mode==1)
		{
			std::vector<unsigned short> index16(sets.index.begin(), sets.index.end());
			return Alloc<unsigned short>(index16.size(), index16.data());
		}
		return Alloc<unsigned int>(sets.index.size(), sets.index.data());
	}

	//! No UPML regions updated by the main kernels (call after changing main_start/main_stop)
	void ResetZSlabs();
	//! Nodes [B, E) of the main kernels: main_start/main_stop and the UPML regions along z they update
	void MainRange(HIP_GridDim& B, HIP_GridDim& E, unsigned int& threads_z);

	//! Wait until all work on the stream is done
	void Flush();

	//! Check for a kernel launch error
	void CheckLaunch(const char* kernel);

	//! Page-lock host memory (once) for fast transfers
	void PinHostMemory(void* ptr, size_t bytes);

	//! Launch geometry: one thread per (i,j,k), i fastest
	static dim3 Block(size_t ni, size_t nj);
	static dim3 Grid(dim3 block, size_t ni, size_t nj, size_t nk);

protected:
	std::vector<void*> m_Allocations;
	std::set<void*> m_Pinned;
};

//! Kernel timing for development (environment variable OPENEMS_HIP_KERNEL_TIMES=1): the
//! time of every launch, per kernel name, printed at exit. Waits for each kernel.
struct HIP_KernelTimes
{
	static bool Enabled();
	static void Begin(hipStream_t stream);
	static void End(hipStream_t stream, const char* name);
};

//! Launch \a kernel with one thread per (i,j,k) on the stream of \a d; the kernel checks its bounds
template <typename Kernel, typename... Args>
void HIP_Launch(GPU_Backend_HIP::Impl* d, const char* name, Kernel kernel, size_t ni, size_t nj, size_t nk, Args... args)
{
	if (ni==0 || nj==0 || nk==0)
		return;
	dim3 block = GPU_Backend_HIP::Impl::Block(ni, nj);
	const bool timed = HIP_KernelTimes::Enabled();
	if (timed)
		HIP_KernelTimes::Begin(d->Stream());
	kernel<<<GPU_Backend_HIP::Impl::Grid(block, ni, nj, nk), block, 0, d->Stream()>>>(args...);
	d->CheckLaunch(name);
	if (timed)
		HIP_KernelTimes::End(d->Stream(), name);
}

//! Factory of a HIP extension: the device implementation of \a eng_ext, or NULL if \a eng_ext is not of its type
typedef GPU_Extension* (*HIP_ExtensionFactory)(GPU_Backend_HIP::Impl* d, Engine_Extension* eng_ext, Engine* eng);
GPU_Extension* HIP_CreateExt_Excitation(GPU_Backend_HIP::Impl* d, Engine_Extension* eng_ext, Engine* eng);
GPU_Extension* HIP_CreateExt_UPML(GPU_Backend_HIP::Impl* d, Engine_Extension* eng_ext, Engine* eng);
GPU_Extension* HIP_CreateExt_Mur_ABC(GPU_Backend_HIP::Impl* d, Engine_Extension* eng_ext, Engine* eng);
GPU_Extension* HIP_CreateExt_LorentzMaterial(GPU_Backend_HIP::Impl* d, Engine_Extension* eng_ext, Engine* eng);
GPU_Extension* HIP_CreateExt_DebyeMaterial(GPU_Backend_HIP::Impl* d, Engine_Extension* eng_ext, Engine* eng);
GPU_Extension* HIP_CreateExt_LumpedRLC(GPU_Backend_HIP::Impl* d, Engine_Extension* eng_ext, Engine* eng);
GPU_Extension* HIP_CreateExt_TFSF(GPU_Backend_HIP::Impl* d, Engine_Extension* eng_ext, Engine* eng);
GPU_Extension* HIP_CreateExt_Absorbing_BC(GPU_Backend_HIP::Impl* d, Engine_Extension* eng_ext, Engine* eng);
GPU_Extension* HIP_CreateExt_SteadyState(GPU_Backend_HIP::Impl* d, Engine_Extension* eng_ext, Engine* eng);
GPU_Extension* HIP_CreateExt_Cylinder(GPU_Backend_HIP::Impl* d, Engine_Extension* eng_ext, Engine* eng);

#endif // HIP_INTERNAL_H
