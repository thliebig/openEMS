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

//! Fields and coefficients use the ArrayNIJK layout of the host:
//! index(n,x,y,z) = n*nx*ny*nz + x*ny*nz + y*nz + z
__device__ __forceinline__ unsigned int nijk(const HIP_GridDim& N, unsigned int n, unsigned int x, unsigned int y, unsigned int z)
{
	return ((n*N.nx + x)*N.ny + y)*N.nz + z;
}

//! Main update coefficients of node i, see GPU_Backend_HIP::Impl::coeff_mode:
//! a[n*s] and b[n*s] are vv/vi (voltages, set_offset 0) or ii/iv (currents, set_offset 6)
//! of direction n. Mode 0: ca/cb are the full arrays, else ca is the table of sets.
struct HIP_MainCoeff { const float* a; const float* b; unsigned int s; };

__device__ __forceinline__ HIP_MainCoeff main_coeff(const void* index, const float* ca, const float* cb,
                                                     unsigned int mode, unsigned int sn, unsigned int i, unsigned int set_offset)
{
	HIP_MainCoeff c;
	if (mode==0)
	{
		c.a = ca + i; c.b = cb + i; c.s = sn;
		return c;
	}
	const unsigned int set = (mode==1) ? (unsigned int)((const unsigned short*)index)[i] : ((const unsigned int*)index)[i];
	c.a = ca + 12*set + set_offset; c.b = c.a + 3; c.s = 1;
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
	std::vector<unsigned int> volt_modified; //!< flat indices of the voltages the extensions register as changed between the half-steps (excitation, lumped RLC)
	HIP_FusedRegions fregions;              //!< the UPML regions in the kernel (see HIP_Ext_UPML::CanFuse())
	unsigned int* fixup;                     //!< main nodes whose currents are recomputed after the voltage extensions
	unsigned int fixup_count;
	//! Decide once whether the fused step is used, and prepare it
	bool DecideFusedStep();

	//! Field snapshots (see GPU_Backend_HIP::SnapshotFields()): the dumped values, evaluated on the device
	//! (gather_dumps) and downloaded on the copy stream into one of two page-locked host slots while the
	//! next timesteps run
	GPU_GatherEntry* snap_entries;   //!< the entries of the voltages, then of the currents
	unsigned int snap_nv, snap_n;     //!< number of voltage entries, of all entries
	float* snap_dev;                  //!< the values of the entries
	float* snap_host[2];              //!< the values of the entries, for each slot
	hipEvent_t snap_evaluated;       //!< snap_dev is written
	hipEvent_t snap_done[2];         //!< the download into each slot is done
	int snap_last;                    //!< slot of the last download from snap_dev, -1: none
	void FreeSnapshots();

	//! Sums of weighted dump entries (frequency domain dumps, see GPU_Backend_HIP::AddFieldDFT())
	struct FieldDFT
	{
		GPU_GatherEntry* entries;
		unsigned int n;       //!< number of entries
		unsigned int count;   //!< number of sums per entry
		bool currents;        //!< entries from the currents, else from the voltages
		float2* sums;         //!< sum k of entry i at [k*n + i]
	};
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
GPU_Extension* HIP_CreateExt_LumpedRLC(GPU_Backend_HIP::Impl* d, Engine_Extension* eng_ext, Engine* eng);
GPU_Extension* HIP_CreateExt_TFSF(GPU_Backend_HIP::Impl* d, Engine_Extension* eng_ext, Engine* eng);
GPU_Extension* HIP_CreateExt_Absorbing_BC(GPU_Backend_HIP::Impl* d, Engine_Extension* eng_ext, Engine* eng);
GPU_Extension* HIP_CreateExt_SteadyState(GPU_Backend_HIP::Impl* d, Engine_Extension* eng_ext, Engine* eng);
GPU_Extension* HIP_CreateExt_Cylinder(GPU_Backend_HIP::Impl* d, Engine_Extension* eng_ext, Engine* eng);

#endif // HIP_INTERNAL_H
