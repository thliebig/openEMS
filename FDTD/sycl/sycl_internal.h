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

// SYCL only: shared by the source files of the SYCL backend

#ifndef SYCL_INTERNAL_H
#define SYCL_INTERNAL_H

#include <sycl/sycl.hpp>
#include "FDTD/kernels/cell_update.h"

#include <algorithm>
#include <memory>
#include <mutex>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "gpu_backend_sycl.h"

//! Grid size (number of mesh lines) as passed to the kernels
struct SYCL_GridDim
{
	unsigned int nx, ny, nz;
};

//! The fields use the ArrayNIJK layout of the host:
//! index(n,x,y,z) = n*nx*ny*nz + x*ny*nz + y*nz + z
template<typename I = size_t>
inline I nijk(const SYCL_GridDim& N, unsigned int n, unsigned int x, unsigned int y, unsigned int z)
{
	return (((I)n*N.nx + x)*N.ny + y)*N.nz + z;
}

//! A complex float in the layout of std::complex<float>
struct SYCL_Float2
{
	float x, y;
};

//! The raw values of a dump on the device, see Engine_GatherRegion: value e of component n is
//! the edge (n, lines[n][0][a], lines[n][1][b], lines[n][2][c]) with e = start[n] + (a*size[n][1] + b)*size[n][2] + c
struct SYCL_GatherRegion
{
	const unsigned int* lines[3][3];
	unsigned int size[3][3];
	size_t start[3];
	size_t total;
};

inline size_t region_field_index(const SYCL_GatherRegion& R, const SYCL_GridDim& N, size_t e)
{
	const int n = (e>=R.start[2]) ? 2 : ((e>=R.start[1]) ? 1 : 0);
	size_t l = e - R.start[n];
	const unsigned int c = l % R.size[n][2];
	l /= R.size[n][2];
	const unsigned int b = l % R.size[n][1];
	const unsigned int a = l / R.size[n][1];
	return nijk(N, n, R.lines[n][0][a], R.lines[n][1][b], R.lines[n][2][c]);
}

//! Run \a f; a sycl::exception becomes a std::runtime_error naming \a what
template <typename F>
auto SYCL_Check(const char* what, F f) -> decltype(f())
{
	try
	{
		return f();
	}
	catch (const sycl::exception& e)
	{
		throw std::runtime_error(std::string("GPU_Backend_SYCL: ") + what + " failed: " + e.what());
	}
}

//! Device and queues, shared by the backends of all grids of a simulation (see NewSubGridBackend())
struct SYCL_Context
{
	sycl::device device;
	sycl::context context;
	sycl::queue queue;        //!< all kernels and the transfers of the work, in order
	sycl::queue copy_queue;   //!< downloads of field snapshots, overlapping the work (see GPU_Backend_SYCL::SnapshotFields())
	std::string name;
	size_t cache_bytes = 0;   //!< the last-level cache, 0 unknown (see GPU_Backend_SYCL::CacheBytes())
	size_t max_alloc = 0;     //!< the largest allocation the device allows
	size_t max_group = 0;     //!< the most work-items of a work-group

	//! Both queues in order, on one context of \a dev; the asynchronous errors are kept (see ThrowAsyncError())
	explicit SYCL_Context(const sycl::device& dev);
	~SYCL_Context();

	//! Throw a std::runtime_error if the device reported an asynchronous error
	void ThrowAsyncError();

protected:
	std::mutex m_ErrorMutex;
	std::string m_AsyncError;
};

//! State of one grid
struct GPU_Backend_SYCL::Impl
{
	std::shared_ptr<SYCL_Context> ctx;

	SYCL_GridDim dim;
	size_t numCells = 0;    //!< nx*ny*nz, the field buffers hold 3*numCells values

	float* volt = NULL;
	float* curr = NULL;
	SYCL_Float2* energy = NULL;          //!< per-line energy sums, see GPU_Backend_SYCL::CalcFastEnergy()
	size_t energy_count = 0;
	std::vector<SYCL_Float2> energy_host;
	//! The main updates on the cell model (see Init(), sycl_cell_model.cpp)
	bool cm_index16 = true;
	void* cm_index[2] = {NULL, NULL};
	CM_Field cm_field[2];
	size_t cm_state_size[2] = {0, 0};   //!< floats of cm_field[f].state
	//! The fused step (see GPU_Backend_SYCL::PrepareFused(), sycl_cell_model.cpp): work-groups of TZ x TY items, through XC lines along x
	struct Fused
	{
		unsigned int TZ = 0, TY = 0, XC = 0;
		float* volt = NULL;      //!< the second field buffers, NULL without the fused step
		float* curr = NULL;
		float* vstate = NULL;    //!< and of the state of the voltages
		uint32_t* late = NULL;   //!< bitmap of the late currents (GPU_Backend::PrepareFused()), NULL none
		uint64_t* list = NULL;   //!< and their list
		size_t count = 0;
	} fused;

	//! Field snapshots (see GPU_Backend_SYCL::SnapshotFields()): the raw values the dumps read, copied on
	//! the device and downloaded on the copy queue into host memory of the runtime while the next timesteps run
	std::vector<SYCL_GatherRegion> snap_regions[2];   //!< of the voltages, of the currents
	std::vector<void*> snap_lines;    //!< device line tables of the regions
	size_t snap_nv = 0, snap_n = 0;   //!< number of raw voltages, of all raw values
	float* snap_dev = NULL;           //!< the raw values
	float* snap_host = NULL;          //!< the downloaded raw values
	sycl::event snap_done;            //!< the download
	bool snap_taken = false;          //!< a download from snap_dev was started
	void FreeSnapshots();

	//! Sums of weighted raw values (frequency domain dumps, see GPU_Backend_SYCL::AddFieldDFT())
	struct FieldDFT
	{
		SYCL_GatherRegion region;
		unsigned int count;    //!< number of sums per raw value
		bool currents;         //!< raw values of the currents, else of the voltages
		SYCL_Float2* sums;     //!< sum k of raw value e at [k*region.total + e]
	};
	//! Device copy of the line tables of \a region (one buffer, see TryAlloc()); NULL if out of memory
	void* UploadRegion(const Engine_GatherRegion& region, SYCL_GatherRegion& dev);
	std::vector<FieldDFT> field_dfts;

	Impl();
	~Impl();

	sycl::queue& Queue() const {return ctx->queue;}

	//! Device buffer of \a count elements, initialized with \a host or zero; freed with this Impl
	template <typename T>
	T* Alloc(size_t count, const T* host=NULL)
	{
		T* ptr = TryAlloc<T>(count);
		if (ptr==NULL)
			throw std::runtime_error("GPU_Backend_SYCL: out of device memory (" + std::to_string(count*sizeof(T)/1048576) + " MiB)");
		const size_t bytes = std::max(count, (size_t)1)*sizeof(T);
		SYCL_Check("upload", [&] {
			if (host && count)
				Queue().memcpy(ptr, host, count*sizeof(T));
			else
				Queue().memset(ptr, 0, bytes);
		});
		Flush();   // the host data may be temporary
		return ptr;
	}

	//! Uninitialized device buffer of \a count elements, freed with this Impl or by Free(); NULL if the
	//! device memory is not sufficient
	template <typename T>
	T* TryAlloc(size_t count)
	{
		const size_t bytes = std::max(count, (size_t)1)*sizeof(T);
		if (bytes>ctx->max_alloc)
			return NULL;
		T* ptr = NULL;
		try
		{
			ptr = sycl::malloc_device<T>(std::max(count, (size_t)1), ctx->device, ctx->context);
		}
		catch (const sycl::exception&)
		{
			ptr = NULL;
		}
		if (ptr)
			m_Allocations.push_back(ptr);
		return ptr;
	}

	//! Free a buffer of Alloc() or TryAlloc() before this Impl, once the work on it is done
	void Free(void* ptr);

	//! Wait until all work on the queue is done, throw its errors
	void Flush();

	//! Make the host memory \a ptr faster to copy to and from (once), if the runtime can
	void PinHostMemory(void* ptr, size_t bytes);

protected:
	std::vector<void*> m_Allocations;
	std::set<void*> m_Pinned;
};

//! Launch \a f(i, j, k) for all i<ni, j<nj, k<nk on the queue of \a d, after \a deps; i is the fastest index
/*!
  Work-groups as the blocks of the HIP backend: 256 items along i, 32 x 8 or 16 x 16 for two dimensions.
  */
template <typename F>
sycl::event SYCL_Launch(GPU_Backend_SYCL::Impl* d, const char* name, size_t ni, size_t nj, size_t nk, F f,
                        const std::vector<sycl::event>& deps = std::vector<sycl::event>())
{
	if (ni==0 || nj==0 || nk==0)
		return sycl::event();
	size_t bi = 256, bj = 1;
	if ((nj>1) && (ni<=16))
		bi = bj = 16;
	else if (nj>1)
	{
		bi = 32;
		bj = 8;
	}
	const sycl::range<3> local(1, bj, bi);
	const sycl::range<3> global(nk, (nj+bj-1)/bj*bj, (ni+bi-1)/bi*bi);
	return SYCL_Check(name, [&] {
		return d->Queue().submit([&](sycl::handler& h) {
			h.depends_on(deps);
			h.parallel_for(sycl::nd_range<3>(global, local), [=](sycl::nd_item<3> it) {
				const size_t i = it.get_global_id(2), j = it.get_global_id(1);
				if ((i<ni) && (j<nj))
					f(i, j, it.get_global_id(0));
			});
		});
	});
}

#endif // SYCL_INTERNAL_H
