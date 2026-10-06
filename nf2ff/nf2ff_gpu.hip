/*
*	Copyright (C) 2026 Hugodonotexit
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

#include "nf2ff_gpu.h"

#include <hip/hip_runtime.h>

#include <algorithm>
#include <cstring>
#include <mutex>
#include <sstream>
#include <vector>

namespace
{

// complex floats are float2 (x real, y imaginary), the layout of std::complex<float>

// ---- device code ----------------------------------------------------------

//! Component ax of the unit vector of angle a (theta index a/n_ph, phi index a%n_ph)
__device__ __forceinline__ float direction(unsigned int a, unsigned int n_ph, const float* theta, const float* phi, int ax)
{
	float sT, cT;
	sincosf(theta[a / n_ph], &sT, &cT);
	if (ax == 2)
		return cT;
	float sP, cP;
	sincosf(phi[a % n_ph], &sP, &cP);
	return ax == 0 ? sT * cP : sT * sP;
}

//! P[a][i] = exp(j k u_ax(a0+a) crd[i]) for one tile of angles
__global__ void phase_kernel(float2* P, const float* crd, unsigned int n_cols, unsigned int n_tile,
                             unsigned int a0, unsigned int n_ph, const float* theta, const float* phi,
                             float k, int ax)
{
	size_t total = (size_t)n_tile * n_cols;
	for (size_t idx = blockIdx.x * (size_t)blockDim.x + threadIdx.x; idx < total; idx += (size_t)gridDim.x * blockDim.x)
	{
		unsigned int a = (unsigned int)(idx / n_cols);
		unsigned int i = (unsigned int)(idx % n_cols);
		float s, c;
		sincosf(k * direction(a0 + a, n_ph, theta, phi, ax) * crd[i], &s, &c);
		P[idx] = make_float2(c, s);
	}
}

// tile of the matrix product: BM x BN outputs per block, BM/TM x BN/TN = 256 threads
#define GEMM_BM 64
#define GEMM_BN 64
#define GEMM_BK 8
#define GEMM_TM 4
#define GEMM_TN 4
#define GEMM_THREADS 256

//! C(M x N) = A(M x K) B(K x N), complex, all row-major
__global__ void __launch_bounds__(GEMM_THREADS) gemm_kernel(const float2* A, const float2* B, float2* C,
                                                            unsigned int M, unsigned int N, unsigned int K)
{
	__shared__ float2 As[GEMM_BK][GEMM_BM];
	__shared__ float2 Bs[GEMM_BK][GEMM_BN];

	const unsigned int tid = threadIdx.x;
	const unsigned int tx = tid % (GEMM_BN / GEMM_TN);
	const unsigned int ty = tid / (GEMM_BN / GEMM_TN);
	const size_t row0 = (size_t)blockIdx.x * GEMM_BM;
	const size_t col0 = (size_t)blockIdx.y * GEMM_BN;

	float2 acc[GEMM_TM][GEMM_TN];
	for (int i = 0; i < GEMM_TM; ++i)
		for (int j = 0; j < GEMM_TN; ++j)
			acc[i][j] = make_float2(0.f, 0.f);

	for (unsigned int k0 = 0; k0 < K; k0 += GEMM_BK)
	{
		for (unsigned int e = tid; e < GEMM_BM * GEMM_BK; e += GEMM_THREADS)
		{
			unsigned int r = e / GEMM_BK, c = e % GEMM_BK;
			As[c][r] = (row0 + r < M && k0 + c < K) ? A[(row0 + r) * K + k0 + c] : make_float2(0.f, 0.f);
		}
		for (unsigned int e = tid; e < GEMM_BK * GEMM_BN; e += GEMM_THREADS)
		{
			unsigned int r = e / GEMM_BN, c = e % GEMM_BN;
			Bs[r][c] = (k0 + r < K && col0 + c < N) ? B[(size_t)(k0 + r) * N + col0 + c] : make_float2(0.f, 0.f);
		}
		__syncthreads();

		for (int kk = 0; kk < GEMM_BK; ++kk)
		{
			float2 a[GEMM_TM], b[GEMM_TN];
			for (int i = 0; i < GEMM_TM; ++i)
				a[i] = As[kk][ty * GEMM_TM + i];
			for (int j = 0; j < GEMM_TN; ++j)
				b[j] = Bs[kk][tx * GEMM_TN + j];
			for (int i = 0; i < GEMM_TM; ++i)
				for (int j = 0; j < GEMM_TN; ++j)
				{
					acc[i][j].x += a[i].x * b[j].x - a[i].y * b[j].y;
					acc[i][j].y += a[i].x * b[j].y + a[i].y * b[j].x;
				}
		}
		__syncthreads();
	}

	for (int i = 0; i < GEMM_TM; ++i)
		for (int j = 0; j < GEMM_TN; ++j)
		{
			size_t r = row0 + ty * GEMM_TM + i, c = col0 + tx * GEMM_TN + j;
			if (r < M && c < N)
				C[r * N + c] = acc[i][j];
		}
}

#define REDUCE_THREADS 128

//! One block per angle: S_c = sum_q T[c][q] Pq[q], then Nt, Np, Lt, Lp
__global__ void __launch_bounds__(REDUCE_THREADS) reduce_kernel(const float2* T, const float2* Pq, unsigned int n_q,
                                                                unsigned int n_tile, unsigned int a0, unsigned int n_ph,
                                                                const float* theta, const float* phi,
                                                                float k, int ax_n, float const_off, float2* out)
{
	// shared memory instead of warp shuffles: the warp size is 32 or 64 depending on the GPU
	__shared__ float2 sh[4][REDUCE_THREADS];

	const unsigned int a = blockIdx.x;
	const float2* Ta = T + (size_t)a * 4 * n_q;
	const float2* Pa = Pq + (size_t)a * n_q;

	float2 acc[4];
	for (int c = 0; c < 4; ++c)
		acc[c] = make_float2(0.f, 0.f);
	for (unsigned int q = threadIdx.x; q < n_q; q += REDUCE_THREADS)
	{
		float2 p = Pa[q];
		for (int c = 0; c < 4; ++c)
		{
			float2 t = Ta[(size_t)c * n_q + q];
			acc[c].x += t.x * p.x - t.y * p.y;
			acc[c].y += t.x * p.y + t.y * p.x;
		}
	}
	for (int c = 0; c < 4; ++c)
		sh[c][threadIdx.x] = acc[c];
	__syncthreads();

	for (unsigned int stride = REDUCE_THREADS / 2; stride > 0; stride >>= 1)
	{
		if (threadIdx.x < stride)
			for (int c = 0; c < 4; ++c)
			{
				sh[c][threadIdx.x].x += sh[c][threadIdx.x + stride].x;
				sh[c][threadIdx.x].y += sh[c][threadIdx.x + stride].y;
			}
		__syncthreads();
	}

	if (threadIdx.x != 0)
		return;

	float2 S[4];
	for (int c = 0; c < 4; ++c)
		S[c] = sh[c][0];

	// the plane position along its normal is one phase factor per angle
	const unsigned int ang = a0 + a;
	float sT, cT, sP, cP;
	sincosf(theta[ang / n_ph], &sT, &cT);
	sincosf(phi[ang % n_ph], &sP, &cP);
	float u[3] = {sT * cP, sT * sP, cT};
	float s, c;
	sincosf(k * u[ax_n] * const_off, &s, &c);
	for (int i = 0; i < 4; ++i)
		S[i] = make_float2(S[i].x * c - S[i].y * s, S[i].x * s + S[i].y * c);

	// the currents have no component along the normal
	float2 N[3] = {make_float2(0, 0), make_float2(0, 0), make_float2(0, 0)};
	float2 L[3] = {make_float2(0, 0), make_float2(0, 0), make_float2(0, 0)};
	int b = (ax_n + 1) % 3, d = (ax_n + 2) % 3;
	N[b] = S[0]; N[d] = S[1];
	L[b] = S[2]; L[d] = S[3];

	float2 Nt, Np, Lt, Lp;
	Nt.x = (N[0].x * cP + N[1].x * sP) * cT - N[2].x * sT;
	Nt.y = (N[0].y * cP + N[1].y * sP) * cT - N[2].y * sT;
	Np.x = -N[0].x * sP + N[1].x * cP;
	Np.y = -N[0].y * sP + N[1].y * cP;
	Lt.x = (L[0].x * cP + L[1].x * sP) * cT - L[2].x * sT;
	Lt.y = (L[0].y * cP + L[1].y * sP) * cT - L[2].y * sT;
	Lp.x = -L[0].x * sP + L[1].x * cP;
	Lp.y = -L[0].y * sP + L[1].y * cP;

	out[0 * (size_t)n_tile + a] = Nt;
	out[1 * (size_t)n_tile + a] = Np;
	out[2 * (size_t)n_tile + a] = Lt;
	out[3 * (size_t)n_tile + a] = Lp;
}

// ---- host code ------------------------------------------------------------

//! A device buffer that only ever grows
struct Buffer
{
	void* ptr = nullptr;
	size_t size = 0;

	bool Reserve(size_t bytes)
	{
		if (bytes <= size)
			return true;
		(void)hipFree(ptr);
		ptr = nullptr;
		size = 0;
		if (hipMalloc(&ptr, bytes) != hipSuccess)
		{
			(void)hipGetLastError();
			return false;
		}
		size = bytes;
		return true;
	}
	template <typename T> T* as() const { return static_cast<T*>(ptr); }
};

struct Context
{
	std::mutex mutex;
	bool ok = false;
	std::string device;
	std::string error;
	Buffer theta, phi, crd_p, crd_q, G, Pp, Pq, T, out;

	Context()
	{
		int n = 0;
		hipError_t err = hipGetDeviceCount(&n);
		if (err != hipSuccess || n < 1)
		{
			(void)hipGetLastError();
			error = err == hipSuccess ? "no GPU found" : hipGetErrorString(err);
			return;
		}
		err = hipSetDevice(0);
		if (err == hipSuccess)
			err = hipFree(0);   // forces the context creation
		if (err != hipSuccess)
		{
			(void)hipGetLastError();
			error = std::string("GPU initialisation failed: ") + hipGetErrorString(err);
			return;
		}
		hipDeviceProp_t prop;
		if (hipGetDeviceProperties(&prop, 0) == hipSuccess)
			device = prop.name;
		ok = true;
	}
};

// Never destroyed: on Windows the GPU runtime is already gone when static
// destructors run, and the OS reclaims the device memory at exit anyway.
Context& ctx()
{
	static Context* c = new Context();
	return *c;
}

bool fail(Context& c, const std::string& what, hipError_t err = hipSuccess)
{
	c.error = what;
	if (err != hipSuccess)
	{
		c.error += std::string(": ") + hipGetErrorString(err);
		(void)hipGetLastError();
	}
	return false;
}

inline unsigned int blocks(size_t n, unsigned int per_block)
{
	return (unsigned int)std::min<size_t>(65535, (n + per_block - 1) / per_block);
}

} // namespace

namespace nf2ff_gpu
{

bool Available()
{
	return ctx().ok;
}

std::string DeviceName()
{
	return ctx().device;
}

std::string LastError()
{
	return ctx().error;
}

bool PlaneSums(float k,
               unsigned int n_th, const float* theta,
               unsigned int n_ph, const float* phi,
               int ax_n, int ax_p, int ax_q,
               unsigned int n_p, const float* crd_p,
               unsigned int n_q, const float* crd_q,
               float const_off,
               const std::complex<float>* G,
               std::complex<float>* out)
{
	Context& c = ctx();
	std::lock_guard<std::mutex> lock(c.mutex);
	if (!c.ok)
		return false;

	const size_t n_ang = (size_t)n_th * n_ph;
	if (n_ang == 0 || n_p == 0 || n_q == 0 || n_ang > 0xffffffffu)
		return fail(c, "empty or too large angle/mesh grid");
	if (n_p > 0x3fffffffu / 4 || n_q > 0x3fffffffu / 4)
		return fail(c, "plane too large for the GPU kernels");

	if (!c.theta.Reserve(sizeof(float) * n_th) || !c.phi.Reserve(sizeof(float) * n_ph)
	    || !c.crd_p.Reserve(sizeof(float) * n_p) || !c.crd_q.Reserve(sizeof(float) * n_q)
	    || !c.G.Reserve(sizeof(float2) * (size_t)n_p * 4 * n_q))
		return fail(c, "out of device memory");

	hipError_t err;
	if ((err = hipMemcpy(c.theta.ptr, theta, sizeof(float) * n_th, hipMemcpyHostToDevice)) != hipSuccess
	    || (err = hipMemcpy(c.phi.ptr, phi, sizeof(float) * n_ph, hipMemcpyHostToDevice)) != hipSuccess
	    || (err = hipMemcpy(c.crd_p.ptr, crd_p, sizeof(float) * n_p, hipMemcpyHostToDevice)) != hipSuccess
	    || (err = hipMemcpy(c.crd_q.ptr, crd_q, sizeof(float) * n_q, hipMemcpyHostToDevice)) != hipSuccess
	    || (err = hipMemcpy(c.G.ptr, G, sizeof(float2) * (size_t)n_p * 4 * n_q, hipMemcpyHostToDevice)) != hipSuccess)
		return fail(c, "upload failed", err);

	// angles per tile: T (4 n_q), Pp (n_p), Pq (n_q) and the result (4) per angle
	size_t free_bytes = 0, total_bytes = 0;
	if (hipMemGetInfo(&free_bytes, &total_bytes) != hipSuccess)
	{
		(void)hipGetLastError();
		free_bytes = (size_t)2 << 30;   // unknown, assume a small card
	}
	size_t budget = std::min<size_t>(free_bytes / 2, (size_t)1 << 30);
	size_t per_angle = sizeof(float2) * ((size_t)4 * n_q + n_p + n_q + 4);
	unsigned int tile = (unsigned int)std::max<size_t>(1, std::min<size_t>(n_ang, budget / per_angle));

	if (!c.Pp.Reserve(sizeof(float2) * (size_t)tile * n_p)
	    || !c.Pq.Reserve(sizeof(float2) * (size_t)tile * n_q)
	    || !c.T.Reserve(sizeof(float2) * (size_t)tile * 4 * n_q)
	    || !c.out.Reserve(sizeof(float2) * (size_t)tile * 4))
		return fail(c, "out of device memory");

	std::vector<std::complex<float> > staging((size_t)tile * 4);

	for (size_t a0 = 0; a0 < n_ang; a0 += tile)
	{
		unsigned int nt = (unsigned int)std::min<size_t>(tile, n_ang - a0);

		phase_kernel<<<blocks((size_t)nt * n_p, 256), 256>>>(
			c.Pp.as<float2>(), c.crd_p.as<float>(), n_p, nt, (unsigned int)a0, n_ph,
			c.theta.as<float>(), c.phi.as<float>(), k, ax_p);
		phase_kernel<<<blocks((size_t)nt * n_q, 256), 256>>>(
			c.Pq.as<float2>(), c.crd_q.as<float>(), n_q, nt, (unsigned int)a0, n_ph,
			c.theta.as<float>(), c.phi.as<float>(), k, ax_q);

		// T(nt x 4n_q) = Pp(nt x n_p) G(n_p x 4n_q)
		dim3 grid((nt + GEMM_BM - 1) / GEMM_BM, (4 * n_q + GEMM_BN - 1) / GEMM_BN);
		gemm_kernel<<<grid, GEMM_THREADS>>>(c.Pp.as<float2>(), c.G.as<float2>(), c.T.as<float2>(),
		                                    nt, 4 * n_q, n_p);

		reduce_kernel<<<nt, REDUCE_THREADS>>>(c.T.as<float2>(), c.Pq.as<float2>(), n_q, nt, (unsigned int)a0, n_ph,
		                                      c.theta.as<float>(), c.phi.as<float>(), k, ax_n, const_off, c.out.as<float2>());

		if ((err = hipGetLastError()) != hipSuccess)
			return fail(c, "kernel launch failed", err);
		if ((err = hipMemcpy(staging.data(), c.out.ptr, sizeof(float2) * (size_t)nt * 4, hipMemcpyDeviceToHost)) != hipSuccess)
			return fail(c, "kernel execution failed", err);

		for (int comp = 0; comp < 4; ++comp)
			std::memcpy(out + (size_t)comp * n_ang + a0, staging.data() + (size_t)comp * nt, sizeof(std::complex<float>) * nt);
	}
	return true;
}

} // namespace nf2ff_gpu
