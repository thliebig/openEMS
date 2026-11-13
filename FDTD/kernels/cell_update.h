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

#ifndef FDTD_KERNELS_CELL_UPDATE_H
#define FDTD_KERNELS_CELL_UPDATE_H

// The update of one node from the cell model (CellModel), per node instead of per run: the
// GPU backends run it with a thread per node, the reference backend on the host threads, all
// with the same data and the same operations. The body is the one of the CPU engine's runs
// (engine_cpu_kernels.inc), operation for operation, except the coefficient (CM_FACTOR): the
// CPU engine multiplies the material part and the mesh factors in double once per run, a
// thread here once per node in float, at a fraction of the double rate on a GPU (HIP on an
// RX 7900 XTX 5-11 % faster) and Metal has no double. So the results are within single
// precision of the CPU engine.
//
// Portable subset, so that HIP, Metal, SYCL and the host compile the same text: plain data,
// no std, the address spaces by macro (CM_GLOBAL: the buffers, CM_PRIVATE: the kernel's own
// memory, which Metal needs on references). The Metal backend compiles this header at run time,
// together with kernels/ext_ops.h and metal/metal_kernels.metal (embedded at build time, see
// FDTD/CMakeLists.txt; the build target metal_check compiles them offline).

#if defined(__METAL_VERSION__)
#define CM_FN inline
#define CM_GLOBAL device
#define CM_PRIVATE thread
#define CM_RESTRICT
typedef ulong uint64_t;
typedef long int64_t;
typedef uint uint32_t;
typedef int int32_t;
typedef ushort uint16_t;
typedef uchar uint8_t;
#else
#if defined(__HIPCC__) || defined(__CUDACC__)
#define CM_FN __host__ __device__ __forceinline__
#else
#define CM_FN inline
#endif
#define CM_GLOBAL
#define CM_PRIVATE
#if defined(_MSC_VER) && !defined(__CUDACC__)
#define CM_RESTRICT __restrict
#else
#define CM_RESTRICT __restrict__
#endif
#include <stdint.h>
#endif

//! The material part times the mesh factors of the line of component \a n, in single precision
#define CM_FACTOR(cb, F, n, line) ((cb).m[n]*(F).gxy[n][line])

//! What a UPML or dispersive component of a node reads beyond the header of CM_Combo: the rest of its CellModel_Set, with the profile offset of its run
struct CM_Comp
{
	float uvv, ufo, ufn;      //!< UPML
	float solve;              //!< DISP: Debye solve
	int32_t p0;               //!< UPML_PROFILE: profile index of z = 0 (set.offset - run.start), DISP: first pole
	uint16_t state;           //!< first float of this component in the node's state slot
	uint8_t nl, nd;           //!< DISP: Drude/Lorentz poles, Debye poles
};

enum {CM_ZERO = 0, CM_LOSSLESS = 1, CM_LOSSY = 2, CM_UPML = 3, CM_DISP = 4, CM_UPML_PROFILE = 5, CM_SKIP = 255};
enum {CM_STATE = 1, CM_LOSSLESS3 = 2};   //!< the bits of CM_Combo::state

//! The three components of a node: what a node index points to
/*!
  The header (class, m, vv per component) is all a ZERO, LOSSLESS or LOSSY
  component reads; the rest of a component (c[n]) only UPML and DISP.
  */
struct CM_Combo
{
	uint8_t cls[3];     //!< the class of each component, CellModel_Set::Class or CM_SKIP: not updated
	uint8_t state;      //!< CM_STATE: a component has state (UPML, DISP), only then the node's slot is read; CM_LOSSLESS3: all three LOSSLESS
	float vv[3];        //!< CellModel_Set::vv
	float m[3];         //!< material part of vi (iv), CellModel_Set::m (CM_FACTOR)
	CM_Comp c[3];
};

//! A pole, as CellModel_Pole
struct CM_Pole
{
	float k_int, k_ext, k_lor;
};

//! What the kernel of one field (voltages or currents) reads besides the fields
struct CM_Field
{
	CM_GLOBAL const CM_Combo* combos;
	CM_GLOBAL const uint32_t* slot;      //!< per node its state slot, 0xffffffff none
	CM_GLOBAL float* state;              //!< float j of slot s at [j*numSlots + s]
	uint64_t numSlots;
	CM_GLOBAL const CM_Pole* poles;
	CM_GLOBAL const float* profile[5];   //!< m, vv, uvv, ufo, ufn
	CM_GLOBAL const float* gxy[3];       //!< per component the mesh factors of the line (x, y) at [x*ny + y]: gx*gy, as the CPU engine's "slow" (CM_FACTOR)
	CM_GLOBAL const float* gz[3];        //!< of each z line (the fast axis)
};

//! The new value of component \a n of a node (line z) from its value \a x and \a curl; its class is neither CM_SKIP nor CM_ZERO
/*!
  Reads the node's state (UPML, dispersive materials) from F.state and, if \a store, writes the
  new one to \a sout: F.state in place, or a second buffer (the fused step). Without \a store
  the same value for a neighbouring tile of the fused step, which reads it, but does not own it.
  */
CM_FN float cm_apply_component(CM_GLOBAL const CM_Combo& cb, CM_PRIVATE const CM_Field& F, int n, uint64_t line, unsigned int z,
                               float x, float curl, uint32_t slot, bool store, CM_GLOBAL float* sout)
{
	const uint8_t cls = cb.cls[n];
	const float a = CM_FACTOR(cb, F, n, line);
	const float gf = F.gz[n][z];
	if (cls==CM_LOSSLESS)
		return x + a*gf*curl;
	if (cls==CM_LOSSY)
		return cb.vv[n]*x + a*gf*curl;
	CM_GLOBAL const CM_Comp& c = cb.c[n];
	switch (cls)
	{
	case CM_UPML:
	{
		const uint64_t o = c.state*F.numSlots + slot;
		const float fl = F.state[o];
		const float Fn = cb.vv[n]*fl + a*gf*curl;
		if (store)
			sout[o] = Fn;
		return c.uvv*x - c.ufo*fl + c.ufn*Fn;
	}
	case CM_UPML_PROFILE:
	{
		const int64_t p = (int64_t)c.p0 + z;
		const uint64_t o = c.state*F.numSlots + slot;
		const float fl = F.state[o];
		const float Fn = F.profile[1][p]*fl + a*F.profile[0][p]*gf*curl;
		if (store)
			sout[o] = Fn;
		return F.profile[2][p]*x - F.profile[3][p]*fl + F.profile[4][p]*Fn;
	}
	case CM_DISP:
	{
		// as UpdateRunDisp(): the ADE pre passes, the main update, the Debye solve and post pass, the apply pass
		CM_GLOBAL const CM_Pole* P = F.poles + c.p0;
		CM_GLOBAL const CM_Pole* D = P + c.nl;
		const uint64_t ns = F.numSlots;
		uint64_t j = c.state;
		float X0 = x;
		for (unsigned int p=0; p<c.nl; ++p)
		{
			const uint64_t oJ = (j++)*ns + slot;
			if (P[p].k_lor!=0)
			{
				const uint64_t oL = (j++)*ns + slot;
				const float l = F.state[oL] + P[p].k_lor*F.state[oJ];
				const float jn = F.state[oJ]*P[p].k_int + P[p].k_ext*(X0 - l);
				if (store)
				{
					sout[oL] = l;
					sout[oJ] = jn;
				}
			}
			else
			{
				const float jn = F.state[oJ]*P[p].k_int + P[p].k_ext*X0;
				if (store)
					sout[oJ] = jn;
			}
		}
		const uint64_t jd = j;
		float pre = 0;
		for (unsigned int p=0; p<c.nd; ++p)
		{
			const uint64_t oW = (jd+p)*ns + slot;
			const float w = D[p].k_int*F.state[oW] + D[p].k_ext*X0;
			pre += w - F.state[oW];
			if (store)
				sout[oW] = w;
		}
		X0 = cb.vv[n]*X0 + a*gf*curl;
		if (c.nd>0)
		{
			X0 = (X0 - pre)*c.solve;
			if (store)
				for (unsigned int p=0; p<c.nd; ++p)
					sout[(jd+p)*ns + slot] += D[p].k_ext*X0;
		}
		// the new currents of the poles: stored above, or the same again from the old state
		j = c.state;
		for (unsigned int p=0; p<c.nl; ++p)
		{
			const uint64_t oJ = j*ns + slot;
			if (store)
				X0 -= sout[oJ];
			else if (P[p].k_lor!=0)
			{
				const float l = F.state[oJ + ns] + P[p].k_lor*F.state[oJ];
				X0 -= F.state[oJ]*P[p].k_int + P[p].k_ext*(x - l);
			}
			else
				X0 -= F.state[oJ]*P[p].k_int + P[p].k_ext*x;
			j += (P[p].k_lor!=0) ? 2 : 1;
		}
		return X0;
	}
	}
	return x;
}

//! The new value of component \a n of a node from its value \a x and \a curl, any class (see cm_apply_component())
CM_FN float cm_component_value(CM_GLOBAL const CM_Combo& cb, CM_PRIVATE const CM_Field& F, int n, uint64_t line, unsigned int z,
                               float x, float curl, uint32_t slot, bool store, CM_GLOBAL float* sout)
{
	const uint8_t cls = cb.cls[n];
	if (cls==CM_SKIP)
		return x;
	if (cls==CM_ZERO)
		return 0;
	return cm_apply_component(cb, F, n, line, z, x, curl, slot, store, sout);
}

//! Update component \a n of field \a X at node \a i (line z), from the curl of \a A (component nPP) and \a B (nP) with neighbour offsets \a dA, \a dB
CM_FN void cm_update_component(CM_GLOBAL const CM_Combo& cb, CM_PRIVATE const CM_Field& F, int n, uint64_t i, uint64_t sN, uint64_t line, unsigned int z,
                               CM_GLOBAL float* CM_RESTRICT X, CM_GLOBAL const float* CM_RESTRICT A, CM_GLOBAL const float* CM_RESTRICT B, int64_t dA, int64_t dB, uint32_t slot)
{
	const uint8_t cls = cb.cls[n];
	if (cls==CM_SKIP)
		return;
	CM_GLOBAL float* x_ = X + n*sN + i;
	if (cls==CM_ZERO)
	{
		*x_ = 0;
		return;
	}
	*x_ = cm_apply_component(cb, F, n, line, z, *x_, A[i] - A[i+dA] - B[i] + B[i+dB], slot, true, F.state);
}

//! The new voltages of node (x, y, z) in \a v0, \a v1, \a v2, from the voltages \a V and the currents \a I at the node and below (at line 0 the
//! neighbour is the node itself); the node's new state to \a sout only if \a store (see cm_apply_component()). The index
//! arithmetic in \a A: 32 bit if all 3*nx*ny*nz values of a field have such an index (fewer instructions on a device)
template <typename Index, typename A = uint64_t>
CM_FN void cm_voltage_values(CM_GLOBAL const float* CM_RESTRICT V, CM_GLOBAL const float* CM_RESTRICT I, CM_GLOBAL const Index* index, CM_PRIVATE const CM_Field& F,
                             unsigned int nx, unsigned int ny, unsigned int nz, unsigned int x, unsigned int y, unsigned int z,
                             CM_PRIVATE float& v0, CM_PRIVATE float& v1, CM_PRIVATE float& v2, bool store, CM_GLOBAL float* sout)
{
	const A sz = 1, sy = nz, sx = (A)ny*nz, sN = (A)nx*sx;
	const A i = x*sx + y*sy + z;
	CM_GLOBAL const CM_Combo& cb = F.combos[index[i]];
	const A line = (A)x*ny + y;
	const A dx = (x>0) ? sx : 0, dy = (y>0) ? sy : 0, dz = (z>0) ? sz : 0;
	CM_GLOBAL const float* Ix = I;
	CM_GLOBAL const float* Iy = I + sN;
	CM_GLOBAL const float* Iz = I + 2*sN;
	// the components written out, a constant n each: no array indexed at run time (a device copies it to local memory)
	const float c0 = Iz[i] - Iz[i-dy] - Iy[i] + Iy[i-dz];
	const float c1 = Ix[i] - Ix[i-dz] - Iz[i] + Iz[i-dx];
	const float c2 = Iy[i] - Iy[i-dx] - Ix[i] + Ix[i-dy];
	if (cb.state & CM_LOSSLESS3)
	{
		// the common node: the operations of cm_apply_component() for LOSSLESS
		const float a0 = CM_FACTOR(cb, F, 0, line), a1 = CM_FACTOR(cb, F, 1, line), a2 = CM_FACTOR(cb, F, 2, line);
		v0 = V[i]      + a0*F.gz[0][z]*c0;
		v1 = V[sN+i]   + a1*F.gz[1][z]*c1;
		v2 = V[2*sN+i] + a2*F.gz[2][z]*c2;
		return;
	}
	const uint32_t slot = (cb.state & CM_STATE) ? F.slot[i] : 0;
	v0 = cm_component_value(cb, F, 0, line, z, V[i], c0, slot, store, sout);
	v1 = cm_component_value(cb, F, 1, line, z, V[sN+i], c1, slot, store, sout);
	v2 = cm_component_value(cb, F, 2, line, z, V[2*sN+i], c2, slot, store, sout);
}

//! The voltages of node (x, y, z), in place
template <typename Index>
CM_FN void cm_update_voltages(CM_GLOBAL float* CM_RESTRICT V, CM_GLOBAL const float* CM_RESTRICT I, CM_GLOBAL const Index* index, CM_PRIVATE const CM_Field& F,
                              unsigned int nx, unsigned int ny, unsigned int nz, unsigned int x, unsigned int y, unsigned int z)
{
	float v0, v1, v2;
	cm_voltage_values<Index>(V, I, index, F, nx, ny, nz, x, y, z, v0, v1, v2, true, F.state);
	const uint64_t sN = (uint64_t)nx*ny*nz, i = ((uint64_t)x*ny + y)*nz + z;
	V[i] = v0;
	V[sN+i] = v1;
	V[2*sN+i] = v2;
}

//! The currents of node (x, y, z): the voltages at it and above (the last lines are CM_SKIP)
template <typename Index>
CM_FN void cm_update_currents(CM_GLOBAL float* CM_RESTRICT I, CM_GLOBAL const float* CM_RESTRICT V, CM_GLOBAL const Index* index, CM_PRIVATE const CM_Field& F,
                              unsigned int nx, unsigned int ny, unsigned int nz, unsigned int x, unsigned int y, unsigned int z)
{
	const uint64_t sz = 1, sy = nz, sx = (uint64_t)ny*nz, sN = (uint64_t)nx*sx;
	const uint64_t i = x*sx + y*sy + z;
	CM_GLOBAL const CM_Combo& cb = F.combos[index[i]];
	const uint64_t line = (uint64_t)x*ny + y;
	if (cb.state & CM_LOSSLESS3)
	{
		// the common node, all loads first (never on the last lines: those are CM_SKIP)
		CM_GLOBAL const float* Vx = V;
		CM_GLOBAL const float* Vy = V + sN;
		CM_GLOBAL const float* Vz = V + 2*sN;
		const float c0 = Vz[i] - Vz[i+sy] - Vy[i] + Vy[i+sz];
		const float c1 = Vx[i] - Vx[i+sz] - Vz[i] + Vz[i+sx];
		const float c2 = Vy[i] - Vy[i+sx] - Vx[i] + Vx[i+sy];
		const float a0 = CM_FACTOR(cb, F, 0, line), a1 = CM_FACTOR(cb, F, 1, line), a2 = CM_FACTOR(cb, F, 2, line);
		I[i]      += a0*F.gz[0][z]*c0;
		I[sN+i]   += a1*F.gz[1][z]*c1;
		I[2*sN+i] += a2*F.gz[2][z]*c2;
		return;
	}
	const uint32_t slot = (cb.state & CM_STATE) ? F.slot[i] : 0;
	const uint64_t s[3] = {sx, sy, sz};
	for (int n=0; n<3; ++n)
	{
		const int nP = (n+1)%3, nPP = (n+2)%3;
		cm_update_component(cb, F, n, i, sN, line, z, I, V + nPP*sN, V + nP*sN, (int64_t)s[nP], (int64_t)s[nPP], slot);
	}
}

#endif // FDTD_KERNELS_CELL_UPDATE_H
