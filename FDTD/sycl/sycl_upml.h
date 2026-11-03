/*
*	Copyright (C) 2026 Sean Mollet (sean@malmoset.com)
*	Copyright (C) 2026 Alexandru Iercosan-Lucaci (alex.lucaci@alacrity.ro)
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

// SYCL only: the fused UPML update of a node, shared by the UPML kernels and the main
// updates (which update the UPML regions along z, see SYCL_ZSlabs)

#ifndef SYCL_UPML_H
#define SYCL_UPML_H

// Uniaxial PML, see Engine_Ext_UPML. One thread per cell of the PML region.
// The flux and coefficient arrays use the local region size.
struct UPMLRegion { unsigned int sx, sy, sz; unsigned int lx, ly, lz; };

template<typename I = size_t>
inline I upml_local(const UPMLRegion& R, unsigned int n, unsigned int x, unsigned int y, unsigned int z)
{
	return (((I)n*R.lx + x)*R.ly + y)*R.lz + z;
}

//! Arguments of the UPML kernels: the flux and UPML coefficients of the region, the main coefficients (fused kernels only)
struct UPMLFusedArgs
{
	float* flux;
	float* flux_out;   //!< the new voltage flux: flux (in place) or the second buffer of the fused step (see update_fused)
	const float* c_old; const float* c_fo; const float* c_fn;   //!< full UPML coefficient arrays (mode 0)
	const void* u_index; const float* u_sets; unsigned int u_mode;
	const void* index; const float* ca; const float* cb; unsigned int mode;   //!< main coefficients, see main_coeff()
};

//! UPML coefficients of the region cell l0: old[n*s], fo[n*s] and fn[n*s] are vv/vvfo/vvfn (set_offset 0)
//! or ii/iifo/iifn (set_offset 9) of direction n. Mode 0: the full arrays, else a 16/32 bit set index and the sets of 18.
template<typename I> struct UPMLCoeffT { const float* old; const float* fo; const float* fn; I s; };
typedef UPMLCoeffT<size_t> UPMLCoeff;

template<typename I = size_t>
inline UPMLCoeffT<I> upml_coeff(const UPMLFusedArgs& A, typename SYCL_NoDeduce<I>::type l0,
                                typename SYCL_NoDeduce<I>::type cells, unsigned int set_offset)
{
	UPMLCoeffT<I> c;
	if (A.u_mode==0)
	{
		c.old = A.c_old + l0; c.fo = A.c_fo + l0; c.fn = A.c_fn + l0; c.s = cells;
		return c;
	}
	const unsigned int set = (A.u_mode==1) ? (unsigned int)((const unsigned short*)A.u_index)[l0] : ((const unsigned int*)A.u_index)[l0];
	c.old = A.u_sets + (size_t)18*set + set_offset; c.fo = c.old + 3; c.fn = c.old + 6; c.s = 1;
	return c;
}

// Fused: upml_pre, the main update (update_voltages) and upml_post of the region node
// (x, y, z), operation by operation the same (see GPU_UPMLFusionBox() for the conditions).
// Reads volt_in, writes volt_out: the same field, or the two buffers of the fused step.
// FACTOR: the main coefficients of modes 3/4 (see main_coeff()).
template<bool FACTOR>
inline void upml_fused_volt_node(float* volt_out, const float* volt_in, const float* curr, const UPMLFusedArgs& A,
                                 const SYCL_GridDim& N, const UPMLRegion& R,
                                 unsigned int x, unsigned int y, unsigned int z)
{
	const unsigned int gx = R.sx + x, gy = R.sy + y, gz = R.sz + z;
	const size_t sn = (size_t)N.nx*N.ny*N.nz;
	const size_t i  = nijk(N, 0, gx, gy, gz);
	const size_t xm = (gx>0) ? (size_t)N.ny*N.nz : 0;
	const size_t ym = (gy>0) ? N.nz : 0;
	const size_t zm = (gz>0) ? 1 : 0;
	const SYCL_MainCoeffT<size_t, FACTOR> C = main_coeff<size_t, FACTOR>(A.index, A.ca, A.cb, A.mode, sn, i, 0, N, gx, gy, gz);
	const size_t cells = (size_t)R.lx*R.ly*R.lz;
	const UPMLCoeff U = upml_coeff(A, (x*R.ly + y)*R.lz + z, cells, 0);
	float curl[3];
	curl[0] = (curr[2*sn+i] - curr[2*sn+i-ym] - curr[sn+i] + curr[sn+i-zm]);
	curl[1] = (curr[i] - curr[i-zm] - curr[2*sn+i] + curr[2*sn+i-xm]);
	curl[2] = (curr[sn+i] - curr[sn+i-xm] - curr[i] + curr[i-ym]);
	for (unsigned int n=0; n<3; ++n)
	{
		const size_t l = upml_local(R, n, x, y, z);
		const size_t g = n*sn + i;
		const float f_help = U.old[n*U.s]*volt_in[g] - U.fo[n*U.s]*A.flux[l];
		float v;
		v  = A.flux[l] * C.a[n*C.s];
		v += C.b[n*C.s] * curl[n];
		A.flux[l] = v;
		volt_out[g] = f_help + U.fn[n*U.s]*v;
	}
}

// same for the currents (update_currents), which are not updated on the last mesh lines
template<bool FACTOR>
inline void upml_fused_curr_node(float* curr_out, const float* curr_in, const float* volt, const UPMLFusedArgs& A,
                                 const SYCL_GridDim& N, const UPMLRegion& R,
                                 unsigned int x, unsigned int y, unsigned int z)
{
	const unsigned int gx = R.sx + x, gy = R.sy + y, gz = R.sz + z;
	const size_t sn = (size_t)N.nx*N.ny*N.nz;
	const size_t i  = nijk(N, 0, gx, gy, gz);
	const bool update = (gx+1<N.nx) && (gy+1<N.ny) && (gz+1<N.nz);
	const size_t cells = (size_t)R.lx*R.ly*R.lz;
	const UPMLCoeff U = upml_coeff(A, (x*R.ly + y)*R.lz + z, cells, 9);
	float curl[3] = {0, 0, 0};
	SYCL_MainCoeffT<size_t, FACTOR> C = {};
	if (update)
	{
		const size_t xp = (size_t)N.ny*N.nz;
		const size_t yp = N.nz;
		C = main_coeff<size_t, FACTOR>(A.index, A.ca, A.cb, A.mode, sn, i, 6, N, gx, gy, gz);
		curl[0] = (volt[2*sn+i] - volt[2*sn+i+yp] - volt[sn+i] + volt[sn+i+1]);
		curl[1] = (volt[i] - volt[i+1] - volt[2*sn+i] + volt[2*sn+i+xp]);
		curl[2] = (volt[sn+i] - volt[sn+i+xp] - volt[i] + volt[i+yp]);
	}
	for (unsigned int n=0; n<3; ++n)
	{
		const size_t l = upml_local(R, n, x, y, z);
		const size_t g = n*sn + i;
		const float f_help = U.old[n*U.s]*curr_in[g] - U.fo[n*U.s]*A.flux[l];
		float c = A.flux[l];
		if (update)
		{
			c  = c * C.a[n*C.s];
			c += C.b[n*C.s] * curl[n];
		}
		A.flux[l] = c;
		curr_out[g] = f_help + U.fn[n*U.s]*c;
	}
}

// The fused step (update_fused) with the UPML z slabs in the kernel: the voltages of a slab node
// (gx, gy, gz) of region R from the current fields and flux (A.flux), the new flux into A.flux_out
// if write (the node's own work-group; others compute the node again as a border voltage), the
// same operations as upml_fused_volt_node. All loads come before the stores, which could alias
// them for the compiler.
template<typename I, bool FACTOR>
inline void upml_volt_value(const float* __restrict__ volt, const float* __restrict__ curr, const UPMLFusedArgs& A,
                            const SYCL_GridDim& N, const UPMLRegion& R,
                            unsigned int gx, unsigned int gy, unsigned int gz, bool write, float out[3])
{
	const I sn = (I)N.nx*N.ny*N.nz;
	const I i  = nijk<I>(N, 0, gx, gy, gz);
	const I xm = (gx>0) ? (I)N.ny*N.nz : 0;
	const I ym = (gy>0) ? N.nz : 0;
	const I zm = (gz>0) ? 1 : 0;
	const SYCL_MainCoeffT<I, FACTOR> C = main_coeff<I, FACTOR>(A.index, A.ca, A.cb, A.mode, sn, i, 0, N, gx, gy, gz);
	const I l0 = ((I)(gx-R.sx)*R.ly + (gy-R.sy))*R.lz + (gz-R.sz);
	const I cells = (I)R.lx*R.ly*R.lz;
	const UPMLCoeffT<I> U = upml_coeff<I>(A, l0, cells, 0);
	float curl[3];
	curl[0] = (curr[2*sn+i] - curr[2*sn+i-ym] - curr[sn+i] + curr[sn+i-zm]);
	curl[1] = (curr[i] - curr[i-zm] - curr[2*sn+i] + curr[2*sn+i-xm]);
	curl[2] = (curr[sn+i] - curr[sn+i-xm] - curr[i] + curr[i-ym]);
	float flux[3], f_help[3];
	for (unsigned int n=0; n<3; ++n)
	{
		flux[n] = A.flux[n*cells + l0];   // upml_local()
		f_help[n] = U.old[n*U.s]*volt[n*sn + i] - U.fo[n*U.s]*flux[n];
	}
	for (unsigned int n=0; n<3; ++n)
	{
		float v;
		v  = flux[n] * C.a[n*C.s];
		v += C.b[n*C.s] * curl[n];
		if (write)
			A.flux_out[n*cells + l0] = v;
		out[n] = f_help[n] + U.fn[n*U.s]*v;
	}
}

// ... and the currents of a slab node (the node's own work-group only, the flux in place), with the
// curl of the new voltages (the operations of upml_fused_curr_node); update: not on the last mesh lines
template<typename I, bool FACTOR>
inline void upml_curr_value(const float* __restrict__ curr_in, float* __restrict__ curr_out, const UPMLFusedArgs& A,
                            const SYCL_GridDim& N, const UPMLRegion& R,
                            unsigned int gx, unsigned int gy, unsigned int gz, bool update, const float curl[3])
{
	const I sn = (I)N.nx*N.ny*N.nz;
	const I i  = nijk<I>(N, 0, gx, gy, gz);
	const I l0 = ((I)(gx-R.sx)*R.ly + (gy-R.sy))*R.lz + (gz-R.sz);
	const I cells = (I)R.lx*R.ly*R.lz;
	const UPMLCoeffT<I> U = upml_coeff<I>(A, l0, cells, 9);
	SYCL_MainCoeffT<I, FACTOR> C = {};
	if (update)
		C = main_coeff<I, FACTOR>(A.index, A.ca, A.cb, A.mode, sn, i, 6, N, gx, gy, gz);
	float c[3], f_help[3];
	for (unsigned int n=0; n<3; ++n)
	{
		c[n] = A.flux[n*cells + l0];
		f_help[n] = U.old[n*U.s]*curr_in[n*sn + i] - U.fo[n*U.s]*c[n];
	}
	for (unsigned int n=0; n<3; ++n)
	{
		if (update)
		{
			c[n] = c[n] * C.a[n*C.s];
			c[n] += C.b[n*C.s] * curl[n];
		}
		A.flux[n*cells + l0] = c[n];
		curr_out[n*sn + i] = f_help[n] + U.fn[n*U.s]*c[n];
	}
}

//! UPML regions at the lower/upper z end covering the x/y range of the main updates: the
//! main kernels update them too (full z lines instead of short rows, see sycl_ext_upml.cpp),
//! the separate ones (update_voltages/update_currents) as the fused step (update_fused)
struct SYCL_ZSlabs
{
	unsigned int z_lo, z_hi;       //!< main update nodes: z_lo <= z < z_hi, below: slab lo, above: slab hi
	unsigned int t_lo, t_int;      //!< separate updates: threads along z of slab lo and of the main update nodes, whole sub-groups (0: no padding)
	UPMLRegion R_lo, R_hi;
	//! the fused step reads the voltage flux from flux and writes the new one into flux_out (two
	//! buffers, the UPML hooks swap them every step), the separate updates in place
	UPMLFusedArgs volt_lo, volt_hi, curr_lo, curr_hi;
};

//! \a lo ? \a a : \a b, field by field: a pick of kernel arguments that stays in registers (a
//! reference to one of two by-value structs makes the work-items copy them to private memory);
//! \a lo should be the same for a work-group, else every work-item holds its own pointers
inline UPMLFusedArgs upml_pick(bool lo, const UPMLFusedArgs& a, const UPMLFusedArgs& b)
{
	UPMLFusedArgs r;
	r.flux = lo ? a.flux : b.flux;
	r.flux_out = lo ? a.flux_out : b.flux_out;
	r.c_old = lo ? a.c_old : b.c_old;
	r.c_fo = lo ? a.c_fo : b.c_fo;
	r.c_fn = lo ? a.c_fn : b.c_fn;
	r.u_index = lo ? a.u_index : b.u_index;
	r.u_sets = lo ? a.u_sets : b.u_sets;
	r.u_mode = lo ? a.u_mode : b.u_mode;
	r.index = a.index; r.ca = a.ca; r.cb = a.cb; r.mode = a.mode;   // the main coefficients, the same in both
	return r;
}

inline UPMLRegion upml_pick(bool lo, const UPMLRegion& a, const UPMLRegion& b)
{
	UPMLRegion r;
	r.sx = lo ? a.sx : b.sx; r.sy = lo ? a.sy : b.sy; r.sz = lo ? a.sz : b.sz;
	r.lx = lo ? a.lx : b.lx; r.ly = lo ? a.ly : b.ly; r.lz = lo ? a.lz : b.lz;
	return r;
}

#endif // SYCL_UPML_H
