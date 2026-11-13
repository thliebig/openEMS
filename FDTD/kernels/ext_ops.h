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

#ifndef FDTD_KERNELS_EXT_OPS_H
#define FDTD_KERNELS_EXT_OPS_H

// The per-point operations of the device extensions (GPU_Ext_*, gpu_ext.h), one struct per
// operation: its data and operator()(n) for point n, run by GPU_Backend::RunPoints() on every
// backend that uses the cell model. Each does what its engine extension does for one point,
// operation for operation, so that the results are bit-identical to it.
//
// The portable subset of kernels/cell_update.h: plain data, no std, address spaces by macro.

#if !defined(__METAL_VERSION__)
#include "cell_update.h"
#endif

//! Soft excitation of the voltages or currents (Engine_Ext_Excitation): field += amp*signal
struct Op_Excitation
{
	CM_GLOBAL float* field;
	CM_GLOBAL const uint64_t* index;     //!< the flat ArrayNIJK index of each point
	CM_GLOBAL const float* amp;
	CM_GLOBAL const uint32_t* delay;
	CM_GLOBAL const float* signal;
	int32_t numTS;
	int32_t period;   //!< numTS+1 without a signal period
	int32_t length;

	CM_FN void operator()(uint64_t n) const
	{
		int exc_pos = numTS - (int)delay[n];
		if (exc_pos<0)   // delay not yet expired
			return;
		exc_pos %= period;
		exc_pos *= (exc_pos<length);
		field[index[n]] = field[index[n]] + amp[n]*signal[exc_pos];
	}
};

//! A boundary sheet of first order Mur (Engine_Ext_Mur_ABC, Engine_Ext_Absorbing_BC): per point of the sheet its two tangential components
/*!
  hook 0 (pre): s = E(shift) - K1*E(line); 1 (post): s += K1*E(shift); 2 (apply): E(line) = s;
  3 (apply, super-absorption of the currents): E(line) = (s*K2 + E(line))/(K2 + 1).
  */
struct Op_Sheet
{
	CM_GLOBAL float* field;
	CM_GLOBAL const uint64_t* line;    //!< per point the node index on the sheet
	CM_GLOBAL const uint64_t* shift;   //!< and on the line next to it
	CM_GLOBAL const float* K1[2];      //!< per component (nyP, nyPP)
	CM_GLOBAL const float* K2[2];
	CM_GLOBAL float* s[2];             //!< the value between the hooks
	uint64_t offset[2];                //!< the component offsets in the field (nyP, nyPP times the node count)
	int32_t hook;

	CM_FN void operator()(uint64_t n) const
	{
		for (int c=0; c<2; ++c)
		{
			const uint64_t il = offset[c] + line[n], is = offset[c] + shift[n];
			switch (hook)
			{
			case 0: s[c][n] = field[is] - K1[c][n]*field[il]; break;
			case 1: s[c][n] += K1[c][n]*field[is]; break;
			case 2: field[il] = s[c][n]; break;
			default: field[il] = (s[c][n]*K2[c][n] + field[il])/(K2[c][n] + 1.0f); break;
			}
		}
	}
};

//! Lumped RLC elements (Engine_Ext_LumpedRLC): mode 0 the parallel inductor current (pre),
//! 1 read the node voltage, 2 the new node voltage and J, 3 write it (the passes of the engine extension)
struct Op_LumpedRLC
{
	CM_GLOBAL float* volt;
	CM_GLOBAL const uint64_t* index;
	CM_GLOBAL const float *ilv, *i2v, *vv2, *vj1, *vj2, *vvd, *ib0, *b1, *b2;
	CM_GLOBAL float* Il;
	CM_GLOBAL float *Vd0, *Vd1, *Vd2;   //!< the rings at this timestep: ts, ts-1, ts-2
	CM_GLOBAL float *J0, *J1, *J2;
	int32_t mode;

	CM_FN void operator()(uint64_t n) const
	{
		switch (mode)
		{
		case 0: Il[n] += (i2v[n])*(ilv[n])*Vd1[n]; break;
		case 1: Vd0[n] = volt[index[n]]; break;
		case 2:
			Vd0[n] = (vvd[n])*(Vd0[n] - Il[n] + (vv2[n])*Vd2[n] + (vj1[n])*J1[n] + (vj2[n])*J2[n]);
			J0[n] = (ib0[n])*(Vd0[n] - Vd2[n]) - ((b1[n])*(ib0[n]))*J1[n] - ((b2[n])*(ib0[n]))*J2[n];
			break;
		default: volt[index[n]] = Vd0[n]; break;
		}
	}
};

//! The planes of a TF/SF box (Engine_Ext_TFSF): field += (1-dd)*amp*signal(delay) + dd*amp*signal(delay+1)
struct Op_TFSF
{
	CM_GLOBAL float* field;
	CM_GLOBAL const uint64_t* index;
	CM_GLOBAL const uint32_t* delay;
	CM_GLOBAL const float* dd;
	CM_GLOBAL const float* amp;
	CM_GLOBAL const float* signal;
	uint32_t ts, length;
	int32_t period;   //!< in timesteps, 0 none

	CM_FN unsigned int Pos(unsigned int d) const
	{
		unsigned int v;
		if (ts<d)
			v = 0;
		else if ((ts-d >= length) && (period==0))
			v = 0;
		else
			v = ts-d;
		if (period>0)
			v = v % period;
		return v;
	}

	CM_FN void operator()(uint64_t n) const
	{
		const float s0 = signal[Pos(delay[n])], s1 = signal[Pos(delay[n]+1)];
		field[index[n]] = field[index[n]] + (1.0f-dd[n])*amp[n]*s0 + dd[n]*amp[n]*s1;
	}
};

//! The probe voltages of the steady-state detection (Engine_Ext_SteadyState), recorded per timestep
struct Op_SteadyState
{
	CM_GLOBAL const float* volt;
	CM_GLOBAL const uint64_t* index;
	CM_GLOBAL float* records;   //!< [probe][length]
	uint32_t length, pos;

	CM_FN void operator()(uint64_t n) const
	{
		records[n*length + pos] = volt[index[n]];
	}
};

//! A closed cylindrical mesh (Engine_Ext_Cylinder); mode 0: the voltage on the r=0 axis from the alpha currents around it, per z;
//! 1: its copy to all alpha and no alpha voltage at r=0, per (alpha, z); 2: the voltages of the first alpha plane from the last,
//! 3: the currents of the last alpha plane from the first, per (r, z)
struct Op_Cylinder
{
	CM_GLOBAL float* volt;
	CM_GLOBAL float* curr;
	CM_GLOBAL const float* vv_R0;
	CM_GLOBAL const float* vi_R0;
	uint32_t ny, nz;      //!< the lines along alpha and z
	uint32_t z0, nzr;     //!< mode 2/3: the first z and the z of a point (n = r*nzr + z - z0)
	uint32_t r0;          //!< mode 2/3: the first r
	uint64_t N;           //!< nodes
	int32_t mode;

	CM_FN uint64_t Node(uint64_t r, uint64_t a, uint64_t z) const {return (r*ny + a)*nz + z;}

	CM_FN void operator()(uint64_t n) const
	{
		switch (mode)
		{
		case 0:
		{
			const uint64_t z = n;
			float v = volt[2*N + Node(0, 0, z)]*vv_R0[z];
			for (uint64_t a=0; a<ny-1; ++a)
				v += vi_R0[z] * curr[N + Node(0, a, z)];
			volt[2*N + Node(0, 0, z)] = v;
			break;
		}
		case 1:
		{
			const uint64_t a = n/nz, z = n%nz;
			volt[N + Node(0, a, z)] = 0;
			if (a>0)   // alpha 0 is the source, which the other points read
				volt[2*N + Node(0, a, z)] = volt[2*N + Node(0, 0, z)];
			break;
		}
		case 2:
		{
			const uint64_t r = r0 + n/nzr, z = z0 + n%nzr;
			volt[Node(r, 0, z)] = volt[Node(r, ny-2, z)];
			volt[2*N + Node(r, 0, z)] = volt[2*N + Node(r, ny-2, z)];
			break;
		}
		default:
		{
			const uint64_t r = r0 + n/nzr, z = z0 + n%nzr;
			curr[Node(r, ny-2, z)] = curr[Node(r, 0, z)];
			curr[2*N + Node(r, ny-2, z)] = curr[2*N + Node(r, 0, z)];
			break;
		}
		}
	}
};

//! The coupling of a cylindrical multi-grid level and its sub-grid (SubGrid_Coupling_Basic); mode 0: the base grid voltages at the
//! r line \a r0 to the sub-grid, per (h, z) (every other base alpha line); 1: the sub-grid field interpolated to the base grid at the
//! r lines from \a r0, per (dr, alpha, z)
struct Op_MultiGrid
{
	CM_GLOBAL float* base;
	CM_GLOBAL float* sub;
	CM_GLOBAL const uint32_t* pos_2p;   //!< [n][alpha], n 0 for r and z, 1 for alpha
	CM_GLOBAL const uint32_t* pos_2pp;
	CM_GLOBAL const float* w_2p;
	CM_GLOBAL const float* w_2pp;
	uint32_t nb[3], ns[3];   //!< the lines of the base and the sub-grid
	uint32_t r0;
	int32_t mode;

	CM_FN uint64_t B(uint32_t n, uint64_t r, uint64_t a, uint64_t z) const {return (((uint64_t)n*nb[0] + r)*nb[1] + a)*nb[2] + z;}
	CM_FN uint64_t S(uint32_t n, uint64_t r, uint64_t a, uint64_t z) const {return (((uint64_t)n*ns[0] + r)*ns[1] + a)*ns[2] + z;}

	CM_FN void operator()(uint64_t i) const
	{
		const uint64_t z = i%nb[2];
		if (mode==0)
		{
			const uint64_t h = i/nb[2], a = 2*h;
			sub[S(0, r0, h, z)] = 0;
			sub[S(2, r0, h, z)] = base[B(2, r0, a, z)];
			sub[S(1, r0, h, z)] = base[B(1, r0, a, z)] + base[B(1, r0, a+1, z)];
			return;
		}
		const uint64_t a = (i/nb[2])%nb[1], r = r0 + i/((uint64_t)nb[1]*nb[2]);
		const uint64_t k1 = nb[1] + a;
		base[B(0, r, a, z)] = w_2p[a]*sub[S(0, r, pos_2p[a], z)] + w_2pp[a]*sub[S(0, r, pos_2pp[a], z)];
		base[B(2, r, a, z)] = w_2p[a]*sub[S(2, r, pos_2p[a], z)] + w_2pp[a]*sub[S(2, r, pos_2pp[a], z)];
		base[B(1, r, a, z)] = w_2p[k1]*sub[S(1, r, pos_2p[k1], z)] + w_2pp[k1]*sub[S(1, r, pos_2pp[k1], z)];
	}
};

//! The dispersive materials the cell model does not fuse (Engine_Ext_LorentzMaterial, Engine_Ext_DebyeMaterial): per point
//! (a cell of one order) its three components, the state and the coefficients of component n at [n*count + p], those of
//! Debye pole k at [(k*3 + n)*count + p]. Mode 0: the Drude/Lorentz ADE currents (pre hook); 1: the Debye pole states as
//! far as V^n fixes them, and their sum (pre); 2: the Debye solve for V^(n+1), its correction, the rest of the pole states
//! (post); 3: the correction off the field (apply)
struct Op_Dispersive
{
	CM_GLOBAL float* field;
	CM_GLOBAL const uint64_t* node;   //!< per point its node (x*ny + y)*nz + z
	CM_GLOBAL float* ade;             //!< the correction of the apply pass: the ADE current, the Debye correction
	CM_GLOBAL float* lor;             //!< Lorentz: its second state
	CM_GLOBAL float* pole;            //!< Debye: the state of each pole
	CM_GLOBAL float* pre;             //!< Debye: the sum of the pre pass
	CM_GLOBAL const float* k_int;     //!< v_int_ADE (i_int_ADE); Debye: v_relax_ADE per pole
	CM_GLOBAL const float* k_ext;     //!< v_ext_ADE (i_ext_ADE); Debye: v_drive_ADE per pole
	CM_GLOBAL const float* k_lor;     //!< v_Lor_ADE (i_Lor_ADE)
	CM_GLOBAL const float* solve;     //!< Debye: v_solve_ADE
	uint64_t sN;                      //!< the nodes of a component
	uint64_t count;                   //!< the points
	uint32_t poles;                   //!< Debye
	int32_t lorentz;                  //!< a Lorentz pole, else Drude
	int32_t mode;

	CM_FN void operator()(uint64_t p) const
	{
		for (uint64_t n=0; n<3; ++n)
		{
			const uint64_t j = n*count + p;
			CM_GLOBAL float& f = field[n*sN + node[p]];
			switch (mode)
			{
			case 0:
				if (lorentz)
				{
					lor[j] = lor[j] + k_lor[j]*ade[j];
					ade[j] = ade[j]*k_int[j];
					ade[j] = ade[j] + k_ext[j]*(f - lor[j]);
				}
				else
				{
					ade[j] = ade[j]*k_int[j];
					ade[j] = ade[j] + k_ext[j]*f;
				}
				break;
			case 1:
			{
				float sum = 0;
				for (uint64_t k=0; k<poles; ++k)
				{
					const uint64_t q = (k*3 + n)*count + p;
					const float W = k_int[q]*pole[q] + k_ext[q]*f;
					sum = sum + (W - pole[q]);
					pole[q] = W;
				}
				pre[j] = sum;
				break;
			}
			case 2:
			{
				const float v = (f - pre[j])*solve[j];
				ade[j] = f - v;
				for (uint64_t k=0; k<poles; ++k)
				{
					const uint64_t q = (k*3 + n)*count + p;
					pole[q] = pole[q] + k_ext[q]*v;
				}
				break;
			}
			default:
				f = f - ade[j];
				break;
			}
		}
	}
};

//! Every point operation as X(kind, struct, kernel): its GPU_PointOp (gpu_backend.h), its struct and the name of its kernels.
//! The backends expand the list to their launches (GPU_Backend::RunPoints()), so a new operation goes here and into
//! GPU_PointOp only
#define CM_POINT_OPS(X) \
	X(GPU_OP_EXCITATION,   Op_Excitation,  op_excitation) \
	X(GPU_OP_SHEET,        Op_Sheet,       op_sheet) \
	X(GPU_OP_LUMPED_RLC,   Op_LumpedRLC,   op_lumped_rlc) \
	X(GPU_OP_TFSF,         Op_TFSF,        op_tfsf) \
	X(GPU_OP_STEADY_STATE, Op_SteadyState, op_steady_state) \
	X(GPU_OP_CYLINDER,     Op_Cylinder,    op_cylinder) \
	X(GPU_OP_MULTIGRID,    Op_MultiGrid,   op_multigrid) \
	X(GPU_OP_DISPERSIVE,   Op_Dispersive,  op_dispersive)

#endif // FDTD_KERNELS_EXT_OPS_H
