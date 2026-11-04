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

// The update kernels of Engine_CPU, included once per instruction set into
// the namespace ENGINE_CPU_NAMESPACE (engine_cpu_kernels*.cpp). Only static
// functions and the two entry points: nothing here may be shared with the
// other instruction sets through the linker (see engine_cpu_kernel.h).
//
// Per line of the fast axis in the box and per component, the runs of the
// line clipped to the box. The curl of a cell reads A and B at k and at
// k+dA, k+dB: the neighbours below for a voltage, above for a current; along
// the fast axis the offset is 1.

namespace ENGINE_CPU_NAMESPACE
{

static inline unsigned int Max(unsigned int a, unsigned int b) {return a>b ? a : b;}
static inline unsigned int Min(unsigned int a, unsigned int b) {return a<b ? a : b;}

//! the cells [k0, k1) of a run, with the coefficients of its set; a = m times the slow geometry factors
static inline void UpdateRun(const Engine_CPU_Run& run, const Engine_CPU_Set& set, float a, unsigned int k0, unsigned int k1,
                             float* __restrict X, const float* __restrict A, const float* __restrict B,
                             ptrdiff_t dA, ptrdiff_t dB, const Engine_CPU_Kernel::Component& c)
{
	const float vv = set.vv;
	const float* __restrict gf = c.gf;
	switch (set.cls)
	{
	case Engine_CPU_Set::ZERO:
		for (unsigned int k=k0; k<k1; ++k)
			X[k] = 0;
		break;
	case Engine_CPU_Set::LOSSLESS:
		for (unsigned int k=k0; k<k1; ++k)
			X[k] += a*gf[k]*(A[k] - A[k+dA] - B[k] + B[k+dB]);
		break;
	case Engine_CPU_Set::LOSSY:
		for (unsigned int k=k0; k<k1; ++k)
			X[k] = vv*X[k] + a*gf[k]*(A[k] - A[k+dA] - B[k] + B[k+dB]);
		break;
	case Engine_CPU_Set::UPML:
	{
		// the extension's pre pass, the main update on the flux and its post pass in one
		const ptrdiff_t j0 = (ptrdiff_t)run.u - (ptrdiff_t)run.start;
		const float uvv = set.uvv, ufo = set.ufo, ufn = set.ufn;
		float* __restrict F = c.state + j0;
		for (unsigned int k=k0; k<k1; ++k)
		{
			const float Fn = vv*F[k] + a*gf[k]*(A[k] - A[k+dA] - B[k] + B[k+dB]);
			X[k] = uvv*X[k] - ufo*F[k] + ufn*Fn;
			F[k] = Fn;
		}
		break;
	}
	case Engine_CPU_Set::DISP:
	{
		// the extensions' passes in one, per chunk of cells: the Drude/Lorentz ADE currents and the
		// Debye pole states from the old field, the main update, the Debye solve for the new field and
		// the rest of its pole states, the Drude/Lorentz apply pass (X -= J)
		const ptrdiff_t len = (ptrdiff_t)run.stop - (ptrdiff_t)run.start;
		const Engine_CPU_Pole* P = c.poles + set.pole;
		const Engine_CPU_Pole* D = P + set.nl;
		float* const S = c.state + ((ptrdiff_t)run.u - (ptrdiff_t)run.start);
		const unsigned int CHUNK = 64;
		float pre[CHUNK];
		for (unsigned int c0=k0; c0<k1; c0+=CHUNK)
		{
			const unsigned int c1 = Min(c0+CHUNK, k1);
			float* s = S;
			for (unsigned int p=0; p<set.nl; ++p)
			{
				const float ki = P[p].k_int, ke = P[p].k_ext, kl = P[p].k_lor;
				float* __restrict J = s;
				s += len;
				if (kl!=0)
				{
					float* __restrict L = s;
					s += len;
					for (unsigned int k=c0; k<c1; ++k)
					{
						L[k] += kl*J[k];
						J[k] = J[k]*ki + ke*(X[k] - L[k]);
					}
				}
				else
					for (unsigned int k=c0; k<c1; ++k)
						J[k] = J[k]*ki + ke*X[k];
			}
			float* const Wd = s;
			if (set.nd>0)
			{
				for (unsigned int k=c0; k<c1; ++k)
					pre[k-c0] = 0;
				for (unsigned int p=0; p<set.nd; ++p)
				{
					const float kr = D[p].k_int, kd = D[p].k_ext;
					float* __restrict W = Wd + p*len;
					for (unsigned int k=c0; k<c1; ++k)
					{
						const float w = kr*W[k] + kd*X[k];
						pre[k-c0] += w - W[k];
						W[k] = w;
					}
				}
			}
			for (unsigned int k=c0; k<c1; ++k)
				X[k] = vv*X[k] + a*gf[k]*(A[k] - A[k+dA] - B[k] + B[k+dB]);
			if (set.nd>0)
			{
				const float solve = set.solve;
				for (unsigned int k=c0; k<c1; ++k)
					X[k] = (X[k] - pre[k-c0])*solve;
				for (unsigned int p=0; p<set.nd; ++p)
				{
					const float kd = D[p].k_ext;
					float* __restrict W = Wd + p*len;
					for (unsigned int k=c0; k<c1; ++k)
						W[k] += kd*X[k];
				}
			}
			s = S;
			for (unsigned int p=0; p<set.nl; ++p)
			{
				const float* __restrict J = s;
				s += (P[p].k_lor!=0) ? 2*len : len;
				for (unsigned int k=c0; k<c1; ++k)
					X[k] -= J[k];
			}
		}
		break;
	}
	}
}

void UpdateVoltages(const Engine_CPU_Kernel& k, const unsigned int start[3], const unsigned int stop[3])
{
	const int F = k.F, S = k.S, M = k.M;
	unsigned int pos[3];
	pos[F] = 0;
	for (pos[S]=start[S]; pos[S]<stop[S]; ++pos[S])
		for (pos[M]=start[M]; pos[M]<stop[M]; ++pos[M])
		{
			const size_t base = pos[S]*k.s[S] + pos[M]*k.s[M];
			const size_t line = (size_t)pos[S]*k.numLines[M] + pos[M];
			for (int n=0; n<3; ++n)
			{
				const int nP = (n+1)%3, nPP = (n+2)%3;
				// V_n += vi * (I_nPP(pos) - I_nPP(pos-e_nP) - I_nP(pos) + I_nP(pos-e_nPP)); at line 0 the neighbour is the cell itself
				const ptrdiff_t dA = -(ptrdiff_t)((nP==F) ? 1 : (pos[nP]>0 ? k.s[nP] : 0));
				const ptrdiff_t dB = -(ptrdiff_t)((nPP==F) ? 1 : (pos[nPP]>0 ? k.s[nPP] : 0));
				float* V = k.field[0] + n*k.sN + base;
				const float* Ia = k.field[1] + nPP*k.sN + base;
				const float* Ib = k.field[1] + nP*k.sN + base;
				const Engine_CPU_Kernel::Component& c = k.comp[0][n];
				const double slow = c.gS[pos[S]] * c.gM[pos[M]];
				for (uint64_t i=c.first[line]; i<c.first[line+1]; ++i)
				{
					const Engine_CPU_Run& run = c.runs[i];
					unsigned int k0 = Max(run.start, start[F]);
					const unsigned int k1 = Min(run.stop, stop[F]);
					if (k0>=k1)
						continue;
					const Engine_CPU_Set& set = c.sets[run.set];
					const float a = (float)(set.m*slow);
					if ((k0==0) && (set.cls!=Engine_CPU_Set::ZERO) && ((nP==F) || (nPP==F)))
					{
						// the first line of the fast axis has no lower neighbour along it
						UpdateRun(run, set, a, 0, 1, V, Ia, Ib, (nP==F) ? 0 : dA, (nPP==F) ? 0 : dB, c);
						k0 = 1;
					}
					UpdateRun(run, set, a, k0, k1, V, Ia, Ib, dA, dB, c);
				}
			}
		}
}

void UpdateCurrents(const Engine_CPU_Kernel& k, const unsigned int start[3], const unsigned int stop[3])
{
	const int F = k.F, S = k.S, M = k.M;
	unsigned int pos[3];
	pos[F] = 0;
	for (pos[S]=start[S]; pos[S]<stop[S]; ++pos[S])
		for (pos[M]=start[M]; pos[M]<stop[M]; ++pos[M])
		{
			const size_t base = pos[S]*k.s[S] + pos[M]*k.s[M];
			const size_t line = (size_t)pos[S]*k.numLines[M] + pos[M];
			for (int n=0; n<3; ++n)
			{
				const int nP = (n+1)%3, nPP = (n+2)%3;
				// I_n += iv * (V_nPP(pos) - V_nPP(pos+e_nP) - V_nP(pos) + V_nP(pos+e_nPP))
				const ptrdiff_t dA = (nP==F) ? 1 : k.s[nP];
				const ptrdiff_t dB = (nPP==F) ? 1 : k.s[nPP];
				float* I = k.field[1] + n*k.sN + base;
				const float* Va = k.field[0] + nPP*k.sN + base;
				const float* Vb = k.field[0] + nP*k.sN + base;
				const Engine_CPU_Kernel::Component& c = k.comp[1][n];
				const double slow = c.gS[pos[S]] * c.gM[pos[M]];
				for (uint64_t i=c.first[line]; i<c.first[line+1]; ++i)
				{
					const Engine_CPU_Run& run = c.runs[i];
					const unsigned int k0 = Max(run.start, start[F]);
					const unsigned int k1 = Min(run.stop, stop[F]);
					if (k0<k1)
					{
						const Engine_CPU_Set& set = c.sets[run.set];
						UpdateRun(run, set, (float)(set.m*slow), k0, k1, I, Va, Vb, dA, dB, c);
					}
				}
			}
		}
}

}
