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

#include "engine_cpu.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

#include "tools/global.h"
#include "tools/useful.h"
#include "extensions/operator_ext_upml.h"
#include "extensions/engine_ext_upml.h"

using std::cout;
using std::endl;

// cells whose material part m differs by less than this (relative) share a run
static const double M_TOLERANCE = 1e-6;

static const char AXIS_LETTER[3] = {'I', 'J', 'K'};

Engine_CPU* Engine_CPU::New(const Operator* op, unsigned int numThreads)
{
	cout << "Create FDTD engine (" << GetTypeName(CPU) << ")" << endl;
	Engine_CPU* e = new Engine_CPU(op);
	e->SetNumThreads(numThreads);
	e->Init();
	return e;
}

Engine_CPU::Engine_CPU(const Operator* op) : Engine(op)
{
	m_type = CPU;
	m_FastAxis = 2;
	m_SplitAxis = 0;
	m_MidAxis = 1;
	m_Layout = "NIJK";
}

Engine_CPU::~Engine_CPU()
{
	Reset();
}

void Engine_CPU::Init()
{
	m_UPML.clear();
	for (size_t n=0; n<Op->GetNumberOfExtentions(); ++n)
	{
		const Operator_Ext_UPML* upml = dynamic_cast<const Operator_Ext_UPML*>(Op->GetExtension(n));
		if (upml)
			m_UPML.push_back(upml);
	}
	CalcGeometry();
	ChooseLayout();
	Engine::Init();   // the fields in m_Layout (GetFieldLayout()), the extensions
	TakeOverUPML();
	BuildRuns();
}

void Engine_CPU::TakeOverUPML()
{
	// the engine updates the UPML cells itself (UpdateVoltages(), UpdateCurrents())
	std::vector<Engine_Extension*> keep;
	for (size_t n=0; n<m_Eng_exts.size(); ++n)
	{
		if (dynamic_cast<Engine_Ext_UPML*>(m_Eng_exts[n]))
			delete m_Eng_exts[n];
		else
			keep.push_back(m_Eng_exts[n]);
	}
	m_Eng_exts = keep;
}

const Operator_Ext_UPML* Engine_CPU::FindUPML(const unsigned int pos[3], unsigned int loc[3]) const
{
	for (const Operator_Ext_UPML* upml : m_UPML)
	{
		bool inside = true;
		for (int a=0; a<3 && inside; ++a)
			inside = (pos[a]>=upml->m_StartPos[a]) && (pos[a]<upml->m_StartPos[a]+upml->m_numLines[a]);
		if (!inside)
			continue;
		for (int a=0; a<3; ++a)
			loc[a] = pos[a] - upml->m_StartPos[a];
		return upml;
	}
	return NULL;
}

void Engine_CPU::Reset()
{
	Engine::Reset();
	for (int f=0; f<2; ++f)
		for (int n=0; n<3; ++n)
		{
			m_Runs[f][n].runs.clear();
			m_Runs[f][n].first.clear();
			m_UPMLCells[f][n].Clear();
		}
}

// The factors g = L/A of a cell, per axis: for a voltage the primary edge
// length along n over the dual area across it, for a current the dual edge
// length over the primary area. Any factorisation is exact (m = vi/g takes
// the rest); this one makes m depend on the material only, on a Cartesian mesh.
void Engine_CPU::CalcGeometry()
{
	for (int a=0; a<3; ++a)
	{
		const unsigned int N = numLines[a];
		std::vector<double> prim(N), dual(N);
		for (unsigned int i=0; i+1<N; ++i)
			prim[i] = Op->GetDiscLine(a, i+1) - Op->GetDiscLine(a, i);
		prim[N-1] = N>1 ? prim[N-2] : 1;
		for (unsigned int i=0; i<N; ++i)
		{
			double lo = (i>0) ? prim[i-1] : prim[0];
			double hi = (i+1<N) ? prim[i] : prim[N-1];
			dual[i] = 0.5*(lo + hi);
		}
		for (int n=0; n<3; ++n)
		{
			std::vector<double>& gv = m_G[0][n][a];
			std::vector<double>& gi = m_G[1][n][a];
			gv.resize(N);
			gi.resize(N);
			for (unsigned int i=0; i<N; ++i)
			{
				gv[i] = (a==n) ? prim[i] : 1.0/dual[i];
				gi[i] = (a==n) ? dual[i] : 1.0/prim[i];
			}
		}
	}
}

double Engine_CPU::CellKey::M() const
{
	if (cls==Engine_CPU_Run::ZERO || sign==0)
		return 0;
	return sign*std::exp(q*M_TOLERANCE);
}

Engine_CPU::CellKey Engine_CPU::GetCellKey(int f, int n, const unsigned int pos[3]) const
{
	CellKey key;
	unsigned int p[3] = {pos[0], pos[1], pos[2]};
	const float s = (f==0) ? Op->GetVV(n, p) : Op->GetII(n, p);
	const float c = (f==0) ? Op->GetVI(n, p) : Op->GetIV(n, p);
	key.vv = s;
	key.q = 0;
	key.sign = 0;
	unsigned int loc[3];
	if (FindUPML(pos, loc))
	{
		key.cls = Engine_CPU_Run::UPML;   // per cell coefficients, one run through the PML
		key.vv = 0;
		return key;
	}
	if ((s==0) && (c==0))
	{
		key.cls = Engine_CPU_Run::ZERO;
		return key;
	}
	key.cls = (s==1.0f) ? Engine_CPU_Run::LOSSLESS : Engine_CPU_Run::LOSSY;
	const double g = m_G[f][n][0][pos[0]] * m_G[f][n][1][pos[1]] * m_G[f][n][2][pos[2]];
	const double m = c/g;
	if (m!=0)
	{
		key.sign = (m>0) ? 1 : -1;
		key.q = std::llround(std::log(std::fabs(m))/M_TOLERANCE);
	}
	return key;
}

void Engine_CPU::ChooseLayout()
{
	// the number of runs along each axis, over all components of both fields
	double runs[3] = {0, 0, 0};
	const size_t N0 = numLines[0], N1 = numLines[1], N2 = numLines[2];
	std::vector<CellKey> keys(N0*N1*N2);
	unsigned int pos[3];
	for (int f=0; f<2; ++f)
		for (int n=0; n<3; ++n)
		{
			for (pos[0]=0; pos[0]<N0; ++pos[0])
				for (pos[1]=0; pos[1]<N1; ++pos[1])
					for (pos[2]=0; pos[2]<N2; ++pos[2])
						keys[(pos[0]*N1 + pos[1])*N2 + pos[2]] = GetCellKey(f, n, pos);
			const size_t stride[3] = {N1*N2, N2, 1};
			for (int a=0; a<3; ++a)
				for (size_t i=0; i<keys.size(); ++i)
				{
					const size_t p_a = (i/stride[a]) % numLines[a];
					if ((p_a==0) || !(keys[i]==keys[i-stride[a]]))
						runs[a] += 1;
				}
		}

	m_FastAxis = 2;
	for (int a=0; a<2; ++a)
		if ((runs[a]<runs[m_FastAxis]) || ((runs[a]==runs[m_FastAxis]) && (numLines[a]>numLines[m_FastAxis])))
			m_FastAxis = a;
	const char* env = std::getenv("OPENEMS_CPU_FAST_AXIS");
	if (env && env[0])
	{
		const char c = env[0];
		if ((c>='x') && (c<='z'))
			m_FastAxis = c-'x';
		else
			throw std::invalid_argument("Engine_CPU: OPENEMS_CPU_FAST_AXIS must be x, y or z");
	}

	// the threads split the slow axis with more lines
	const int a1 = (m_FastAxis+1)%3, a2 = (m_FastAxis+2)%3;
	m_SplitAxis = (numLines[a1]>=numLines[a2]) ? a1 : a2;
	m_MidAxis = 3 - m_FastAxis - m_SplitAxis;
	m_Layout = std::string("N") + AXIS_LETTER[m_SplitAxis] + AXIS_LETTER[m_MidAxis] + AXIS_LETTER[m_FastAxis];

	if (g_settings.GetVerboseLevel()>0)
	{
		const double cells = 6.0*N0*N1*N2;
		cout << "Engine_CPU: mean run length along x, y, z: " << cells/runs[0] << ", " << cells/runs[1]
		     << ", " << cells/runs[2] << "; fast axis " << (char)('x'+m_FastAxis) << ", layout " << m_Layout << endl;
	}
}

void Engine_CPU::BuildRuns()
{
	const int F = m_FastAxis, S = m_SplitAxis, M = m_MidAxis;
	size_t count[4] = {0, 0, 0, 0};   // cells per class
	for (int f=0; f<2; ++f)
		for (int n=0; n<3; ++n)
		{
			Engine_CPU_Runs& r = m_Runs[f][n];
			Engine_CPU_UPML& u = m_UPMLCells[f][n];
			r.runs.clear();
			u.Clear();
			r.first.assign((size_t)numLines[S]*numLines[M] + 1, 0);
			m_Gf[f][n].assign(m_G[f][n][F].begin(), m_G[f][n][F].end());
			// the currents are updated on [0, numLines-1) per axis
			const unsigned int last = (f==0) ? 0 : 1;
			unsigned int pos[3];
			for (pos[S]=0; pos[S]<numLines[S]; ++pos[S])
				for (pos[M]=0; pos[M]<numLines[M]; ++pos[M])
				{
					const size_t line = (size_t)pos[S]*numLines[M] + pos[M];
					r.first[line] = r.runs.size();
					if ((pos[S]+last>=numLines[S]) || (pos[M]+last>=numLines[M]))
						continue;
					const double slow = m_G[f][n][S][pos[S]] * m_G[f][n][M][pos[M]];
					CellKey prev = CellKey();
					for (pos[F]=0; pos[F]+last<numLines[F]; ++pos[F])
					{
						CellKey key = GetCellKey(f, n, pos);
						++count[key.cls];
						if (key.cls==Engine_CPU_Run::UPML)
						{
							unsigned int loc[3];
							const Operator_Ext_UPML* upml = FindUPML(pos, loc);
							Operator_Ext_UPML* w = const_cast<Operator_Ext_UPML*>(upml);   // its Get*() are not const
							unsigned int p[3] = {pos[0], pos[1], pos[2]};
							u.vv.push_back((f==0) ? Op->GetVV(n, p) : Op->GetII(n, p));
							u.vi.push_back((f==0) ? Op->GetVI(n, p) : Op->GetIV(n, p));
							u.uvv.push_back((f==0) ? w->GetVV(n, loc)   : w->GetII(n, loc));
							u.ufo.push_back((f==0) ? w->GetVVFO(n, loc) : w->GetIIFO(n, loc));
							u.ufn.push_back((f==0) ? w->GetVVFN(n, loc) : w->GetIIFN(n, loc));
							u.flux.push_back(0);
						}
						if ((pos[F]>0) && (key==prev))
						{
							++r.runs.back().stop;
							continue;
						}
						Engine_CPU_Run run;
						run.start = pos[F];
						run.stop = pos[F]+1;
						run.cls = key.cls;
						run.vv = (key.cls==Engine_CPU_Run::LOSSLESS) ? 1.0f : key.vv;
						run.a = (float)(key.M()*slow);
						run.u = (key.cls==Engine_CPU_Run::UPML) ? u.flux.size()-1 : 0;
						r.runs.push_back(run);
						prev = key;
					}
				}
			r.first.back() = r.runs.size();
		}

	if (g_settings.GetVerboseLevel()>0)
	{
		size_t runs = 0;
		for (int f=0; f<2; ++f)
			for (int n=0; n<3; ++n)
				runs += m_Runs[f][n].runs.size();
		const double cells = count[0] + count[1] + count[2] + count[3];
		cout << "Engine_CPU: " << runs << " runs, mean length " << cells/runs << "; cells zero "
		     << 100.0*count[0]/cells << " %, lossless " << 100.0*count[1]/cells << " %, lossy "
		     << 100.0*count[2]/cells << " %, upml " << 100.0*count[3]/cells << " %" << endl;
	}
}

// ----------------------------------------------------------------------------
// the kernels: per line of the fast axis in the box, per component, its runs
// clipped to the box. The curl of a cell is c(k) = A[k] - A[k-oA] - B[k] + B[k-oB]
// (voltages, the neighbours below) or A[k] - A[k+oA] - B[k] + B[k+oB] (currents,
// above); along the fast axis the offset is 1.

namespace
{
//! the cells [k0, k1) of a run; the curl at k reads A, B at k and k+dA, k+dB
inline void UpdateRun(const Engine_CPU_Run& run, const Engine_CPU_UPML& u, unsigned int k0, unsigned int k1,
                      FDTD_FLOAT* __restrict X, const FDTD_FLOAT* __restrict A, const FDTD_FLOAT* __restrict B,
                      ptrdiff_t dA, ptrdiff_t dB, const float* __restrict gf, float* __restrict flux)
{
	const float a = run.a, vv = run.vv;
	switch (run.cls)
	{
	case Engine_CPU_Run::ZERO:
		for (unsigned int k=k0; k<k1; ++k)
			X[k] = 0;
		break;
	case Engine_CPU_Run::LOSSLESS:
		for (unsigned int k=k0; k<k1; ++k)
			X[k] += a*gf[k]*(A[k] - A[k+dA] - B[k] + B[k+dB]);
		break;
	case Engine_CPU_Run::LOSSY:
		for (unsigned int k=k0; k<k1; ++k)
			X[k] = vv*X[k] + a*gf[k]*(A[k] - A[k+dA] - B[k] + B[k+dB]);
		break;
	case Engine_CPU_Run::UPML:
	{
		// the extension's pre pass, the main update on the flux and its post pass in one
		const ptrdiff_t j0 = (ptrdiff_t)run.u - (ptrdiff_t)run.start;
		const float* __restrict mvv = u.vv.data() + j0;
		const float* __restrict mvi = u.vi.data() + j0;
		const float* __restrict uvv = u.uvv.data() + j0;
		const float* __restrict ufo = u.ufo.data() + j0;
		const float* __restrict ufn = u.ufn.data() + j0;
		float* __restrict F = flux + j0;
		for (unsigned int k=k0; k<k1; ++k)
		{
			const float Fn = mvv[k]*F[k] + mvi[k]*(A[k] - A[k+dA] - B[k] + B[k+dB]);
			X[k] = uvv[k]*X[k] - ufo[k]*F[k] + ufn[k]*Fn;
			F[k] = Fn;
		}
		break;
	}
	}
}
}

void Engine_CPU::UpdateVoltages(const Engine_Box& box)
{
	const int F = m_FastAxis, S = m_SplitAxis, M = m_MidAxis;
	ArrayLib::ArrayNIJK<FDTD_FLOAT>& volt = *volt_ptr;
	ArrayLib::ArrayNIJK<FDTD_FLOAT>& curr = *curr_ptr;
	const size_t sN = volt.stride(0);
	const size_t s[3] = {volt.stride(1), volt.stride(2), volt.stride(3)};
	unsigned int pos[3];
	pos[F] = 0;
	for (pos[S]=box.start[S]; pos[S]<box.stop[S]; ++pos[S])
		for (pos[M]=box.start[M]; pos[M]<box.stop[M]; ++pos[M])
		{
			const size_t base = pos[S]*s[S] + pos[M]*s[M];
			const size_t line = (size_t)pos[S]*numLines[M] + pos[M];
			for (int n=0; n<3; ++n)
			{
				const int nP = (n+1)%3, nPP = (n+2)%3;
				// V_n += vi * (I_nPP(pos) - I_nPP(pos-e_nP) - I_nP(pos) + I_nP(pos-e_nPP)); at line 0 the neighbour is the cell itself
				const ptrdiff_t dA = -(ptrdiff_t)((nP==F) ? 1 : (pos[nP]>0 ? s[nP] : 0));
				const ptrdiff_t dB = -(ptrdiff_t)((nPP==F) ? 1 : (pos[nPP]>0 ? s[nPP] : 0));
				FDTD_FLOAT* V = volt.data() + n*sN + base;
				const FDTD_FLOAT* Ia = curr.data() + nPP*sN + base;
				const FDTD_FLOAT* Ib = curr.data() + nP*sN + base;
				const Engine_CPU_Runs& r = m_Runs[0][n];
				Engine_CPU_UPML& u = m_UPMLCells[0][n];
				for (uint64_t i=r.first[line]; i<r.first[line+1]; ++i)
				{
					const Engine_CPU_Run& run = r.runs[i];
					unsigned int k0 = std::max(run.start, box.start[F]);
					const unsigned int k1 = std::min(run.stop, box.stop[F]);
					if (k0>=k1)
						continue;
					if ((k0==0) && (run.cls!=Engine_CPU_Run::ZERO) && ((nP==F) || (nPP==F)))
					{
						// the first line of the fast axis has no lower neighbour along it
						UpdateRun(run, u, 0, 1, V, Ia, Ib, (nP==F) ? 0 : dA, (nPP==F) ? 0 : dB, m_Gf[0][n].data(), u.flux.data());
						k0 = 1;
					}
					UpdateRun(run, u, k0, k1, V, Ia, Ib, dA, dB, m_Gf[0][n].data(), u.flux.data());
				}
			}
		}
}

void Engine_CPU::UpdateCurrents(const Engine_Box& box)
{
	const int F = m_FastAxis, S = m_SplitAxis, M = m_MidAxis;
	ArrayLib::ArrayNIJK<FDTD_FLOAT>& volt = *volt_ptr;
	ArrayLib::ArrayNIJK<FDTD_FLOAT>& curr = *curr_ptr;
	const size_t sN = curr.stride(0);
	const size_t s[3] = {curr.stride(1), curr.stride(2), curr.stride(3)};
	unsigned int pos[3];
	for (pos[S]=box.start[S]; pos[S]<box.stop[S]; ++pos[S])
		for (pos[M]=box.start[M]; pos[M]<box.stop[M]; ++pos[M])
		{
			const size_t base = pos[S]*s[S] + pos[M]*s[M];
			const size_t line = (size_t)pos[S]*numLines[M] + pos[M];
			for (int n=0; n<3; ++n)
			{
				const int nP = (n+1)%3, nPP = (n+2)%3;
				// I_n += iv * (V_nPP(pos) - V_nPP(pos+e_nP) - V_nP(pos) + V_nP(pos+e_nPP))
				const ptrdiff_t dA = (nP==F) ? 1 : s[nP];
				const ptrdiff_t dB = (nPP==F) ? 1 : s[nPP];
				FDTD_FLOAT* I = curr.data() + n*sN + base;
				const FDTD_FLOAT* Va = volt.data() + nPP*sN + base;
				const FDTD_FLOAT* Vb = volt.data() + nP*sN + base;
				const Engine_CPU_Runs& r = m_Runs[1][n];
				Engine_CPU_UPML& u = m_UPMLCells[1][n];
				for (uint64_t i=r.first[line]; i<r.first[line+1]; ++i)
				{
					const Engine_CPU_Run& run = r.runs[i];
					const unsigned int k0 = std::max(run.start, box.start[F]);
					const unsigned int k1 = std::min(run.stop, box.stop[F]);
					if (k0<k1)
						UpdateRun(run, u, k0, k1, I, Va, Vb, dA, dB, m_Gf[1][n].data(), u.flux.data());
				}
			}
		}
}

void Engine_CPU::UpdateVoltages(unsigned int startX, unsigned int numX)
{
	Engine_Box box = {{startX, 0, 0}, {startX+numX, numLines[1], numLines[2]}};
	UpdateVoltages(box);
}

void Engine_CPU::UpdateCurrents(unsigned int startX, unsigned int numX)
{
	Engine_Box box = {{startX, 0, 0}, {startX+numX, numLines[1]-1, numLines[2]-1}};
	UpdateCurrents(box);
}
