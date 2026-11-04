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
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <map>
#include <tuple>
#if defined(__unix__) || defined(__APPLE__)
#include <unistd.h>
#endif
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
	// the UPML cells hold their coefficients now: free the extension's (24 bytes per node of a component pair)
	for (const Operator_Ext_UPML* upml : m_UPML)
	{
		Operator_Ext_UPML* w = const_cast<Operator_Ext_UPML*>(upml);
		w->vv.Reset();
		w->vvfo.Reset();
		w->vvfn.Reset();
		w->ii.Reset();
		w->iifo.Reset();
		w->iifn.Reset();
	}
	InitKernels();
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
	if (cls==Engine_CPU_Set::ZERO || sign==0)
		return 0;
	return sign*std::exp(q*M_TOLERANCE);
}

uint32_t Engine_CPU::CellKey::Hash() const
{
	// FNV-1a over the fields; a collision only merges two runs in the statistics
	uint32_t h = 2166136261u;
	auto mix = [&h](const void* p, size_t n) {
		for (size_t i=0; i<n; ++i)
			h = (h ^ ((const unsigned char*)p)[i]) * 16777619u;
	};
	mix(&cls, sizeof(cls));
	mix(&vv, sizeof(vv));
	mix(&q, sizeof(q));
	mix(&sign, sizeof(sign));
	mix(&uvv, sizeof(uvv));
	mix(&ufo, sizeof(ufo));
	mix(&ufn, sizeof(ufn));
	return h;
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
	key.uvv = key.ufo = key.ufn = 0;
	unsigned int loc[3];
	const Operator_Ext_UPML* upml = FindUPML(pos, loc);
	if (upml)
	{
		// the main coefficients (rewritten by the UPML for its flux) and the UPML ones, graded with the depth
		Operator_Ext_UPML* w = const_cast<Operator_Ext_UPML*>(upml);   // its Get*() are not const
		key.cls = Engine_CPU_Set::UPML;
		key.uvv = (f==0) ? w->GetVV(n, loc)   : w->GetII(n, loc);
		key.ufo = (f==0) ? w->GetVVFO(n, loc) : w->GetIIFO(n, loc);
		key.ufn = (f==0) ? w->GetVVFN(n, loc) : w->GetIIFN(n, loc);
	}
	else if ((s==0) && (c==0))
	{
		key.cls = Engine_CPU_Set::ZERO;
		return key;
	}
	else
		key.cls = (s==1.0f) ? Engine_CPU_Set::LOSSLESS : Engine_CPU_Set::LOSSY;
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
	std::vector<uint32_t> keys(N0*N1*N2);   // 4 bytes per node, the peak of the setup
	unsigned int pos[3];
	for (int f=0; f<2; ++f)
		for (int n=0; n<3; ++n)
		{
			for (pos[0]=0; pos[0]<N0; ++pos[0])
				for (pos[1]=0; pos[1]<N1; ++pos[1])
					for (pos[2]=0; pos[2]<N2; ++pos[2])
						keys[(pos[0]*N1 + pos[1])*N2 + pos[2]] = GetCellKey(f, n, pos).Hash();
			const size_t stride[3] = {N1*N2, N2, 1};
			for (int a=0; a<3; ++a)
				for (size_t i=0; i<keys.size(); ++i)
				{
					const size_t p_a = (i/stride[a]) % numLines[a];
					if ((p_a==0) || (keys[i]!=keys[i-stride[a]]))
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
			r.sets.clear();
			u.Clear();
			// the index of each distinct key in r.sets
			typedef std::tuple<int, float, int64_t, int, float, float, float> SetKey;
			std::map<SetKey, uint32_t> setIndex;
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
					CellKey prev = CellKey();
					for (pos[F]=0; pos[F]+last<numLines[F]; ++pos[F])
					{
						CellKey key = GetCellKey(f, n, pos);
						++count[key.cls];
						if (key.cls==Engine_CPU_Set::UPML)
							u.flux.push_back(0);
						if ((pos[F]>0) && (key==prev))
						{
							++r.runs.back().stop;
							continue;
						}
						const SetKey sk(key.cls, key.vv, key.q, key.sign, key.uvv, key.ufo, key.ufn);
						auto it = setIndex.find(sk);
						if (it==setIndex.end())
						{
							Engine_CPU_Set set;
							set.cls = key.cls;
							set.m = key.M();
							set.vv = (key.cls==Engine_CPU_Set::LOSSLESS) ? 1.0f : key.vv;
							set.uvv = key.uvv;
							set.ufo = key.ufo;
							set.ufn = key.ufn;
							it = setIndex.insert(std::make_pair(sk, (uint32_t)r.sets.size())).first;
							r.sets.push_back(set);
						}
						Engine_CPU_Run run;
						run.start = pos[F];
						run.stop = pos[F]+1;
						run.set = it->second;
						if (u.flux.size() > 0xffffffffu)
							throw std::length_error("Engine_CPU: more than 2^32 UPML cells in a component");
						run.u = (key.cls==Engine_CPU_Set::UPML) ? (uint32_t)(u.flux.size()-1) : 0;
						r.runs.push_back(run);
						prev = key;
					}
				}
			r.first.back() = r.runs.size();
		}

	if (g_settings.GetVerboseLevel()>0)
	{
		size_t runs = 0, sets = 0;
		for (int f=0; f<2; ++f)
			for (int n=0; n<3; ++n)
			{
				runs += m_Runs[f][n].runs.size();
				sets += m_Runs[f][n].sets.size();
			}
		const double cells = count[0] + count[1] + count[2] + count[3];
		cout << "Engine_CPU: " << runs << " runs (" << runs*sizeof(Engine_CPU_Run)/1048576.0 << " MiB), " << sets
		     << " coefficient sets, mean length " << cells/runs << "; cells zero "
		     << 100.0*count[0]/cells << " %, lossless " << 100.0*count[1]/cells << " %, lossy "
		     << 100.0*count[2]/cells << " %, upml " << 100.0*count[3]/cells << " %" << endl;
	}
}

// ----------------------------------------------------------------------------
// the kernels are in engine_cpu_kernels.inc, compiled per instruction set

void Engine_CPU::InitKernels()
{
	Engine_CPU_Kernel& k = m_Kernel;
	k.F = m_FastAxis;
	k.S = m_SplitAxis;
	k.M = m_MidAxis;
	for (int a=0; a<3; ++a)
	{
		k.numLines[a] = numLines[a];
		k.s[a] = volt_ptr->stride(1+a);
	}
	k.sN = volt_ptr->stride(0);
	k.field[0] = volt_ptr->data();
	k.field[1] = curr_ptr->data();
	for (int f=0; f<2; ++f)
		for (int n=0; n<3; ++n)
		{
			Engine_CPU_Kernel::Component& c = k.comp[f][n];
			c.runs = m_Runs[f][n].runs.data();
			c.first = m_Runs[f][n].first.data();
			c.sets = m_Runs[f][n].sets.data();
			c.flux = m_UPMLCells[f][n].flux.data();
			c.gf = m_Gf[f][n].data();
			c.gS = m_G[f][n][m_SplitAxis].data();
			c.gM = m_G[f][n][m_MidAxis].data();
		}

	// the instruction set: AVX2 and FMA if the CPU has both (OPENEMS_CPU_ISA=generic overrides)
	m_ISA = "generic";
	m_UpdateVoltages = Engine_CPU_Generic::UpdateVoltages;
	m_UpdateCurrents = Engine_CPU_Generic::UpdateCurrents;
	const char* env = std::getenv("OPENEMS_CPU_ISA");
	const std::string want = env ? env : "";
	if ((want!="") && (want!="generic") && (want!="avx2"))
		throw std::invalid_argument("Engine_CPU: OPENEMS_CPU_ISA must be generic or avx2");
#ifdef OPENEMS_CPU_AVX2
	if (want!="generic")
	{
		__builtin_cpu_init();
		if (__builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma"))
		{
			m_ISA = "avx2";
			m_UpdateVoltages = Engine_CPU_AVX2::UpdateVoltages;
			m_UpdateCurrents = Engine_CPU_AVX2::UpdateCurrents;
		}
		else if (want=="avx2")
			throw std::runtime_error("Engine_CPU: OPENEMS_CPU_ISA=avx2, but the CPU lacks AVX2 or FMA");
	}
#else
	if (want=="avx2")
		throw std::runtime_error("Engine_CPU: OPENEMS_CPU_ISA=avx2, but this build has no AVX2 kernels");
#endif
	if (g_settings.GetVerboseLevel()>0)
		cout << "Engine_CPU: kernels for " << m_ISA << endl;
}

// the bytes of the last level cache, 32 MiB if the system does not say
static size_t CacheBytes()
{
#ifdef _SC_LEVEL3_CACHE_SIZE
	long l3 = sysconf(_SC_LEVEL3_CACHE_SIZE);
	if (l3>0)
		return l3;
#endif
	return 32u<<20;
}

bool Engine_CPU::BuildBlockedPlans(unsigned int numThreads, std::vector<Engine_BlockedPlan>& plans)
{
	const char* env = std::getenv("OPENEMS_CPU_TBLOCK");
	std::string why;
	unsigned int k = 0, W = 0;
	if (env && (std::string(env)=="0"))
		why = "OPENEMS_CPU_TBLOCK=0";
	else if (m_SubGrid)
		why = "a cylindrical multi-grid";
	for (size_t n=0; n<m_Eng_exts.size() && why.empty(); ++n)
		if (!m_Eng_exts[n]->SupportsBoxes())
			why = "the " + m_Eng_exts[n]->GetExtensionName() + " has no box hooks";
	if (why.empty() && env && env[0])
	{
		if (sscanf(env, "%u:%u", &k, &W)<1 || k==0)
			throw std::invalid_argument("Engine_CPU: OPENEMS_CPU_TBLOCK must be 0, k or k:W");
		if (W==0)
			W = 2*k+1;
	}
	else if (why.empty())
	{
		// the tile: the cache for the fields, the UPML flux and the runs of its planes (measured: fuller tiles are faster than half the cache)
		size_t runs = 0, flux = 0;
		for (int f=0; f<2; ++f)
			for (int n=0; n<3; ++n)
			{
				runs += m_Runs[f][n].runs.size();
				flux += m_UPMLCells[f][n].flux.size();
			}
		const double plane = (double)numLines[m_MidAxis]*numLines[m_FastAxis];
		const double bytes = plane*24.0 + (4.0*flux + sizeof(Engine_CPU_Run)*runs)/numLines[m_SplitAxis];
		W = (unsigned int)(CacheBytes()/bytes);
		k = std::min(8u, (W>0) ? (W-1)/2 : 0);
		if (k<2)
			why = "a plane is too large for the cache";
		else if (numLines[m_SplitAxis] < 2*W)
			why = "the model fits the cache";
	}
	if (!why.empty())
	{
		if (g_settings.GetVerboseLevel()>0)
			cout << "Engine_CPU: no temporal blocking: " << why << endl;
		return false;
	}
	// k, then halves of it down to 2: the processing between two calls of IterateTS() may leave fewer timesteps
	for (unsigned int kk=k; kk>=2; kk/=2)
		plans.push_back(Engine_BlockedPlan::Trapezoids(numLines, numThreads, m_SplitAxis, m_MidAxis, W, kk));
	if (g_settings.GetVerboseLevel()>0)
		cout << "Engine_CPU: temporal blocking, up to " << k << " timesteps per block, tiles of " << W << " lines along "
		     << (char)('x'+m_SplitAxis) << ", " << plans[0].phases.size() << " phases" << endl;
	return true;
}

void Engine_CPU::UpdateVoltages(const Engine_Box& box)
{
	m_UpdateVoltages(m_Kernel, box.start, box.stop);
}

void Engine_CPU::UpdateCurrents(const Engine_Box& box)
{
	m_UpdateCurrents(m_Kernel, box.start, box.stop);
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
