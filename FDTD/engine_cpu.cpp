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
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>
#include <thread>
#include <typeinfo>
#include <utility>
#ifdef OPENEMS_CPU_AVX2
#if defined(_MSC_VER) && !defined(__clang__)
#include <intrin.h>   // __cpuid(), __cpuidex(), _xgetbv()
#else
#include <cpuid.h>    // GCC and Clang, clang-cl too
#endif
#endif
#if defined(__unix__) || defined(__APPLE__)
#include <unistd.h>
#endif
#if defined(__linux__)
#include <sched.h>    // sched_getaffinity()
#endif
#if defined(__APPLE__)
#include <sys/sysctl.h>
#endif
#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#endif

#include "tools/global.h"
#include "tools/useful.h"
#include "extensions/engine_extension.h"
#include "extensions/engine_ext_upml.h"
#include "extensions/engine_ext_lorentzmaterial.h"
#include "extensions/engine_ext_debyematerial.h"
#include "extensions/engine_ext_absorbing_bc.h"
#include "extensions/operator_ext_absorbing_bc.h"
#include "extensions/operator_ext_cylinder.h"

using std::cout;
using std::endl;

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
	m_GatherBudget = 0;
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
	// OPENEMS_CPU_INERT=0: update the inert currents too
	const char* inert = std::getenv("OPENEMS_CPU_INERT");
	m_Model.reset(new CellModel(Op));
	m_Model->Analyse(!(inert && (std::string(inert)=="0")));
	ChooseLayout();
	Engine::Init();   // the fields in m_Layout (GetFieldLayout()), the extensions
	DropFusedExtensions();
	m_Model->Build(m_FastAxis, m_SplitAxis, m_MidAxis);
	if (Op->GetDebugFlags() & Operator::debugOperator)
		m_Model->Dump("cellmodel_dump");
	m_Model->FreeFused();
	for (int f=0; f<2; ++f)
		for (int n=0; n<3; ++n)
		{
			m_State[f][n].assign(m_Model->Component(f, n).state, 0);
			const std::vector<double>& g = m_Model->G(f, n, m_FastAxis);
			m_Gf[f][n].assign(g.begin(), g.end());
		}
	m_GatherBudget = volt_ptr->bytes() + curr_ptr->bytes();
	for (int f=0; f<2; ++f)
		for (int n=0; n<3; ++n)
			m_GatherBudget += m_State[f][n].size()*sizeof(float);
	m_GatherBudget = std::max(m_GatherBudget, (size_t)256<<20);
	InitKernels();
}

void Engine_CPU::DropFusedExtensions()
{
	// the engine updates their cells itself (the classes UPML, UPML_PROFILE, DISP)
	std::vector<Engine_Extension*> keep;
	for (size_t n=0; n<m_Eng_exts.size(); ++n)
	{
		Engine_Extension* e = m_Eng_exts[n];
		if (!m_Model->IsFused(e->GetOperatorExtension()))
		{
			keep.push_back(e);
			continue;
		}
		const std::type_info& type = typeid(*e);
		if ((type!=typeid(Engine_Ext_UPML)) && (type!=typeid(Engine_Ext_LorentzMaterial)) && (type!=typeid(Engine_Ext_DebyeMaterial)))
			throw std::logic_error("Engine_CPU: an unknown engine extension of the fused " + e->GetExtensionName());
		delete e;
	}
	m_Eng_exts = keep;
}

void Engine_CPU::Reset()
{
	Engine::Reset();
	m_Model.reset();
	for (int f=0; f<2; ++f)
		for (int n=0; n<3; ++n)
		{
			m_State[f][n].clear();
			m_Gf[f][n].clear();
		}
}

void Engine_CPU::ChooseLayout()
{
	m_FastAxis = m_Model->RecommendFastAxis();
	// the box plans: a closed alpha mesh couples the first and the last alpha plane (Engine_Ext_Cylinder),
	// a local absorber its sheet and the next line along its normal within a half step; as the fast
	// axis, which every box spans, the blocked plans can run them (BuildBlockedPlans())
	bool alpha = false;
	int normal = -1;
	for (size_t n=0; n<Op->GetNumberOfExtentions(); ++n)
	{
		const Operator_Extension* e = Op->GetExtension(n);
		const Operator_Ext_Cylinder* cyl = dynamic_cast<const Operator_Ext_Cylinder*>(e);
		if (cyl && cyl->ClosedAlpha())
			alpha = true;
		const Operator_Ext_Absorbing_BC* abc = dynamic_cast<const Operator_Ext_Absorbing_BC*>(e);
		if (abc && (normal<0))
			normal = abc->GetNormalAxis();
		else if (abc && (normal!=abc->GetNormalAxis()))
			normal = 3;   // several: no blocking
	}
	if (alpha)
		m_FastAxis = 1;
	else if ((normal>=0) && (normal<3))   // the cylinder's alpha first: its box hooks need it
		m_FastAxis = normal;
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
		cout << "Engine_CPU: fast axis " << (char)('x'+m_FastAxis) << ", layout " << m_Layout << endl;
}

// ----------------------------------------------------------------------------
// the kernels are in engine_cpu_kernels.inc, compiled per instruction set

#ifdef OPENEMS_CPU_AVX2
// AVX2 and FMA, and the AVX registers saved by the system: from CPUID and XGETBV directly.
// __builtin_cpu_supports() needs __cpu_model of the compiler runtime, which e.g. clang
// on Windows does not link (x86-64 only, see CMakeLists.txt)
static bool CpuHasAVX2FMA()
{
	unsigned int b, c;   // CPUID leaf 1 ECX, leaf 7 EBX
	unsigned long long xcr0;
#if defined(_MSC_VER) && !defined(__clang__)
	int r[4];
	__cpuid(r, 0);
	if (r[0]<7)
		return false;
	__cpuid(r, 1);
	c = (unsigned int)r[2];
	if (!((c>>27) & 1))   // OSXSAVE: XGETBV available
		return false;
	xcr0 = _xgetbv(0);
	__cpuidex(r, 7, 0);
	b = (unsigned int)r[1];
#else
	unsigned int a, d, c7;
	if (!__get_cpuid(1, &a, &b, &c, &d))
		return false;
	if (!((c>>27) & 1))   // OSXSAVE: XGETBV available
		return false;
	unsigned int lo, hi;
	__asm__ volatile ("xgetbv" : "=a"(lo), "=d"(hi) : "c"(0));
	xcr0 = ((unsigned long long)hi<<32) | lo;
	if (!__get_cpuid_count(7, 0, &a, &b, &c7, &d))   // c keeps leaf 1
		return false;
#endif
	const bool fma = (c>>12) & 1, avx = (c>>28) & 1;
	return fma && avx && ((xcr0 & 6)==6)   // the SSE and the AVX state saved
	       && ((b>>5) & 1);                // AVX2
}
#endif

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
			const CellModel_Component& mc = m_Model->Component(f, n);
			c.runs = mc.runs.data();
			c.first = mc.first.data();
			c.sets = mc.sets.data();
			c.state = m_State[f][n].data();
			c.poles = mc.poles.data();
			for (int q=0; q<5; ++q)
				c.profile[q] = mc.profile[q].data();
			c.gf = m_Gf[f][n].data();
			c.gS = m_Model->G(f, n, m_SplitAxis).data();
			c.gM = m_Model->G(f, n, m_MidAxis).data();
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
		if (CpuHasAVX2FMA())
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

#if defined(__linux__)
// the distinct caches of \a level of the CPUs this process may run on (shared_cpu_list in sysfs), 0 if it does not say
// (other systems: 1)
static unsigned int CacheInstancesSysfs(int level)
{
	cpu_set_t set;
	if (sched_getaffinity(0, sizeof(set), &set)!=0)
		return 1;
	std::set<std::string> distinct;
	for (int cpu=0; cpu<CPU_SETSIZE; ++cpu)
	{
		if (!CPU_ISSET(cpu, &set))
			continue;
		for (int index=0; index<8; ++index)
		{
			const std::string dir = "/sys/devices/system/cpu/cpu" + std::to_string(cpu) + "/cache/index" + std::to_string(index) + "/";
			std::ifstream fl(dir + "level"), ft(dir + "type"), fs(dir + "shared_cpu_list");
			int l = 0;
			std::string type, cpus;
			if ((fl >> l) && (ft >> type) && (fs >> cpus) && (l==level) && (type!="Instruction"))
			{
				distinct.insert(cpus);
				break;
			}
		}
	}
	return distinct.size();
}
#else
static unsigned int CacheInstancesSysfs(int) {return 1;}
#endif

// the bytes of one last level cache, where they come from, and how many such caches the machine has (\a instances,
// e.g. one per CCD of a Zen CPU or per socket); 32 MiB and one if the system does not say
static size_t CacheBytes(std::string& source, unsigned int& instances)
{
	instances = 1;
#if defined(_WIN32)
	// the largest cache of the highest level
	DWORD len = 0;
	GetLogicalProcessorInformationEx(RelationCache, NULL, &len);
	if (len>0)
	{
		std::vector<char> buf(len);
		SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX* info = (SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX*)buf.data();
		if (GetLogicalProcessorInformationEx(RelationCache, info, &len))
		{
			int level = 0;
			size_t bytes = 0;
			for (DWORD off=0; off<len; )
			{
				const SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX* e = (const SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX*)(buf.data()+off);
				if ((e->Relationship==RelationCache) && (e->Cache.Type!=CacheInstruction))
				{
					const int l = e->Cache.Level;
					if ((l>level) || ((l==level) && (e->Cache.CacheSize>bytes)))
					{
						level = l;
						bytes = e->Cache.CacheSize;
					}
				}
				off += e->Size;
			}
			if (bytes>0)
			{
				unsigned int count = 0;   // the caches of that level, each an entry with its processors
				for (DWORD off=0; off<len; )
				{
					const SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX* e = (const SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX*)(buf.data()+off);
					if ((e->Relationship==RelationCache) && (e->Cache.Type!=CacheInstruction) && (e->Cache.Level==level))
						++count;
					off += e->Size;
				}
				instances = std::max(count, 1u);
				source = "L" + std::to_string(level) + " (Windows)";
				return bytes;
			}
		}
	}
#elif defined(__APPLE__)
	// an Intel Mac has an L3; Apple silicon has none (its system cache is not the cores'): the L2 of the performance cores
	const char* names[] = {"hw.l3cachesize", "hw.perflevel0.l2cachesize", "hw.l2cachesize"};
	for (const char* name : names)
	{
		int64_t bytes = 0;
		size_t size = sizeof(bytes);
		if ((sysctlbyname(name, &bytes, &size, NULL, 0)==0) && (bytes>0))
		{
			source = name;
			// the performance cores of a larger Apple chip are several clusters, an L2 each
			int cores = 0, perL2 = 0;
			size_t s1 = sizeof(cores), s2 = sizeof(perL2);
			if ((std::string(name)=="hw.perflevel0.l2cachesize")
			    && (sysctlbyname("hw.perflevel0.physicalcpu", &cores, &s1, NULL, 0)==0)
			    && (sysctlbyname("hw.perflevel0.cpusperl2", &perL2, &s2, NULL, 0)==0) && (cores>0) && (perL2>0))
				instances = std::max(1, cores/perL2);
			return (size_t)bytes;
		}
	}
#else
#ifdef _SC_LEVEL3_CACHE_SIZE
	const long l3 = sysconf(_SC_LEVEL3_CACHE_SIZE);
	if (l3>0)
	{
		source = "L3 (sysconf)";
		instances = std::max(CacheInstancesSysfs(3), 1u);
		return l3;
	}
#endif
	// some systems (e.g. ARM) do not say through sysconf: the cache descriptions of cpu0, the highest level
	size_t best = 0;
	int bestLevel = 0;
	for (int index=0; index<8; ++index)
	{
		const std::string dir = "/sys/devices/system/cpu/cpu0/cache/index" + std::to_string(index) + "/";
		std::ifstream fl(dir + "level"), ft(dir + "type"), fs(dir + "size");
		int level = 0;
		std::string type, size;
		if (!(fl >> level) || !(ft >> type) || !(fs >> size) || (type=="Instruction") || size.empty())
			continue;
		size_t bytes = std::strtoul(size.c_str(), NULL, 10);
		const char unit = size.back();
		if ((unit=='K') || (unit=='k'))
			bytes <<= 10;
		else if ((unit=='M') || (unit=='m'))
			bytes <<= 20;
		if ((level>bestLevel) && (bytes>0))
		{
			best = bytes;
			bestLevel = level;
		}
	}
	if (best>0)
	{
		source = "L" + std::to_string(bestLevel) + " (sysfs)";
		instances = std::max(CacheInstancesSysfs(bestLevel), 1u);
		return best;
	}
#endif
	source = "assumed, the system does not say";
	return 32u<<20;
}

bool Engine_CPU::BuildBlockedPlans(unsigned int numThreads, Engine_BlockedChoice& choice)
{
	const char* env = std::getenv("OPENEMS_CPU_TBLOCK");
	std::string why;
	unsigned int k = 0, W = 0;
	bool explicitTile = false;
	std::vector< std::pair<unsigned int, unsigned int> > tiles;   // OPENEMS_CPU_TBLOCK=k:W,k:W,...
	if (env && (std::string(env)=="0"))
		why = "OPENEMS_CPU_TBLOCK=0";
	else if (m_SubGrid)
		why = "a cylindrical multi-grid";
	for (size_t n=0; n<m_Eng_exts.size() && why.empty(); ++n)
	{
		const Engine_Ext_Absorbing_BC* abc = dynamic_cast<const Engine_Ext_Absorbing_BC*>(m_Eng_exts[n]);
		if (!m_Eng_exts[n]->SupportsBoxes())
			why = "the " + m_Eng_exts[n]->GetExtensionName() + " has no box hooks";
		else if (abc && (abc->GetNormalAxis()!=m_FastAxis))
			why = "a local absorber not normal to the fast axis";
	}
	if (why.empty() && env && env[0])
	{
		// k or k:W, several separated by commas
		const char* p = env;
		while (true)
		{
			char* end = NULL;
			unsigned int kk = (unsigned int)std::strtoul(p, &end, 10), WW = 0;
			if ((end!=NULL) && (*end==':'))
				WW = (unsigned int)std::strtoul(end+1, &end, 10);
			if ((kk==0) || (end==NULL) || ((*end!='\0') && (*end!=',')))
				throw std::invalid_argument("Engine_CPU: OPENEMS_CPU_TBLOCK must be 0, or k or k:W, several separated by commas");
			tiles.push_back(std::make_pair(kk, std::max(WW, 2*kk+1)));
			if (*end=='\0')
				break;
			p = end+1;
		}
		k = tiles[0].first;
		W = tiles[0].second;
		explicitTile = true;
	}
	else if (why.empty())
	{
		// the tile: the fields, the UPML and ADE state and the runs of its planes in the cache. The timesteps
		// per block from the whole cache, the width then 3/4 of it as long as it holds them (measured: a tile
		// that fills the whole cache is much slower, the rest of the working set evicts it; but fewer
		// timesteps per block cost more than that)
		size_t runs = 0, state = 0;
		for (int f=0; f<2; ++f)
			for (int n=0; n<3; ++n)
			{
				runs += m_Model->Component(f, n).runs.size();
				state += m_State[f][n].size();
			}
		double ext = 0;   // the per-cell state of the extensions (dispersive materials)
		for (size_t n=0; n<m_Eng_exts.size(); ++n)
			ext += m_Eng_exts[n]->BoxBytes();
		const double plane = (double)numLines[m_MidAxis]*numLines[m_FastAxis];
		const double fieldBytes = 6*sizeof(FDTD_FLOAT);   // the voltages and currents of a node
		const double bytes = plane*fieldBytes + (sizeof(float)*state + sizeof(CellModel_Run)*runs + ext)/numLines[m_SplitAxis];
		// the threads keep their lines of a tile through its timesteps, so each cache holds the part of its cores:
		// the caches of the fewest the threads span (all at the default threads, all physical cores; with fewer
		// threads the scheduler may spread them wider, the trials try larger tiles)
		std::string source;
		unsigned int instances = 1;
		size_t cache = CacheBytes(source, instances);
		const unsigned int coresPer = std::max(1u, PhysicalCores()/instances);
		const unsigned int spanned = std::min(instances, (numThreads+coresPer-1)/coresPer);
		cache *= spanned;
		source += ", " + std::to_string(spanned) + " of " + std::to_string(instances) + " such caches for " + std::to_string(numThreads) + " threads";
		if (const char* mib = std::getenv("OPENEMS_CPU_CACHE"))   // to try other cache sizes
		{
			char* end = NULL;
			const double v = std::strtod(mib, &end);
			if ((end==mib) || (*end!='\0') || (v<=0))
				throw std::invalid_argument("Engine_CPU: OPENEMS_CPU_CACHE must be the cache size in MiB");
			cache = (size_t)(v*1048576.0);
			source = "OPENEMS_CPU_CACHE";
		}
		if (g_settings.GetVerboseLevel()>0)
			cout << "Engine_CPU: cache " << cache/1048576.0 << " MiB, " << source << endl;
		const unsigned int full = (unsigned int)(cache/bytes);
		k = std::min(8u, (full>0) ? (full-1)/2 : 0);
		W = std::max(2*k+1, (unsigned int)(0.75*full));
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
	for (size_t n=0; n<m_Eng_exts.size(); ++n)
		m_Eng_exts[n]->PrepareBoxes(numLines);
	// a set of plans per tile, every k down to 2: the processing between two calls of IterateTS() may leave fewer timesteps
	auto add = [&](unsigned int kmax, unsigned int width) {
		std::vector<Engine_BlockedPlan> plans;
		for (unsigned int kk=kmax; kk>=2; --kk)
			plans.push_back(Engine_BlockedPlan::Trapezoids(numLines, numThreads, m_SplitAxis, m_MidAxis, width, kk));
		choice.candidates.push_back(plans);
	};
	// the seconds per trial: OPENEMS_CPU_TRIAL, default 0.5 with the automatic tiles, none with given ones
	choice.trialSeconds = explicitTile ? 0 : 0.5;
	if (const char* trial = std::getenv("OPENEMS_CPU_TRIAL"))
	{
		char* end = NULL;
		choice.trialSeconds = std::strtod(trial, &end);
		if ((end==trial) || (*end!='\0') || (choice.trialSeconds<0))
			throw std::invalid_argument("Engine_CPU: OPENEMS_CPU_TRIAL must be the seconds per trial, 0 for none");
	}
	if (!explicitTile)
	{
		// the tile from the estimate of the cache, and as candidates tiles 4x and 2x, half and a quarter as wide:
		// the estimate is off by a cache shared or split other than it says, a working set larger than counted.
		// The thread pool searches them from the estimate (Engine_BlockedChoice; measured on one CPU: with the
		// cache assumed 4x too large or too small the search ends at the speed of the right tile)
		tiles.clear();
		const unsigned int Wmax = numLines[m_SplitAxis]/2;   // two tiles at least
		for (unsigned int m=4; m>=2; m/=2)
			if (m*W<=Wmax)
				tiles.push_back(std::make_pair(std::min(8u, (m*W-1)/2), m*W));
		choice.start = tiles.size();
		tiles.push_back(std::make_pair(k, W));
		for (unsigned int d=2; d<=4; d*=2)
		{
			const unsigned int Wd = W/d, kd = std::min(k, (Wd>0) ? (Wd-1)/2 : 0);
			if ((kd>=2) && (Wd<tiles.back().second))
				tiles.push_back(std::make_pair(kd, Wd));
		}
	}
	for (size_t t=0; t<tiles.size(); ++t)
		add(tiles[t].first, tiles[t].second);
	if (g_settings.GetVerboseLevel()>0)
	{
		cout << "Engine_CPU: temporal blocking, up to " << k << " timesteps per block, tiles of " << W << " lines along "
		     << (char)('x'+m_SplitAxis) << ", " << choice.candidates[choice.start][0].phases.size() << " phases" << endl;
		if ((choice.trialSeconds>0) && (tiles.size()>1))
		{
			cout << "Engine_CPU: tiles to try (timesteps:lines)";
			for (size_t t=0; t<tiles.size(); ++t)
				cout << " " << t+1 << ") " << tiles[t].first << ":" << tiles[t].second;
			cout << ", and the flat sweep, " << choice.trialSeconds << " s each" << endl;
		}
		else if (choice.trialSeconds>0)
			cout << "Engine_CPU: trying the tile and the flat sweep, " << choice.trialSeconds << " s each" << endl;
	}
	return true;
}

unsigned int Engine_CPU::DefaultThreads() const
{
	// all physical cores: the blocked plans keep them busy (no ramp as for the SSE engine); SMT siblings
	// share the core and its caches, measured slower (pml M: 8 cores 1200 MC/s, 16 threads 1050)
	return std::min(PhysicalCores(), AvailableThreads());
}

bool Engine_CPU::CalcFastEnergy(const unsigned int numNodes[3], double& E_energy, double& H_energy) const
{
	const unsigned int nS = numNodes[m_SplitAxis], nM = numNodes[m_MidAxis], nF = numNodes[m_FastAxis];
	const size_t sN = volt_ptr->stride(0), sS = volt_ptr->stride(1+m_SplitAxis), sM = volt_ptr->stride(1+m_MidAxis);
	const float* const V = volt_ptr->data();
	const float* const I = curr_ptr->data();
	// a sum per line of the split axis, added in line order: the same result for any number of threads
	std::vector<double> e(nS, 0), h(nS, 0);
	auto lines = [&](unsigned int s0, unsigned int s1) {
		for (unsigned int s=s0; s<s1; ++s)
		{
			double es = 0, hs = 0;
			for (int n=0; n<3; ++n)
				for (unsigned int m=0; m<nM; ++m)
				{
					const float* v = V + n*sN + s*sS + m*sM;
					const float* c = I + n*sN + s*sS + m*sM;
					for (unsigned int f=0; f<nF; ++f)
					{
						es += v[f]*v[f];   // the products in float, as the interface did
						hs += c[f]*c[f];
					}
				}
			e[s] = es;
			h[s] = hs;
		}
	};
	const unsigned int numThreads = std::max(1u, std::min(AvailableThreads(), nS));
	std::vector<std::thread> threads;
	for (unsigned int t=1; t<numThreads; ++t)
		threads.push_back(std::thread(lines, (unsigned int)((size_t)nS*t/numThreads), (unsigned int)((size_t)nS*(t+1)/numThreads)));
	lines(0, nS/numThreads);
	for (std::thread& th : threads)
		th.join();
	E_energy = H_energy = 0;
	for (unsigned int s=0; s<nS; ++s)
	{
		E_energy += e[s];
		H_energy += h[s];
	}
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
